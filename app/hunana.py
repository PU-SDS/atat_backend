import math
import time
import numpy
import random
import sys, getopt
from collections import OrderedDict, Counter
from Bio import SeqIO
from operator import itemgetter
from scipy import stats
from decimal import Decimal
import json

class Hunana:
    resultList = []
    entropy_iteration = 10
    enable_entropy = False
    random.seed()

    def processRecord(self,input="",maxIteration="",kmer=9):
        index = 0
        limit=index+kmer
        result = []
        while index <= maxIteration:
            count=0
            positionRecord = {'position':index, 'count':count}
            limit=index+kmer
            unwantedchar = set('-XBJZOU')
            for record in input:
                sequence = record.seq[index:limit]
                if any((c in unwantedchar) for c in sequence):
                    continue
                else:
                    if str(sequence) in positionRecord:
                        positionRecord[str(sequence)]=positionRecord[str(sequence)]+1
                    else:
                        positionRecord[str(sequence)]=1
                    positionRecord['count']=positionRecord['count']+1
            if self.enable_entropy:
                #normalized_entropy=self.calcNormalizedEntropy(positionRecord,self.entropy_iteration,positionRecord['count'])
                normalized_entropy=self.calcNormalizedEntropy(positionRecord,self.entropy_iteration, 10000)
                positionRecord['entropy']=normalized_entropy
            else:
                positionRecord['entropy']=0
            result.append(positionRecord)
            index+=1
        return result

    def run(self,id="",inputfile="",outputfile="",alignmentFormat="fasta",kmer=9,protein="",source="",entropy=False):
        self.enable_entropy=entropy
        filestream=open(inputfile)
        alignedSeqData=list(SeqIO.parse(filestream, alignmentFormat))
        filestream.close()
        dataLength = len(str(alignedSeqData[0].seq))
        maxIteration = dataLength - kmer
        finalResultList=self.processRecord(alignedSeqData,maxIteration)
        resultListSorted = sorted(finalResultList,key=itemgetter("position"))
        return self.printOutput(id,resultListSorted,protein,outputfile,source,kmer)

    def flattenDataSet(self,ListOfData=""):
        random.seed()
        flatDataSet=[]
        for sequence in list(ListOfData.keys()):
          if sequence == 'position' or sequence == 'count':
            continue
          else:
            for i in range(0,int(ListOfData[sequence])):
              flatDataSet.append(sequence)
        random.shuffle(flatDataSet)
        return flatDataSet

    def calcNormalizedEntropy(self,ListOfData="",iteration="",maxSample=""):
        normalizedEntropy = 0
        totalEntropy = []
        inverse = []
        flatDataSet=self.flattenDataSet(ListOfData)
        for iteration in range(iteration):
            #sample = random.choice(range(maxSample))+1
            sample = random.randint(1,maxSample)
            if sample > maxSample:
                sample=maxSample
            totalSample = 0
            sampleDic = {}
            for sampleID in range(sample):
                sequence = ''
                while True:
                    #sequence = random.choice(list(ListOfData.keys()))
                    sequence = random.choice(flatDataSet)
                    if sequence == 'position' or sequence == 'count':
                        continue
                    else:
                        totalSample += 1
                        if sequence in sampleDic:
                            sampleDic[sequence]=sampleDic[sequence]+1
                        else:
                            sampleDic[sequence]=1
                        break
            tTotalEntropy = 0
            for record in sampleDic:
                entropy = (float(sampleDic[record])/float(sample)) * (math.log2(float(sampleDic[record])/float(sample)))
                tTotalEntropy = tTotalEntropy + entropy
            totalEntropy.append(tTotalEntropy*-1)
            inverse.append(float(1)/float(sample))
            '''
            totalEntropy.append(float(round(Decimal(tTotalEntropy),2)))
            inverse.append(float(round(Decimal(1)/Decimal(sample),2)))
            '''
        try:
           slope, intercept, r_value, p_value, std_err = stats.linregress(inverse,totalEntropy)
           normalizedEntropy = intercept
        except:
            pass
        return normalizedEntropy

    def printOutput(self,id="",resultListSorted="",protein="",outputfile="",source="",kmer=9):
        index='hunana_'+source+"_"+str(id)
        result=[]
        for recordDic in resultListSorted:
            recordForTheFollowingPosition={}
            #oRecordDic = OrderedDict(recordDic.items)
            if recordDic:
                #position = recordDic['position']+int((numpy.median(range(1,kmer+1))))
                position = recordDic['position']
                count = recordDic['count']
                entropy = recordDic['entropy']
                if count > 0:
                    recordDicSorted = sorted(recordDic.items(), key=itemgetter(1),reverse=True)
                    firstRecord=True
                    sequences={}
                    for key,value in recordDicSorted:
                        if key == "position" or key == "count" or key == "entropy":
                            continue
                        percentage_conservation=float(value)*100/float(count)
                        if firstRecord:
                            recordForTheFollowingPosition['position']=position
                            #recordForTheFollowingPosition['variants']=len(recordDic)-4
                            recordForTheFollowingPosition['variants']=len(recordDic)
                            recordForTheFollowingPosition['supports']=count
                            recordForTheFollowingPosition['entropy']=entropy
                            sequences[key]=[str(value),str(percentage_conservation)]
                            firstRecord=False
                        else:
                            sequences[key]=[str(value),str(percentage_conservation)]
                    recordForTheFollowingPosition['sequences']=sequences
                else:
                    recordForTheFollowingPosition['sequences']="NOT_ANALYZED"
                result.append(recordForTheFollowingPosition)
        return result
