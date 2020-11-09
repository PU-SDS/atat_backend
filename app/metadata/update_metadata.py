import requests
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import Iterable, Tuple, Union
from io import StringIO

from app.metadata import UpdateMetaConsts


class UpdateMetadata(object):
    def __init__(self, sequences: Iterable[SeqRecord]):
        """
        Updates the headers of all the sequences provided and returns an iterable of type SeqRecord

        Args:
            sequences: An iterable of type SeqRecord
        """

        self.sequences = sequences

    def update(self) -> Iterable[SeqRecord]:
        for seq in self.sequences:
            strain = ''.join([x for x in seq.description.split('|') if '/' in x])
            strain = re.search('A/.*', strain).group()

            results = self._try_get_seqs_from_strain(strain=strain)

            if not results:
                print(f'Skipping: {strain}')
                continue

            correct_meta = self._choose_correct_strain(seqs=results, seq_len=len(seq.seq))
            seq.description = correct_meta

            yield seq

    def _try_get_seqs_from_strain(self, strain: str) -> Union[Iterable[SeqRecord], None]:
        max_retries = len(strain.split('/'))
        results = None

        for retry in range(max_retries - 1):
            results = self._get_seqs_from_strain(strain=strain, retry=retry)

            if results:
                break

        return results

    def _get_seqs_from_strain(self, strain: str, retry: int = 0) -> Union[Iterable[SeqRecord], None]:

        if retry != 0:
            strain = self._strain_prog_trim(strain=strain, retry=retry)

        query = {
            'datatype': 'protein',
            'output': UpdateMetaConsts.FLUDB_OUTPUT_FORMAT,
            'strainname': strain,
            'metadata': ','.join(UpdateMetaConsts.META_ORDER)
        }

        r = requests.get(
            url=f'{UpdateMetaConsts.FLUDB_API_ENDPOINT}',
            params=query
        )

        if not r.status_code == 200:
            return None

        result = r.text

        if not result or len(result) < 40:
            return None

        file_handle = StringIO(r.text)

        seqs_iterable = SeqIO.parse(
            handle=file_handle,
            format='fasta'
        )

        return seqs_iterable

    @classmethod
    def _choose_correct_strain(cls, seqs: Iterable[SeqRecord], seq_len: int) -> str:
        # match = ''.join([x.description for x in seqs if len(x) == seq_len and x])
        match = ''.join([None if not (len(x) == seq_len or x) else x.description for x in seqs][0])

        if not match:
            print('no')

        # print(f'Match is {match}')
        return match

    @classmethod
    def _strain_prog_trim(cls, strain: str, retry: int) -> str:
        stripped_strain_data = strain.rsplit('/', retry)

        return stripped_strain_data[0]
