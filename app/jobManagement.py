import mysql.connector
import random
import string


# Generate analysisId
def getAnalysisId(stringLength=10):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

# Register Atat Analysis Task
def registerAtatAnalysis(app, analysisId,proteinName):
	mydb = mysql.connector.connect(host=app.config["DB_HOSTNAME"],user=app.config["DB_USERNAME"],passwd=app.config["DB_PASSWORD"],database=app.config["DB_DBNAME"])
	mycursor = mydb.cursor()
	sql = "INSERT INTO proteinAnalysis (analysisId,proteinName,date) VALUES (%s,%s,CURDATE())"
	val = (analysisId,proteinName)
	mycursor.execute(sql, val)
	try:
		mydb.commit()
		return True
	except:
		return False
