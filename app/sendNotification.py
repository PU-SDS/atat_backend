import os
import smtplib
from email.message import EmailMessage
from email.utils import make_msgid
import mimetypes


def send_email_notification(app,user_email,analysisId):
	email_username=app.config['EMAIL_USERNAME']
	email_password=app.config['EMAIL_PASSWORD']

	atatURL=app.config['ATAT_URL']
	atatResultURL=app.config['ATAT_RESULT_URL']+"?analysisId="+analysisId

	mail_template_file = os.path.join(app.config['ATAT_TEMPLATE_DIR'],'email.html')

	msg = EmailMessage()
	msg['Subject'] = app.config['EMAIL_SUBJECT']
	msg['From'] = app.config['EMAIL_FROM']
	msg['To'] = user_email

	msg.set_content('Dear ATAT User,\n\nThe antigenic transmissibility analysis for your sequences with analysis id: '+analysisId+' is completed. To view the results, you can revisit ATAT website ('+atatURL+') and select "View Result" from the menu, then enter your analysis id.\n\nYou could also view the results directly by visiting this URL: '+atatResultURL+'.\n\nThank You,\nViral Diversity Team@PU-ScDS')

	temp_file = open(mail_template_file, "r")
	if temp_file.mode == 'r':
		contents = temp_file.read()
		msg.add_alternative(contents.format(analysisId=analysisId, atatURL=atatURL, atatResultURL=atatResultURL), subtype='html')

	with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
		smtp.login(email_username, email_password)
		smtp.send_message(msg)
