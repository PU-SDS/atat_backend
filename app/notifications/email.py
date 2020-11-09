from jinja2 import Template
import urllib.parse
import requests


class EmailNotification(object):
    SENDER_ADDRESS = 'notifications@perdanauniversity.com'
    SENDER_NAME = 'Perdana University Bioinformatics Team'

    URL_PREFIX = 'https://atat.perdanauniversity/results?'

    def __init__(self, email: str, job_id: str):
        self.email = email
        self.job_id = job_id

    def send(self):
        raw_template = ''
        job_url = urllib.parse.urljoin(self.URL_PREFIX, self.job_id)

        with open('C:\\Reps\\atat-single\\app\\notifications\\template.html', 'r') as f:
            raw_template = f.read()

        custom_template = Template(raw_template).render(job_url=job_url)

        return requests.post(
            "https://api.mailgun.net/v3/mg.shanetw.tech/messages",
            auth=("api", "adadcdd75fbd715ec0bb1c2f71353141-ea44b6dc-5869f1d2"),
            data={"from": "Perdana University Bioinformatics Team <no-reply@mg.shanetw.tech>",
                  "to": self.email,
                  "subject": f"Your job {self.job_id} has completed!",
                  "html": custom_template})
