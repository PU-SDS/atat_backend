from typing import Tuple

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from app.warehousing.mongodb.constants import JobStatuses
from app.warehousing.mongodb.read import MongoDBRead

from dashboard.eventhandlers import position_dropdown_update
from dashboard.events.constants import PrecheckStatusStyles


class PrecheckPageEvents(object):
    def __init__(self, dash_app):
        self.dash_app = dash_app

    def set(self):
        @self.dash_app.callback(
            [
                Output('precheck-page-interval', 'disabled'),
                Output('results-page-container', 'className'),
                Output('precheck-page-container', 'className'),
                Output('status-box', 'className'),
                Output('precheck-status', 'children'),
                Output('position-dropdown', 'options')
            ],
            [Input('precheck-page-interval', 'n_intervals'), Input('url', 'pathname')]
        )
        def job_status_check(interval: int, job_id: str) -> Tuple:
            if job_id is None or job_id =='/':
                raise PreventUpdate

            job_id = job_id.split('/')[1]
            job_status = MongoDBRead(job_id).get_job_status()

            if job_status == JobStatuses.JOB_NONEXIST:
                status = JobStatuses.JOB_NONEXIST

                return True, \
                       'col-sm-12 d-none', \
                       'col-sm-12 d-block', \
                       PrecheckStatusStyles().set(status), \
                       status, \
                       []
            elif job_status == JobStatuses.JOB_FAILED:
                status = JobStatuses.JOB_FAILED

                return True, \
                       'col-sm-12 d-none', \
                       'col-sm-12 d-block', \
                       PrecheckStatusStyles().set(status), \
                       status, \
                       []
            elif job_status == JobStatuses.JOB_PROCESSING:
                status = JobStatuses.JOB_PROCESSING

                return False, \
                       'col-sm-12 d-none', \
                       'col-sm-12 d-block', \
                       PrecheckStatusStyles().set(status), \
                       status, \
                       []

            status = JobStatuses.JOB_COMPLETED

            return True, \
                   'col-sm-12 d-block', \
                   'col-sm-12 d-none', \
                   PrecheckStatusStyles().set(status), \
                   status, \
                   position_dropdown_update(job_id=job_id)
