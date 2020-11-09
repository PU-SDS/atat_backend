import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

from dashboard.componants import SequenceUploader, InputRow


class ResultsPreCheck(object):
    def __init__(self, precheck_display: str):
        self.precheck_display = precheck_display

    def set(self):
        submit_job_page = html.Div(id='precheck-page-container', className=f'col-sm-12 d-none', children=[
            html.Div(className='row w-100 d-flex justify-content-center align-items-center m-0', children=[
                html.Div(className='col-4 d-flex justify-content-center p-0', children=[
                    html.Div(id='status-box', className='alert alert-info d-flex justify-content-center', **{'role': 'alert'}, children=[
                        html.H3('PROCESSING', id='precheck-status', className='pr-3'),
                        dbc.Spinner(color='info')
                    ])
                ]),
                dcc.Interval(id='precheck-page-interval', interval=5000, disabled=False)
            ])
        ])

        return submit_job_page
