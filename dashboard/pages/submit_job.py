import dash_html_components as html
import dash_core_components as dcc

from dashboard.componants import SequenceUploader, InputRow


class SubmitJobPage(object):
    def __init__(self, submit_job_display: str):
        self.submit_job_display = submit_job_display

    def set(self):
        submit_job_page = html.Div(id='submit-job-container', className=f'col-sm-12 d-none', children=[
            html.Div(className='d-flex justify-content-center align-items-center', children=[
                html.Div(className='row border rounded w-100', children=[
                    html.Div(className='col-12 mt-2 mb-2', children=[
                        SequenceUploader(context='host').set(),
                        SequenceUploader(context='reservoir').set(),
                        InputRow(label='K-mer Length', content=[
                            dcc.Input(id='kmer-length',
                                          className='form-control border-top-0 border-left-0 border-right-0 '
                                                    'border-bottom-1 shadow-none rounded-0',
                                          placeholder='Default: 9')
                        ]).set(),
                        InputRow(label='E-mail', content=[
                            dcc.Input(id='email',
                                          className='form-control border-top-0 border-left-0 border-right-0 '
                                                    'border-bottom-1 shadow-none rounded-0',
                                          placeholder='Ex: example@example.com')
                        ]).set(),
                        html.Div(className='form-group row align-items-center mb-1', children=[
                            html.Div(className='col-12 d-flex justify-content-center', children=[
                                dcc.Loading(id='lalalala', fullscreen=True, type='graph', children=[html.Button('Submit', id='submit-job', className='btn btn-info')],
                                            style={'content': 'shalalalala'})
                            ])
                        ])
                    ])
                ])
            ]),
            dcc.Location(id='submitjob-page-location', refresh=True)
        ])

        return submit_job_page
