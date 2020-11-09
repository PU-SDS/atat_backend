from typing import Tuple

import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc

from dash.dependencies import Input, Output

from dashboard.componants import Header
from dashboard.pages import SubmitJobPage, ResultsPage, ResultsPreCheck


from flask import Flask


server = Flask(__name__)

dash_app = dash.Dash(server=server, external_stylesheets=[
    dbc.themes.BOOTSTRAP,
    'https://use.fontawesome.com/releases/v5.7.0/css/all.css',
    'https://drive.google.com/uc?export=view&id=15jYfXXzrcjjg-5c-31XICQo58uNBjvAf'
],
                     external_scripts=[
                         'https://code.jquery.com/jquery-3.5.1.js',
                         'https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js',
                         'https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js'
                     ])

dash_app.server.config.from_object('config.MongoDBs')
dash_app.server.config.from_object('config.LocalPaths')

kk = dash_app.server.config

with dash_app.server.app_context():
    from dashboard.events import ResultsPageEvents, SubmitJobPageEvents, PrecheckPageEvents

dash_app.layout = html.Div(children=[
    Header().set(),
    html.Div(id='container', className='container pt-1 container-center', style={'max-width': '90%'}, children=[
        SubmitJobPage('d-none').set(), ResultsPage('d-none').set(), ResultsPreCheck('d-none').set()
    ]),
    dcc.Location(
        id='url',
        refresh=True
    )
])


@dash_app.callback(
    Output('submit-job-container', 'className'),
    [Input('url', 'pathname')]
)
def result_page_load_event(job_id: str):
    if job_id == None or job_id == '/':
        return 'col-sm-12 d-block'

    return 'col-sm-12 d-none'


ResultsPageEvents(dash_app).set()
SubmitJobPageEvents(dash_app).set()
PrecheckPageEvents(dash_app).set()

if __name__ == "__main__":
    dash_app.run_server(port=8057)
