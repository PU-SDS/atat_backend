from typing import Tuple

import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

from dashboard.componants import Collapse, Header
from dashboard.tabs import Statistics, Variants, VariantOrigins, VariantsSources, Position

from flask import Flask

import plotly.express as px

server = Flask(__name__)

dash_app = dash.Dash(server=server, external_stylesheets=[
    dbc.themes.BOOTSTRAP,
    'https://use.fontawesome.com/releases/v5.7.0/css/all.css'
],
                     external_scripts=[
                         'https://code.jquery.com/jquery-3.5.1.js',
                         'https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js',
                         'https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js'
                     ])

dash_app.server.config.from_object('config.MongoDBs')

kk = dash_app.server.config

with dash_app.server.app_context():
    from dashboard.eventhandlers import result_page_load_eventhandler, position_update_variants_eventhandler, \
        host_variant_update_eventhandler, position_update_statistics_eventhandler

dash_app.layout = html.Div(children=[
    Header().set(),
    html.Div(className='container pt-1', style={'max-width': '100%'}, children=[
        html.Div(className='row', children=[
            Position().set()
        ]),
        html.Div(className='row', children=[
            # Container row for host and reservoir columns
            html.Div(className='col-sm-6', children=[
                # Host container
                Statistics(context='host').set(),
                Variants(context='host').set(),
                VariantOrigins(context='host').set()
            ]),
            html.Div(className='col-sm-6', children=[
                # Reservoir container
                Statistics(context='reservoir').set(),
                Variants(context='reservoir').set(),
                VariantOrigins(context='reservoir').set(),
                VariantsSources(context='reservoir').set()
            ])
        ])
    ]),
    dcc.Location(
        id='url',
        refresh=False
    )
])


@dash_app.callback(
    Output('position-dropdown', 'options'),
    [Input('url', 'pathname')]
)
def result_page_load_event(job_id: str) -> Tuple:
    return result_page_load_eventhandler(job_id=job_id)


@dash_app.callback(
    [
        Output('host-variants', 'data'),
        Output('reservoir-variants', 'data'),
        Output('host-position-entropy', 'children'),
        Output('host-position-variants', 'children'),
        Output('host-position-support', 'children'),
        Output('reservoir-position-entropy', 'children'),
        Output('reservoir-position-variants', 'children'),
        Output('reservoir-position-support', 'children'),
    ],
    [Input('position-dropdown', 'value'), Input('url', 'pathname')]
)
def position_update_event(position: int, job_id: str) -> Tuple:
    if position is None:
        raise PreventUpdate

    host, reservoir = position_update_statistics_eventhandler(position=position, job_id=job_id)
    return *position_update_variants_eventhandler(position=position, job_id=job_id), \
           host.get('entropy'), host.get('variants'), host.get('supports'), \
           reservoir.get('entropy'), reservoir.get('variants'), reservoir.get('supports')


@dash_app.callback(
    Output('host-variant-details', 'data'),
    [
        Input('host-variants', 'selected_rows'),
        Input('host-variants', 'derived_virtual_data')
    ]
)
def host_variant_update_event(host_selected_row: list, host_row_data: list):
    return host_variant_update_eventhandler(host_selected_row=host_selected_row, host_row_data=host_row_data)


if __name__ == "__main__":
    dash_app.run_server(port=8057)
