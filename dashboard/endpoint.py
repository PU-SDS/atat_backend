import itertools

import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc

from dash.dependencies import Input, Output
from dash_table import DataTable

from dashboard.componants import INTRODUCTION_COMPONENT, SEARCH_JOB_COMPONENT, Collapse
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
    from app.warehousing.mongodb.read import MongoDBRead


dash_app.layout = html.Div(children=[
    html.Nav(className='navbar nav-bar-light bg-light', children=[
        html.Div(className='col-sm-4 p-0 pb-1 d-flex justify-content-center', children=[
            html.A(className='navbar-brand text-dark', href='#', children=[
                html.Img(className='d-inline-block align-top', width=40, height=30, src="./assets/images/logo.png"),
                'Perdana University'
            ])
        ]),
        html.Div(className='col-sm-4 p-0 pb-1 d-flex justify-content-center', children=[
            html.A('ATAT: Antigenic Transmissibility Analysis Tool',
                   className='navbar-brand d-inline-block align-top text-wrap text-justify text-center text-dark', href='#')
        ]),
        html.Div(className='col-sm-4 p-0 pb-1 d-flex justify-content-center', children=[
            html.Div(className='btn-group', role='group', children=[
                html.Button('Check Results', className='btn btn-info navbar-toggle collapsed',
                            **{'data-toggle': 'collapse', 'data-target': '#check-results'}),
                html.Button('Introduction', className='btn btn-outline-secondary navbar-toggle collapsed',
                            **{'data-toggle': 'collapse', 'data-target': '#introduction'})
            ])
        ]),
        SEARCH_JOB_COMPONENT,
        INTRODUCTION_COMPONENT
    ]),
    html.Div(className='container pt-1', style={'max-width': '100%'}, children=[
        html.Div(className='row', children=[
            html.Div(className='col-sm-12', children=[
                Collapse(
                    title='Entropy',
                    content=dcc.Graph(
                        id='entropy-graph'
                    ),
                    ele_id='entropy-collapse'
                ).set()
            ])
        ]),
        html.Div(className='row', children=[
            # Container row for host and reservoir columns
            html.Div(className='col-sm-6', children=[
                # Host container
                Collapse(
                    title='Variants',
                    content=DataTable(
                        id='host-variants',
                        data=[],
                        style_cell={'textAlign': 'left', 'padding': '5px'},
                        style_header={'backgroundColor':  'white', 'fontWeight': 'bold'},
                        style_as_list_view=True,
                        row_selectable='single',
                        selected_rows=[0],
                        page_size=10,
                        sort_action="native",
                        sort_mode="single"
                    ),
                    ele_id='host-variants-collapse'
                ).set(),
                Collapse(
                    title='Variant Origins',
                    content=DataTable(
                        id='host-variant-details',
                        data=[],
                        style_cell={'textAlign': 'left', 'padding': '5px'},
                        style_header={'backgroundColor':  'white', 'fontWeight': 'bold'},
                        style_as_list_view=True,
                        page_size=10,
                        filter_action="native",
                    ),
                    ele_id='host-variant-details-collapse'
                ).set(),
                Collapse(
                    title='Variant Sources',
                    content=dcc.Graph(
                        id='host-sources-distribution'
                    ),
                    ele_id='host-variant-sources-collapse'
                ).set()
            ]),
            html.Div(className='col-sm-6', children=[
                # Reservoir container
                Collapse(
                    title='Variants',
                    content=DataTable(
                        id='reservoir-variants',
                        data=[],
                        style_cell={'textAlign': 'left', 'padding': '5px'},
                        style_header={'backgroundColor':  'white', 'fontWeight': 'bold'},
                        style_as_list_view=True,
                        row_selectable='single',
                        selected_rows=[0],
                        page_size=10,
                        sort_action="native",
                        sort_mode="single"
                    ),
                    ele_id='reservoir-variants-collapse'
                ).set(),
                Collapse(
                    title='Variant Origins',
                    content=DataTable(
                        id='reservoir-variant-details',
                        data=[],
                        style_cell={'textAlign': 'left', 'padding': '5px'},
                        style_header={'backgroundColor':  'white', 'fontWeight': 'bold'},
                        style_as_list_view=True,
                        page_size=10,
                        filter_action="native",
                    ),
                    ele_id='reservoir-variant-details-collapse'
                ).set(),
                Collapse(
                    title='Variant Sources',
                    content=dcc.Graph(
                        id='reservoir-sources-distribution'
                    ),
                    ele_id='reservoir-variant-sources-collapse'
                ).set()
            ])
        ])
    ]),
    dcc.Location(
        id='url',
        refresh=False
    ),
    html.Textarea(id='lala')
])

@dash_app.callback(
    Output('entropy-graph', 'figure'),
    [Input('url', 'pathname')]
)
def _update_entropy_graph(job_id: str):
    if not job_id:
        return []

    db_reader = MongoDBRead(
        job_id=job_id.replace('/', '')
    )
    aaaa = db_reader.get_variants(2)

    job_data = db_reader.get_all_positions(
        columns={
            'position': 1,
            'entropy': 1,
            'sequences': 1,
            '_id': 0
        }
    )


    figg = job_data.plot(x='position', y=['entropy-host', 'entropy-reservoir'])

    figg.update_layout(clickmode='event+select')
    return figg

@dash_app.callback(
    [Output('host-variants', 'data'), Output('host-variants', 'columns'), Output('reservoir-variants', 'data'),
     Output('reservoir-variants', 'columns')],
    [Input('entropy-graph', 'clickData'), Input('url', 'pathname')]
)
def _lala(click_data, job_id):
    if not click_data:
        return [], [], [], []

    db_reader = MongoDBRead(
        job_id=job_id.replace('/', '')
    )

    position = click_data['points'][0]['x']
    variants_host, variants_reservoir = db_reader.get_variants(position=position)

    return variants_host, \
           [{"id": "position", "name": "Position"}, {"id": "sequence", "name": "Variant"}, {"id": "conservation", "name": "Conservation"}, {"id": "count", "name": "Supports"}, {"id": "motif_long", "name": "Motif"}],\
           variants_reservoir, \
           [{"id": "position", "name": "Position"}, {"id": "sequence", "name": "Variant"}, {"id": "conservation", "name": "Conservation"}, {"id": "count", "name": "Supports"}, {"id": "motif_long", "name": "Motif"}]

@dash_app.callback(
    [
        Output('host-variant-details', 'data'),
        Output('host-variant-details', 'columns'),
        Output('host-sources-distribution', 'figure'),
        Output('reservoir-variant-details', 'data'),
        Output('reservoir-variant-details', 'columns'),
        Output('reservoir-sources-distribution', 'figure')
    ],
    [
        Input('host-variants', 'selected_rows'),
        Input('host-variants', 'derived_virtual_data'),
        Input('reservoir-variants', 'selected_rows'),
        Input('reservoir-variants', 'derived_virtual_data')
    ]
)
def _lala2(a, b, c, d):
    if not a:
        return [], []

    kk = b[a[0]]

    sequences = []
    for id, strain, host in zip(kk.get('id'), kk.get('strain'), kk.get('host')):
        sequences.append({'id': id, 'strain': strain, 'host': host})

    figg = px.pie(sequences, names='host')

    kk2 = b[a[0]]

    sequences2 = []
    for id, strain, host in zip(kk2.get('id'), kk2.get('strain'), kk2.get('host')):
        sequences2.append({'id': id, 'strain': strain, 'host': host})

    figg2 = px.pie(sequences, names='host')

    return sequences, \
           [{"id": "id", "name": "NCBI ID"}, {"id": "strain", "name": "Strain"}, {"id": "host", "name": "Host"}], \
           figg, \
           sequences2, \
           [{"id": "id", "name": "NCBI ID"}, {"id": "strain", "name": "Strain"}, {"id": "host", "name": "Host"}], \
           figg2

if __name__ == "__main__":
    dash_app.run_server(port=8057)