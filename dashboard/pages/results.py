import dash_html_components as html
import dash_core_components as dcc
from dash_extensions import Download

from dashboard.tabs import Statistics, Variants, VariantOrigins, VariantsSources, VariantsCountry, Position, MotifAnalysis, Controls


class ResultsPage(object):
    def __init__(self, results_display: str):
        self.results_display = results_display

    def set(self):
        results_page = html.Div(id='results-page-container', className=f'col-sm-12 d-none', children=[
            html.Div(className='row', children=[
                Position().set(),
                Controls().set()
            ]),
            html.Div(className='row', children=[
                # Container row for host and reservoir columns
                html.Div(className='col-sm-6', children=[
                    # Host container
                    Statistics(context='host').set(),
                    Variants(context='host').set(),
                    VariantOrigins(context='host').set(),
                    VariantsCountry(context='host').set()
                ]),
                html.Div(className='col-sm-6', children=[
                    # Reservoir container
                    Statistics(context='reservoir').set(),
                    Variants(context='reservoir').set(),
                    VariantOrigins(context='reservoir').set(),
                    VariantsSources(context='reservoir').set()
                ])
            ]),
            html.Div(className='row', children=[
                MotifAnalysis().set()
            ]),
            dcc.Location(id='results-page-location'),
            Download(id='download-results-driver')
        ])

        return results_page
