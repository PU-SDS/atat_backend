from typing import Tuple

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash

from dashboard.eventhandlers import position_update_statistics_eventhandler, position_update_variants_eventhandler, \
    host_variant_update_eventhandler, reservoir_variant_update_eventhandler, motif_search_eventhandler, \
    download_results_json_eventhandler, download_results_csv_eventhandler


class ResultsPageEvents(object):
    def __init__(self, dash_app):
        self.dash_app = dash_app

    def set(self):
        @self.dash_app.callback(
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
            [Input('position-dropdown', 'value'), Input('results-page-location', 'pathname')]
        )
        def position_update_event(position: int, job_id: str) -> Tuple:
            if position is None or job_id is None or job_id == '/':
                raise PreventUpdate

            host, reservoir = position_update_statistics_eventhandler(position=position, job_id=job_id)

            return *position_update_variants_eventhandler(position=position, job_id=job_id), \
                   host.get('entropy'), host.get('variants'), host.get('supports'), \
                   reservoir.get('entropy'), reservoir.get('variants'), reservoir.get('supports')

        @self.dash_app.callback(
            [
                Output('host-variant-details', 'data'),
                Output('host-variants-countries', 'figure')
            ],
            [
                Input('host-variants', 'selected_rows'),
                Input('host-variants', 'derived_virtual_data')
            ]
        )
        def host_variant_update_event(host_selected_row: list, host_row_data: list):

            return host_variant_update_eventhandler(host_selected_row, host_row_data)

        @self.dash_app.callback(
            [
                Output('reservoir-variant-details', 'data'),
                Output('reservoir-sources-distribution', 'figure')
            ],
            [
                Input('reservoir-variants', 'selected_rows'),
                Input('reservoir-variants', 'derived_virtual_data')
            ]
        )
        def reservoir_variant_update_event(reservoir_selected_row: list, reservoir_row_data: list):
            return reservoir_variant_update_eventhandler(reservoir_selected_row, reservoir_row_data)

        @self.dash_app.callback(
            Output('motif-analysis-table', 'data'),
            [
                Input('search-motif', 'n_clicks'),
                Input('results-page-location', 'pathname')
            ],
            [
                State('source-motif-dropdown', 'value'),
                State('reservoir-motif-dropdown', 'value')
            ]
        )
        def motif_search(clicks: int, job_id: str, source_motif: str, reservoir_motif: str):
            if job_id is None or job_id == '/' or clicks is None:
                raise PreventUpdate

            return motif_search_eventhandler(job_id, source_motif, reservoir_motif)

        @self.dash_app.callback(
            Output('download-results-driver', 'data'),
            [
                Input('results-page-location', 'pathname'),
                Input('download-json', 'n_clicks'),
                Input('download-csv', 'n_clicks')
            ], prevent_initial_call=True
        )
        def download_results(job_id: str, click_json: int, click_csv: int) -> dict:
            if click_json is None and click_csv is None:
                raise PreventUpdate

            clicked: dict = dash.callback_context.triggered[0]

            if 'download-json' in clicked.get('prop_id'):
                return download_results_json_eventhandler(job_id)
            else:
                return download_results_csv_eventhandler(job_id)
