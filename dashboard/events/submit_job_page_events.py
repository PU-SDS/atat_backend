import base64

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app.core import SaveUploads
from app.tasks.run_atat import run_atat


class SubmitJobPageEvents(object):
    def __init__(self, dash_app):
        self.dash_app = dash_app

    def set(self):
        @self.dash_app.callback(
            Output('host-filename', 'children'),
            Input('host-seq', 'filename')
        )
        def host_seq_updated(filename: str):
            if filename is None:
                raise PreventUpdate

            return filename

        @self.dash_app.callback(
            Output('reservoir-filename', 'children'),
            Input('reservoir-seq', 'filename')
        )
        def reservoir_seq_updated(filename: str):
            if filename is None:
                raise PreventUpdate

            return filename

        @self.dash_app.callback(
            Output('submitjob-page-location', 'pathname'),
            Input('submit-job', 'n_clicks'),
            [State('host-seq', 'contents'), State('reservoir-seq', 'contents'), State('email', 'value'),
             State('kmer-length', 'value')]
        )
        def reservoir_seq_updated(n_clicks: int, host_seq, reservoir_seq, email, kmer_length):
            if host_seq is None or reservoir_seq is None or n_clicks is None:
                raise PreventUpdate

            host_file_data = host_seq.split(',')[1]
            reservoir_file_data = reservoir_seq.split(',')[1]

            host_seq_bytes = base64.b64decode(host_file_data)
            reservoir_seq_bytes = base64.b64decode(reservoir_file_data)

            job_id = SaveUploads(host_file=host_seq_bytes, reservoir_file=reservoir_seq_bytes).save2()
            run_atat.delay(job_id, email, int(kmer_length))

            return f'/{job_id}'
