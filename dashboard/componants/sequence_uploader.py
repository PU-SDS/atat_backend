import dash_html_components as html
import dash_core_components as dcc


class SequenceUploader(object):
    def __init__(self, context: str):
        self.context = context

    def set(self):
        sequence_upload_block = html.Div(className='form-group row align-items-center', children=[
                                    html.Label(className='col-4 mb-0', children=f'{self.context.capitalize()} '
                                                                                f'Sequences:'),
                                    html.Div(className='col-4', children=[
                                        html.Div(className='btn-group', children=[
                                            html.Span('sequences.fasta', id=f'{self.context}-filename',
                                                      className='input-group-text rounded-0'),
                                            dcc.Upload(id=f'{self.context}-seq', children=[
                                                html.Button('Upload', className='btn btn-secondary rounded-0')
                                            ])
                                        ])
                                    ])
                                ])

        return sequence_upload_block
