import dash_html_components as html


class InputRow(object):
    def __init__(self, label: str, content):
        self.label = label
        self.content = content

    def set(self):
        input_row = html.Div(className='form-group row align-items-center', children=[
                                    html.Label(className='col-4 mb-0', children=self.label),
                                    html.Div(className='col-4', children=self.content)
                                ])

        return input_row
