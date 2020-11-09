import dash_html_components as html


class Badge(object):
    def __init__(self, idx: str, text_color: str, border_color: str, title: str, value: float):
        self.idx = idx
        self.text_color = text_color
        self.border_color = border_color
        self.title = title
        self.value = value

    def set(self):
        card = html.Div(className=f'card {self.text_color} {self.border_color} w-100 d-flex text-align-center',
                        children=[
                            html.Div(className='card-header pb-1 pt-1 d-flex justify-content-center',
                                     children=self.title),
                            html.Div(className='card-body pb-1 pt-1 d-flex justify-content-center '
                                               'align-items-center',
                                     children=html.H5(self.value, id=self.idx, className='card-title m-1'))
                        ])

        return card
