import dash_html_components as html


class Collapse(object):
    def __init__(self, ele_id, title, content):
        self.ele_id = ele_id
        self.title = title
        self.content = content

    def set(self):
        collapse = html.Div(className='card', children=[
            html.Div(self.title, className='card-header bg-secondary text-white pt-1 pb-1', **{'data-toggle': 'collapse',
                                                                                'data-target': f'#{self.ele_id}'}),

            html.Div(className='panel-collapse collapse show', id=self.ele_id, children=[
                html.Div(className='card-body overflow-auto mh-50', children=[
                    self.content
                ])
            ])
        ])

        return collapse
