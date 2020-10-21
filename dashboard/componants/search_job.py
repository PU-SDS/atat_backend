import dash_html_components as html

SEARCH_JOB_COMPONENT =  html.Div(className='pos-f-t col-sm-12', children=[
    html.Div(className='collapse navbar-collapse', id='check-results', children=[
        html.Div(className='col-sm-12 d-flex justify-content-center bg-light p-4', children=[
            html.Form(className='form-inline', children=[
                html.Div(className='input-group', children=[
                    html.Div(className='input-group-prepend', children=[
                        html.Span('Job ID:', className='input-group-text', id='job-id-desc')
                    ]),
                    html.Textarea(cols=25, style={'resize': 'none'}, rows=1, className='form-control',
                                  placeholder='Ex: fyysdc', **{'aria-label': 'Job ID'}),
                    html.Div(className='input-group-append', children=[
                        html.Button(className='btn btn-link btn-danger', children=[
                            html.Span(className='fas fa-search text-white')
                        ])
                    ])
                ])
            ])
        ])
    ])
])