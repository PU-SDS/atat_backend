import dash_html_components as html
from dashboard.componants import INTRODUCTION_COMPONENT, SEARCH_JOB_COMPONENT


class Header(object):
    def set(self):
        header = html.Nav(className='navbar nav-bar-light bg-light', children=[
            html.Div(className='col-sm-4 p-0 pb-1 d-flex justify-content-center', children=[
                html.A(className='navbar-brand text-dark', href='#', children=[
                    html.Img(className='d-inline-block align-top', width=40, height=30, src="./assets/images/logo.png"),
                    'Perdana University'
                ])
            ]),
            html.Div(className='col-sm-4 p-0 pb-1 d-flex justify-content-center', children=[
                html.A('ATAT: Antigenic Transmissibility Analysis Tool',
                       className='navbar-brand d-inline-block align-top text-wrap text-justify text-center text-dark',
                       href='#')
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
        ])

        return header
