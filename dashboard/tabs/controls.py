import dash_html_components as html
import dash_core_components as dcc
from dashboard.componants import Collapse


class Controls(object):
    def set(self):
        position_tab = html.Div(className='col-sm-8', children=[
            Collapse(
                ele_id='page-controls',
                title='Controls',
                content=[
                    html.Div(
                        className='btn-group mr-2',
                        children=[
                            html.Button(
                                id='print-results',
                                className='btn btn-info dropdown-toggle',
                                children=[
                                    html.Span(className='fas fa-download text-white'),
                                    ' Download'
                                ],
                                **{'data-toggle': 'dropdown'}

                            ),
                            html.Div(
                                className='dropdown-menu',
                                children=[
                                    html.A(
                                        id='download-json',
                                        className='dropdown-item',
                                        href='#',
                                        children=[
                                            html.Span(className='fa fa-file-alt text-black'),
                                            ' JSON'
                                        ]
                                    ),
                                    html.A(
                                        id='download-csv',
                                        className='dropdown-item',
                                        href='#',
                                        children=[
                                            html.Span(className='fa fa-file-excel text-black'),
                                            ' CSV'
                                        ]
                                    )
                                ]
                            ),
                            html.Button(
                                id='download-results',
                                className='btn btn-outline-dark',
                                children=[
                                    html.Span(className='fas fa-print text-black'),
                                    ' Print'
                                ]

                            )
                        ]
                    ),
                    html.Div(
                        className='btn-group mr-2',
                        children=[
                            html.Button(
                                id='motif-control',
                                className='btn btn-outline-dark',
                                children=[
                                    html.Span(className='fas fa-dna text-black'),
                                    ' Motif Analysis'
                                ]
                            ),
                            html.Button(
                                id='statistics-control',
                                className='btn btn-outline-dark',
                                children=[
                                    html.Span(className='fas fa-cogs text-black'),
                                    ' Statistics'
                                ]

                            )
                        ]
                    )
                ]
            ).set()
        ])

        return position_tab
