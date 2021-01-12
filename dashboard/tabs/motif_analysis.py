import dash_html_components as html
import dash_core_components as dcc
from dashboard.componants import Collapse, Badge
from dashboard.constants import Constants
from dash_table import DataTable


class MotifAnalysis(object):
    def set(self):
        motif_analysis_tab = html.Div(className='col-sm-12',
                                      children=[
                                          Collapse(
                                              ele_id='motif_analysis',
                                              title='Motif Analysis',
                                              content=[
                                                  html.Div(
                                                      className='row mb-2',
                                                      children=[
                                                          html.Div(
                                                              className='col-sm-4',
                                                              children=[
                                                                  # Source motif dropdown
                                                                  dcc.Dropdown(
                                                                      id='source-motif-dropdown',
                                                                      clearable=False,
                                                                      options=Constants().MOTIF_CODES,
                                                                      placeholder='Select Source Motif'
                                                                  )
                                                              ]
                                                          ),
                                                          html.Div(
                                                              className='col-sm-2 d-flex align-content-center flex-wrap justify-content-center',
                                                              children=[
                                                                  # Switch icon

                                                                  html.Span(className='fa fa-arrow-right text-black')

                                                              ]
                                                          ),
                                                          html.Div(
                                                              className='col-sm-4',
                                                              children=[
                                                                  # Reservoir motif dropdown
                                                                  dcc.Dropdown(
                                                                      id='reservoir-motif-dropdown',
                                                                      clearable=False,
                                                                      options=Constants().MOTIF_CODES,
                                                                      placeholder='Select Reservoir Motif'
                                                                  )
                                                              ]
                                                          ),
                                                          html.Div(
                                                              className='col-sm-2',
                                                              children=[
                                                                  # Reservoir motif dropdown
                                                                  html.Button(
                                                                      id='search-motif',
                                                                      className='btn btn-block btn-danger',
                                                                      children=[
                                                                          html.Span(
                                                                              className='fas fa-search text-white'
                                                                          ),
                                                                          ' Search'
                                                                      ])
                                                              ]
                                                          )
                                                      ]
                                                  ),
                                                  DataTable(
                                                      id='motif-analysis-table',
                                                      data=[],
                                                      style_cell={'textAlign': 'left', 'padding': '5px'},
                                                      style_header={'backgroundColor': 'white',
                                                                    'fontWeight': 'bold'},
                                                      style_as_list_view=True,
                                                      row_selectable='single',
                                                      selected_rows=[0],
                                                      page_size=10,
                                                      sort_action="native",
                                                      sort_mode="single",
                                                      columns=Constants.MOTIF_ANALYSIS_COLUMNS
                                                  )
                                              ]).set()
                                      ])

        return motif_analysis_tab
