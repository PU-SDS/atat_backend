import dash_html_components as html
import dash_core_components as dcc
from dashboard.componants import Collapse


class Position(object):
    def set(self):
        position_tab = html.Div(className='col-sm-4', children=[
            Collapse(
                ele_id='select-position',
                content=dcc.Dropdown('position-dropdown', clearable=False),
                title='Position'
            ).set()
        ])

        return position_tab
