from dash_table import DataTable
from dashboard.componants import Collapse
from dashboard.constants import Constants


class Variants(object):
    def __init__(self, context: str):
        self.context = context

    def set(self):
        variants_tab = Collapse(
            title='Variants',
            content=DataTable(
                id=f'{self.context}-variants',
                data=[],
                style_cell={'textAlign': 'left', 'padding': '5px'},
                style_header={'backgroundColor': 'white', 'fontWeight': 'bold'},
                style_as_list_view=True,
                row_selectable='single',
                selected_rows=[0],
                page_size=10,
                sort_action="native",
                sort_mode="single",
                columns=Constants.VARIANTS_COLUMN_LIST
            ),
            ele_id=f'{self.context}-variants-collapse'
        ).set()

        return variants_tab
