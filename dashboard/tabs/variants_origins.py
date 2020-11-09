from dash_table import DataTable
from dashboard.componants import Collapse
from dashboard.constants import Constants


class VariantOrigins(object):
    def __init__(self, context: str):
        self.context = context

    def set(self):
        variants_origins_tab = Collapse(
                    title='Variant Origins',
                    content=DataTable(
                        id=f'{self.context}-variant-details',
                        data=[],
                        style_cell={'textAlign': 'left', 'padding': '5px'},
                        style_header={'backgroundColor': 'white', 'fontWeight': 'bold'},
                        style_as_list_view=True,
                        page_size=10,
                        filter_action="native",
                        columns=Constants.VARIANT_ORIGINS_COLUMN_LIST
                    ),
                    ele_id=f'{self.context}-variant-details-collapse'
                ).set()

        return variants_origins_tab
