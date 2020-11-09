import dash_core_components as dcc
from dashboard.componants import Collapse


class VariantsSources(object):
    def __init__(self, context: str):
        self.context = context

    def set(self):
        variants_sources_tab = Collapse(
            title='Variant Sources',
            content=dcc.Graph(
                id=f'{self.context}-sources-distribution'
            ),
            ele_id=f'{self.context}-variant-sources-collapse'
        ).set()

        return variants_sources_tab
