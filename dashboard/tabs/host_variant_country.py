import dash_core_components as dcc
from dashboard.componants import Collapse


class VariantsCountry(object):
    def __init__(self, context: str):
        self.context = context

    def set(self):
        variants_country_tab = Collapse(
            title='Variant Countries',
            content=dcc.Graph(
                id=f'{self.context}-variants-countries'
            ),
            ele_id=f'{self.context}-variant-countries-collapse'
        ).set()

        return variants_country_tab
