import dash_html_components as html
from dashboard.componants import Collapse, Badge


class Statistics(object):
    def __init__(self, context: str):
        self.context = context

    def set(self):
        statistics_tab = Collapse(
            ele_id=f'{self.context}-statistics',
            content=[
                html.Div(
                    className='col-sm-4 border-right d-flex justify-content-center',
                    children=Badge(
                        idx=f'{self.context}-position-entropy',
                        text_color='text-secondary',
                        border_color='border-secondary',
                        title='Entropy',
                        value=0.842254
                    ).set()
                ),
                html.Div(
                    className='col-sm-4 border-right d-flex justify-content-center',
                    children=Badge(
                        idx=f'{self.context}-position-variants',
                        text_color='text-secondary',
                        border_color='border-secondary',
                        title='Variants',
                        value=20
                    ).set()
                ),
                html.Div(
                    className='col-sm-4 d-flex justify-content-center',
                    children=Badge(
                        idx=f'{self.context}-position-support',
                        text_color='text-secondary',
                        border_color='border-secondary',
                        title='Supports',
                        value=586
                    ).set()
                )
            ],
            title=self.context.capitalize()
        ).set()

        return statistics_tab
