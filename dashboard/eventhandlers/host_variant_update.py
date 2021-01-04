from dashboard.eventhandlers import variant_update_eventhandler
import plotly.express as px


def host_variant_update_eventhandler(host_selected_row: list, host_row_data: list):
    eventhandler = variant_update_eventhandler(host_selected_row, host_row_data)()
    host_variant_country = px.pie(eventhandler, names='country')

    return eventhandler, host_variant_country
