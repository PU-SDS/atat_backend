from dashboard.eventhandlers import variant_update_eventhandler
import plotly.express as px

def reservoir_variant_update_eventhandler(reservoir_selected_row: list, reservoir_row_data: list):
    eventhandler = variant_update_eventhandler(reservoir_selected_row, reservoir_row_data)()
    host_distribution = px.pie(eventhandler, names='host')

    return eventhandler, host_distribution
