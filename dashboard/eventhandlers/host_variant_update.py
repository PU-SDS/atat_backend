from dashboard.eventhandlers import variant_update_eventhandler


def host_variant_update_eventhandler(host_selected_row: list, host_row_data: list):
    eventhandler = variant_update_eventhandler(host_selected_row, host_row_data)

    return eventhandler()
