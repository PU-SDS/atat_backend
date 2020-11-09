from dash.exceptions import PreventUpdate


def variant_update_eventhandler(host_selected_row: list, host_row_data: list):
    def _eventhandler() -> list:
        if not host_row_data:
            raise PreventUpdate

        selected_row = host_selected_row[0]
        row_data = host_row_data[selected_row]
        arranged_variant_data = zip(row_data.get('id'), row_data.get('strain'), row_data.get('host'))

        variant_data = [
            {'id': id_, 'strain': strain, 'host': host}
            for id_, strain, host in arranged_variant_data]

        return variant_data

    return _eventhandler
