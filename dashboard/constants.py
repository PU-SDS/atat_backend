class Constants(object):
    VARIANTS_COLUMN_LIST = [{"id": "sequence", "name": "Variant"},
                            {"id": "conservation", "name": "Conservation"}, {"id": "count", "name": "Supports"},
                            {"id": "motif_long", "name": "Motif"}]

    VARIANT_ORIGINS_COLUMN_LIST = [{"id": "id", "name": "NCBI ID"}, {"id": "strain", "name": "Strain"},
                                   {"id": "host", "name": "Host"}]

    MOTIF_CODES = [
        {'label': 'Index', 'value': 'I'},
        {'label': 'Major', 'value': 'Ma'},
        {'label': 'Minor', 'value': 'Mi'},
        {'label': 'Unique', 'value': 'U'}
    ]

    MOTIF_ANALYSIS_COLUMNS = [{"id": "position", "name": "Position"},
                            {"id": "sequence", "name": "Sequence"}]
