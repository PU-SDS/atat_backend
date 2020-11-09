from typing import Tuple

from dash.exceptions import PreventUpdate

from app.warehousing.mongodb.read import MongoDBRead


def position_update_variants_eventhandler(position: int, job_id: str) -> Tuple:
    if position is None:
        raise PreventUpdate

    db_reader = MongoDBRead(
        job_id=job_id.replace('/', '')
    )

    variants_host, variants_reservoir = db_reader.get_variants(position=position)

    return variants_host, variants_reservoir
