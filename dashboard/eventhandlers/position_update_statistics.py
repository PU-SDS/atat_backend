from typing import Tuple

from dash.exceptions import PreventUpdate

from app.warehousing.mongodb.read import MongoDBRead


def position_update_statistics_eventhandler(position: int, job_id: str) -> Tuple:
    if job_id is None:
        raise PreventUpdate

    db_reader = MongoDBRead(
        job_id=job_id.replace('/', '')
    )

    statistics_host, statistics_reservoir = db_reader.get_position_statistics(position=position)

    return statistics_host[0], statistics_reservoir[0]
