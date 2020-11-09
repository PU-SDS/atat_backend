from app.warehousing.mongodb.read import MongoDBRead


def position_dropdown_update(job_id: str):
    db_reader = MongoDBRead(
        job_id=job_id
    )

    max_position = db_reader.get_max_positions()
    dropdown_items = [{'label': x, 'value': x} for x in range(max_position)]

    return dropdown_items
