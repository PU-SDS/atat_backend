from app.warehousing.mongodb.read import MongoDBRead


def motif_search_eventhandler(job_id: str, source_motif: str, reservoir_motif: str) -> list:
    db_reader = MongoDBRead(
        job_id=job_id.replace('/', '')
    )

    source_motifs, reservoir_motifs = db_reader.get_motifs(source_motif, reservoir_motif)

    # Convert the  list to a dictionary so we can access the position numbers efficiently
    source_motifs_dict = _mongo_list_to_dict(source_motifs)
    reservoir_motifs_dict = _mongo_list_to_dict(reservoir_motifs)

    motif_switches = []

    for position, source_variants in source_motifs_dict.items():
        reservoir_variants = reservoir_motifs_dict.get(position)

        if reservoir_variants is None:
            continue

        for source_variant in source_variants:
            # Get the variant sequence from the source variant dictionary
            source_variant_seq = source_variant.get('sequences').get('sequence')

            for reservoir_variant in reservoir_variants:
                reservoir_variant_seq = reservoir_variant.get('sequences').get('sequence')

                if source_variant_seq != reservoir_variant_seq:
                    continue

                motif_switches.append(
                    {
                        'position': position,
                        'sequence': reservoir_variant_seq
                    }
                )
    return motif_switches


def _mongo_list_to_dict(motif_list: list) -> dict:
    result_dict = {position.get('_id'): position.get('variants') for position in motif_list}

    return result_dict
