from app import ATATData, Variant


class ATAT(object):
    def __init__(self, host, reservoir):
        """
            ATAT Analysis Module

            TODO: Write better description.

            :param host: Iterable of host kmer positions
            :param reservoir: Iterable of reservoir kmer positions.

            :type host: iterable
            :type reservoir: iterable
        """

        self.host = host
        self.reservoir = reservoir

    def run(self) -> list:
        """
            Creates a list of ATATData objects for each variant at each position within the alignment.

            :return: A list of ATATData objects
        """

        positions = zip(self.host, self.reservoir)
        data = []

        for p_host, p_res in positions:
            for variant in p_host.sequences:
                match_idx = [str(index) for index, _ in enumerate(p_res.sequences) if _.sequence == variant.sequence]
                if len(match_idx) != 0:
                    match_idx = ''.join(match_idx)
                    atat_data = ATATData(
                        sequence=variant.sequence,
                        position=variant.position,
                        reservoir=p_res.sequences[int(match_idx)],
                        host=variant)
                else:
                    atat_data = ATATData(
                        position=variant.position,
                        sequence=variant.sequence,
                        reservoir=Variant(
                            position=variant.position,
                            sequence='',
                            count=0,
                            conservation=0,
                            motif_short='',
                            motif_long='',
                            idx=[],
                            strain=[],
                            country=[],
                            host=[]
                        )
                        ,
                        host=variant)

                data.append(atat_data)

            for variant in p_res.sequences:
                match_idx = [str(index) for index, _ in enumerate(p_host.sequences) if
                             _.sequence == variant.sequence]
                if len(match_idx) != 0:
                    match_idx = ''.join(match_idx)
                    atat_data = ATATData(
                        sequence=variant.sequence,
                        position=variant.position,
                        reservoir=variant,
                        host=p_host.sequences[int(match_idx)])
                else:
                    atat_data = ATATData(
                        position=variant.position,
                        sequence=variant.sequence,
                        reservoir=variant,
                        host=Variant(
                            position=variant.position,
                            sequence='',
                            count=0,
                            conservation=0,
                            motif_short='',
                            motif_long='',
                            idx=[],
                            strain=[],
                            country=[],
                            host=[]
                        )
                    )

                data.append(atat_data)

        return data
