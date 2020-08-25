from app import ATATData

class ATAT(object):
    def __init__(self, host, reservoir):
        self.host = host
        self.reservoir = reservoir

    def run(self):
        lala = []
        for position in self.reservoir:
            position_num = position.position
            for idx, variant in enumerate(position.sequences):
                a = ATATData(position_num, 'index', 'master', variant.sequence, variant,
                             self.host[position_num].sequences[idx])
                lala.append(a)
        return lala

    def lala(self):
        a = list(zip(self.host, self.reservoir.sequences))
        for position in a:
