import itertools
from operator import itemgetter


class MotifClasses(object):
    I = 'Index'
    Ma = 'Major'
    Mi = 'Minor'
    U = 'Unique'


class PList(list):
    def get_all_variants(self):
        all_variants = (x.sequences for x in self)
        all_variants = itertools.chain(*all_variants)
        return all_variants

    def get_total_support(self):
        total = [y.count for y in itertools.chain(*(x.sequences for x in self))]
        total = sum(total)
        return total


class Position(object):
    def __init__(self, position, sequences, variants_flattened, entropy=None):
        self.position = position
        self.entropy = entropy
        self.variants_flattened = variants_flattened
        self.supports = len(self.variants_flattened)
        self.sequences = sequences
        self.variants = len(self.sequences)

    def __setattr__(self, key, value):
        if key == 'sequences':
            value = sorted(value, key=itemgetter('count'), reverse=True)

            for idx, variant in enumerate(value):
                if idx == 0:
                    value[idx].motif_short = 'I'
                    value[idx].motif_long = MotifClasses.I
                elif idx == 1:
                    value[idx].motif_short = 'Ma'
                    value[idx].motif_long = MotifClasses.Ma
                elif value[idx].count > 1:
                    value[idx].motif_short = 'Mi'
                    value[idx].motif_long = MotifClasses.Mi
                else:
                    value[idx].motif_short = 'U'
                    value[idx].motif_long = MotifClasses.U

                value[idx].support_percentage = round(
                    (value[idx].count * 100) / self.supports
                )

        super(Position, self).__setattr__(key, value)


class Variant(object):
    def __init__(self, position, sequence, count, conservation, motif_short=None,
                 motif_long=None, support_percentage=None):
        self.position = position
        self.sequence = sequence
        self.count = count
        self.conservation = conservation
        self.motif_short = motif_short
        self.motif_long = motif_long
        self.support_percentage = support_percentage

    def __getitem__(self, item):
        return getattr(self, item)
