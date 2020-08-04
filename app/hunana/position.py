class Position(object):
    POSITION = None
    ENTROPY = None
    VARIANTS = None
    VARIANTS_FLATTENED = None

    @property
    def get_variant_count(self):
        return len(self.VARIANTS)

    @property
    def get_support_count(self):
        return len(self.VARIANTS_FLATTENED)

class Variant(object):
    SEQUENCE = None
    COUNT = None
    CONSERVATION = None

    def __init__(self, sequence, count, conservation):
        self.SEQUENCE = sequence
        self.COUNT = count
        self.CONSERVATION = conservation


