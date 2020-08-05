class Position(object):

    def __init__(self, position, variants, variants_flattened, entropy=None):
        self.position = position
        self.entropy = entropy
        self.variants = variants
        self.variants_flattened = variants_flattened

    @property
    def get_variant_count(self):
        return len(self.variants)

    @property
    def get_support_count(self):
        return len(self.variants_flattened)
    
    def __repr__(self):
        return f"<Position: {self.position} {{{len(self.variants)} | {self.entropy:.2f}}}>"


class Variant(object):

    def __init__(self, sequence, count, conservation):
        self.sequence = sequence
        self.count = count
        self.conservation = conservation
