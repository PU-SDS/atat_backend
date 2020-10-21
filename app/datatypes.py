import itertools
from collections import Counter, defaultdict
from operator import itemgetter

"""
    Datatypes Module
    
    Describes multiple datatypes used in the Hunana and ATAT packages
"""

class MotifClasses(object):
    """ Defines the short and long name of Motifs. """

    I = 'Index'
    Ma = 'Major'
    Mi = 'Minor'
    U = 'Unique'


class ATATData(object):
    def __init__(self, reservoir, host):
        """
            ATAT Result Data Structure

            Describes the data structure for storing ATAT results. It stores the host and reservoir variants at each
            position.

            :param position: The zero based position number within the sequences
            :param sequence:
            :param reservoir: The
            :param host:
        """

        self.reservoir = reservoir
        self.host = host


class VariantDict(defaultdict):
    """
        Data structure for storing sequence indexes

        Describes the data structure for storing the indexes of sequences where a specific variant originates from

        Methods:
            get_counter: Returns a Counter object for all the variants of the current position
    """

    def get_counter(self) -> Counter:
        """
            Get a counter object for the current position's variants

            :return: Returns a counter object for variants at current position
        """

        position = Counter()
        for x in self:
            position[x] = len(self[x])
        return position


class PList(list):
    def get_all_variants(self):
        all_variants = (x.sequences for x in self)
        all_variants = itertools.chain(*all_variants)
        return all_variants

    def get_total_support(self):
        total = [y.count for y in itertools.chain(*(x.sequences for x in self))]
        total = sum(total)
        return total


class Position(dict):
    SEQ_DESCRIPTIONS = None

    def __init__(self, position, sequences, variants_flattened, variant_dict, entropy=None):
        """
            Data structure for kmer positions

            Describes the data structure for kmer positions within the alignment

            :param position: A zero-based index for position number
            :param sequences: A list of type Variant containing variants seen at current position
            :param variants_flattened: A flattened list of all variants. Note: Only used for processing
            :param variant_dict: A default dictionary of lists to store the sequence idx for later getting description data
            :param entropy: The Shannon entropy at the current kmer position

            :type position: int
            :type sequences: list
            :type variants_flattened: list
            :type variant_dict: defaultdict
            :type entropy: float
        """

        self.position = position
        self.entropy = entropy
        self.variants_flattened = variants_flattened
        self.supports = len(self.variants_flattened)
        self.sequences = self._motif_classify(sequences)
        self.variants = len(self.sequences)

        self._set_desc_data(variant_dict)

    def _set_desc_data(self, variant_dict):
        """
            Sets the description data for each variant within the Position object

            :param variant_dict: A defaultdict of list containing idx of sequences
            :type variant_dict: defaultdict
        """

        for dict_variant, variant in itertools.product(variant_dict, self.sequences):
            if dict_variant != variant.sequence:
                continue
            idx = variant_dict[dict_variant]
            variant.id = self._get_id(idx)
            variant.strain = self._get_strain(idx)
            variant.country = self._get_country(idx)
            variant.host = self._get_host(idx)

    def __setattr__(self, key, value):
        # This really makes me sad, find a way around this mess TODO: fix this please, also in Variant
        dict.__setitem__(self, key, value)

    def __getattr__(self, item):
        return dict.__getitem__(self, item)

    @classmethod
    def _motif_classify(cls, variants) -> list:
        """
        Given a list of Variant objects sort by the count of each, and  classify each variant into motif classes

        :param variants: A list of Variant objects
        :type variants: list

        :return: A list containing Variant objects sorted by count and each variant classified into motif classes
        """

        variants = sorted(variants, key=lambda x: x.count, reverse=True)

        for idx, variant in enumerate(variants):
            if idx == 0:
                variants[idx].motif_short = 'I'
                variants[idx].motif_long = MotifClasses.I
            elif idx == 1:
                variants[idx].motif_short = 'Ma'
                variants[idx].motif_long = MotifClasses.Ma
            elif variants[idx].count > 1:
                variants[idx].motif_short = 'Mi'
                variants[idx].motif_long = MotifClasses.Mi
            else:
                variants[idx].motif_short = 'U'
                variants[idx].motif_long = MotifClasses.U

        return variants

    def _get_id(self, indexs) -> list:
        """
            Given a list of sequence indexes return a list of sequence ids

            :param indexs: A list of sequence indexs
            :type indexs: list

            :return: A list containing sequence ids
        """

        return [self.SEQ_DESCRIPTIONS[x].split('|')[0] for x in indexs]

    def _get_strain(self, indexs) -> list:
        """
            Given a list of sequence indexes return a list of sequence strains

            :param indexs: A list of sequence indexs
            :type indexs: list

            :return: A list containing sequence ids
        """

        return [self.SEQ_DESCRIPTIONS[x].split('|')[1] for x in indexs]

    def _get_country(self, indexs) -> list:
        """
            Given a list of sequence indexes return a list of sequence countries

            :param indexs: A list of sequence indexs
            :type indexs: list

            :return: A list containing sequence ids
        """

        return [self.SEQ_DESCRIPTIONS[x].split('|')[2] for x in indexs]

    def _get_host(self, indexs) -> list:
        """
            Given a list of sequence indexes return a list of sequence hosts

            :param indexs: A list of sequence indexs
            :type indexs: list

            :return: A list containing sequence ids
        """

        return [self.SEQ_DESCRIPTIONS[x].split('|')[3] for x in indexs]


class Variant(dict):
    def __init__(self, position, sequence, count, conservation, motif_short=None, motif_long=None, idx=None,
                 strain=None, country=None, host=None):
        """
            Data Structure for kmer variants

            Describes the data structure for kmer variants

            :param position: The position kmer position within the sequence
            :param sequence: The sequence of the given variant
            :param count: The number of times the variant is seen in the alignment at this position
            :param conservation: The percentage conservation of this variant for this position
            :param motif_short: The short name motif classification for this variant
            :param motif_long: The long name motif classification for this variant
            :param idx: A list containing the idx of the sequences where this variant was observed
            :param strain: A list containing the strain names of the sequences where this variant was observed
            :param country: A list containing the country names of the sequences where this variant was observed
            :param host: A list containing the host names of the sequences where this variant was observed

            :type position: int
            :type sequence: str
            :type count: int
            :type conservation: float
            :type motif_short: str
            :type motif_long: str
            :type idx: list
            :type strain: list
            :type country: list
            :type host: list
        """

        self.position = position
        self.sequence = sequence
        self.count = count
        self.conservation = conservation
        self.motif_short = motif_short
        self.motif_long = motif_long
        self.id = idx
        self.strain = strain
        self.country = country
        self.host = host

    def __setattr__(self, key, value):
        # This really makes me sad, find a way around this mess TODO: fix this please, also in Variant
        dict.__setitem__(self, key, value)

    def __getattr__(self, item):
        return dict.__getitem__(self, item)
