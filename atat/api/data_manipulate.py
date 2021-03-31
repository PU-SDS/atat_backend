from collections import defaultdict

from mongoengine.base import BaseList


class DataManipulate(object):
    @staticmethod
    def baselist_to_dict(base_list: BaseList, context: str) -> dict:
        """
            When a MongoDB query returns a BaseList we need to convert it into a dict before we can serialize it.
            This method does exactly that.

            :param base_list: The base list returned by Mongo.
            :param context: The context which is used as the key of the dictionary

            :type base_list: BaseList
            :type context: str

            :returns: A dictionary that contains the data within the BaseList in a list under the given context as key.
        """

        variant_dict = defaultdict(list)

        for variant in base_list:
            variant_dict[context].append(variant)

        return variant_dict
