import json
from collections import defaultdict

from flask import make_response
from bson import ObjectId
from mongoengine.base import EmbeddedDocumentList


class JSONSerializer(object):
    def __init__(self, api):
        """
            We need this because Flask-RESTful does not provide JSON serialization by default. But Flask does.

            :param api: The instance of Flask-RESTful API.
        """

        self.api = api

    def serializer(self):
        @self.api.representation('application/json')
        def _(data, code, headers=None):
            if isinstance(data, EmbeddedDocumentList):
                data = self._embedded_document_to_dict(data)

            resp = make_response(data, code)
            resp.headers['Access-Control-Allow-Origin'] = '*'
            resp.headers.extend(headers or {})
            return resp

        return _

    @classmethod
    def _embedded_document_to_dict(cls, document: EmbeddedDocumentList):
        """
            This just converts the EmbeddedDocuments inside a EmbeddedDocumentList into a dictionary.
        """

        items = defaultdict(list)

        for count, item in enumerate(document):
            items[document._name].append(item.to_mongo().to_dict())

        return items


class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)

        return json.JSONEncoder.default(self, o)
