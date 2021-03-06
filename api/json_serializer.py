from flask import make_response
import json
from bson import ObjectId


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
            resp = make_response(data, code)
            resp.headers.extend(headers or {})
            return resp

        return _


class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return json.JSONEncoder.default(self, o)
