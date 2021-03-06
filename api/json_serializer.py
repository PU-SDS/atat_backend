from flask import make_response


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
