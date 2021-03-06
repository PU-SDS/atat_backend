from flask import Flask
from flask_restful import Resource, Api, abort
from mongoengine import DoesNotExist

from atat_single.models import Job, Result
from atat_single.api.json_serializer import JSONSerializer, JSONEncoder

from atat_single.api.job_queries import JobQueries

app = Flask(__name__)
app.json_encoder = JSONEncoder
api = Api(app)

JSONSerializer(api).serializer()


class GetJob(Resource):
    def get(self, jobid: str):
        return JobQueries.get_job(jobid).to_json(indent=2), 200


class GetGroupedPosition(Resource):
    def get(self, jobid: str, position: int):
        result_id = JobQueries.get_result(jobid).id

        try:
            return Result.objects.get_grouped_position(result_id, position), 200
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not have position {position}.')


class GetSourcePosition(Resource):
    def get(self, jobid: str, position: int):
        try:
            return JobQueries.get_result(jobid).source.get(position=position).to_json(indent=2), 200
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not have position {position}.')


class GetReservoirPosition(Resource):
    def get(self, jobid: str, position: int):
        try:
            return JobQueries.get_result(jobid).source.get(position=position).to_json(indent=2), 200
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not have position {position}.')


class GetResult(Resource):
    def get(self, jobid: str):
        return JobQueries.get_result(jobid).to_json(indent=2), 200


api.add_resource(GetJob, '/info/<string:jobid>')
api.add_resource(GetResult, '/results/<string:jobid>')
api.add_resource(GetGroupedPosition, '/results/<string:jobid>/position/grouped/<int:position>')
api.add_resource(GetSourcePosition, '/results/<string:jobid>/position/source/<int:position>')
api.add_resource(GetReservoirPosition, '/results/<string:jobid>/position/reservoir/<int:position>')

if __name__ == '__main__':
    app.run(debug=True)
