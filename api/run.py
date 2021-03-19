from collections import defaultdict

from flask import Flask
from flask_restful import Resource, Api, abort
from mongoengine import DoesNotExist

from atat_single.api.data_manipulate import DataManipulate
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
            return JobQueries.get_result(jobid).reservoir.get(position=position).to_json(indent=2), 200
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not have position {position}.')


class GetResult(Resource):
    def get(self, jobid: str):
        return JobQueries.get_result(jobid).to_json(indent=2), 200


class GetPositionCount(Resource):
    def get(self, jobid: str):
        return str(JobQueries.get_result(jobid).source.count()), 200


class GetAllMotifSwitches(Resource):
    def get(self, jobid: str):
        return JobQueries.get_result(jobid).switches, 200


class GetPositionMotifSwitches(Resource):
    def get(self, jobid: str, position: int):
        try:
            return JobQueries.get_result(jobid).switches.get(position=position), 200
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not have any motif switches at position {position}.')


class GetSourcePositionVariants(Resource):
    def get(self, jobid: str, position: int):
        try:
            variants_list = JobQueries.get_result(jobid).source.get(position=position).variants
            variants_dict = DataManipulate.baselist_to_dict(variants_list, 'variants')

            return variants_dict, 200
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not have position {position}.')


class GetReservoirPositionVariants(Resource):
    def get(self, jobid: str, position: int):
        try:
            variants_list = JobQueries.get_result(jobid).reservoir.get(position=position).variants
            variants_dict = DataManipulate.baselist_to_dict(variants_list, 'variants')

            return variants_dict, 200
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not have position {position}.')


api.add_resource(GetJob, '/info/<string:jobid>')
api.add_resource(GetResult, '/results/<string:jobid>')
api.add_resource(GetGroupedPosition, '/results/<string:jobid>/positions/<int:position>/grouped')
api.add_resource(GetSourcePosition, '/results/<string:jobid>/positions/<int:position>/source')
api.add_resource(GetReservoirPosition, '/results/<string:jobid>/positions/<int:position>/reservoir')
api.add_resource(GetPositionCount, '/results/<string:jobid>/positions/count')
api.add_resource(GetAllMotifSwitches, '/results/<string:jobid>/switches')
api.add_resource(GetPositionMotifSwitches, '/results/<string:jobid>/positions/<int:position>/switches')
api.add_resource(GetSourcePositionVariants, '/results/<string:jobid>/positions/<int:position>/source/variants')
api.add_resource(GetReservoirPositionVariants, '/results/<string:jobid>/positions/<int:position>/reservoir/variants')

if __name__ == '__main__':
    app.run(debug=True)
