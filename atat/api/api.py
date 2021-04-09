from flask import Flask
from flask_restful import Resource, Api, abort, reqparse
from mongoengine import DoesNotExist

from .json_serializer import JSONSerializer, JSONEncoder
from ..models import Result
from .data_manipulate import DataManipulate
from .job_queries import JobQueries
from ..tasks.run_job import run_job

app = Flask(__name__)
app.json_encoder = JSONEncoder
api = Api(app)

argparser = reqparse.RequestParser()
argparser.add_argument('email', type=str, required=True, help='Please provide your email address so we can let you '
                                                              'know once we are done.')
argparser.add_argument('kmer_length', type=int, required=True, help='Please provide a k-mer length.')
argparser.add_argument('source', type=str, required=True, help='Please provide a source dataset in FASTA format '
                                                               '(co-aligned with the reservoir dataset).')
argparser.add_argument('reservoir', type=str, required=True, help='Please provide a reservoir dataset in FASTA format '
                                                                  '(co-aligned with the source dataset).')

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


class SubmitJob(Resource):
    def post(self):
        args = argparser.parse_args()
        jobid = DataManipulate.get_random_jobid()

        run_job.delay(
            source_seqs=args.get('source'),
            reservoir_seqs=args.get('reservoir'),
            jobid=jobid,
            kmer_len=args.get('kmer_length'),
            header_decode=True,
            header_format='(accession)|(strain)|(host)|(country)'
        )

        return jobid, 200


class GetJobStatus(Resource):
    def get(self, jobid: str):
        return JobQueries.get_job(jobid).status.to_mongo().to_dict(), 200


class GetJobLog(Resource):
    def get(self, jobid: str):
        return JobQueries.get_job(jobid).log.to_mongo().to_dict(), 200


api.add_resource(GetJob, '/info/<string:jobid>')
api.add_resource(GetResult, '/results/<string:jobid>')
api.add_resource(GetJobStatus, '/status/<string:jobid>')
api.add_resource(GetJobLog, '/log/<string:jobid>')
api.add_resource(GetGroupedPosition, '/results/<string:jobid>/positions/<int:position>/grouped')
api.add_resource(GetSourcePosition, '/results/<string:jobid>/positions/<int:position>/source')
api.add_resource(GetReservoirPosition, '/results/<string:jobid>/positions/<int:position>/reservoir')
api.add_resource(GetPositionCount, '/results/<string:jobid>/positions/count')
api.add_resource(GetAllMotifSwitches, '/results/<string:jobid>/switches')
api.add_resource(GetPositionMotifSwitches, '/results/<string:jobid>/positions/<int:position>/switches')
api.add_resource(GetSourcePositionVariants, '/results/<string:jobid>/positions/<int:position>/source/variants')
api.add_resource(GetReservoirPositionVariants, '/results/<string:jobid>/positions/<int:position>/reservoir/variants')
api.add_resource(SubmitJob, '/submit')
