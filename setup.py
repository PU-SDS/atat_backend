"""
This script is used to set up the environment, and resources on first run.
"""

import time

from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError, OperationFailure, CollectionInvalid

from atat_backend.settings import ResourceSettings


class Setup(object):
    def __init__(self):
        self.settings = ResourceSettings()
        self._con_retries = 0

        self.client = MongoClient(
            host=f'{self.settings.mongo_host}',
            port=27017,
            username=self.settings.mongo_root_username,
            password=self.settings.mongo_root_password,
        )

    def run(self):
        self._check_if_connected()
        self._create_db()
        self._create_user()

    def _check_if_connected(self):
        try:
            print(self.client.server_info())
        except ServerSelectionTimeoutError as e:
            if self._con_retries > 10:
                print(e)
                exit(1)

            time.sleep(5)
            self._con_retries += 1

            print(f'Retrying Database Connection x{self._con_retries}/10')

            self._check_if_connected()

    def _create_db(self):
        dbs = self.client.list_database_names()
        atat_db = self.settings.mongo_atat_database
        rabbitmq_db = self.settings.mongo_rabbitmq_database

        if atat_db in dbs or rabbitmq_db in dbs:
            print('Databases already present')
            return

        try:
            db = self.client[atat_db]
            result = db.create_collection('delete_me')
            print(result)
        except CollectionInvalid as e:
            print(e)
            exit(2)

        try:
            db = self.client[rabbitmq_db]
            result = db.create_collection('delete_me')
            print(result)
        except CollectionInvalid as e:
            print(e)
            exit(2)

    def _create_user(self):
        atat_db = self.settings.mongo_atat_database
        rabbitmq_db = self.settings.mongo_rabbitmq_database

        try:
            db = self.client[atat_db]
            result = db.command(
                'createUser',
                self.settings.mongo_atat_username,
                pwd=self.settings.mongo_atat_password,
                roles=[{'role': 'readWrite', 'db': atat_db}],
            )

            print(result)
        except OperationFailure as e:
            print(e.details)

        try:
            db = self.client[rabbitmq_db]
            result = db.command(
                'createUser',
                self.settings.mongo_rabbitmq_username,
                pwd=self.settings.mongo_rabbitmq_password,
                roles=[{'role': 'readWrite', 'db': rabbitmq_db}],
            )

            print(result)
        except OperationFailure as e:
            print(e.details)


def main():
    Setup().run()


if __name__ == "__main__":
    main()
