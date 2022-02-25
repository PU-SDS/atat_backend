from typing import Literal

from pydantic import BaseSettings


class ResourceSettings(BaseSettings):
    """
    All constants and settings will be here.
    For more info: https://pydantic-docs.helpmanual.io/usage/settings/
    """

    env_state: Literal["prod", "dev"] = "dev"  # Change this when going production

    mongo_host: str
    mongo_port: int = 27017
    mongo_replicaset_name: str = "rs0"
    mongo_root_username: str = None  # Only needed in DEV to set up the database so we default to None in Prod
    mongo_root_password: str = None  # Only needed in DEV to set up the database so we default to None in Prod

    mongo_atat_username: str
    mongo_atat_password: str

    mongo_rabbitmq_username: str
    mongo_rabbitmq_password: str

    mongo_atat_database: str
    mongo_rabbitmq_database: str

    rabbitmq_host: str
    rabbitmq_port: int
    rabbitmq_username: str
    rabbitmq_password: str

    class Config:
        env_file: str = ".env"
        validate_assignment: bool = True

        # Environment variables will always take priority over values loaded from a dotenv file.
        fields = {
            "env_state": {"env": ["ENVIRONMENT_STATE"]},
            "mongo_host": {"env": ["MONGO_ATAT_HOST"]},
            "mongo_port": {"env": ["MONGO_ATAT_PORT"]},
            "mongo_replicaset_name": {"env": ["MONGO_REPLICASET_NAME"]},
            "mongo_root_username": {"env": ["MONGO_INITDB_ROOT_USERNAME"]},
            "mongo_root_password": {"env": ["MONGO_INITDB_ROOT_PASSWORD"]},
            "mongo_atat_username": {"env": ["MONGODB_ATAT_USERNAME"]},
            "mongo_atat_password": {"env": ["MONGODB_ATAT_PASSWORD"]},
            "mongo_rabbitmq_username": {"env": ["MONGODB_RABBITMQ_USERNAME"]},  # For results backend
            "mongo_rabbitmq_password": {"env": ["MONGODB_RABBITMQ_PASSWORD"]},  # For results backend
            "mongo_atat_database": {"env": ["MONGO_ATAT_DATABASE"]},
            "mongo_rabbitmq_database": {"env": ["MONGO_RABBITMQ_DATABASE"]},
            "rabbitmq_host": {"env": ["RABBITMQ_SERVICE_HOST"]},
            "rabbitmq_port": {"env": ["RABBITMQ_SERVICE_PORT"]},
            "rabbitmq_username": {"env": ["RABBITMQ_BROKER_USERNAME"]},
            "rabbitmq_password": {"env": ["RABBITMQ_BROKER_PASSWORD"]},
        }
