from flask import Flask
import os
from app import views
from flask_bootstrap import Bootstrap


def create_app(debug=True):
    app = Flask(__name__)

    if debug is True:
        app.config.from_object('config.DevelopmentConfig')
    else:
        app.config.from_object('config.ProductionConfig')

    # Load paths of commandline tools
    app.config.from_object('config.BioApps')
    app.config.from_object('config.LocalPaths')
    app.config.from_object('config.FileNames')
    app.config.from_object('config.Containers')
    app.config.from_object('config.MongoDBs')

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    app.register_blueprint(views.blueprint)
    Bootstrap(app)

    return app
