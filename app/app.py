from flask import Flask
import os
from flask_bootstrap import Bootstrap


def create_app(debug=True):
    app = Flask(__name__)

    if debug is True:
        app.config.from_object('config.DevelopmentConfig')
    else:
        app.config.from_object('config.ProductionConfig')

    bootstrap = Bootstrap(app)

    # Load paths of commandline tools
    app.config.from_object('config.BioApps')
    app.config.from_object('config.LocalPaths')
    app.config.from_object('config.MongoDBs')
    app.config.from_object('config.CeleryConfig')

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    with app.app_context():
        from app import views
        app.register_blueprint(views.blueprint)

    return app
