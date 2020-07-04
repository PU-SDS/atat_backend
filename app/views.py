from flask import Blueprint, url_for

blueprint = Blueprint('root', __name__)

# - Endpoint: root.index
# - URL: /index2

@blueprint.route('/index')
def index():
    return url_for('root.index')

@blueprint.route('/home')
def home():
    return "<a href='{}'> Go back to index</a>".format(url_for('root.index'))