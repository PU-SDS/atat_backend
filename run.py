from app.app import create_app

app = create_app()
app.run('0.0.0.0', 8080, debug=True, use_debugger=False, use_reloader=False, passthrough_errors=True)
