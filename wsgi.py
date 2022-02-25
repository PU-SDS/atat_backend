import uvicorn

from atat_backend.v1.api import app


def main():
    uvicorn.run(app, port=8090)


if __name__ == "__main__":
    main()
