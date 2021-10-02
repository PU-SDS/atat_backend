import uvicorn

from viva_atat.v1.api import app


def main():
    uvicorn.run(app, port=8090)


if __name__ == "__main__":
    main()
