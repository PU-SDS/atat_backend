from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .endpoints import job, results

app = FastAPI(
    title="Antigenic Transmissibility Analysis Tool RESTful API",
    description="This is the RESTful API of the Antigenic Transmissibility Analysis Tool.",
    contact={'name': 'Shan Tharanga', 'email': 'stwm2@student.london.ac.uk'},
    version="1.0.1",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(job.router)
app.include_router(results.router)
