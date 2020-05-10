# Statement for enabling the development environment
DEBUG = True

# Define the application directory
import os
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Application threads. A common general assumption is
# using 2 per available processor cores - to handle
# incoming requests using one and performing background
# operations using the other.
THREADS_PER_PAGE = 2

# Enable protection agains *Cross-site Request Forgery (CSRF)*
CSRF_ENABLED     = True

# Use a secure, unique and absolutely secret key for
# signing the data.
CSRF_SESSION_KEY = "LxH1AYbL3aUTTLbmyuxn1A0WiMhJLVhZ"

# Secret key for signing cookies
SECRET_KEY = "LxH1AYbL3aUTTLbmyuxn1A0WiMhJLVhZ"

UPLOADS_DEFAULT_DEST = '/home/centos/atat-single/app/input/'
UPLOADED_FILES_ALLOW = ('fasta')
