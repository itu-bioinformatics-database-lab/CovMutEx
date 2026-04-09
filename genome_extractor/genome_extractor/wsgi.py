"""
WSGI config for genome_extractor project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/5.1/howto/deployment/wsgi/
"""

import os

from .runtime_env import configure_runtime_tempdir
from django.core.wsgi import get_wsgi_application

configure_runtime_tempdir()
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'genome_extractor.settings')

application = get_wsgi_application()
