"""
ASGI config for genome_extractor project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/5.1/howto/deployment/asgi/
"""

import os

from .runtime_env import configure_runtime_tempdir
from django.core.asgi import get_asgi_application

configure_runtime_tempdir()
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'genome_extractor.settings')

application = get_asgi_application()
