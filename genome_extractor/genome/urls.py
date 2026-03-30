from django.urls import path
from .viewsUpdated import predict_genome, home, generate_weblogo, get_models, get_model_parameters

urlpatterns = [
    path('api/predict/', predict_genome, name='predict_genome'),
    path('api/models/', get_models, name='get_models'),
    path('api/model-parameters/', get_model_parameters, name='get_model_parameters'),
    path('generate-weblogo/', generate_weblogo, name='generate_weblogo'),
    path('', home, name='home'),
]
