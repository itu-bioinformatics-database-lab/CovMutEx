from django.urls import path
from . import viewsUpdated as views
from .viewsUpdated import predict_genome, home, get_model_parameters

urlpatterns = [
    path('api/predict/', predict_genome, name='predict_genome'),
    path('api/model-parameters/', get_model_parameters, name='get_model_parameters'),
    path('api/models/', get_models, name='get_models'),
    path('generate-weblogo/', views.generate_weblogo, name='generate_weblogo'),
    path('', home, name='home'), 
]



# urlpatterns = [
#     path('api/predict/', views.predict_genome, name='predict_genome'),
# ]
