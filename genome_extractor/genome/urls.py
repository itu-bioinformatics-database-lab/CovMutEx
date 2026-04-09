from django.urls import path

from . import views
from .views import known_hotspot_case_study, home, predict_genome

urlpatterns = [
    path('api/predict/', predict_genome, name='predict_genome'),
    path(
        'api/case-studies/known-hotspot/',
        known_hotspot_case_study,
        name='known_hotspot_case_study',
    ),
    path('generate-weblogo/', views.generate_weblogo, name='generate_weblogo'),
    path('', home, name='home'), 
]



# urlpatterns = [
#     path('api/predict/', views.predict_genome, name='predict_genome'),
# ]
