from django.urls import path
from . import views


from django.urls import path
from .views import predict_genome, home  

urlpatterns = [
    path('api/predict/', predict_genome, name='predict_genome'),
    path('generate-weblogo/', views.generate_weblogo, name='generate_weblogo'),
    path('', home, name='home'), 
]



# urlpatterns = [
#     path('api/predict/', views.predict_genome, name='predict_genome'),
# ]
