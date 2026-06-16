from django.urls import path
from . import viewsUpdated as views
from .viewsUpdated import predict_genome, home, get_model_parameters, get_models
# , delete_model
from .benchmark_views import (
    run_benchmark_view, run_dataset_benchmark_view,
    get_benchmark_results, list_benchmarks,
    benchmark_datasets_view, export_benchmark_view,
    search_cov_spectrum_variants_view,
)

urlpatterns = [
    # Existing endpoints
    path('api/predict/', predict_genome, name='predict_genome'),
    path('api/models/', get_models, name='get_models'),
    path('api/model-parameters/', get_model_parameters, name='get_model_parameters'),
    # path('api/model_parameters/', get_model_parameters, name='get_model_parameters_alt'),
    # path('api/models/delete/', delete_model, name='delete_model'),
    path('generate-weblogo/', views.generate_weblogo, name='generate_weblogo'),
    
    # Benchmark endpoints
    path('api/benchmark/run/', run_benchmark_view, name='benchmark_run'),
    path('api/benchmark/run-dataset/', run_dataset_benchmark_view, name='benchmark_run_dataset'),
    path('api/benchmark/results/', get_benchmark_results, name='benchmark_results'),
    path('api/benchmark/list/', list_benchmarks, name='benchmark_list'),
    path('api/benchmark/datasets/', benchmark_datasets_view, name='benchmark_datasets'),
    path('api/benchmark/export/', export_benchmark_view, name='benchmark_export'),
    path('api/benchmark/cov-spectrum/search-variants/', search_cov_spectrum_variants_view, name='benchmark_cov_spectrum_search'),
    
    path('', home, name='home'),
]
