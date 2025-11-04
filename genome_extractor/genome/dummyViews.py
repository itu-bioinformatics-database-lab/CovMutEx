import time
import traceback
from django.http import JsonResponse
from rest_framework.decorators import api_view
from .runtime import get_model 
from .tasks import run_model_prediction
from celery.exceptions import TimeoutError

@api_view(["GET", "POST"])
def predict_genome(request):
    if request.method in ['POST', 'GET']:
        return handle_prediction_sync(request, timeout=5.0) # NFR-1: 5s latency
    return JsonResponse({"error": "Invalid request method"}, status=400)

def handle_prediction_sync(request, timeout):
    try:
        start_time = time.time()
        data = request.data if request.method == 'POST' else request.query_params
        
        inputs = {
            "nodeId": data.get('nodeId'),
            "elapsedDay": int(data.get('elapsedDay', 0)),
            "selectedModel": data.get('selectedModel'),
            "selectedProteinRegion": data.get('selectedProteinRegion')
        }
        
        model_name = inputs.get("selectedModel", "legacy_v1") 
        adapter = get_model(model_name)
        if not adapter:
            return JsonResponse({"error": f"Model '{model_name}' bulunamadı."}, status=404)
        
        preprocessed_data = adapter.preprocess(inputs)
        
        task_result = run_model_prediction.delay(
            model_name=model_name, 
            preprocessed_data=preprocessed_data
        )
        
        try:
            final_output = task_result.get(timeout=timeout) # FR-1.4
        except TimeoutError:
            task_result.revoke(terminate=True) 
            return JsonResponse({
                "error": f"Tahmin {timeout} saniyeden uzun sürdü (NFR-1)."
            }, status=504)
        except Exception as e:
            return JsonResponse({"error": f"Model çalıştırılırken hata oluştu: {str(e)}"}, status=500)

        response_data = {
            "nodeId": inputs["nodeId"],
            "elapsedDay": inputs["elapsedDay"],
            "selectedModel": model_name,
            "selectedProteinRegion": inputs["selectedProteinRegion"],
            "model_metadata": adapter.metadata(),
            "results": final_output 
        }
        
        end_time = time.time()
        print(f"Total runtime (sync call): {end_time - start_time:.2f} seconds")

        return JsonResponse(response_data)

    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": f"Genel bir hata oluştu: {str(e)}"}, status=500)