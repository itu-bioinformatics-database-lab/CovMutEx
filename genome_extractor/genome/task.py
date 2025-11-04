from celery import shared_task
from .runtime import get_model
import time

@shared_task(bind=True, time_limit=60, soft_time_limit=50) # FR-1.4: Timeout
def run_model_prediction(self, model_name: str, preprocessed_data: tuple):
    try:
        adapter = get_model(model_name)
        if not adapter:
            raise ValueError(f"Model '{model_name}' bulunamadı.")
        
        prediction_data = adapter.predict(preprocessed_data)
        final_output = adapter.postprocess(prediction_data)
        
        return final_output
    except Exception as e:
        self.retry(exc=e, countdown=5, max_retries=1)
        return {"error": str(e)}