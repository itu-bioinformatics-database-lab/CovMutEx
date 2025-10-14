# /genome_extractor/genome/views.py dosyasının içinde

from .legacy_adapter import LegacyModelAdapter # Yeni adaptör sınıfımızı import ediyoruz

def handle_prediction(request):
    try:
        start_time = time.time()
        print("[ADAPTER-DRIVEN] Handling prediction request")

        adapter = LegacyModelAdapter()
        data = request.data if request.method == 'POST' else request.GET

        preprocessed_bundle = adapter.preprocess(data)
        prediction_bundle = adapter.predict(preprocessed_bundle)
        results = adapter.postprocess(prediction_bundle) # Artık içinde her şey var

        final_probabilities_dict = results["genome_data"]
        protein_mutation_probs = results["protein_summary"]
        
        response_data = {
            "nodeId": data.get('nodeId'),
            "elapsedDay": data.get('elapsedDay'),
            "selectedModel": data.get('selectedModel'),
            "selectedProteinRegion": data.get('selectedProteinRegion'),
            "genomeSequence": preprocessed_bundle[1],
            "genomeData": list(final_probabilities_dict.values()),
            "protein_mutation_probs": protein_mutation_probs,
            "proteinRegionPossibilities": adapter.protein_regions,
        }
        
        measure_time("total_request_handling_with_adapter", start_time)
        return JsonResponse(response_data)

    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": f"Bir hata oluştu: {str(e)}"}, status=500)