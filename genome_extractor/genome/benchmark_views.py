"""
CovMutEx-X Benchmark API Views - Complete FR-3

Endpoints:
- POST /api/benchmark/run/          Run single-variant benchmark
- POST /api/benchmark/run-dataset/  Run multi-variant dataset benchmark
- GET  /api/benchmark/results/      Get stored benchmark results
- GET  /api/benchmark/list/         List all past benchmark runs
- GET  /api/benchmark/datasets/     List available benchmark datasets
- POST /api/benchmark/datasets/     Create custom dataset
- GET  /api/benchmark/export/       Export results as JSON/CSV/HTML
"""
import os, json, time, traceback
import numpy as np
from django.http import JsonResponse, HttpResponse
from rest_framework.decorators import api_view
from .feature_extractor_updated import parse_mutations
from .configs import configs
from .covmutex_models import load_model as load_covmutex_model
from .covmutex_feature_extractors import load_feature_extractor
from .helpers import predict_mutations, read_genome_sequence
from .benchmark_engine import run_benchmark, build_ground_truth, compute_total_mutation_probability
from .benchmark_datasets import (get_available_datasets, get_dataset, save_custom_dataset,
    set_reproducibility_seed, capture_environment, create_reproducibility_pack,
    aggregate_multi_variant_results, DEFAULT_SEED)
from .benchmark_reports import (export_json, export_csv, export_csv_multi_variant,
    export_html, export_html_multi_variant)

codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')
cache_path = os.path.join(os.path.dirname(__file__), 'node_features.h5')
UPLOADED_MODELS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'uploaded_models')
BENCHMARK_RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'benchmark_results')
PROTEIN_REGIONS = {
    "ORF1ab": [266, 21555], "S": [21563, 25384], "ORF3a": [25393, 26220],
    "E": [26245, 26472], "M": [26523, 27191], "ORF6": [27202, 27387],
    "ORF7a": [27394, 27759], "ORF7b": [27756, 27887], "ORF8": [27894, 28259],
    "N": [28274, 29533], "ORF10": [29558, 29674],
}

def resolve_model_path(model_identifier):
    if model_identifier.startswith('uploaded:'):
        folder_name = model_identifier.replace('uploaded:', '')
        folder_path = os.path.join(UPLOADED_MODELS_DIR, folder_name)
        if not os.path.exists(folder_path):
            raise FileNotFoundError(f"Uploaded model folder not found: {folder_name}")
        model_path = None
        for f in os.listdir(folder_path):
            if f.startswith('model.') and f.endswith(('.keras', '.h5', '.pt', '.pth')):
                model_path = os.path.join(folder_path, f)
                break
        if not model_path:
            for f in os.listdir(folder_path):
                if f.endswith(('.keras', '.h5', '.pt', '.pth')):
                    model_path = os.path.join(folder_path, f)
                    break
        if not model_path:
            raise FileNotFoundError(f"No model file in: {folder_name}")
        params = {}
        params_file = os.path.join(folder_path, 'custom_parameters.json')
        if os.path.exists(params_file):
            with open(params_file, 'r') as f:
                raw = json.load(f)
                for k, v in raw.items():
                    params[k] = v['value'] if isinstance(v, dict) and 'value' in v else v
        return {'path': model_path, 'name': folder_name, 'source': 'uploaded', 'custom_parameters': params}
    else:
        model_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'covid19_models', 'models')
        model_path = os.path.join(model_dir, f"{model_identifier}.keras")
        if not os.path.exists(model_path):
            for ext in ['.h5', '.pt', '.pth']:
                alt = os.path.join(model_dir, f"{model_identifier}{ext}")
                if os.path.exists(alt):
                    model_path = alt
                    break
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Server model not found: {model_identifier}")
        return {'path': model_path, 'name': model_identifier, 'source': 'server', 'custom_parameters': {}}

def run_single_prediction(model_path, model_name, source, node_id, elapsed_day,
    mutations, genome_sequence, protein_regions, selected_protein_region=None, custom_parameters=None):
    model_wrapper = load_covmutex_model(model_path=model_path, model_name=model_name, source=source)
    extractor = load_feature_extractor(extractor_type='default')
    pr_dict = {}
    if selected_protein_region and selected_protein_region in protein_regions:
        pr_dict = {selected_protein_region: protein_regions[selected_protein_region]}
    params = {
        'cache_path': cache_path, 'genome_seq': genome_sequence, 'mutations': mutations,
        'codon_mapper': codon_mapping_path, 'config_file': configs(), 'node_ids': node_id,
        'elapsed_day': elapsed_day, 'protein_regions': pr_dict, 'selectedModel': model_name,
        'model_wrapper': model_wrapper, 'feature_extractor': extractor,
    }
    if custom_parameters:
        params.update(custom_parameters)
    return predict_mutations(**params)

def _clean_nan(obj):
    """Recursively replace NaN/Inf with None for JSON compatibility."""
    if isinstance(obj, dict):
        return {k: _clean_nan(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_clean_nan(v) for v in obj]
    elif isinstance(obj, float):
        if obj != obj or obj == float('inf') or obj == float('-inf'):  # NaN check
            return None
        return obj
    elif isinstance(obj, (np.floating,)):
        val = float(obj)
        if val != val or val == float('inf') or val == float('-inf'):
            return None
        return val
    elif isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, np.ndarray):
        return _clean_nan(obj.tolist())
    return obj

def _save(results, prefix=""):
    os.makedirs(BENCHMARK_RESULTS_DIR, exist_ok=True)
    bid = results.get('benchmark_id', f"bench_{int(time.time())}")
    fn = f"{prefix}{bid}.json" if prefix else f"{bid}.json"
    path = os.path.join(BENCHMARK_RESULTS_DIR, fn)
    cleaned = _clean_nan(results)
    with open(path, 'w') as f:
        json.dump(cleaned, f, indent=2)
    return path

def _serialize(results):
    return _clean_nan(results)

@api_view(['POST'])
def run_benchmark_view(request):
    try:
        data = request.data
        model_ids = data.get('models', [])
        node_id = data.get('nodeId')
        elapsed_day = int(data.get('elapsedDay', 60))
        region = data.get('selectedProteinRegion')
        seed = int(data.get('seed', DEFAULT_SEED))
        if not model_ids: return JsonResponse({'error': 'At least 1 model required'}, status=400)
        if not node_id: return JsonResponse({'error': 'nodeId required'}, status=400)
        set_reproducibility_seed(seed)
        print(f"\n{'='*60}\n[BENCHMARK] Single | Seed:{seed} | Models:{model_ids}\n{'='*60}")
        configs_list = []
        for mid in model_ids:
            try:
                c = resolve_model_path(mid)
                c['model_path'] = c.pop('path')
                configs_list.append(c)
            except FileNotFoundError as e:
                return JsonResponse({'error': f"Model not found: {mid}"}, status=404)
        base_dir = os.path.dirname(os.path.abspath(__file__))
        genome_seq = read_genome_sequence(os.path.join(base_dir, 'genome.txt'))
        mutations = parse_mutations(node_id)
        results = run_benchmark(configs_list, node_id, elapsed_day, mutations, genome_seq, PROTEIN_REGIONS, run_single_prediction, region)
        results['reproducibility'] = create_reproducibility_pack(results, seed)
        results['environment'] = capture_environment()
        _save(results)
        return JsonResponse(_serialize(results), status=200)
    except Exception as e:
        print(f"[BENCHMARK ERROR] {traceback.format_exc()}")
        return JsonResponse({'error': f'Benchmark failed: {str(e)}'}, status=500)

@api_view(['POST'])
def run_dataset_benchmark_view(request):
    try:
        data = request.data
        model_ids = data.get('models', [])
        dataset_id = data.get('datasetId')
        region = data.get('selectedProteinRegion')
        seed = int(data.get('seed', DEFAULT_SEED))
        if not model_ids: return JsonResponse({'error': 'At least 1 model required'}, status=400)
        if not dataset_id: return JsonResponse({'error': 'datasetId required'}, status=400)
        dataset = get_dataset(dataset_id)
        if not dataset:
            return JsonResponse({'error': f'Dataset "{dataset_id}" not found', 'available': list(get_available_datasets().keys())}, status=404)
        variants = dataset.get('variants', [])
        if not variants: return JsonResponse({'error': 'Dataset has no variants'}, status=400)
        set_reproducibility_seed(seed)
        print(f"\n{'='*60}\n[DS BENCH] {dataset.get('name')} | {len(variants)} variants | Seed:{seed}\n{'='*60}")
        configs_list = []
        for mid in model_ids:
            try:
                c = resolve_model_path(mid)
                c['model_path'] = c.pop('path')
                configs_list.append(c)
            except FileNotFoundError as e:
                return JsonResponse({'error': f"Model not found: {mid}"}, status=404)
        base_dir = os.path.dirname(os.path.abspath(__file__))
        genome_seq = read_genome_sequence(os.path.join(base_dir, 'genome.txt'))
        per_variant = []
        for i, v in enumerate(variants):
            print(f"\n[DS BENCH] Variant {i+1}/{len(variants)}: {v.get('label','')}")
            mutations = parse_mutations(v['nodeId'])
            vr = run_benchmark(configs_list, v['nodeId'], v['elapsedDay'], mutations, genome_seq, PROTEIN_REGIONS, run_single_prediction, region)
            vr['variant_label'] = v.get('label', f"Variant {i+1}")
            per_variant.append(vr)
        aggregated = aggregate_multi_variant_results(per_variant)
        final = {
            'benchmark_id': f"ds_bench_{int(time.time())}", 'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'type': 'dataset', 'dataset': {'id': dataset_id, 'name': dataset.get('name'), 'description': dataset.get('description'), 'num_variants': len(variants)},
            'aggregated': aggregated, 'per_variant': per_variant,
            'reproducibility': create_reproducibility_pack({'parameters': {'dataset_id': dataset_id}, 'models': {c['name']: {'source': c['source']} for c in configs_list}}, seed),
            'environment': capture_environment(),
        }
        _save(final, "ds_")
        return JsonResponse(_serialize(final), status=200)
    except Exception as e:
        print(f"[DS BENCH ERROR] {traceback.format_exc()}")
        return JsonResponse({'error': str(e)}, status=500)

@api_view(['GET', 'POST'])
def benchmark_datasets_view(request):
    if request.method == 'GET':
        return JsonResponse({'datasets': get_available_datasets()}, status=200)
    else:
        try:
            data = request.data
            did = data.get('id')
            if not did: return JsonResponse({'error': 'id required'}, status=400)
            ds = {'name': data.get('name', did), 'description': data.get('description', ''), 'category': data.get('category', 'custom'), 'variants': data.get('variants', [])}
            if not ds['variants']: return JsonResponse({'error': 'variants required'}, status=400)
            save_custom_dataset(did, ds)
            return JsonResponse({'message': f'Dataset "{did}" created', 'dataset': ds}, status=201)
        except Exception as e:
            return JsonResponse({'error': str(e)}, status=500)

@api_view(['GET', 'POST'])
def export_benchmark_view(request):
    try:
        # Accept params from GET query string or POST body
        if request.method == 'POST':
            bid = request.data.get('id')
            fmt = request.data.get('format', 'json').lower()
        else:
            bid = request.GET.get('id')
            fmt = request.GET.get('format', 'json').lower()
        if not bid: return JsonResponse({'error': 'id required'}, status=400)
        
        print(f"[EXPORT] Requested: id={bid}, format={fmt}")
        print(f"[EXPORT] Results dir: {BENCHMARK_RESULTS_DIR}")
        
        # List files in results dir for debug
        os.makedirs(BENCHMARK_RESULTS_DIR, exist_ok=True)
        files_in_dir = os.listdir(BENCHMARK_RESULTS_DIR)
        print(f"[EXPORT] Files in dir: {files_in_dir}")
        
        # Try exact match first
        path = os.path.join(BENCHMARK_RESULTS_DIR, f"{bid}.json")
        if not os.path.exists(path):
            # Try with ds_ prefix for dataset benchmarks
            path = os.path.join(BENCHMARK_RESULTS_DIR, f"ds_{bid}.json")
        if not os.path.exists(path):
            # Try searching for file containing the benchmark_id
            for fn in files_in_dir:
                if bid in fn and fn.endswith('.json'):
                    path = os.path.join(BENCHMARK_RESULTS_DIR, fn)
                    print(f"[EXPORT] Found via search: {fn}")
                    break
        if not os.path.exists(path):
            print(f"[EXPORT] NOT FOUND: {bid}")
            return JsonResponse({'error': f'{bid} not found', 'available_files': files_in_dir}, status=404)
        
        print(f"[EXPORT] Loading: {path}")
        with open(path, 'r') as f:
            results = json.load(f)
        
        is_ds = results.get('type') == 'dataset'
        
        if fmt == 'json':
            content = export_json(results)
            r = HttpResponse(content, content_type='application/json')
            r['Content-Disposition'] = f'attachment; filename="{bid}.json"'
        elif fmt == 'csv':
            content = export_csv_multi_variant(results.get('aggregated', {})) if is_ds else export_csv(results)
            r = HttpResponse(content, content_type='text/csv; charset=utf-8')
            r['Content-Disposition'] = f'attachment; filename="{bid}.csv"'
        elif fmt == 'html':
            if is_ds:
                content = export_html_multi_variant(results.get('aggregated', {}), results.get('per_variant', []))
            else:
                content = export_html(results)
            r = HttpResponse(content, content_type='text/html; charset=utf-8')
            r['Content-Disposition'] = f'attachment; filename="{bid}.html"'
        else:
            return JsonResponse({'error': f'Unknown format: {fmt}', 'supported': ['json','csv','html']}, status=400)
        
        # Add CORS headers
        r['Access-Control-Allow-Origin'] = '*'
        print(f"[EXPORT] Success: {fmt}, size={len(content)} bytes")
        return r
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

@api_view(['GET'])
def get_benchmark_results(request):
    try:
        bid = request.GET.get('id')
        if not bid: return JsonResponse({'error': 'id required'}, status=400)
        path = os.path.join(BENCHMARK_RESULTS_DIR, f"{bid}.json")
        if not os.path.exists(path):
            path = os.path.join(BENCHMARK_RESULTS_DIR, f"ds_{bid}.json")
        if not os.path.exists(path):
            return JsonResponse({'error': f'{bid} not found'}, status=404)
        with open(path, 'r') as f:
            return JsonResponse(json.load(f), status=200)
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

@api_view(['GET'])
def list_benchmarks(request):
    try:
        os.makedirs(BENCHMARK_RESULTS_DIR, exist_ok=True)
        benchmarks = []
        for fn in sorted(os.listdir(BENCHMARK_RESULTS_DIR), reverse=True):
            if not fn.endswith('.json'): continue
            try:
                with open(os.path.join(BENCHMARK_RESULTS_DIR, fn), 'r') as f:
                    d = json.load(f)
                e = {'benchmark_id': d.get('benchmark_id'), 'timestamp': d.get('timestamp'), 'type': d.get('type', 'single')}
                if d.get('type') == 'dataset':
                    e['dataset_name'] = d.get('dataset', {}).get('name')
                    e['num_variants'] = d.get('dataset', {}).get('num_variants')
                    e['model_names'] = list(d.get('aggregated', {}).get('models', {}).keys())
                else:
                    e['parameters'] = d.get('parameters', {})
                    e['model_names'] = list(d.get('models', {}).keys())
                    e['summary'] = d.get('summary', {})
                e['num_models'] = len(e.get('model_names', []))
                benchmarks.append(e)
            except: continue
        return JsonResponse({'benchmarks': benchmarks, 'total': len(benchmarks)}, status=200)
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)