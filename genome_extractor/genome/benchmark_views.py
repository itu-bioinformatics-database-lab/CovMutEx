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
from .cache_paths import CACHE_DIR, NODE_FEATURES_CACHE_PATH
from .feature_extractor_updated import parse_mutations
from .configs import configs
from .covmutex_models import load_model as load_covmutex_model
from .covmutex_adapters import load_model_adapter
from .covmutex_feature_extractors import load_feature_extractor
from .helpers import predict_mutations, read_genome_sequence
from .plugin_runtime import (
    BundleResolution,
    resolve_uploaded_bundle,
    resolve_server_model,
)
from .organism_registry import (
    read_builtin_organism,
    read_custom_organism_from_bundle,
    find_variant_subtype,
    read_variant,
)
from .benchmark_engine import (
    run_benchmark, build_ground_truth, compute_total_mutation_probability,
    build_ground_truth_from_csv, build_ground_truth_from_api_data,
    fetch_cov_spectrum_mutations, search_cov_spectrum_variants,
)
from .benchmark_datasets import (get_available_datasets, get_dataset, save_custom_dataset,
    set_reproducibility_seed, capture_environment, create_reproducibility_pack,
    aggregate_multi_variant_results, DEFAULT_SEED)
from .benchmark_reports import (export_json, export_csv, export_csv_multi_variant,
    export_html, export_html_multi_variant)

codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')
cache_path = NODE_FEATURES_CACHE_PATH
UPLOADED_MODELS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'uploaded_models')
BENCHMARK_RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'benchmark_results')
PROTEIN_REGIONS = {
    "ORF1ab": [266, 21555], "S": [21563, 25384], "ORF3a": [25393, 26220],
    "E": [26245, 26472], "M": [26523, 27191], "ORF6": [27202, 27387],
    "ORF7a": [27394, 27759], "ORF7b": [27756, 27887], "ORF8": [27894, 28259],
    "N": [28274, 29533], "ORF10": [29558, 29674],
}

def _resolve_bundle(model_identifier: str, request_custom_parameters: dict = None) -> BundleResolution:
    """Resolve model identifier to a BundleResolution via plugin_runtime."""
    if model_identifier.startswith('uploaded:'):
        return resolve_uploaded_bundle(UPLOADED_MODELS_DIR, model_identifier, request_custom_parameters or {})
    model_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'covid19_models', 'models')
    return resolve_server_model(model_dir, model_identifier)


def _bundle_to_config(bundle: BundleResolution) -> dict:
    """Convert a BundleResolution to the flat dict format consumed by run_benchmark."""
    return {
        'model_path': bundle.model_path,
        'name': bundle.bundle_name or bundle.model_name,
        'source': bundle.source,
        'custom_parameters': bundle.custom_parameters,
        'extractor_path': bundle.extractor_path,
        'adapter_path': bundle.adapter_path,
        'bundle_dir': bundle.bundle_dir,
        'organism': bundle.organism,
        'genome_file': bundle.genome_file,
        'protein_regions_file': bundle.protein_regions_file,
        'variant': bundle.variant,
    }


def _resolve_benchmark_genome(configs_list: list) -> tuple:
    """Return (genome_seq, protein_regions) for the benchmark run.

    Uses the first uploaded bundle's organism declaration. Falls back to the
    built-in COVID genome when all models are server-side or organism='covid'.
    Influenza bundles without a variant use the H1N1 reference by default.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    covid_genome = lambda: read_genome_sequence(os.path.join(base_dir, 'genome.txt'))

    uploaded = next((c for c in configs_list if c.get('source') == 'uploaded'), None)
    if uploaded is None:
        return covid_genome(), PROTEIN_REGIONS

    organism = uploaded.get('organism', 'covid')
    bundle_dir = uploaded.get('bundle_dir') or ''
    variant = uploaded.get('variant')

    if organism == 'custom':
        try:
            return read_custom_organism_from_bundle(
                bundle_dir,
                uploaded.get('genome_file') or '',
                uploaded.get('protein_regions_file'),
            )
        except Exception as exc:
            print(f"[BENCHMARK] custom organism load failed ({exc}), falling back to COVID")
            return covid_genome(), PROTEIN_REGIONS

    if organism in ('influenza', 'influenza_h1n1', 'influenza_h3n2', 'influenza_h5n1'):
        subtype = organism if organism != 'influenza' else 'influenza_h1n1'
        if variant:
            subtype_from_variant = find_variant_subtype(variant)
            if subtype_from_variant:
                subtype = subtype_from_variant
                try:
                    genome_seq, protein_regions = read_builtin_organism(subtype)
                    variant_seq = read_variant(subtype, variant)
                    return variant_seq, protein_regions
                except Exception:
                    pass
        try:
            return read_builtin_organism(subtype)
        except Exception as exc:
            print(f"[BENCHMARK] influenza organism load failed ({exc}), falling back to COVID")
            return covid_genome(), PROTEIN_REGIONS

    return covid_genome(), PROTEIN_REGIONS

def run_single_prediction(model_path, model_name, source, node_id, elapsed_day,
    mutations, genome_sequence, protein_regions, selected_protein_region=None,
    custom_parameters=None, extractor_path=None, adapter_path=None, base_dir=None):
    """Run prediction using the plugin architecture — honours model_adapter.py when present."""
    # Use model adapter when the bundle ships one; fall back to default Keras wrapper.
    if adapter_path and os.path.exists(adapter_path):
        print(f"[BENCHMARK] Using model adapter: {adapter_path}")
        model_wrapper = load_model_adapter(adapter_path, model_path, model_name=model_name, source=source)
    else:
        model_wrapper = load_covmutex_model(model_path=model_path, model_name=model_name, source=source)

    # Use uploaded extractor if available, otherwise default
    if extractor_path and os.path.exists(extractor_path):
        print(f"[BENCHMARK] Using uploaded extractor: {extractor_path}")
        extractor = load_feature_extractor(extractor_type='uploaded', module_path=extractor_path)
    else:
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
    payload = predict_mutations(**params)
    # predict_mutations returns a PredictionPayload dict (v2.0). Extract the
    # raw values array so benchmark_engine can compute scalar mutation probs.
    return np.asarray(payload["predictions"]["values"])

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

def _parse_models_list(data):
    """Accept models either as a JSON list or as a comma-separated string (multipart forms)."""
    raw = data.get('models', [])
    if isinstance(raw, list):
        return raw
    if isinstance(raw, str):
        # Try JSON first (frontend may send JSON-encoded string with multipart)
        try:
            parsed = json.loads(raw)
            if isinstance(parsed, list):
                return parsed
        except (ValueError, TypeError):
            pass
        return [m.strip() for m in raw.split(',') if m.strip()]
    return []


@api_view(['POST'])
def run_benchmark_view(request):
    try:
        os.makedirs(CACHE_DIR, exist_ok=True)
        data = request.data
        model_ids = _parse_models_list(data)
        node_id = data.get('nodeId')
        elapsed_day = int(data.get('elapsedDay', 60))
        region = data.get('selectedProteinRegion') or None
        seed = int(data.get('seed', DEFAULT_SEED))
        gt_source = (data.get('groundTruthSource') or 'mutations_txt').strip()
        if not model_ids: return JsonResponse({'error': 'At least 1 model required'}, status=400)
        if not node_id: return JsonResponse({'error': 'nodeId required'}, status=400)
        set_reproducibility_seed(seed)
        print(f"\n{'='*60}\n[BENCHMARK] Single | Seed:{seed} | Models:{model_ids} | GT:{gt_source}\n{'='*60}")
        configs_list = []
        for mid in model_ids:
            try:
                configs_list.append(_bundle_to_config(_resolve_bundle(mid)))
            except FileNotFoundError:
                return JsonResponse({'error': f"Model not found: {mid}"}, status=404)
        genome_seq, protein_regions_map = _resolve_benchmark_genome(configs_list)
        mutations = parse_mutations(node_id)

        # --- Optional custom ground truth (cov-spectrum.org CSV upload or API fetch) ---
        custom_gt = None
        custom_gt_full = None
        gt_meta = {}
        genome_length = len(genome_seq) if genome_seq else 29904
        region_tuple = tuple(protein_regions_map[region]) if (region and region in protein_regions_map) else None

        def _apply_gt(full_result, source_label, extra_meta):
            """Populate custom_gt / custom_gt_full / gt_meta from a build_ground_truth_* result dict."""
            nonlocal custom_gt, custom_gt_full, gt_meta
            custom_gt_full = full_result['ground_truth']
            if region_tuple is not None:
                # Slice to region (we built full-genome above; slice here)
                start, end = region_tuple
                custom_gt = custom_gt_full[start:end + 1]
            else:
                custom_gt = custom_gt_full
            gt_meta = {
                'source': source_label,
                'num_mutations_parsed': full_result['num_mutations_parsed'],
                'num_mutations_skipped': full_result['num_mutations_skipped'],
                'skipped_examples': full_result['skipped_examples'],
                **extra_meta,
            }

        if gt_source == 'cov_spectrum_csv':
            csv_file = request.FILES.get('groundTruthCsv') if hasattr(request, 'FILES') else None
            if not csv_file:
                return JsonResponse({'error': 'groundTruthCsv file required when groundTruthSource=cov_spectrum_csv'}, status=400)
            try:
                csv_content = csv_file.read().decode('utf-8', errors='replace')
            except Exception as e:
                return JsonResponse({'error': f'Failed to read CSV: {e}'}, status=400)

            full_result = build_ground_truth_from_csv(csv_content, genome_length=genome_length)
            _apply_gt(full_result, 'cov_spectrum_csv', {
                'filename': getattr(csv_file, 'name', 'uploaded.csv'),
            })

        elif gt_source == 'cov_spectrum_api':
            # Required: variant (Pango lineage). Dates are optional — empty = all time.
            lineage = (data.get('covSpectrumLineage') or '').strip() or None
            date_from = (data.get('covSpectrumDateFrom') or '').strip() or None
            date_to = (data.get('covSpectrumDateTo') or '').strip() or None

            if not lineage:
                return JsonResponse({
                    'error': 'A variant (Pango lineage) is required. Use the search to pick one.'
                }, status=400)

            try:
                api_result = fetch_cov_spectrum_mutations(
                    lineage=lineage, date_from=date_from, date_to=date_to,
                )
            except RuntimeError as e:
                return JsonResponse({'error': f'cov-spectrum API: {e}'}, status=502)

            if not api_result['data']:
                return JsonResponse({
                    'error': f'cov-spectrum returned 0 mutations for variant "{lineage}"'
                             + (f' between {date_from} and {date_to}' if (date_from or date_to) else '')
                             + '. Try a different variant or widen the date range.',
                    'query': api_result['query_params'],
                }, status=404)

            full_result = build_ground_truth_from_api_data(api_result['data'], genome_length=genome_length)
            _apply_gt(full_result, 'cov_spectrum_api', {
                'lineage': lineage,
                'date_from': date_from,
                'date_to': date_to,
                'query_params': api_result['query_params'],
                'api_url': api_result['url'],
                'num_api_rows': len(api_result['data']),
            })

        results = run_benchmark(
            configs_list, node_id, elapsed_day, mutations, genome_seq,
            protein_regions_map, run_single_prediction, region,
            custom_ground_truth=custom_gt,
            custom_ground_truth_full=custom_gt_full,
            ground_truth_source=gt_source,
            ground_truth_meta=gt_meta,
        )
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
        os.makedirs(CACHE_DIR, exist_ok=True)
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
                configs_list.append(_bundle_to_config(_resolve_bundle(mid)))
            except FileNotFoundError:
                return JsonResponse({'error': f"Model not found: {mid}"}, status=404)
        genome_seq, protein_regions_map = _resolve_benchmark_genome(configs_list)
        per_variant = []
        for i, v in enumerate(variants):
            print(f"\n[DS BENCH] Variant {i+1}/{len(variants)}: {v.get('label','')}")
            mutations = parse_mutations(v['nodeId'])
            vr = run_benchmark(configs_list, v['nodeId'], v['elapsedDay'], mutations, genome_seq, protein_regions_map, run_single_prediction, region)
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


@api_view(['GET'])
def search_cov_spectrum_variants_view(request):
    """
    Proxy to cov-spectrum.org LAPIS — search Pango lineages by substring.
    Query params:
        q: search string (empty = top-N most common)
        limit: max results (default 50)
    """
    query = request.GET.get('q', '').strip()
    try:
        limit = max(1, min(10000, int(request.GET.get('limit', 50))))
    except (TypeError, ValueError):
        limit = 50
    try:
        variants = search_cov_spectrum_variants(query, limit=limit)
    except RuntimeError as e:
        return JsonResponse({'error': str(e)}, status=502)
    return JsonResponse({'variants': variants, 'query': query, 'count': len(variants)}, status=200)


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
