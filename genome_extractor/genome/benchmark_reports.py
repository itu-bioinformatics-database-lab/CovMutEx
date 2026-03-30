"""
CovMutEx-X Benchmark Report Exporter

FR-3.4: Reports -- Static JSON + HTML; exportable CSV.
"""
import os
import json
import csv
import io
import time
from typing import Dict, Any


def _fmt(val, decimals=4):
    """Safely format a numeric value, returning '—' for None/NaN."""
    if val is None:
        return "—"
    try:
        return f"{float(val):.{decimals}f}"
    except (ValueError, TypeError):
        return str(val)


def _fmt_int(val):
    """Safely format an integer value with comma separator."""
    if val is None:
        return "N/A"
    try:
        return f"{int(val):,}"
    except (ValueError, TypeError):
        return str(val)


def export_json(results: Dict) -> str:
    return json.dumps(results, indent=2, default=str)


def export_csv(results: Dict) -> str:
    output = io.StringIO()
    models = results.get("models", {})
    if not models:
        output.write("No models found in benchmark results\n")
        return output.getvalue()

    non_scalar = {'roc_curve', 'pr_curve', 'calibration_curve'}
    metric_keys = set()
    for model_data in models.values():
        if model_data.get("metrics"):
            metric_keys.update(model_data["metrics"].keys())
    metric_keys = sorted(metric_keys - non_scalar)

    writer = csv.writer(output)
    writer.writerow(['model_name', 'source', 'status', 'runtime_seconds'] + metric_keys)

    for model_name, model_data in models.items():
        row = [
            model_name,
            model_data.get('source', ''),
            model_data.get('status', ''),
            _fmt(model_data.get('runtime_seconds'), 3),
        ]
        for mk in metric_keys:
            val = model_data.get('metrics', {}).get(mk)
            row.append(_fmt(val, 6) if isinstance(val, (int, float)) else str(val or ''))
        writer.writerow(row)

    # Per-protein section
    first_model = next(iter(models.values()), {})
    if first_model.get('per_protein'):
        writer.writerow([])
        writer.writerow(['--- Per-Protein Region Breakdown ---'])
        writer.writerow(['protein_region', 'model_name', 'mean_prediction', 'max_prediction',
                        'num_positions', 'num_mutations', 'mutation_rate', 'auroc', 'auprc', 'brier_score'])
        for model_name, model_data in models.items():
            for protein, pp in model_data.get('per_protein', {}).items():
                writer.writerow([
                    protein, model_name,
                    _fmt(pp.get('mean_prediction'), 6),
                    _fmt(pp.get('max_prediction'), 6),
                    pp.get('num_positions', ''),
                    pp.get('num_mutations', ''),
                    _fmt(pp.get('mutation_rate'), 6),
                    _fmt(pp.get('auroc'), 4),
                    _fmt(pp.get('auprc'), 4),
                    _fmt(pp.get('brier_score'), 6),
                ])

    return output.getvalue()


def export_csv_multi_variant(aggregated_results: Dict) -> str:
    output = io.StringIO()
    writer = csv.writer(output)
    models = aggregated_results.get("models", {})
    if not models:
        output.write("No aggregated results found\n")
        return output.getvalue()

    metric_keys = set()
    for model_data in models.values():
        if model_data.get("metrics"):
            metric_keys.update(model_data["metrics"].keys())
    metric_keys = sorted(metric_keys)

    header = ['model_name', 'success_rate', 'num_variants']
    for mk in metric_keys:
        header.extend([f'{mk}_mean', f'{mk}_std', f'{mk}_min', f'{mk}_max'])
    writer.writerow(header)

    for model_name, model_data in models.items():
        row = [
            model_name,
            _fmt(model_data.get('success_rate', 0), 2),
            model_data.get('num_success', 0),
        ]
        for mk in metric_keys:
            md = model_data.get('metrics', {}).get(mk)
            if md and isinstance(md, dict):
                row.extend([_fmt(md.get('mean'), 6), _fmt(md.get('std'), 6),
                           _fmt(md.get('min'), 6), _fmt(md.get('max'), 6)])
            else:
                row.extend(['', '', '', ''])
        writer.writerow(row)

    rankings = aggregated_results.get("rankings", {})
    if rankings:
        writer.writerow([])
        writer.writerow(['--- Rankings ---'])
        writer.writerow(['metric', 'rank', 'model', 'value'])
        for metric, ranked in rankings.items():
            for entry in ranked:
                writer.writerow([metric, entry['rank'], entry['model'], _fmt(entry['value'], 6)])

    return output.getvalue()


def export_html(results: Dict, title: str = "CovMutEx-X Benchmark Report") -> str:
    models = results.get("models", {})
    params = results.get("parameters", {})
    agreement = results.get("model_agreement", {})
    model_names = [k for k, v in models.items() if v.get("status") == "success"]
    colors = ["#3B82F6", "#EF4444", "#10B981", "#F59E0B", "#8B5CF6"]

    gt_total = params.get('ground_truth_total')
    gt_positives = params.get('ground_truth_positives')
    node_id_str = str(params.get('node_id', 'N/A'))
    if len(node_id_str) > 40:
        node_id_str = node_id_str[:40] + "..."

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         background: #f8fafc; color: #1e293b; padding: 2rem; }}
  .container {{ max-width: 1100px; margin: 0 auto; }}
  h1 {{ font-size: 1.8rem; color: #1e40af; margin-bottom: 0.5rem; }}
  h2 {{ font-size: 1.3rem; color: #334155; margin: 1.5rem 0 0.8rem; padding-bottom: 0.4rem;
        border-bottom: 2px solid #e2e8f0; }}
  .meta {{ color: #64748b; font-size: 0.85rem; margin-bottom: 1.5rem; }}
  .meta span {{ margin-right: 1.5rem; }}
  .card {{ background: white; border-radius: 12px; padding: 1.5rem; margin-bottom: 1.5rem;
           box-shadow: 0 1px 3px rgba(0,0,0,0.08); border: 1px solid #e2e8f0; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
  th {{ background: #f1f5f9; padding: 0.6rem 0.8rem; text-align: left; font-weight: 600;
       color: #475569; text-transform: uppercase; font-size: 0.75rem; }}
  td {{ padding: 0.6rem 0.8rem; border-bottom: 1px solid #f1f5f9; }}
  tr:hover td {{ background: #f8fafc; }}
  .best {{ color: #059669; font-weight: 700; }}
  .badge {{ display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 0.7rem;
            font-weight: 600; color: white; }}
  .metric-note {{ font-size: 0.75rem; color: #94a3b8; margin-top: 0.5rem; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 1rem; margin-bottom: 1.5rem; }}
  .summary-card {{ background: white; border-radius: 10px; padding: 1rem; text-align: center;
                   border: 1px solid #e2e8f0; }}
  .summary-card .value {{ font-size: 1.5rem; font-weight: 700; color: #1e40af; }}
  .summary-card .label {{ font-size: 0.75rem; color: #64748b; text-transform: uppercase; margin-top: 0.3rem; }}
  .footer {{ text-align: center; color: #94a3b8; font-size: 0.75rem; margin-top: 2rem; padding-top: 1rem;
             border-top: 1px solid #e2e8f0; }}
  .corr-high {{ background: #ecfdf5; color: #059669; }}
  .corr-mid {{ background: #fffbeb; color: #d97706; }}
  .corr-low {{ background: #fef2f2; color: #dc2626; }}
  @media print {{
    body {{ padding: 1rem; }}
    .card {{ box-shadow: none; border: 1px solid #ccc; }}
  }}
</style>
</head>
<body>
<div class="container">

<h1>{title}</h1>
<div class="meta">
  <span>ID: {results.get('benchmark_id', 'N/A')}</span>
  <span>Date: {results.get('timestamp', 'N/A')}</span>
  <span>Node: {node_id_str}</span>
  <span>Elapsed: {params.get('elapsed_day', 'N/A')} days</span>
  <span>Mutations: {params.get('num_mutations', 'N/A')}</span>
</div>

<div class="summary-grid">
  <div class="summary-card">
    <div class="value">{len(model_names)}</div>
    <div class="label">Models Compared</div>
  </div>
  <div class="summary-card">
    <div class="value">{_fmt_int(gt_positives)}</div>
    <div class="label">Mutation Positions</div>
  </div>
  <div class="summary-card">
    <div class="value">{_fmt_int(gt_total)}</div>
    <div class="label">Total Positions</div>
  </div>
  <div class="summary-card">
    <div class="value">{params.get('selected_protein_region') or 'Full Genome'}</div>
    <div class="label">Region</div>
  </div>
</div>
"""

    # Metrics Table
    html += '<div class="card">\n<h2>Metrics Comparison</h2>\n<table>\n<thead><tr><th>Metric</th>\n'
    for i, name in enumerate(model_names):
        c = colors[i % len(colors)]
        html += f'<th><span class="badge" style="background:{c}">{name}</span></th>\n'
    html += '<th>Best</th>\n</tr></thead>\n<tbody>\n'

    metrics_info = [
        ("AUROC", "auroc", True), ("AUPRC", "auprc", True),
        ("Brier Score", "brier_score", False), ("ECE", "ece", False),
        ("Runtime (s)", "runtime_seconds", False),
        ("Mean Prediction", "mean_prediction", None), ("Std Deviation", "std_prediction", None),
        ("Num Positions", "num_positions", None), ("Num Mutations", "num_mutations", None),
    ]

    for label, key, higher_better in metrics_info:
        values = [models[n].get("metrics", {}).get(key) for n in model_names]
        valid = [(i, v) for i, v in enumerate(values) if v is not None and isinstance(v, (int, float))]
        best_idx = -1
        if higher_better is not None and valid:
            best_idx = (max if higher_better else min)(valid, key=lambda x: x[1])[0]

        html += f'<tr><td><strong>{label}</strong></td>\n'
        for i, val in enumerate(values):
            cls = ' class="best"' if i == best_idx else ''
            star = " ★" if i == best_idx else ""
            html += f'<td{cls}>{_fmt(val)}{star}</td>\n'
        direction = "Higher" if higher_better else ("Lower" if higher_better is False else "—")
        html += f'<td style="color:#94a3b8;font-size:0.75rem">{direction}</td>\n</tr>\n'

    html += '</tbody></table>\n'
    html += '<p class="metric-note">★ = best performer. Higher/Lower indicates preferred direction.</p>\n</div>\n'

    # Per-Protein Table
    first_model = model_names[0] if model_names else None
    if first_model and models[first_model].get("per_protein"):
        html += '<div class="card">\n<h2>Per-Protein Region Breakdown</h2>\n<table>\n<thead><tr>\n'
        html += '<th>Region</th><th>Positions</th><th>Mutations</th>\n'
        for i, name in enumerate(model_names):
            c = colors[i % len(colors)]
            html += f'<th colspan="2"><span class="badge" style="background:{c}">{name}</span> Mean / AUROC</th>\n'
        html += '</tr></thead>\n<tbody>\n'
        for protein in models[first_model]["per_protein"]:
            pp0 = models[first_model]["per_protein"][protein]
            html += f'<tr><td><strong>{protein}</strong></td>'
            html += f'<td>{pp0.get("num_positions", "")}</td>'
            html += f'<td>{pp0.get("num_mutations", "")}</td>'
            for name in model_names:
                pp = models[name].get("per_protein", {}).get(protein, {})
                html += f'<td>{_fmt(pp.get("mean_prediction"))}</td>'
                html += f'<td>{_fmt(pp.get("auroc"), 3)}</td>'
            html += '</tr>\n'
        html += '</tbody></table>\n</div>\n'

    # Model Agreement
    if agreement and len(model_names) > 1:
        html += '<div class="card">\n<h2>Model Agreement (Pearson Correlation)</h2>\n<table>\n<thead><tr><th></th>\n'
        for name in model_names:
            html += f'<th>{name}</th>\n'
        html += '</tr></thead>\n<tbody>\n'
        for rn in model_names:
            html += f'<tr><td><strong>{rn}</strong></td>\n'
            for cn in model_names:
                val = agreement.get(rn, {}).get(cn)
                if rn == cn:
                    html += '<td style="background:#f1f5f9;color:#94a3b8">1.000</td>\n'
                elif val is not None:
                    cls = "corr-high" if val > 0.8 else ("corr-mid" if val > 0.5 else "corr-low")
                    html += f'<td class="{cls}">{_fmt(val, 3)}</td>\n'
                else:
                    html += '<td>—</td>\n'
            html += '</tr>\n'
        html += '</tbody></table>\n</div>\n'

    # Reproducibility info
    repro = results.get('reproducibility', {})
    env = results.get('environment', {})
    if repro or env:
        html += '<div class="card">\n<h2>Reproducibility</h2>\n'
        html += f'<p><strong>Seed:</strong> {repro.get("seed", "N/A")}</p>\n'
        if env:
            html += f'<p><strong>Python:</strong> {env.get("python_version", "N/A")} | '
            html += f'<strong>OS:</strong> {env.get("os", "N/A")} | '
            pkgs = env.get("packages", {})
            html += f'<strong>TensorFlow:</strong> {pkgs.get("tensorflow", "N/A")} | '
            html += f'<strong>NumPy:</strong> {pkgs.get("numpy", "N/A")}</p>\n'
        html += '</div>\n'

    html += f"""
<div class="footer">
  Generated by CovMutEx-X Benchmark Suite | {time.strftime('%Y-%m-%d %H:%M:%S')}
</div>
</div></body></html>"""

    return html


def export_html_multi_variant(aggregated: Dict, per_variant: list,
                               title: str = "Multi-Variant Benchmark Report") -> str:
    models_agg = aggregated.get("models", {})
    rankings = aggregated.get("rankings", {})
    model_names = list(models_agg.keys())
    colors = ["#3B82F6", "#EF4444", "#10B981", "#F59E0B", "#8B5CF6"]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>{title}</title>
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ font-family: -apple-system, sans-serif; background:#f8fafc; color:#1e293b; padding:2rem; }}
  .container {{ max-width:1100px; margin:0 auto; }}
  h1 {{ font-size:1.8rem; color:#1e40af; margin-bottom:1rem; }}
  h2 {{ font-size:1.3rem; color:#334155; margin:1.5rem 0 0.8rem; border-bottom:2px solid #e2e8f0; padding-bottom:0.4rem; }}
  .card {{ background:white; border-radius:12px; padding:1.5rem; margin-bottom:1.5rem; box-shadow:0 1px 3px rgba(0,0,0,0.08); border:1px solid #e2e8f0; }}
  table {{ width:100%; border-collapse:collapse; font-size:0.85rem; }}
  th {{ background:#f1f5f9; padding:0.6rem 0.8rem; text-align:left; font-weight:600; color:#475569; text-transform:uppercase; font-size:0.75rem; }}
  td {{ padding:0.6rem 0.8rem; border-bottom:1px solid #f1f5f9; }}
  .best {{ color:#059669; font-weight:700; }}
  .badge {{ display:inline-block; padding:2px 8px; border-radius:12px; font-size:0.7rem; font-weight:600; color:white; }}
  .rank-1 {{ background:#fef3c7; font-weight:700; }}
  .footer {{ text-align:center; color:#94a3b8; font-size:0.75rem; margin-top:2rem; padding-top:1rem; border-top:1px solid #e2e8f0; }}
</style>
</head>
<body>
<div class="container">
<h1>{title}</h1>
<p style="color:#64748b;margin-bottom:1.5rem;">
  {aggregated.get('num_variants', 0)} variants | {aggregated.get('total_models', 0)} models
</p>

<div class="card">
<h2>Aggregated Metrics (Mean +/- Std across Variants)</h2>
<table><thead><tr><th>Metric</th>
"""
    for i, name in enumerate(model_names):
        c = colors[i % len(colors)]
        html += f'<th><span class="badge" style="background:{c}">{name}</span></th>\n'
    html += '</tr></thead>\n<tbody>\n'

    for mk in ['auroc', 'auprc', 'brier_score', 'ece', 'runtime_seconds']:
        html += f'<tr><td><strong>{mk.upper()}</strong></td>\n'
        for name in model_names:
            md = models_agg.get(name, {}).get("metrics", {}).get(mk)
            if md and isinstance(md, dict) and md.get("mean") is not None:
                html += f'<td>{_fmt(md["mean"])} +/- {_fmt(md["std"])}</td>\n'
            else:
                html += '<td>—</td>\n'
        html += '</tr>\n'
    html += '</tbody></table>\n</div>\n'

    # Rankings
    if rankings:
        html += '<div class="card">\n<h2>Rankings</h2>\n<table>\n'
        html += '<thead><tr><th>Metric</th><th>Rank</th><th>Model</th><th>Score</th></tr></thead>\n<tbody>\n'
        for metric, ranked in rankings.items():
            for entry in ranked:
                cls = ' class="rank-1"' if entry['rank'] == 1 else ''
                medal = ["", "1st", "2nd", "3rd"]
                m = medal[entry['rank']] if entry['rank'] <= 3 else f"#{entry['rank']}"
                html += f'<tr{cls}><td>{metric.upper()}</td><td>{m}</td><td>{entry["model"]}</td><td>{_fmt(entry["value"], 6)}</td></tr>\n'
        html += '</tbody></table>\n</div>\n'

    html += f"""
<div class="footer">
  Generated by CovMutEx-X Benchmark Suite | {time.strftime('%Y-%m-%d %H:%M:%S')}
</div>
</div></body></html>"""

    return html