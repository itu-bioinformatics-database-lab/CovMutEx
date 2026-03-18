# Future Work: Plugin Architecture

This document lists the planned follow-up work for evolving the current
CovMutEx plugin runtime into a more general mutation-visualization platform.

Current priority remains:

- make the COVID-focused plugin upload and prediction flow stable
- keep frontend compatibility with genome mutation probability visualizations
- let users upload their own model, feature extractor, helper files, and custom parameters

Items below are intentionally deferred.

## 1. Generalize Core Runtime Responsibilities

- Reduce core runtime responsibility to:
  - request parsing
  - upload bundle handling
  - plugin loading
  - plugin execution
  - frontend-compatible response building
- Move more data-domain logic out of core runtime and into plugins.
- Make genome loading, mutation parsing, annotation loading, and variant construction plugin-owned when needed.

## 2. Decouple COVID-Specific Defaults

- Keep SARS-CoV-2 defaults as optional convenience behavior, not hard requirements.
- Gradually isolate:
  - default `genome.txt`
  - default `mutations.txt`
  - default `depth_date.json`
  - default protein region boundaries
  - default mutation parsing logic
- Let plugins explicitly override these defaults without depending on core internals.

## 3. Dynamic Genome Length Support

- Remove implicit assumptions around fixed SARS-CoV-2 genome length.
- Allow plugins to return predictions for arbitrary genome lengths.
- Ensure backend aggregation and response shaping work for non-29904 genomes.
- Review frontend assumptions for chart sizing, zooming, and region mapping.

## 4. Dynamic Annotation Support

- Replace hardcoded `PROTEIN_REGIONS` dependency with plugin-provided annotations.
- Support optional annotation/config files uploaded with the plugin bundle.
- Allow different organisms or model families to define their own regions/features.
- Support future protein-level models via genome-to-protein mapping configs.

## 5. Stronger Output Contract

- Keep the current mandatory frontend contract centered on position-wise
  `A/T/G/C x position` mutation probabilities.
- Define optional output sections for:
  - genome sequence
  - annotations
  - protein boundaries
  - plugin-specific metadata
  - model card metadata
- Version the output contract separately from implementation details.

## 6. Input Contract Generalization

- Make the canonical input contract less COVID-specific.
- Keep COVID-oriented fields available as optional convenience fields.
- Support alternative sample identifiers, sequence inputs, and future non-COVID datasets.
- Define which inputs are runtime-level and which are plugin-level.

## 7. Validation and Smoke Testing

- Add a validation phase for uploaded plugin bundles.
- Validate:
  - required files
  - importability
  - contract compliance
  - output shape compatibility
- Add a small smoke-test prediction before accepting a bundle as valid.

## 8. Registry and Metadata Layer

- Add model registry APIs and persistent metadata management.
- Support:
  - model versions
  - owners/authors
  - tags
  - provenance fields
  - shareable model metadata
- Store richer plugin metadata beyond local folder conventions.

## 9. Provenance and Reproducibility

- Record which model, extractor, helper files, and parameters were used for each run.
- Track:
  - plugin contract version
  - model version
  - preprocessing configuration
  - runtime environment details
  - seed and reproducibility metadata

## 10. Benchmarking Layer

- Add benchmarking APIs and result storage.
- Support side-by-side evaluation of multiple models on common datasets.
- Include predictive and system metrics such as:
  - AUROC
  - AUPRC
  - calibration
  - runtime
  - memory
- Keep this separate from the real-time prediction endpoint.

## 11. Security Hardening

- Add stronger sandboxing for uploaded plugins.
- Introduce limits for:
  - execution time
  - memory usage
  - file size
  - per-user quotas
- Add a clearer approval model for risky plugin behaviors if needed later.

## 12. Bundle Format Evolution

- Keep the current convention-based bundle format for now.
- Consider adding a manifest file later for:
  - declared entrypoints
  - contract version
  - helper file inventory
  - annotation files
  - model capabilities
- Only introduce a manifest when it provides clear value over conventions.

## 13. Frontend Generalization

- Keep current frontend focus on genome mutation probability views.
- Later support:
  - dynamic annotation rendering
  - non-COVID region naming
  - plugin-provided overlays
  - richer model metadata display
- Preserve compatibility with existing genome visualization UX.

## 14. Multi-Domain Expansion

- Extend the architecture to support influenza and other sequence domains.
- Reuse the same plugin contract where possible.
- Add domain-specific adapters only where necessary.
- Avoid baking organism-specific rules into core runtime logic.

## 15. Developer-Facing SDK and Examples

- Publish a clearer SDK-style reference for plugin authors.
- Provide example bundles for:
  - Keras
  - PyTorch
  - sklearn
  - future non-COVID cases
- Document how helper files should be resolved from `base_dir`.
