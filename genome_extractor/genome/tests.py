"""Unit tests for the v2.0 PredictionPayload contract.

These cover the two public helpers exported from covmutex_models:

    * wrap_predictions_as_payload(...) — builds payloads
    * validate_prediction_payload(...) — enforces the schema

The bundle adapters and viewsUpdated.run_prediction() both rely on these
invariants holding; if a test here breaks, those paths will too.
"""

import unittest

import numpy as np
from django.test import SimpleTestCase

from .covmutex_models import (
    MODEL_CONTRACT_VERSION,
    validate_prediction_payload,
    wrap_predictions_as_payload,
)


def _base_payload(overrides=None):
    """Return a fresh, valid binary_per_position payload, then apply dotted-path
    overrides (e.g. {"task.kind": "scalar_per_position"})."""
    payload = {
        "contract_version": "2.0",
        "task": {"kind": "binary_per_position"},
        "domain": {"total_length": 100, "region": None},
        "predictions": {
            "indexing": "absolute",
            "values": [0.0] * 100,
            "value_kind": "probability",
            "value_range": [0.0, 1.0],
        },
        "annotations": {},
    }
    for path, value in (overrides or {}).items():
        keys = path.split(".")
        target = payload
        for key in keys[:-1]:
            target = target[key]
        target[keys[-1]] = value
    return payload


class WrapPayloadHappyPathTests(SimpleTestCase):
    """wrap_predictions_as_payload should always produce schema-valid output."""

    def test_default_is_covid_atgc(self):
        payload = wrap_predictions_as_payload(np.zeros((29903, 4)))
        validate_prediction_payload(payload)
        self.assertEqual(payload["task"]["kind"], "categorical_per_position")
        self.assertEqual(payload["task"]["labels"], ["A", "T", "G", "C"])
        self.assertEqual(payload["domain"]["total_length"], 29903)
        self.assertIsNone(payload["domain"]["region"])

    def test_binary_with_region(self):
        payload = wrap_predictions_as_payload(
            np.zeros(1260),
            kind="binary_per_position",
            total_length=29903,
            region={"start": 28273, "end": 29533},
        )
        validate_prediction_payload(payload)
        self.assertEqual(payload["task"]["kind"], "binary_per_position")
        self.assertNotIn("labels", payload["task"])
        self.assertEqual(payload["domain"]["region"], {"start": 28273, "end": 29533})

    def test_scalar_arbitrary_length(self):
        # An influenza-shaped genome length, no COVID assumptions anywhere.
        payload = wrap_predictions_as_payload(
            np.zeros(13500),
            kind="scalar_per_position",
            total_length=13500,
            value_kind="score",
            value_range=None,
        )
        validate_prediction_payload(payload)
        self.assertEqual(payload["task"]["kind"], "scalar_per_position")
        self.assertEqual(payload["domain"]["total_length"], 13500)

    def test_contract_version_is_2_0(self):
        payload = wrap_predictions_as_payload(np.zeros((1, 4)))
        self.assertEqual(payload["contract_version"], MODEL_CONTRACT_VERSION)
        self.assertEqual(payload["contract_version"], "2.0")

    def test_categorical_passes_through_custom_labels(self):
        payload = wrap_predictions_as_payload(
            np.zeros((10, 3)),
            kind="categorical_per_position",
            labels=["X", "Y", "Z"],
            total_length=10,
        )
        validate_prediction_payload(payload)
        self.assertEqual(payload["task"]["labels"], ["X", "Y", "Z"])


class ValidatorAcceptsTests(SimpleTestCase):
    """Round-trip: every shape declared as valid in the contract must pass."""

    def test_accepts_categorical_with_region(self):
        # A 4-channel categorical track restricted to a sub-region.
        payload = wrap_predictions_as_payload(
            np.zeros((500, 4)),
            kind="categorical_per_position",
            labels=["A", "T", "G", "C"],
            total_length=10_000,
            region={"start": 4000, "end": 4500},
        )
        validate_prediction_payload(payload)

    def test_accepts_relative_to_region_indexing(self):
        payload = wrap_predictions_as_payload(
            np.zeros(200),
            kind="binary_per_position",
            total_length=29903,
            region={"start": 100, "end": 300},
            indexing="relative_to_region",
        )
        validate_prediction_payload(payload)


class ValidatorRejectsTests(SimpleTestCase):
    """Each invariant the contract documents should be enforced."""

    def test_rejects_non_dict_root(self):
        with self.assertRaises(ValueError):
            validate_prediction_payload([1, 2, 3])

    def test_rejects_unknown_task_kind(self):
        with self.assertRaisesRegex(ValueError, "task.kind"):
            validate_prediction_payload(_base_payload({"task.kind": "regression_per_position"}))

    def test_rejects_categorical_without_labels(self):
        payload = _base_payload({
            "task.kind": "categorical_per_position",
            "predictions.values": [[0.25, 0.25, 0.25, 0.25]] * 100,
        })
        # task.labels intentionally absent
        with self.assertRaisesRegex(ValueError, "labels"):
            validate_prediction_payload(payload)

    def test_rejects_categorical_shape_mismatch(self):
        payload = _base_payload({
            "task.kind": "categorical_per_position",
            "predictions.values": [[0.5, 0.5]] * 100,  # 2 cols, 4 declared labels
        })
        payload["task"]["labels"] = ["A", "T", "G", "C"]
        with self.assertRaisesRegex(ValueError, "categorical values must have shape"):
            validate_prediction_payload(payload)

    def test_rejects_binary_shape_mismatch(self):
        # total_length=100 but only 50 values.
        payload = _base_payload({"predictions.values": [0.0] * 50})
        with self.assertRaisesRegex(ValueError, "binary_per_position values"):
            validate_prediction_payload(payload)

    def test_rejects_scalar_2d_values(self):
        payload = _base_payload({
            "task.kind": "scalar_per_position",
            "predictions.values": [[0.0, 0.0]] * 100,
        })
        with self.assertRaisesRegex(ValueError, "scalar_per_position values"):
            validate_prediction_payload(payload)

    def test_rejects_region_out_of_bounds(self):
        payload = _base_payload({
            "domain.region": {"start": 50, "end": 200},   # end > total_length
            "predictions.values": [0.0] * 150,
        })
        with self.assertRaisesRegex(ValueError, "region"):
            validate_prediction_payload(payload)

    def test_rejects_region_start_greater_than_end(self):
        payload = _base_payload({
            "domain.region": {"start": 80, "end": 50},
            "predictions.values": [0.0] * 0,
        })
        with self.assertRaisesRegex(ValueError, "region"):
            validate_prediction_payload(payload)

    def test_rejects_zero_total_length(self):
        payload = _base_payload({"domain.total_length": 0, "predictions.values": []})
        with self.assertRaisesRegex(ValueError, "total_length"):
            validate_prediction_payload(payload)

    def test_rejects_bad_value_kind(self):
        payload = _base_payload({"predictions.value_kind": "magic"})
        with self.assertRaisesRegex(ValueError, "value_kind"):
            validate_prediction_payload(payload)

    def test_rejects_bad_indexing(self):
        payload = _base_payload({"predictions.indexing": "yolo"})
        with self.assertRaisesRegex(ValueError, "indexing"):
            validate_prediction_payload(payload)

    def test_rejects_missing_values(self):
        payload = _base_payload()
        del payload["predictions"]["values"]
        with self.assertRaisesRegex(ValueError, "values"):
            validate_prediction_payload(payload)


class CustomOrganismEndToEndTests(SimpleTestCase):
    """Drive ``run_prediction`` against ``example_custom_organism_bundle`` to
    prove the custom-organism path works end-to-end:

    * genome is read from a bundle helper file (10800 bp synthetic FASTA)
    * protein_regions are parsed from a bundle helper CSV (5 regions)
    * the payload declares the right ``total_length`` and ``protein_regions``
      annotation, and the validator accepts it

    This is the smallest realistic case for the user's "yeni bir virüs
    yükle" vision, and serves as a working contract example for new bundles.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # Import lazily — these modules pull in TF/PyTorch and are heavy.
        from .plugin_runtime import resolve_uploaded_bundle
        from .viewsUpdated import (
            UPLOADED_MODELS_DIR,
            PredictionContext,
            _resolve_organism_data,
            run_prediction,
        )

        bundle = resolve_uploaded_bundle(
            UPLOADED_MODELS_DIR,
            "uploaded:example_custom_organism_bundle",
            {},
        )
        genome_seq, organism_protein_regions = _resolve_organism_data(bundle)
        ctx = PredictionContext(
            node_id=None,
            node_ids=[],
            elapsed_day=0,
            selected_protein_region=None,
            selected_protein_regions=None,
            depth=0,
            bundle=bundle,
            custom_parameters=bundle.custom_parameters,
            genome_sequence=genome_seq,
            organism_protein_regions=organism_protein_regions,
        )

        class _StubRequest:  # Stand-in for a Django request object.
            session = {}

        cls.bundle = bundle
        cls.response = run_prediction(ctx, _StubRequest())
        cls.payload = cls.response["predictionPayload"]

    def test_bundle_resolved_as_custom_organism(self):
        self.assertEqual(self.bundle.organism, "custom")
        self.assertEqual(self.bundle.genome_file, "synthetic_genome.fasta")
        self.assertEqual(self.bundle.protein_regions_file, "synthetic_proteins.csv")

    def test_payload_declares_organism_length(self):
        # The 10800 bp synthetic genome's length must reach the payload —
        # this is the COVID-free path of the platform.
        self.assertEqual(self.payload["domain"]["total_length"], 10800)
        self.assertIsNone(self.payload["domain"]["region"])

    def test_payload_is_binary_per_position(self):
        self.assertEqual(self.payload["task"]["kind"], "binary_per_position")
        values = np.asarray(self.payload["predictions"]["values"])
        self.assertEqual(values.shape, (10800,))
        # Adapter seeded the RNG, so we can pin a sanity range.
        self.assertTrue(0.0 <= values.min() and values.max() <= 1.0)

    def test_annotations_carry_all_five_protein_regions(self):
        regions = self.payload.get("annotations", {}).get("protein_regions", {})
        self.assertEqual(
            set(regions.keys()),
            {"Capsid", "prM", "Env", "NS1", "NS5"},
        )
        # Boundaries should match the CSV verbatim (1-based inclusive).
        self.assertEqual(regions["NS5"], [9000, 10800])

    def test_validator_accepts_payload(self):
        # Already validated inside run_prediction, but assert explicitly.
        validate_prediction_payload(self.payload)

    def test_no_legacy_covid_fields(self):
        # binary_per_position payloads don't carry the ATGC legacy block.
        self.assertNotIn("genomeData", self.response)
        self.assertNotIn("protein_mutation_probs", self.response)


class OrganismRegistryTests(SimpleTestCase):
    """The built-in organism registry plus the variant catalog drive every
    dispatch path the upload UI can take. Each test pins one invariant a
    future refactor must not break.
    """

    def test_four_builtin_organisms_registered(self):
        from .organism_registry import list_builtin_organisms

        self.assertEqual(
            set(list_builtin_organisms()),
            {"covid", "influenza_h1n1", "influenza_h3n2", "influenza_h5n1"},
        )

    def test_h3n2_reference_loads(self):
        from .organism_registry import read_builtin_organism

        seq, regions = read_builtin_organism("influenza_h3n2")
        # A/Hong Kong/1/1968 HA — 1736 bp; H3 HA CDS starts at position 1
        # (NCBI submissions for modern HA segments typically drop the 5' UTR).
        self.assertEqual(len(seq), 1736)
        self.assertEqual(regions, {"HA": [1, 1701]})

    def test_h3n2_variants_catalog(self):
        # Faz 1.6: references became regular catalog entries so the
        # consolidated 9-strain dropdown can offer them next to the variants.
        # HK/1/1968 is now first, marked as the H3N2 reference.
        from .organism_registry import list_variants

        catalog = list_variants("influenza_h3n2")
        names = [v["name"] for v in catalog]
        self.assertEqual(names, ["hk_1_1968", "hk_4801_2014", "darwin_9_2021"])
        # Darwin/9/2021 is the current WHO H3N2 vaccine reference — pin its
        # accession so a silent catalog edit can't change which strain we ship.
        darwin = next(v for v in catalog if v["name"] == "darwin_9_2021")
        self.assertEqual(darwin["accession"], "PX230218.1")
        self.assertEqual(darwin["subtype"], "H3N2")
        # Reference flag must mark the founder strain so the UI can badge it.
        hk_1_1968 = next(v for v in catalog if v["name"] == "hk_1_1968")
        self.assertEqual(hk_1_1968["is_reference"], "true")

    def test_h1n1_reference_loads(self):
        from .organism_registry import read_builtin_organism

        seq, regions = read_builtin_organism("influenza_h1n1")
        # A/PR/8/34 segment 4 — 1778 bp; HA CDS 33..1733.
        self.assertEqual(len(seq), 1778)
        self.assertEqual(regions, {"HA": [33, 1733]})

    def test_h5n1_reference_loads(self):
        from .organism_registry import read_builtin_organism

        seq, regions = read_builtin_organism("influenza_h5n1")
        # A/goose/Guangdong/1/96 HA — 1760 bp; H5 HA CDS 22..1728.
        self.assertEqual(len(seq), 1760)
        self.assertEqual(regions, {"HA": [22, 1728]})

    def test_h1n1_variants_catalog(self):
        # Faz 1.6: the H1N1 reference (A/PR/8/34) is now a catalog entry too,
        # so the consolidated 9-strain influenza dropdown can list it next to
        # Cal/07 and Mich/45.
        from .organism_registry import list_variants

        catalog = list_variants("influenza_h1n1")
        names = [v["name"] for v in catalog]
        self.assertEqual(names, ["pr_8_34", "cal_07_2009", "mich_45_2015"])
        # Spot-check that NCBI accession + CDS coordinates survive parsing.
        cal = next(v for v in catalog if v["name"] == "cal_07_2009")
        self.assertEqual(cal["accession"], "CY121680.1")
        self.assertEqual((cal["cds_start"], cal["cds_end"]), ("21", "1721"))
        # PR/8/34 carries the reference flag for the UI badge.
        pr8 = next(v for v in catalog if v["name"] == "pr_8_34")
        self.assertEqual(pr8["is_reference"], "true")

    def test_h5n1_variants_include_2024_bovine(self):
        from .organism_registry import list_variants

        catalog = list_variants("influenza_h5n1")
        names = [v["name"] for v in catalog]
        self.assertIn("cattle_texas_2024", names)
        bovine = next(v for v in catalog if v["name"] == "cattle_texas_2024")
        self.assertEqual(bovine["accession"], "PX718764.1")
        self.assertEqual(bovine["subtype"], "H5N1")

    def test_read_variant_returns_sequence(self):
        from .organism_registry import read_variant

        seq = read_variant("influenza_h5n1", "cattle_texas_2024")
        self.assertEqual(len(seq), 1704)  # PX718764.1 length

    def test_read_unknown_variant_raises(self):
        from .organism_registry import read_variant

        with self.assertRaisesRegex(FileNotFoundError, "not in"):
            read_variant("influenza_h5n1", "not_a_real_variant")


class VariantDispatchTests(SimpleTestCase):
    """``_resolve_organism_data`` must combine the built-in organism's protein
    regions with the selected variant's genome — that's the whole point of
    the variant catalog. Custom organisms must never accept a variant.
    """

    def _bundle(self, **fields):
        from .plugin_runtime import BundleResolution

        return BundleResolution(
            source="test",
            model_path="",
            extractor_path=None,
            adapter_path=None,
            bundle_dir=None,
            bundle_name=None,
            model_name="x",
            selected_model="x",
            custom_parameters={},
            saved_folder=None,
            **fields,
        )

    def test_reference_path_no_variant(self):
        from .viewsUpdated import _resolve_organism_data

        seq, regions = _resolve_organism_data(self._bundle(organism="influenza_h1n1"))
        self.assertEqual(len(seq), 1778)
        self.assertEqual(regions, {"HA": [33, 1733]})

    def test_variant_replaces_genome_keeps_protein_regions(self):
        # The H5N1 cattle 2024 variant has its own length (1704 bp), but the
        # frontend's protein overlay should still draw the reference HA region
        # — that's the biological contract for HA subtype variants.
        from .viewsUpdated import _resolve_organism_data

        seq, regions = _resolve_organism_data(
            self._bundle(organism="influenza_h5n1", variant="cattle_texas_2024")
        )
        self.assertEqual(len(seq), 1704)
        self.assertEqual(regions, {"HA": [22, 1728]})  # reference's, not variant's

    def test_variant_dispatch_for_each_cataloged_strain(self):
        from .organism_registry import list_variants
        from .viewsUpdated import _resolve_organism_data

        for organism in ("influenza_h1n1", "influenza_h3n2", "influenza_h5n1"):
            for entry in list_variants(organism):
                seq, _ = _resolve_organism_data(
                    self._bundle(organism=organism, variant=entry["name"])
                )
                self.assertGreater(len(seq), 1000, f"{organism}/{entry['name']}")

    def test_h3n2_darwin_variant_keeps_reference_protein_regions(self):
        # Variant FASTAs may have different lengths from the reference, but
        # the HA region overlay must stay the organism's authoritative one.
        from .viewsUpdated import _resolve_organism_data

        seq, regions = _resolve_organism_data(
            self._bundle(organism="influenza_h3n2", variant="darwin_9_2021")
        )
        self.assertEqual(len(seq), 1718)  # PX230218.1 length
        self.assertEqual(regions, {"HA": [1, 1701]})  # H3N2 reference's, not variant's

    def test_normalize_rejects_custom_with_variant(self):
        from .plugin_runtime import _normalize_organism_fields

        with self.assertRaisesRegex(ValueError, "Custom organisms cannot declare"):
            _normalize_organism_fields(
                {"organism": "custom", "genome_file": "g.fasta", "variant": "x"}
            )

    def test_normalize_accepts_generic_influenza(self):
        # Faz 1.6: "influenza" is now a first-class organism (the generic
        # subtype-agnostic bundle declaration). Variant required at predict
        # time, not normalization time.
        from .plugin_runtime import _normalize_organism_fields

        fields = _normalize_organism_fields({"organism": "influenza"})
        self.assertEqual(fields["organism"], "influenza")
        self.assertIsNone(fields["variant"])


class GenericInfluenzaDispatchTests(SimpleTestCase):
    """Faz 1.6: the prediction-screen variant picker sends ``organism="influenza"``
    plus a strain name (one of the 9 cataloged). The runtime must route that
    to the right subtype's protein_regions without the bundle declaring a
    subtype, and reject the no-variant case before it reaches the model.
    """

    def _bundle(self, **fields):
        from .plugin_runtime import BundleResolution

        return BundleResolution(
            source="test",
            model_path="",
            extractor_path=None,
            adapter_path=None,
            bundle_dir=None,
            bundle_name=None,
            model_name="x",
            selected_model="x",
            custom_parameters={},
            saved_folder=None,
            **fields,
        )

    def test_consolidated_catalog_has_9_strains(self):
        from .organism_registry import list_all_influenza_variants

        catalog = list_all_influenza_variants()
        self.assertEqual(len(catalog), 9)

        # 3 references must be present (one per subtype) and marked as such.
        refs = [v for v in catalog if v.get("is_reference") == "true"]
        self.assertEqual(
            {v["name"] for v in refs},
            {"pr_8_34", "hk_1_1968", "goose_guangdong_1_96"},
        )

    def test_consolidated_catalog_has_2024_bovine_h5n1(self):
        from .organism_registry import list_all_influenza_variants

        bovine = next(
            v
            for v in list_all_influenza_variants()
            if v["name"] == "cattle_texas_2024"
        )
        self.assertEqual(bovine["subtype"], "H5N1")
        self.assertEqual(bovine["organism"], "influenza_h5n1")
        self.assertEqual(bovine["accession"], "PX718764.1")

    def test_find_variant_subtype_routes_each_strain(self):
        # Every cataloged strain must be discoverable by name; this is what
        # makes the generic ``organism="influenza"`` dispatch possible.
        from .organism_registry import find_variant_subtype, list_all_influenza_variants

        for entry in list_all_influenza_variants():
            self.assertEqual(
                find_variant_subtype(entry["name"]),
                entry["organism"],
                f"variant {entry['name']} should route to its catalog organism",
            )

    def test_find_variant_subtype_returns_none_for_unknown(self):
        from .organism_registry import find_variant_subtype

        self.assertIsNone(find_variant_subtype("not_a_real_strain"))

    def test_generic_dispatch_routes_to_right_subtype(self):
        # organism="influenza" + cattle_texas_2024 → H5N1 protein_regions
        # without the bundle declaring "influenza_h5n1" explicitly.
        from .viewsUpdated import _resolve_organism_data

        seq, regions = _resolve_organism_data(
            self._bundle(organism="influenza", variant="cattle_texas_2024")
        )
        self.assertEqual(len(seq), 1704)
        self.assertEqual(regions, {"HA": [22, 1728]})  # H5N1 reference's regions

    def test_generic_dispatch_h3n2_uses_h3n2_regions(self):
        from .viewsUpdated import _resolve_organism_data

        seq, regions = _resolve_organism_data(
            self._bundle(organism="influenza", variant="darwin_9_2021")
        )
        self.assertEqual(len(seq), 1718)
        self.assertEqual(regions, {"HA": [1, 1701]})  # H3N2 reference's regions

    def test_generic_dispatch_reference_strain_works_too(self):
        # References (PR/8/34, HK/1/68, goose/Guangdong/1/96) are first-class
        # catalog entries now, so the prediction picker can offer them.
        from .viewsUpdated import _resolve_organism_data

        seq, regions = _resolve_organism_data(
            self._bundle(organism="influenza", variant="pr_8_34")
        )
        self.assertEqual(len(seq), 1778)  # A/PR/8/34's full segment
        self.assertEqual(regions, {"HA": [33, 1733]})  # H1N1 reference's regions

    def test_generic_dispatch_without_variant_raises(self):
        # The prediction screen must pick a strain; the backend won't guess.
        from .viewsUpdated import _resolve_organism_data

        with self.assertRaisesRegex(ValueError, "require a variant"):
            _resolve_organism_data(self._bundle(organism="influenza"))

    def test_generic_dispatch_unknown_variant_raises(self):
        from .viewsUpdated import _resolve_organism_data

        with self.assertRaisesRegex(ValueError, "Unknown influenza variant"):
            _resolve_organism_data(
                self._bundle(organism="influenza", variant="bogus_strain")
            )


class GetModelsOrganismFieldTests(SimpleTestCase):
    """Faz 1.6: the prediction-screen model dropdown renders organism badges
    by reading the ``models`` field on /api/models/ — make sure that field
    is populated and stays in sync with bundle metadata.
    """

    def test_response_includes_models_field(self):
        import json

        from django.test import RequestFactory

        from .viewsUpdated import get_models

        resp = get_models(RequestFactory().get("/api/models/"))
        data = json.loads(resp.content)
        self.assertIn("models", data)
        self.assertIn("available_models", data)  # legacy back-compat field

    def test_every_model_carries_organism_tag(self):
        import json

        from django.test import RequestFactory

        from .viewsUpdated import get_models

        resp = get_models(RequestFactory().get("/api/models/"))
        data = json.loads(resp.content)
        for entry in data["models"]:
            self.assertIn(entry["organism"], ("covid", "influenza", "custom"))
            self.assertIn("uploaded", entry)

    def test_custom_bundle_surfaces_custom_organism(self):
        # example_custom_organism_bundle is the smoke-test bundle that
        # declares organism="custom" in its bundle_metadata.json.
        import json

        from django.test import RequestFactory

        from .viewsUpdated import get_models

        resp = get_models(RequestFactory().get("/api/models/"))
        data = json.loads(resp.content)
        bundle = next(
            (m for m in data["models"] if m["name"] == "example_custom_organism_bundle"),
            None,
        )
        self.assertIsNotNone(bundle, "expected the example custom-organism bundle in catalog")
        self.assertEqual(bundle["organism"], "custom")
        self.assertTrue(bundle["uploaded"])


if __name__ == "__main__":  # pragma: no cover - convenience for stand-alone runs
    unittest.main()
