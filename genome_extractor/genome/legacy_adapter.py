# /genome_extractor/genome/legacy_adapter.py

import os
import json
import tensorflow as tf
import numpy as np
from typing import Dict, Any, Tuple, List

from typing import Protocol
class CovMutExModel(Protocol):
    def metadata(self) -> Dict[str, Any]: ...
    def input_schema(self) -> Dict[str, Any]: ...
    def preprocess(self, inputs: Dict[str, Any]) -> Any: ...
    def predict(self, batch: Any) -> Any: ...
    def postprocess(self, raw: Any) -> Dict[str, Any]: ...

from .feature_extractor import parse_mutations, cache_precomputed_features, construct_variant_genome
from .configs import configs


def read_genome_sequence(file_path: str) -> str:
    with open(file_path, 'r') as f:
        next(f)
        return ''.join(line.strip() for line in f)

def _get_biologically_distributed_probs(p_no_mutation: float, p_mutation: float, ref_nuc: str, ti_tv_ratio: float = 2.0) -> Dict[str, float]:
    probabilities = {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}
    probabilities[ref_nuc] = p_no_mutation
    
    purines, pyrimidines = {'A', 'G'}, {'T', 'C'}
    
    transition_target = ""
    if ref_nuc in purines:
        transition_target = list(purines - {ref_nuc})[0]
    elif ref_nuc in pyrimidines:
        transition_target = list(pyrimidines - {ref_nuc})[0]
        
    if p_mutation > 0 and transition_target:
        prob_tv = p_mutation / (ti_tv_ratio + 2)
        prob_ti = prob_tv * ti_tv_ratio
    else:
        prob_tv, prob_ti = 0.0, 0.0
    
    for nuc in probabilities:
        if nuc == ref_nuc: continue
        elif nuc == transition_target: probabilities[nuc] = prob_ti
        else: probabilities[nuc] = prob_tv

    total_sum = sum(probabilities.values())
    if total_sum > 1e-9:
        for nuc in probabilities: probabilities[nuc] /= total_sum
    else:
        for nuc in probabilities: probabilities[nuc] = 1.0 if nuc == ref_nuc else 0.0
            
    return probabilities


class LegacyModelAdapter(CovMutExModel):
    
    def __init__(self):
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        model_path = os.path.join(base_dir, 'covid19_models', 'models', 'balanced_data_model.keras')
        
        print(f"LegacyAdapter: Model yükleniyor: {model_path}")
        self.model = tf.keras.models.load_model(model_path)
        print("LegacyAdapter: Model başarıyla yüklendi.")

        genome_file_path = os.path.join(os.path.dirname(__file__), "genome.txt")
        self.reference_genome = read_genome_sequence(genome_file_path)

        self.protein_regions = configs().get('protein regions', {})

    def metadata(self) -> dict:
        return {
            "name": "Balanced Data Model (Legacy)",
            "version": "1.0.0",
            "author": "Original Research Team",
            "description": "Mevcut feature_extractor mantığını kullanan temel model."
        }

    def input_schema(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "nodeId": {"type": "string", "description": "Variant ID (örn: EGY/CCHE57357...)"},
                "elapsedDay": {"type": "integer", "default": 0},
                "selectedProteinRegion": {
                    "type": ["string", "null"],
                    "description": "İsteğe bağlı protein bölgesi (örn: 'S', 'N')"
                },
            },
            "required": ["nodeId"],
        }

    def preprocess(self, inputs: dict) -> tuple:
        print("LegacyAdapter: Ön işleme (preprocess) başlıyor...")
        
        nodeId = inputs.get('nodeId')
        elapsedDay = inputs.get('elapsedDay', 0)
        selectedProteinRegion = inputs.get('selectedProteinRegion')
        
        mutations = parse_mutations(nodeId)
        variant_genome_sequence = construct_variant_genome(self.reference_genome, mutations)

        codon_mapping_path = os.path.join(os.path.dirname(__file__), 'codon_aa_mapping.json')
        cache_path = os.path.join(os.path.dirname(__file__), 'node_features.h5')
        
        protein_region_to_process = {}
        if selectedProteinRegion and selectedProteinRegion in self.protein_regions:
            protein_region_to_process = {selectedProteinRegion: self.protein_regions[selectedProteinRegion]}

        features = cache_precomputed_features(
            cache_path=cache_path, 
            genome_seq=variant_genome_sequence, 
            mutations=mutations,
            codon_mapper=json.load(open(codon_mapping_path)), 
            config_file=configs(), 
            node_ids=[nodeId],
            elapsed_day=elapsedDay, 
            protein_regions=protein_region_to_process or None
        )
        
        features = features.squeeze()
        print(f"LegacyAdapter: Ön işleme bitti. Özellik (features) şekli: {features.shape}")
        
        return (features, variant_genome_sequence, protein_region_to_process)

    def predict(self, batch: tuple) -> tuple:
        print("LegacyAdapter: Tahmin (predict) yapılıyor...")
        
        features, variant_genome_sequence, protein_region = batch
        
        raw_predictions = self.model.predict(features, verbose=0)
            
        print(f"LegacyAdapter: Tahmin bitti. Ham tahmin (raw_predictions) şekli: {raw_predictions.shape}")
        
        return (raw_predictions, variant_genome_sequence, protein_region)

    def postprocess(self, prediction_bundle: tuple) -> dict:
        print("LegacyAdapter: Son işleme (postprocess) başlıyor...")
        raw_predictions, variant_genome_sequence, protein_region = prediction_bundle

        if protein_region:
            region_name = list(protein_region.keys())[0]
            start, end = self.protein_regions[region_name]
            end = min(end + 1, len(variant_genome_sequence))
            start = max(0, start)
            
            processed_genome_sequence = variant_genome_sequence[start:end]
            
            if len(raw_predictions) != len(processed_genome_sequence):
                print(f"UYARI: Tahmin uzunluğu ({len(raw_predictions)}) ile bölge uzunluğu ({len(processed_genome_sequence)}) eşleşmiyor!")
                min_len = min(len(raw_predictions), len(processed_genome_sequence))
                raw_predictions = raw_predictions[:min_len]
                processed_genome_sequence = processed_genome_sequence[:min_len]
        else:
            processed_genome_sequence = variant_genome_sequence
            
            min_len = min(len(raw_predictions), len(processed_genome_sequence))
            raw_predictions = raw_predictions[:min_len]
            processed_genome_sequence = processed_genome_sequence[:min_len]


        genome_data_lists = {'A': [], 'T': [], 'G': [], 'C': []}
        
        for idx in range(len(raw_predictions)):
            p_mutation = raw_predictions[idx][0] 
            p_no_mutation = 1.0 - p_mutation
            
            current_nucleotide = processed_genome_sequence[idx]
            
            position_probs_dict = _get_biologically_distributed_probs(
                p_no_mutation, p_mutation, current_nucleotide
            )
            
            genome_data_lists['A'].append(position_probs_dict['A'])
            genome_data_lists['T'].append(position_probs_dict['T'])
            genome_data_lists['G'].append(position_probs_dict['G'])
            genome_data_lists['C'].append(position_probs_dict['C'])

        
        protein_mutation_probs = {}
        
        if not protein_region:
            mutation_only_predictions = raw_predictions.flatten()
            
            for protein, (start, end) in self.protein_regions.items():
                start_idx = max(0, start)
                end_idx = min(len(mutation_only_predictions), end + 1)
                
                if start_idx >= end_idx:
                    protein_mutation_probs[protein] = 0.0
                    continue

                relevant_predictions = mutation_only_predictions[start_idx:end_idx]
                avg_prob = float(np.mean(relevant_predictions))
                protein_mutation_probs[protein] = avg_prob
        else:
            region_name = list(protein_region.keys())[0]
            avg_prob = float(np.mean(raw_predictions.flatten()))
            protein_mutation_probs[region_name] = avg_prob


        print("LegacyAdapter: Son işleme bitti.")
        
        return {
            "genome_data": genome_data_lists,
            "protein_summary": protein_mutation_probs,
            "metadata": {
                "processed_length": len(processed_genome_sequence),
                "is_partial": bool(protein_region)
            }
        }