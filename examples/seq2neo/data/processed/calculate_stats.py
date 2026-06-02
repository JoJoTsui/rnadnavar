import json
import os
from collections import defaultdict

def calculate_stats(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)

    dna_normal_size = 0
    dna_tumor_size = 0
    rna_tumor_size = 0
    disease_counts = defaultdict(int)

    # To ensure each unique fastq file is counted exactly once for size
    seen_files = set()
    # To track file usage across samples
    file_to_samples = defaultdict(list)
    
    # Track samples from the sets to count as per references
    referenced_samples = set()
    for set_file in ['set1_samples.tsv', 'set2_samples.tsv', 'set3_samples.tsv', 'set4_samples.tsv']:
        set_path = os.path.join(os.path.dirname(json_path), set_file)
        if os.path.exists(set_path):
            with open(set_path, 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        proj, pat, dis, status = parts[0], parts[1], parts[2], parts[3]
                        # Only count standard/extra samples mentioned in the sets
                        referenced_samples.add((proj, pat))
                        disease_counts[dis] += 1

    for project in data.get('projects', []):
        project_id = project.get('project_id')
        for sample in project.get('samples', []):
            patient_id = sample.get('patient_id')
            sample_id = f"{project_id}_{patient_id}"
            
            modalities = sample.get('modalities', {})
            
            def process_modality(mod_name, size_accumulator):
                nonlocal dna_normal_size, dna_tumor_size, rna_tumor_size
                mod_data = modalities.get(mod_name, {})
                mod_size = 0
                for pair in mod_data.get('pairs', []):
                    for key in ['r1', 'r2']:
                        path = pair.get(key)
                        if path:
                            file_to_samples[path].append(f"{sample_id}({mod_name})")
                            if path not in seen_files:
                                if os.path.exists(path):
                                    size = os.path.getsize(path)
                                    mod_size += size
                                    seen_files.add(path)
                return mod_size

            dna_normal_size += process_modality('DN', 'dna_normal_size')
            dna_tumor_size += process_modality('DT', 'dna_tumor_size')
            rna_tumor_size += process_modality('RT', 'rna_tumor_size')

    def to_gb(bytes):
        return bytes / (1024**3)

    print(f"Total DNA Normal file size: {to_gb(dna_normal_size):.2f} GB")
    print(f"Total DNA Tumor file size: {to_gb(dna_tumor_size):.2f} GB")
    print(f"Combined DNA size: {to_gb(dna_normal_size + dna_tumor_size):.2f} GB")
    print(f"Total RNA Tumor file size: {to_gb(rna_tumor_size):.2f} GB")
    
    print("\nSample Disease Count (from Set Tables):")
    sorted_diseases = sorted(disease_counts.items(), key=lambda x: x[1], reverse=True)
    for disease, count in sorted_diseases:
        print(f"- {disease}: {count}")

    print("\nFile Reuse Check:")
    reused_files = {path: samples for path, samples in file_to_samples.items() if len(set(samples)) > 1}
    if not reused_files:
        print("No files reused across multiple modality/sample contexts.")
    else:
        print(f"Found {len(reused_files)} files reused across multiple samples:")
        for path, samples in list(reused_files.items())[:5]: # Show first 5
            print(f"- {path} used in: {', '.join(set(samples))}")
        if len(reused_files) > 5:
            print(f"... and {len(reused_files)-5} more.")

if __name__ == "__main__":
    json_file = "/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/examples/seq2neo/data/processed/merged.json"
    calculate_stats(json_file)
