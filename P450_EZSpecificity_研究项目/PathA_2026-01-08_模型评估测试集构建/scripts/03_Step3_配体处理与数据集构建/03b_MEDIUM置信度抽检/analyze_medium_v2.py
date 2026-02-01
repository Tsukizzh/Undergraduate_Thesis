import json
import random
from collections import Counter

# Already reviewed in 03a (needs_human_review=true)
already_reviewed = {'REC_0195', 'REC_0205', 'REC_0206', 'REC_0207', 'REC_0208',
                    'REC_0209', 'REC_0210', 'REC_0213', 'REC_0214', 'REC_0216',
                    'REC_0220', 'REC_0354', 'REC_0500', 'REC_0515', 'REC_0590'}

# Read MEDIUM confidence records
with open('131条MEDIUM置信度记录.jsonl', 'r', encoding='utf-8') as f:
    all_records = [json.loads(line) for line in f]

# Filter out already reviewed records
records = [r for r in all_records if r['record_id'] not in already_reviewed]

print(f"Total MEDIUM confidence records: {len(all_records)}")
print(f"After excluding 03a reviewed: {len(records)}\n")

# Analyze by enzyme
enzymes = Counter()
for r in records:
    enzyme = r['enzyme_name'].split('(')[0].strip()
    enzymes[enzyme] += 1

print("=== Distribution by Enzyme (after exclusion) ===")
for enzyme, count in enzymes.most_common(10):
    print(f"{enzyme}: {count}")

# Analyze by classification
classes = Counter()
for r in records:
    classes[r['verified_class']] += 1

print("\n=== Distribution by Classification ===")
for cls, count in classes.most_common():
    print(f"{cls}: {count}")

# Analyze by trap detection
traps = Counter()
for r in records:
    traps[r['trap_detected']] += 1

print("\n=== Distribution by Trap Detection ===")
for trap, count in traps.most_common():
    print(f"{trap}: {count}")

# Stratified sampling strategy
print("\n=== Sampling Strategy ===")

# 1. All records with traps (priority)
with_traps = [r for r in records if r['trap_detected'] != 'none']
print(f"Records with traps: {len(with_traps)}")

# 2. Random sample from remaining records
without_traps = [r for r in records if r['trap_detected'] == 'none']
target_sample_size = 75  # Total target
remaining_needed = target_sample_size - len(with_traps)

random.seed(42)  # For reproducibility
sampled_without_traps = random.sample(without_traps, min(remaining_needed, len(without_traps)))

print(f"Random sample from non-trap records: {len(sampled_without_traps)}")
print(f"Total sampled: {len(with_traps) + len(sampled_without_traps)}")

# Combine and save
sampled_records = with_traps + sampled_without_traps
with open('75条抽检样本.jsonl', 'w', encoding='utf-8') as f:
    for r in sampled_records:
        f.write(json.dumps(r, ensure_ascii=False) + '\n')

print(f"\nSampled records saved to: 75条抽检样本.jsonl")

# Print sampled record IDs for verification
print("\n=== Sampled Record IDs ===")
for i, r in enumerate(sampled_records):
    if i < 20:
        print(f"{r['record_id']}: {r['enzyme_name'][:30]} | {r['ligand_name'][:30]} | {r['verified_class']}")
print(f"... and {len(sampled_records) - 20} more")
