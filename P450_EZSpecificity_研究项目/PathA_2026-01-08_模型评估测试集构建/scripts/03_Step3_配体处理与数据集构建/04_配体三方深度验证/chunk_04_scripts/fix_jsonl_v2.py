import json
import os
import re

base_dir = r"c:\Users\Administrator\Desktop\EZSpecificity_Project\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\data\03_Step3_配体处理与数据集构建\04_配体三方深度验证\02_输出数据_results\chunk_04_results"
input_file = os.path.join(base_dir, 'verified_results.jsonl')
output_file = os.path.join(base_dir, 'verified_results_fixed.jsonl')

# Read as binary, clean up invalid bytes
with open(input_file, 'rb') as f:
    raw_data = f.read()

# Replace invalid UTF-8 sequences with placeholders
# First, decode with errors='replace' to find problem areas
text = raw_data.decode('utf-8', errors='replace')

# The replacement character is U+FFFD which shows as "�"
# We need to find patterns like "三方共�?PRODUCT�?" and fix them
# Common patterns: ending with "�?" which is a truncated string

# Fix truncated endings that break JSON
# Pattern: a string value that ends with truncated bytes followed by }
# e.g., "reasoning": "some text�?"} - the truncated bytes break the closing quote

print(f"Total length: {len(text)} chars")
print(f"Contains replacement chars: {chr(0xFFFD) in text}")

# Split into lines
lines = text.split('\n')
print(f"Number of lines: {len(lines)}")

fixed_records = []
skipped_lines = []

for line_num, line in enumerate(lines, 1):
    line = line.strip()
    if not line:
        continue

    # Try to fix common issues:
    # 1. Replace lone replacement characters with empty string
    # 2. Ensure strings are properly closed

    # Replace sequences of replacement characters followed by ? with just ...
    fixed_line = re.sub(r'\ufffd+\??', '...', line)

    try:
        record = json.loads(fixed_line)
        fixed_records.append((line_num, record))
    except json.JSONDecodeError as e:
        # Try more aggressive fixing
        # Find all strings and ensure they're properly closed
        try:
            # Last resort: find the JSON structure and fix it
            # Look for the pattern where a string ends with truncated bytes
            fixed_line2 = re.sub(r'([^\\])"[^"]*\ufffd[^"]*"', r'\1"[text truncated]"', fixed_line)
            record = json.loads(fixed_line2)
            fixed_records.append((line_num, record))
        except:
            skipped_lines.append((line_num, str(e)[:80]))
            print(f"  Line {line_num} still failed after fixes")

print(f"\nSuccessfully parsed: {len(fixed_records)} records")
print(f"Failed: {len(skipped_lines)} lines")

if skipped_lines:
    print("\nSkipped lines:")
    for ln, err in skipped_lines[:5]:
        print(f"  Line {ln}: {err}")

# Process records - remove duplicates and the incorrect REC_0180
seen_ids = set()
final_records = []

for line_num, record in fixed_records:
    record_id = record.get('record_id')

    # Skip the incorrect REC_0180 (classification_changed=true and final_class=SUBSTRATE)
    if record_id == 'REC_0180':
        if record.get('classification_changed') == True and record.get('final_class') == 'SUBSTRATE':
            print(f"Skipping incorrect REC_0180 at line {line_num}")
            continue

    # Check for duplicates
    if record_id in seen_ids:
        print(f"Skipping duplicate {record_id} at line {line_num}")
        continue

    seen_ids.add(record_id)
    final_records.append(record)

# Sort by record_id
final_records.sort(key=lambda x: x['record_id'])

# Write output
with open(output_file, 'w', encoding='utf-8') as f:
    for record in final_records:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')

print(f'\n=== Final Summary ===')
print(f'Total unique records: {len(final_records)}')
print(f'Expected: 57 records')
print(f'Output: {output_file}')

# Statistics
if final_records:
    classification_counts = {}
    changed_count = 0
    for r in final_records:
        fc = r.get('final_class', 'UNKNOWN')
        classification_counts[fc] = classification_counts.get(fc, 0) + 1
        if r.get('classification_changed'):
            changed_count += 1

    print(f'\n=== Classification Statistics ===')
    for cls, count in sorted(classification_counts.items()):
        print(f'{cls}: {count}')
    print(f'classification_changed=true: {changed_count}')

    # List records with classification_changed
    print(f'\n=== Records with classification_changed=true ===')
    for r in final_records:
        if r.get('classification_changed'):
            print(f"  {r['record_id']}: {r.get('verified_class')} -> {r.get('final_class')}")
