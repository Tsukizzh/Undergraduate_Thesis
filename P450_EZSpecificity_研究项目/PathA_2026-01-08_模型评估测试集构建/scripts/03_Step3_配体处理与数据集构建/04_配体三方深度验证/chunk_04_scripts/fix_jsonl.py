import json
import os
import chardet

base_dir = r"c:\Users\Administrator\Desktop\EZSpecificity_Project\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\data\03_Step3_配体处理与数据集构建\04_配体三方深度验证\02_输出数据_results\chunk_04_results"
input_file = os.path.join(base_dir, 'verified_results.jsonl')
output_file = os.path.join(base_dir, 'verified_results_fixed.jsonl')

# First, detect the encoding
with open(input_file, 'rb') as f:
    raw_data = f.read()
    detected = chardet.detect(raw_data)
    print(f'Detected encoding: {detected}')

# Try multiple encodings
encodings = ['utf-8', 'gbk', 'gb18030', 'utf-16', 'latin-1']
successful_encoding = None

for encoding in encodings:
    try:
        print(f'\nTrying encoding: {encoding}')
        with open(input_file, 'r', encoding=encoding) as f:
            lines = f.readlines()
        print(f'  Successfully read {len(lines)} lines')

        # Try to parse first line as JSON
        if lines and lines[0].strip():
            json.loads(lines[0])
            print(f'  Successfully parsed JSON')
            successful_encoding = encoding
            break
    except Exception as e:
        print(f'  Failed: {e}')

if not successful_encoding:
    print('\nNo suitable encoding found. Using utf-8 with ignore errors.')
    successful_encoding = 'utf-8'

print(f'\nUsing encoding: {successful_encoding}')

seen_ids = set()
records = []
duplicates = []
parse_errors = []

with open(input_file, 'r', encoding=successful_encoding, errors='ignore') as f:
    for line_num, line in enumerate(f, 1):
        line = line.strip()
        if not line:
            continue

        try:
            record = json.loads(line)
            record_id = record.get('record_id')

            # Skip the incorrect REC_0180 (line 3 with classification_changed=true)
            if record_id == 'REC_0180' and record.get('classification_changed') == True:
                print(f'Skipping incorrect REC_0180 at line {line_num} (classification_changed=true)')
                continue

            # Check for duplicates
            if record_id in seen_ids:
                print(f'Warning: Duplicate {record_id} at line {line_num}')
                duplicates.append(record_id)
                continue

            seen_ids.add(record_id)
            records.append(record)

        except json.JSONDecodeError as e:
            parse_errors.append(f'Line {line_num}: {str(e)[:100]}')

if parse_errors:
    print(f'\n=== Parse Errors ({len(parse_errors)}) ===')
    for err in parse_errors[:5]:  # Show first 5
        print(err)
    if len(parse_errors) > 5:
        print(f'... and {len(parse_errors) - 5} more')

# Sort by record_id to ensure consistent ordering
records.sort(key=lambda x: x['record_id'])

# Write fixed file
with open(output_file, 'w', encoding='utf-8') as f:
    for record in records:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')

print(f'\n=== Summary ===')
print(f'Total unique records: {len(records)}')
print(f'Expected: 57 records')
print(f'Duplicates removed: {duplicates}')
print(f'Output written to: {output_file}')

# Also show classification statistics
if records:
    classification_counts = {}
    changed_count = 0
    for r in records:
        fc = r.get('final_class', 'UNKNOWN')
        classification_counts[fc] = classification_counts.get(fc, 0) + 1
        if r.get('classification_changed'):
            changed_count += 1

    print(f'\n=== Classification Statistics ===')
    for cls, count in sorted(classification_counts.items()):
        print(f'{cls}: {count}')
    print(f'classification_changed=true: {changed_count}')

