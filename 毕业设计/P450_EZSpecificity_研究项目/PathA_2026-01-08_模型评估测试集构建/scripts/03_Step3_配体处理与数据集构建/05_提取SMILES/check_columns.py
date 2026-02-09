import csv
import sys

csv_path = "source_data/01_核心数据/修复后最终版/独立测试集_682条.csv"
with open(csv_path, encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    print("Columns:", reader.fieldnames)
