"""
阶段A recon：查看样本 .pt 里的字段，重点看有没有残基粒度的信息
"""
import torch
from pathlib import Path

sample_path = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_allfix_unified/random/test/samples/000/sample_000000.pt")
s = torch.load(sample_path, map_location="cpu", weights_only=False)

print("=" * 60)
print(f"sample: {sample_path.name}")
print("=" * 60)
print()
print("=== 全部字段 ===")
for k, v in s.items():
    if hasattr(v, "shape"):
        print(f"  {k:30s} shape={tuple(v.shape)!s:20s} dtype={v.dtype}")
    else:
        print(f"  {k:30s} {type(v).__name__} = {str(v)[:80]}")

print()
print("=== 蛋白原子信息 ===")
print(f"protein_pos 原子数: {s['protein_pos'].shape[0]}")
print(f"protein_aa_type (前20): {s['protein_aa_type'][:20].tolist()}")
print(f"protein_is_backbone (前20): {s['protein_is_backbone'][:20].tolist()}")
print(f"protein_element (前20): {s['protein_element'][:20].tolist()}")

# 尝试推断残基数：通过 aa_type 变化点
aa = s['protein_aa_type'].tolist()
# 不能简单按连续相同来数，因为相邻残基氨基酸可能相同
# 看看有没有别的线索
print()
print("=== 检查残基粒度信息 ===")
found = []
for k in s.keys():
    kl = k.lower()
    if "res" in kl or "idx" in kl or "pocket" in kl:
        found.append(k)
if found:
    print(f"  找到可能字段: {found}")
else:
    print("  没有找到 residue/idx/pocket 相关字段")
    print("  => 残基编号信息在 .pt 里没保留，必须从原始 PDB 回头读")

print()
print("=== ligand_index_raw 是什么 ===")
if "ligand_index_raw" in s:
    print(f"  shape: {s['ligand_index_raw'].shape}")
    print(f"  前20: {s['ligand_index_raw'][:20].tolist()}")
    print(f"  max: {s['ligand_index_raw'].max().item()}")
    print(f"  min: {s['ligand_index_raw'].min().item()}")
