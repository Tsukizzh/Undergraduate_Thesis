import sys, time, torch
sys.path.insert(0, '/root/rivermind-data/EZSpecificity/src')
sys.path.insert(0, '/root/rivermind-data/EZSpecificity/scripts/10_Step10_pt训练管线')
from pt_dataset import PtCacheDataset

ds = PtCacheDataset(
    cache_dir='/root/rivermind-data/EZSpecificity/data/10_Step10_pt训练/ezspec_pt_v1',
    split='train', edge_mode='legacy_bug'
)

# Warmup
for i in range(5):
    _ = ds[i]

# Profile 200 __getitem__ calls
times = []
for i in range(200):
    t0 = time.perf_counter()
    _ = ds[i]
    times.append((time.perf_counter() - t0) * 1000)

print(f"__getitem__ over 200 samples:")
print(f"  mean:  {sum(times)/len(times):.2f} ms")
print(f"  min:   {min(times):.2f} ms")
print(f"  max:   {max(times):.2f} ms")
print(f"  p50:   {sorted(times)[100]:.2f} ms")
print(f"  p95:   {sorted(times)[190]:.2f} ms")
print(f"  single-thread: {200000/sum(times):.0f} samples/s")
print(f"  7 workers theoretical: {7 * 200000/sum(times):.0f} samples/s")
print(f"  needed for 2.86 it/s * bs48 * 2gpu = {2.86*48*2:.0f} samples/s")

# Breakdown per function
import time as _time

def profile_one(idx):
    result = {}

    t0 = _time.perf_counter()
    sample_id = int(ds._index["sample_ids"][idx])
    sub_dir = ds.cache_dir / ds.split / "samples" / f"{sample_id // 1000:03d}"
    sample_path = sub_dir / f"sample_{sample_id:06d}.pt"
    s = torch.load(sample_path, map_location="cpu", weights_only=False)
    result["torch_load"] = (_time.perf_counter() - t0) * 1000

    t0 = _time.perf_counter()
    from pt_dataset import rebuild_protein_x, rebuild_ligand_x, rebuild_edge_features, GaussianSmearing
    n_lig = s["ligand_pos"].shape[0]
    n_prot = s["protein_pos"].shape[0]
    _ = rebuild_protein_x(s["protein_element"].long(), s["protein_aa_type"].long(), s["protein_is_backbone"].long(), n_lig)
    _ = rebuild_ligand_x(s["ligand_element"].long(), s["ligand_atom_aux"].long(), n_prot)
    result["rebuild_features"] = (_time.perf_counter() - t0) * 1000

    t0 = _time.perf_counter()
    smear = GaussianSmearing(stop=10.0, num_gaussians=32)
    _ = rebuild_edge_features(
        ligand_pos=s["ligand_pos"].float(), protein_pos=s["protein_pos"].float(),
        bond_index=s["bond_index"].long(), bond_type=s["bond_type"].long(),
        knn_edge_index=s["knn_edge_index"].long(), n_lig=n_lig,
        edge_mode="legacy_bug", dist_noise=False, smear=smear
    )
    result["rebuild_edges"] = (_time.perf_counter() - t0) * 1000

    t0 = _time.perf_counter()
    ds._ensure_file_handles()
    _ = ds._load_enzyme_embedding(int(s["enzyme_id"]))
    result["load_enzyme"] = (_time.perf_counter() - t0) * 1000

    t0 = _time.perf_counter()
    _ = ds._load_substrate_features(int(s["substrate_id"]))
    result["load_substrate"] = (_time.perf_counter() - t0) * 1000

    return result

# Run breakdown on 20 samples
breakdowns = [profile_one(i) for i in range(20)]
print("\nBreakdown (mean of 20 samples):")
for key in breakdowns[0]:
    vals = [b[key] for b in breakdowns]
    print(f"  {key:20s}: {sum(vals)/len(vals):6.2f} ms ({sum(vals)/sum(sum(b.values()) for b in breakdowns)*100:5.1f}%)")
