"""Inspect index.pt in per-sample layout and shard layout to find a common ID scheme."""
import torch, glob, os

def show(name, obj, depth=0):
    prefix = "  " * depth
    if isinstance(obj, dict):
        for k, v in list(obj.items())[:8]:
            if hasattr(v, 'shape'):
                print(f'{prefix}{k}: tensor {tuple(v.shape)} {v.dtype}')
                if v.numel() <= 10:
                    print(f'{prefix}  values: {v.tolist()}')
                else:
                    print(f'{prefix}  first 5: {v[:5].tolist()}')
            elif isinstance(v, (list, tuple)):
                print(f'{prefix}{k}: list/tuple len={len(v)}')
                if len(v) <= 10: print(f'{prefix}  {v}')
                else: print(f'{prefix}  first 5: {v[:5]}')
            elif isinstance(v, dict):
                print(f'{prefix}{k}: dict keys={list(v.keys())[:8]}')
            else:
                print(f'{prefix}{k}: {type(v).__name__} = {v if len(str(v))<100 else str(v)[:100]+"..."}')

BASE = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'

print('=== pt_cache_fixed/random/train/index.pt (per-sample) ===')
idx = torch.load(f'{BASE}/pt_cache_fixed/random/train/index.pt', map_location='cpu', weights_only=False)
print(f'top-level type: {type(idx).__name__}')
if isinstance(idx, dict):
    print(f'keys: {list(idx.keys())}')
    show('idx', idx)

print()
print('=== pt_cache_heme_fixed/random/train/graph_0000.pt (shard) ===')
sh = torch.load(f'{BASE}/pt_cache_heme_fixed/random/train/graph_0000.pt', map_location='cpu', weights_only=False)
print(f'type: {type(sh).__name__}')
show('sh', sh)

print()
print('=== pt_cache_geom_fixed/train/graph_0000.pt (shard) ===')
sh2 = torch.load(f'{BASE}/pt_cache_geom_fixed/train/graph_0000.pt', map_location='cpu', weights_only=False)
show('sh2', sh2)

# Also check samples dir
samples_dir = f'{BASE}/pt_cache_fixed/random/train/samples'
print()
print(f'=== {samples_dir} ===')
if os.path.islink(samples_dir):
    print(f'symlink → {os.readlink(samples_dir)}')
    samples_dir = os.path.realpath(samples_dir)
names = os.listdir(samples_dir)[:5]
print(f'first 5 files: {names}')

if names:
    s = torch.load(os.path.join(samples_dir, names[0]), map_location='cpu', weights_only=False)
    print(f'sample 0 type: {type(s).__name__}, keys: {list(s.keys()) if isinstance(s,dict) else ""}')
    if isinstance(s, dict):
        for k in ['enzyme_id','substrate_id','sample_id','dock_id','global_id','enzyme_idx','substrate_idx']:
            if k in s: print(f'  {k}: {s[k]}')
