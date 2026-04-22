"""Count sample IDs in each pt_cache, compute intersections, report overlaps."""
import torch, glob, os

def load_per_sample_ids(cache_dir, split):
    """pt_cache (per-sample layout): read index.pt + samples/ dir names."""
    idx_path = os.path.join(cache_dir, split, 'index.pt')
    if not os.path.exists(idx_path):
        return None
    idx = torch.load(idx_path, map_location='cpu', weights_only=False)
    # Try common fields
    for k in ['sample_ids', 'global_ids', 'dock_ids', 'keys']:
        if k in idx:
            ids = idx[k]
            if hasattr(ids, 'tolist'):
                return set(int(x) for x in ids.tolist())
            return set(int(x) for x in ids)
    # Fall back: enumerate samples dir
    samples_dir = os.path.join(cache_dir, split, 'samples')
    if os.path.isdir(samples_dir):
        names = os.listdir(samples_dir)
        return set(int(n.replace('.pt','')) for n in names if n.endswith('.pt'))
    print(f'  ! unknown index format: keys={list(idx.keys())}')
    return None

def load_shard_ids(cache_dir, split):
    """shard layout: load every shard, collect sample_ids (dataset_id * 10M + local_id)."""
    shards = sorted(glob.glob(os.path.join(cache_dir, split, 'graph_*.pt')))
    if not shards:
        return None
    ids = set()
    for sp in shards:
        s = torch.load(sp, map_location='cpu', weights_only=False)
        # Identify ID field
        if 'dock_ids' in s:
            for x in s['dock_ids'].tolist(): ids.add(int(x))
        elif 'sample_ids' in s:
            for x in s['sample_ids'].tolist(): ids.add(int(x))
        elif 'enzyme_ids' in s and 'substrate_ids' in s:
            ez = s['enzyme_ids'].tolist()
            sb = s['substrate_ids'].tolist()
            for e, b in zip(ez, sb):
                ids.add(int(e) * 100000 + int(b))   # synthetic ID
        else:
            print(f'  ! shard {sp} no id field: keys={list(s.keys())[:8]}')
            return None
    return ids

BASE = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'

caches = [
    ('bare',  f'{BASE}/pt_cache_fixed/random',        'per_sample'),
    ('heme',  f'{BASE}/pt_cache_heme_fixed/random',   'shard'),
    ('geom',  f'{BASE}/pt_cache_geom_fixed',          'shard'),
]

all_ids = {}
for split in ['train', 'val', 'test']:
    print(f'=== {split} ===')
    sets = {}
    for name, path, layout in caches:
        if layout == 'per_sample':
            s = load_per_sample_ids(path, split)
        else:
            s = load_shard_ids(path, split)
        if s is None:
            print(f'  {name}: FAIL')
        else:
            sets[name] = s
            print(f'  {name}: {len(s)} samples')
    if len(sets) == 3:
        inter = sets['bare'] & sets['heme'] & sets['geom']
        union = sets['bare'] | sets['heme'] | sets['geom']
        print(f'  INTERSECTION (all 3): {len(inter)}')
        print(f'  UNION (any):          {len(union)}')
        for name, s in sets.items():
            dropped = len(s) - len(inter)
            print(f'    {name} loses: {dropped} ({100*dropped/len(s):.1f}%)')
    print()
