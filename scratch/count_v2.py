"""Use (enzyme_id, substrate_id) as common ID across per-sample / shard layouts."""
import torch, glob, os

def load_ids_bare(split):
    """pt_cache_fixed uses index.pt with enzyme_ids + substrate_ids tensors."""
    idx = torch.load(f'/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_fixed/random/{split}/index.pt',
                     map_location='cpu', weights_only=False)
    ez = idx['enzyme_ids'].tolist()
    sb = idx['substrate_ids'].tolist()
    return set(zip(ez, sb))

def load_ids_shard(base, split):
    """heme/geom: iterate shards, collect enzyme_ids/substrate_ids pairs."""
    shards = sorted(glob.glob(f'{base}/{split}/graph_*.pt'))
    ids = set()
    for sp in shards:
        s = torch.load(sp, map_location='cpu', weights_only=False)
        ez = s['enzyme_ids'].tolist()
        sb = s['substrate_ids'].tolist()
        ids.update(zip(ez, sb))
    return ids

for split in ['train', 'val', 'test']:
    print(f'=== {split} ===')
    bare = load_ids_bare(split)
    heme = load_ids_shard('/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_heme_fixed/random', split)
    geom = load_ids_shard('/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_geom_fixed', split)
    print(f'  bare: {len(bare):6d}')
    print(f'  heme: {len(heme):6d}')
    print(f'  geom: {len(geom):6d}')
    inter_all = bare & heme & geom
    print(f'  intersection (bare ∩ heme ∩ geom): {len(inter_all):6d}')
    print(f'  bare - heme: {len(bare - heme)}, heme - bare: {len(heme - bare)}')
    print(f'  bare - geom: {len(bare - geom)}, geom - bare: {len(geom - bare)}')
    print(f'  heme - geom: {len(heme - geom)}, geom - heme: {len(geom - heme)}')
    print()
    # Report what each cache would lose
    for name, s in [('bare', bare), ('heme', heme), ('geom', geom)]:
        lose = len(s) - len(inter_all)
        pct = 100 * lose / len(s) if s else 0
        print(f'    if unified to intersection: {name} loses {lose} ({pct:.2f}%)')
    print()
