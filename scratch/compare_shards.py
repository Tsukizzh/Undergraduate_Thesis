import torch, glob, os

for name, path in [('heme', '/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_heme_fixed/random/train'),
                   ('geom', '/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_geom_fixed/train')]:
    shards = sorted(glob.glob(os.path.join(path, 'graph_*.pt')))
    print(f'=== {name} ({len(shards)} shards) ===')
    s0 = torch.load(shards[0], weights_only=False)
    print(f'  shard keys: {list(s0.keys())}')
    for k in list(s0.keys())[:30]:
        v = s0[k]
        shape = getattr(v, 'shape', None)
        dtype = getattr(v, 'dtype', None)
        print(f'    {k}: shape={shape} dtype={dtype}')
    print()
