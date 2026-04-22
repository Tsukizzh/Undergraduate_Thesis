import re, sys

def parse_epochs(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        txt = f.read()
    # Match patterns like: Epoch NN: 100%|...| 64/64 [MM:SS<...,
    pat = re.compile(r"Epoch (\d+): 100%\|[^|]+\|\s+(\d+)/(\d+)\s+\[(\d+):(\d+)<")
    seen = {}
    for m in pat.finditer(txt):
        ep = int(m.group(1))
        cur = int(m.group(2))
        tot = int(m.group(3))
        mm = int(m.group(4))
        ss = int(m.group(5))
        if cur == tot:
            seen[ep] = mm * 60 + ss
    return seen

a = parse_epochs("/root/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP003_fixed/logs/train.log")
b = parse_epochs("/root/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP002a_fixed/logs/nohup.out")

def report(name, d):
    print(name + ":")
    print(f"  epochs logged: {len(d)}")
    if not d:
        return None
    vals = sorted(d.values())
    mean = sum(vals) / len(vals)
    med = vals[len(vals) // 2]
    print(f"  min/max:  {min(vals)}s / {max(vals)}s")
    print(f"  mean:     {mean:.1f}s")
    print(f"  median:   {med}s")
    print(f"  first 5:  {[d[k] for k in sorted(d)[:5]]}")
    print(f"  last 5:   {[d[k] for k in sorted(d)[-5:]]}")
    return mean

ma = report("EXP003_fixed", a)
print()
mb = report("EXP002a_fixed", b)

if ma and mb:
    print()
    print(f"EXP002a_fixed / EXP003_fixed mean ratio: {mb/ma:.3f}x")
    print(f"abs diff per epoch: +{mb-ma:.1f}s")
