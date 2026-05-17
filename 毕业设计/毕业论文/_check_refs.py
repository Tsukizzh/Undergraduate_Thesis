import re, sys, os
sys.stdout.reconfigure(encoding='utf-8')

bib_keys = set()
with open('reference.bib', encoding='utf-8') as f:
    for line in f:
        m = re.match(r'@\w+\{([^,]+),', line)
        if m:
            bib_keys.add(m.group(1).strip())

cite_re = re.compile(r'\\cite[a-z]*\{([^}]+)\}')
cited = {}
files = ['docs/chap01.tex','docs/chap02.tex','docs/chap03.tex','docs/chap04.tex','docs/chap05.tex','docs/chap06.tex','docs/abstract.tex']
for fn in files:
    if not os.path.exists(fn): continue
    with open(fn, encoding='utf-8') as f:
        for line in f:
            for m in cite_re.findall(line):
                for k in m.split(','):
                    k = k.strip()
                    cited.setdefault(k, []).append(fn)

print(f'.bib 中定义 : {len(bib_keys)} 条')
print(f'正文引用    : {len(cited)} 条')
print()
print('--- 已定义但正文未引用（孤儿） ---')
orphan = bib_keys - set(cited.keys())
for k in sorted(orphan): print(f'  • {k}')
if not orphan: print('  无')
print()
print('--- 正文引用但 .bib 未定义（缺失） ---')
missing = set(cited.keys()) - bib_keys
for k in sorted(missing): print(f'  • {k}')
if not missing: print('  无')
print()
print('--- 各文献被引用的章节 ---')
for k in sorted(cited):
    chaps = sorted(set([f.replace('docs/','').replace('.tex','') for f in cited[k]]))
    print(f'  {k:<22} -> {",".join(chaps)}')
