#!/usr/bin/env python3
"""Section 95: pre-registered tests H1-H3 of lifetime against start angle from
lifetime_vs_start_angle.cpp output (permutation p-values).

Usage: tools/lifetime_vs_start_angle.py <dir with step2_[0-3].out>
"""
import sys
import re, math, random
from collections import Counter
SP=sys.argv[1] if len(sys.argv)>1 else '.'
rows=[]
for w in range(4):
    try: txt=open(f'{SP}/step2_{w}.out').read()
    except FileNotFoundError: continue
    for m in re.finditer(r'RESULT seed (\d+) cos0 (\S+) t_ps (\S+) outcome (\d+) photons (\d+)',txt):
        rows.append((int(m.group(1)),float(m.group(2)),float(m.group(3)),int(m.group(4)),int(m.group(5))))
print(f'ukonczonych wynikow: {len(rows)}')
good=[r for r in rows if r[3]==0]
print(f'bez cenzury (outcome 0): {len(good)}, ocenzurowane lub awaria: {len(rows)-len(good)}')
if not good: raise SystemExit
print('czasy zycia zaokraglone do 0.1 ps:', sorted(Counter(round(r[2],1) for r in good).items()))
print('liczba fotonow:', sorted(Counter(r[4] for r in good).items()))
def spearman(x,y):
    def rank(v):
        o=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0
        while i<len(o):
            j=i
            while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
            for k in range(i,j+1): r[o[k]]=(i+j)/2
            i=j+1
        return r
    a,b=rank(x),rank(y); ma,mb=sum(a)/len(a),sum(b)/len(b)
    num=sum((p-ma)*(q-mb) for p,q in zip(a,b))
    den=math.sqrt(sum((p-ma)**2 for p in a)*sum((q-mb)**2 for q in b))
    return num/den if den>0 else float('nan')
def perm_p(x,y,n=20000,seed=12345):
    rho=spearman(x,y); rng=random.Random(seed); yy=list(y); hit=0
    for _ in range(n):
        rng.shuffle(yy)
        if abs(spearman(x,yy))>=abs(rho)-1e-12: hit+=1
    return rho,(hit+1)/(n+1)
long_=[r for r in good if r[2]>208]; short=[r for r in good if r[2]<=208]
print(f'\ngalaz dluga (>208 ps): {len(long_)}, krotka: {len(short)}')
print('\n=== H1: cos0 < -0.4 ===')
low=[r for r in good if r[1]<-0.4]
lowlong=sum(1 for r in low if r[2]>208)
rest=[r for r in good if r[1]>=-0.4]
f=sum(1 for r in rest if r[2]>208)/len(rest) if rest else float('nan')
print(f'  trajektorii z cos0 < -0.4: {len(low)}, w galezi dlugiej: {lowlong}')
print(f'  udzial galezi dlugiej przy cos0 >= -0.4: {f:.3f}')
if lowlong==0 and low:
    print(f'  szansa zera przy niezaleznosci, (1-f)^k: {(1-f)**len(low):.2e}')
print('\n=== H2: galaz ~ cos0 ===')
rho,p=perm_p([r[1] for r in good],[1.0 if r[2]>208 else 0.0 for r in good])
print(f'  Spearman {rho:+.3f}, p (permutacje, dwustronnie) {p:.4f}, n={len(good)}')
print('\n=== H3: czas zycia ~ cos0 wewnatrz galezi ===')
for name,g in (('krotka',short),('dluga',long_)):
    if len(g)<5: print(f'  {name}: za malo ({len(g)})'); continue
    rho,p=perm_p([r[1] for r in g],[r[2] for r in g])
    ts=[r[2] for r in g]
    print(f'  {name}: n={len(g)}, czas {min(ts):.4f}..{max(ts):.4f} ps, Spearman {rho:+.3f}, p {p:.4f}')
print('\n=== opisowo: galaz dluga w przedzialach cos0 ===')
e=[-1,-0.6,-0.2,0.2,0.6,1.0001]
for i in range(5):
    g=[r for r in good if e[i]<=r[1]<e[i+1]]
    print(f'  [{e[i]:+.1f},{min(e[i+1],1):+.1f})  n={len(g):2d}  dluga {sum(1 for r in g if r[2]>208):2d}')
