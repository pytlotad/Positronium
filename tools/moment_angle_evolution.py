#!/usr/bin/env python3
"""Section 94: flatness of the mutual-angle distribution, drift, libration
width and lifetime branches from moment_angle_evolution.cpp output.

Usage: tools/moment_angle_evolution.py <dir with cosevo_[0-3].out/.err>
"""
import sys
import gzip, os
def _open_err(SP,w):
    plain=f'{SP}/cosevo_{w}.err'
    return open(plain) if os.path.exists(plain) else gzip.open(plain+'.gz','rt')
import re, glob, math, sys
SP=sys.argv[1] if len(sys.argv)>1 else '.'
cre=re.compile(r'cos=([-0-9.eE+]+)')
trajs=[]
for w in range(4):
    try: out=open(f'{SP}/cosevo_{w}.out').read().splitlines()
    except FileNotFoundError: continue
    res={}
    for l in out:
        m=re.match(r'RESULT seed (\d+) t_ps (\S+) outcome (\d+) photons (\d+)',l)
        if m: res[int(m.group(1))]=(float(m.group(2)),int(m.group(3)),int(m.group(4)))
    cur=None; samples={}
    for l in _open_err(SP,w):
        if l.startswith('TRAJ'):
            cur=int(l.split()[3]); samples[cur]=[]
        elif l.startswith('ALIGN') and cur is not None:
            m=cre.search(l)
            if m: samples[cur].append(float(m.group(1)))
    for seed,c in samples.items():
        if seed not in res: continue          # still running: skip, never guess
        trajs.append((seed,c,res[seed]))
def mean(x): return sum(x)/len(x) if x else float('nan')
def frac(x,f): return sum(1 for v in x if f(v))/len(x) if x else float('nan')
done=[t for t in trajs if len(t[1])>=20]
print(f'trajektorii ukonczonych z probkami: {len(done)}  (odrzucone bez probek: {len(trajs)-len(done)})')
cens=[t for t in done if t[2][1]!=0]
print(f'ocenzurowane lub awaria (outcome != 0): {len(cens)}')
if not done: sys.exit()
print()
print(' seed    probek  cos0     sr.pierwsze10%  sr.ostatnie10%  min     max    |cos|>0.9 pocz/kon  t_ps      fot outc')
rows=[]
for seed,c,(t,o,p) in sorted(done,key=lambda t:t[1][0]):
    n=len(c); k=max(1,n//10); a=c[:k]; b=c[-k:]
    rows.append((c[0],mean(a),mean(b),min(c),max(c),frac(a,lambda v:abs(v)>0.9),frac(b,lambda v:abs(v)>0.9),t,o,p,a,b))
    print(f' {seed:<7d} {n:<7d} {c[0]:+.4f}  {mean(a):+.4f}         {mean(b):+.4f}         {min(c):+.3f}  {max(c):+.3f}  {frac(a,lambda v:abs(v)>0.9):.2f}/{frac(b,lambda v:abs(v)>0.9):.2f}          {t:<9.4f} {p:<3d} {o}')
bins=[-1.0,-0.6,-0.2,0.2,0.6,1.0001]
def hist(vals):
    h=[0]*5
    for v in vals:
        for i in range(5):
            if bins[i]<=v<bins[i+1]: h[i]+=1; break
    s=sum(h); return [x/s for x in h] if s else h
print()
print('rozklad (5 przedzialow od -1 do +1)      [-1,-.6) [-.6,-.2) [-.2,.2) [.2,.6) [.6,1]')
print('  cos na starcie, po trajektoriach       ', '  '.join(f'{x:.3f}' for x in hist([r[0] for r in rows])))
print('  probki z pierwszych 10%, zbiorczo      ', '  '.join(f'{x:.3f}' for x in hist([v for r in rows for v in r[10]])))
print('  probki z ostatnich 10%, zbiorczo       ', '  '.join(f'{x:.3f}' for x in hist([v for r in rows for v in r[11]])))
print('  rozklad jednostajny, odniesienie        0.200  0.200  0.200  0.200  0.200')
drift=[r[2]-r[1] for r in rows]
print()
print(f'dryf sredniego cos (ostatnie - pierwsze 10%): srednia {mean(drift):+.4f}, '
      f'odchylenie {math.sqrt(mean([(d-mean(drift))**2 for d in drift])):.4f}, n={len(drift)}')
print(f'udzial |cos|>0.9: pierwsze 10% {mean([r[5] for r in rows]):.3f}, ostatnie 10% {mean([r[6] for r in rows]):.3f} (jednostajny: 0.100)')

# --- next part ---
import sys
import re, math
SP=sys.argv[1] if len(sys.argv)>1 else '.'
cre=re.compile(r'cos=([-0-9.eE+]+)')
T=[]
for w in range(4):
    res={}
    for l in open(f'{SP}/cosevo_{w}.out'):
        m=re.match(r'RESULT seed (\d+) t_ps (\S+) outcome (\d+) photons (\d+)',l)
        if m: res[int(m.group(1))]=(float(m.group(2)),int(m.group(4)))
    cur=None; S={}
    for l in _open_err(SP,w):
        if l.startswith('TRAJ'): cur=int(l.split()[3]); S[cur]=[]
        elif l.startswith('ALIGN'):
            m=cre.search(l)
            if m: S[cur].append(float(m.group(1)))
    for s,c in S.items(): T.append((c[0],min(c),max(c),res[s][0],c))
T.sort()
print('=== 1. szerokosc pasma libracji wedlug kata startowego ===')
edges=[-1,-0.6,-0.2,0.2,0.6,1.0001]
for i in range(5):
    g=[t for t in T if edges[i]<=t[0]<edges[i+1]]
    if not g: continue
    w=[t[2]-t[1] for t in g]
    print(f'  cos0 w [{edges[i]:+.1f},{min(edges[i+1],1):+.1f})  n={len(g):2d}  szerokosc srednia {sum(w)/len(w):.3f}  min {min(w):.3f}  max {max(w):.3f}')
print()
print('=== 2. czy start lezy na brzegu pasma (w granicach 0.02) ===')
edge=sum(1 for t in T if abs(t[0]-t[1])<0.02 or abs(t[0]-t[2])<0.02)
print(f'  na brzegu: {edge} z {len(T)}')
print()
print('=== 3. galaz czasu zycia wedlug kata startowego ===')
for i in range(5):
    g=[t for t in T if edges[i]<=t[0]<edges[i+1]]
    if not g: continue
    long_=sum(1 for t in g if t[3]>210)
    print(f'  cos0 w [{edges[i]:+.1f},{min(edges[i+1],1):+.1f})  n={len(g):2d}  galaz 216.7 ps: {long_}  galaz 199.4 ps: {len(g)-long_}')
longAll=sum(1 for t in T if t[3]>210)
print(f'  razem: 216.7 ps {longAll}, 199.4 ps {len(T)-longAll}')
vals=sorted(set(round(t[3],1) for t in T)); print('  wartosci czasu zycia wystepujace:',vals)
# Spearman: branch indicator vs cos0
def rank(x):
    o=sorted(range(len(x)),key=lambda i:x[i]); r=[0]*len(x)
    i=0
    while i<len(o):
        j=i
        while j+1<len(o) and x[o[j+1]]==x[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2
        i=j+1
    return r
a=rank([t[0] for t in T]); b=rank([1.0 if t[3]>210 else 0.0 for t in T])
ma=sum(a)/len(a); mb=sum(b)/len(b)
num=sum((x-ma)*(y-mb) for x,y in zip(a,b)); den=math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
print(f'  korelacja rangowa Spearmana galaz~cos0: {num/den:+.3f} (n=40)')
print()
print('=== 4. |cos|>0.9, per trajektoria, a nie zbiorczo ===')
up=down=same=0
for t in T:
    c=t[4]; k=max(1,len(c)//10)
    f1=sum(1 for v in c[:k] if abs(v)>0.9)/k; f2=sum(1 for v in c[-k:] if abs(v)>0.9)/k
    if f2>f1: up+=1
    elif f2<f1: down+=1
    else: same+=1
print(f'  wzrost {up}, spadek {down}, bez zmiany {same}')
