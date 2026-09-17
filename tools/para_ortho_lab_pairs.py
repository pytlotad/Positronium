#!/usr/bin/env python3
"""Sections 98 and 100b: paired ortho - para lab lifetimes (sign test, exact
Wilcoxon), the free-draw correlation, and the tilt correlation when tilt.txt
is present.

Usage: tools/para_ortho_lab_pairs.py <dir with poA_*, poB_*, poC_*.out [tilt.txt]>
"""
import sys
import re, glob, math, random
SP=sys.argv[1] if len(sys.argv)>1 else '.'
def pairs(prefix):
    out={}
    for f in glob.glob(f'{SP}/{prefix}_*.out'):
        for m in re.finditer(r'PAIR seed (\d+) para_lab_ps (\S+) orto_lab_ps (\S+) para_out (\d+) orto_out (\d+) para_ph (\d+) orto_ph (\d+)',open(f).read()):
            out[int(m.group(1))]=(float(m.group(2)),float(m.group(3)),int(m.group(4)),int(m.group(5)),int(m.group(6)),int(m.group(7)))
    return out
def sign_p(pos,n):
    k=min(pos,n-pos); return min(1.0,2*sum(math.comb(n,i) for i in range(k+1))/2**n)
def wilcoxon_p(d):
    d=[x for x in d if x!=0]; n=len(d)
    if n==0: return float('nan')
    a=sorted(range(n),key=lambda i:abs(d[i])); r=[0.0]*n; i=0
    while i<n:
        j=i
        while j+1<n and abs(d[a[j+1]])==abs(d[a[i]]): j+=1
        for q in range(i,j+1): r[a[q]]=(i+j)/2+1
        i=j+1
    W=sum(r[i] for i in range(n) if d[i]>0)
    # exact distribution with integer-doubled ranks
    rr=[int(round(2*x)) for x in r]; tot=sum(rr)
    dist={0:1}
    for v in rr:
        nd=dict(dist)
        for s,c in dist.items(): nd[s+v]=nd.get(s+v,0)+c
        dist=nd
    w2=int(round(2*W)); mean=tot/2
    ext=sum(c for s,c in dist.items() if abs(s-mean)>=abs(w2-mean)-1e-9)
    return ext/2**n
def report(tag,P):
    good={s:v for s,v in P.items() if v[2]==0 and v[3]==0}
    print(f'=== {tag}: par {len(P)}, bez cenzury {len(good)} ===')
    if not good: return None
    d=[(v[1]-v[0])*1000 for v in good.values()]  # fs
    n=len(d); m=sum(d)/n; sd=math.sqrt(sum((x-m)**2 for x in d)/(n-1)) if n>1 else float('nan')
    pos=sum(1 for x in d if x>0); neg=sum(1 for x in d if x<0)
    print(f'  orto - para: srednia {m:+.3f} fs, odchylenie {sd:.3f} fs, blad sredniej {sd/math.sqrt(n):.3f} fs')
    print(f'  orto dluzsze w {pos}, krotsze w {neg}, rowne {n-pos-neg}; test znakow p = {sign_p(pos,pos+neg):.4f}; Wilcoxon p = {wilcoxon_p(d):.4f}')
    print(f'  para srednio {sum(v[0] for v in good.values())/n:.6f} ps, orto {sum(v[1] for v in good.values())/n:.6f} ps')
    phd=sum(1 for v in good.values() if v[4]!=v[5]); print(f'  liczba fotonow rozna w {phd} parach')
    for s in sorted(good): v=good[s]; print(f'    {s}  para {v[0]:.6f}  orto {v[1]:.6f}  roznica {(v[1]-v[0])*1000:+8.3f} fs  fotony {v[4]}/{v[5]}')
    return dict(zip(sorted(good),[ (good[s][1]-good[s][0])*1000 for s in sorted(good)]))
A=report('A: kwantyzacja, s_max 0.30',pairs('poA'))
B=report('B: kwantyzacja, s_max 0.15',pairs('poB'))
if A and B:
    common=sorted(set(A)&set(B))
    if common:
        print(f'\n=== A kontra B na wspolnych ziarnach ({len(common)}) ===')
        ma=sum(A[s] for s in common)/len(common); mb=sum(B[s] for s in common)/len(common)
        print(f'  srednia roznica: A {ma:+.3f} fs, B {mb:+.3f} fs')
        same=sum(1 for s in common if A[s]*B[s]>0); print(f'  zgodny znak w {same} z {len(common)} par')
C=[]
for f in glob.glob(f'{SP}/poC_*.out'):
    for m in re.finditer(r'FREE seed (\d+) cos0 (\S+) lab_ps (\S+) out (\d+) ph (\d+)',open(f).read()):
        C.append((float(m.group(2)),float(m.group(3)),int(m.group(4))))
if C:
    g=[c for c in C if c[2]==0]
    print(f'\n=== C: losowanie swobodne, s_max 0.30: {len(C)} trajektorii, bez cenzury {len(g)} ===')
    def rank(v):
        o=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0
        while i<len(o):
            j=i
            while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
            for q in range(i,j+1): r[o[q]]=(i+j)/2
            i=j+1
        return r
    def sp(x,y):
        a,b=rank(x),rank(y); ma,mb=sum(a)/len(a),sum(b)/len(b)
        return sum((p-ma)*(q-mb) for p,q in zip(a,b))/math.sqrt(sum((p-ma)**2 for p in a)*sum((q-mb)**2 for q in b))
    x=[c[0] for c in g]; y=[c[1] for c in g]; rho=sp(x,y)
    rng=random.Random(7); yy=list(y); hit=0; N=20000
    for _ in range(N):
        rng.shuffle(yy)
        if abs(sp(x,yy))>=abs(rho)-1e-12: hit+=1
    print(f'  Spearman czas~cos0 {rho:+.3f}, p {(hit+1)/(N+1):.4f}')
    print(f'  czasy: {min(y):.6f} .. {max(y):.6f} ps (rozrzut {(max(y)-min(y))*1000:.2f} fs)')
    pa=[c[1] for c in g if c[0]>=0.5]; po=[c[1] for c in g if c[0]<0.5]
    for name,v in (('etykieta para',pa),('etykieta orto',po)):
        if len(v)>1:
            mv=sum(v)/len(v); sd=math.sqrt(sum((t-mv)**2 for t in v)/(len(v)-1))
            print(f'  {name}: n={len(v)}, srednia {mv:.6f} ps, blad {sd/math.sqrt(len(v))*1000:.3f} fs')

# Section 100b: correlation with the first moment's tilt, when tilt.txt exists.
import os
tp=os.path.join(SP,'tilt.txt')
if os.path.exists(tp) and A:
    cz={int(m.group(1)):float(m.group(2)) for m in re.finditer(r'TILT (\d+) cosz (\S+)',open(tp).read())}
    S=sorted(s for s in A if s in cz)
    def rank(v):
        o=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0
        while i<len(o):
            j=i
            while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
            for q in range(i,j+1): r[o[q]]=(i+j)/2
            i=j+1
        return r
    def spr(x,y):
        a,b=rank(x),rank(y); ma,mb=sum(a)/len(a),sum(b)/len(b)
        return sum((p-ma)*(q-mb) for p,q in zip(a,b))/math.sqrt(sum((p-ma)**2 for p in a)*sum((q-mb)**2 for q in b))
    x=[cz[s] for s in S]; y=[A[s] for s in S]; rho=spr(x,y)
    rng=random.Random(3); yy=list(y); hit=0
    for _ in range(20000):
        rng.shuffle(yy)
        if abs(spr(x,yy))>=abs(rho)-1e-12: hit+=1
    up=[A[s] for s in S if cz[s]>0]; dn=[A[s] for s in S if cz[s]<0]
    print(f'\n=== tilt (section 100b): Spearman ortho-para ~ cos z {rho:+.3f}, p {(hit+1)/20001:.4f}')
    print(f'  mean ortho-para: cos z > 0 {sum(up)/len(up):+.1f} fs (n={len(up)}), cos z < 0 {sum(dn)/len(dn):+.1f} fs (n={len(dn)})')
