#!/usr/bin/env python3
"""Section 96: which checkpoint each photon falls in, for pairs of trajectories
from different lifetime branches.

Usage: tools/checkpoint_photon_trace.py <dir with pairs.txt, trace_<seed>.out/.err>
"""
import sys, gzip, os
import re, sys
SP=sys.argv[1] if len(sys.argv)>1 else '.'
def load(seed):
    cps=[]; photons=[]
    for l in (open(f'{SP}/trace_{seed}.err') if os.path.exists(f'{SP}/trace_{seed}.err') else gzip.open(f'{SP}/trace_{seed}.err.gz','rt')):
        if l.startswith('CREM_SKIP'):
            g=lambda k: float(re.search(k+r'=(\S+)',l).group(1))
            cps.append(dict(t=g('t'),skip=int(g('skip')),req=g('requested'),
                            E=g('Emag'),L=g('Lorb'),mm=g('mm'),period=g('period')))
        elif l.startswith('PHOTON_LADDER'):
            photons.append((len(cps),float(re.search(r'E_eV=(\S+)',l).group(1)),
                            float(re.search(r'a/a_pair=(\S+)',l).group(1))))
    res=open(f'{SP}/trace_{seed}.out').read().strip()
    return cps,photons,res
pairs=[l.split()[:2] for l in open(f'{SP}/pairs.txt')]
for a,b in pairs:
    print('='*100)
    for s in (a,b):
        cps,ph,res=load(s)
        print(f'--- ziarno {s}: {res}')
        print(f'    punktow kontrolnych: {len(cps)}')
        for k,(idx,e,ap) in enumerate(ph):
            t=cps[idx-1]['t']*1e12 if idx>0 else 0.0
            tn=cps[idx]['t']*1e12 if idx<len(cps) else float('nan')
            print(f'    foton {k+1}: po punkcie {idx} (t {t:.4f} ps, nastepny {tn:.4f} ps), E {e:.3f} eV, a/a_pair {ap:.4f}')
        print('    pierwsze 8 punktow: t[ps] skip requested E L mm')
        for c in cps[:8]:
            print(f"      {c['t']*1e12:10.4f} {c['skip']:7d} {c['req']:14.4f} {c['E']:.6e} {c['L']:.6f} {c['mm']:.4f}")
