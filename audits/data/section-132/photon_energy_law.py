# R132b/R132c: is the photon energy the Kepler correspondence quantum 2R/n^3,
# and is the landing exact once the MEASURED energy is used (audit 132)?
import re,sys,glob
R=float(sys.argv[1])
print("R = %.9f eV (the pair's ground binding, from the model)"%R)
print("%-16s %10s %14s %14s %8s %12s %12s %8s"%("run","n_before",
      "E measured eV","2R/n^3 eV","ratio","landing meas","landing law","ratio"))
for f in sorted(glob.glob(sys.argv[2])):
    ph=[];fin=None
    for line in open(f):
        if line.startswith('PHOTON_LADDER'):
            d=dict(re.findall(r'(\w+)=([0-9.e+-]+)',line))
            ph.append((float(d['n_E']),float(d['E_eV'])))
        elif line.startswith('ACTION'):
            fin=float(dict(re.findall(r'(\w+)=([0-9.e+-]+)',line))['n'])
    landings=[p[0] for p in ph[1:]]+[fin]
    for (n,e),land in zip(ph,landings):
        kepler=2.0*R/n**3
        exact=(R/(R/n**2+e))**0.5
        print("%-16s %10.6f %14.6f %14.6f %8.4f %12.6f %12.6f %8.4f"%(
            f.replace('_err','').replace('.txt',''),n,e,kepler,e/kepler,
            land,exact,land/exact))
