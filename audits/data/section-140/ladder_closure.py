# The per-photon angular-momentum closure, read from the PHOTON_LADDER trace,
# which prints Jphi_over_h = L/hbar at each emission (audit 140).  CREM_LUPD
# cannot be used here: under CREM_SPIN_MAGNITUDE it still prints the CLASSICAL
# candidate as Lnew, not the value the switch installs.
import re,sys,glob
R=6.802846579
print("%-16s %6s %12s %12s %12s %12s %10s"%("run","photon","E [eV]",
      "L/hbar","surrendered","2R/n^3 [eV]","E/quantum"))
sur=[]
for f in sorted(glob.glob(sys.argv[1])):
    rows=[]
    for line in open(f):
        if not line.startswith('PHOTON_LADDER'): continue
        d=dict(re.findall(r'(\w+)=([-\d.e+]+)',line))
        rows.append((float(d['E_eV']),float(d['n_E']),
                     float(d['Jphi_over_h'])))
    for i,(e,n,l) in enumerate(rows):
        nxt=rows[i+1][2] if i+1<len(rows) else None
        s=(l-nxt) if nxt is not None else None
        if s is not None: sur.append(s)
        quantum=2.0*R/n**3
        print("%-16s %6d %12.6f %12.6f %12s %12.6f %10.4f"%(
            f.replace('spin_err_','seed ').replace('.txt',''),i+1,e,l,
            "%.6f"%s if s is not None else "-",quantum,e/quantum))
if sur:
    sur.sort()
    print("\nemissions with a successor %d; surrendered mean %.4f, "
          "range %.4f..%.4f"%(len(sur),sum(sur)/len(sur),sur[0],sur[-1]))
