# The TRUE per-photon angular-momentum closure: Lold from CREM_LUPD and the
# value actually installed from CREM_LUSED, paired line by line (audit 140).
import re,sys,glob
R=6.802846579
for mode in ('classical','spin','axial'):
    rows=[]
    for f in sorted(glob.glob(mode+'_short_err_*.txt')):
        old=[];used=[];energies=[];levels=[]
        for line in open(f):
            if line.startswith('CREM_LUPD'):
                d=dict(re.findall(r'(\w+)=([-\d.e+]+)',line))
                old.append((float(d['Lold']),float(d['e2used'])))
            elif line.startswith('CREM_LUSED'):
                used.append(float(line.split('=')[1]))
            elif line.startswith('PHOTON_LADDER'):
                d=dict(re.findall(r'(\w+)=([-\d.e+]+)',line))
                energies.append(float(d['E_eV'])); levels.append(float(d['n_E']))
        for i in range(min(len(old),len(used))):
            quantum=2.0*R/levels[i]**3 if i<len(levels) else float('nan')
            rows.append((f.split('_')[-1].replace('.txt',''),i+1,old[i][0],
                         used[i],old[i][0]-used[i],old[i][1],
                         energies[i] if i<len(energies) else float('nan'),
                         quantum))
    print("=== %s ==="%mode)
    print("%6s %6s %11s %11s %12s %10s %11s %10s"%("seed","photon","L before",
          "L after","surrendered","e^2","E [eV]","E/quantum"))
    for r in rows:
        print("%6s %6d %11.6f %11.6f %12.6f %10.6f %11.6f %10.4f"%r)
    if rows:
        s=sorted(r[4] for r in rows)
        neg=sum(1 for r in rows if r[5]<0)
        eq=[r[6]/r[7] for r in rows if r[7]==r[7] and r[7]>0]
        print("emissions %d; surrendered mean %.4f range %.4f..%.4f; "
              "e^2<0 in %d; E/quantum mean %.4f\n"%(
              len(rows),sum(s)/len(s),s[0],s[-1],neg,
              sum(eq)/len(eq) if eq else float('nan')))
