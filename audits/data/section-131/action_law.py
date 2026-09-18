# The action law behind the cascade (audit 131).  The emission rule is
# E_photon = hbar omega_orb = 2R/n^3 with R the pair's Rydberg, so
# E' = -R/n^2 - 2R/n^3 gives n' = n/sqrt(1+2/n) and therefore
#   Delta n = n [1 - (1+2/n)^(-1/2)],
# which tends to 1 only for n >> 1.  Checked against both measured ladders.
import re,sys,math
def law(n): return n*(1.0-(1.0+2.0/n)**-0.5)
print("closed form Delta n = n[1-(1+2/n)^-1/2]")
for n in (1.0,2.0,3.0,10.0,100.0,1000.0):
    print("  n = %8.1f  Delta n = %.6f  landing %.6f"%(n,law(n),n-law(n)))
print()
def ladder(name,seq):
    print(name)
    for a,b in zip(seq,seq[1:]):
        print("  %.6f -> %.6f   measured Delta %.6f   law %.6f   ratio %.4f"
              %(a,b,a-b,law(a),(a-b)/law(a)))
    steps=[a-b for a,b in zip(seq,seq[1:])]
    steps.sort()
    print("  median |Delta n| %.4f, range %.4f..%.4f"
          %(steps[len(steps)//2],steps[0],steps[-1]))
ladder("this run, seed 42 (prepared at n=1)",
       [1.000000,0.577269,0.272406,0.096794])
ladder("audit 29b, prepared at 4 a_pair (n=2)",
       [1.999995,1.414210,0.910168,0.508943,0.228165])
print()
# How much of the action leaves in photon steps and how much drifts away
# continuously between them.
rows=[]
for line in open(sys.argv[1] if len(sys.argv)>1 else 'err_42.txt'):
    f=dict(re.findall(r'(\w+)=(\S+)',line))
    if line.startswith('ACTION'): rows.append(('cp',float(f['n']),
        float(f['Jphi_over_h']),float(f['Jr_over_h'])))
    elif line.startswith('PHOTON_LADDER'): rows.append(('ph',float(f['n_E']),
        float(f['Jphi_over_h']),None))
first=[r for r in rows if r[0]=='cp'][0]; last=[r for r in rows if r[0]=='cp'][-1]
print("prepared  n %.9f  Jphi/h %.9f  Jr/h %.9f"%(first[1],first[2],first[3]))
print("stopped   n %.9f  Jphi/h %.9f  Jr/h %.9f"%(last[1],last[2],last[3]))
print("total action drop %.6f h"%(first[1]-last[1]))
cp=[r for r in rows if r[0]=='cp']
drift=sum(max(0.0,a[1]-b[1]) for a,b in zip(cp,cp[1:])
          if abs(a[1]-b[1])<0.05)
print("continuous drift between checkpoints (steps < 0.05 h) %.6f h"%drift)
jumps=sum(a[1]-b[1] for a,b in zip(cp,cp[1:]) if a[1]-b[1]>=0.05)
print("checkpoint-to-checkpoint jumps >= 0.05 h (the photons) %.6f h"%jumps)
jr=[abs(r[3]) for r in cp]
print("|Jr/h| along the whole collapse: max %.6f, mean %.6f"
      %(max(jr),sum(jr)/len(jr)))
