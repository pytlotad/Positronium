# The averages of audit 153, recomputed on the TERMINAL moments (section 154).
# Each trajectory's block on stderr starts with a TRAJ marker; its last
# SPINPAIR line carries mu1, mu2 and the orbital vector L at annihilation.
import re,sys,glob,math
def vec(text):
    return [float(x) for x in text.strip('()').split(',')]
def norm(v): return math.sqrt(sum(x*x for x in v))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
mu=9.2740100657e-24
for phenomenon in (1,2):
    rows=[]
    for f in sorted(glob.glob(sys.argv[1]+'/ph%d_w*_err.txt'%phenomenon)):
        current=None
        for line in open(f):
            if line.startswith('TRAJ'):
                if current: rows.append(current)
                current=None
            elif line.startswith('SPINPAIR'):
                d=dict(re.findall(r'(\w+)=(\S+)',line))
                current=(vec(d['mu1']),vec(d['mu2']),vec(d['L']))
        if current: rows.append(current)
    if not rows:
        print("phenomenon %d: no trajectories"%phenomenon); continue
    cos=[];coh=[];cohSq=[];iso=[];ten=[]
    for first,second,orbital in rows:
        n1,n2=norm(first),norm(second)
        if n1<=0 or n2<=0: continue
        cos.append(dot(first,second)/(n1*n2))
        m=[first[i]+second[i] for i in range(3)]
        coh.append(norm(m)/mu); cohSq.append((norm(m)/mu)**2)
        ln=norm(orbital)
        nhat=[x/ln for x in orbital] if ln>0 else [0,0,1]
        iso.append(dot(first,second)/(mu*mu))
        ten.append(dot(first,nhat)*dot(second,nhat)/(mu*mu))
    k=len(cos); mean=lambda v:sum(v)/len(v)
    print("phenomenon %d, %d trajectories"%(phenomenon,k))
    print("  <cos> %+.4f   <|m|/mu> %.4f   <|m|^2/mu^2> %.4f"%(
        mean(cos),mean(coh),mean(cohSq)))
    print("  <mu1.mu2>/mu^2 %+.4f   <(mu1.n)(mu2.n)>/mu^2 %+.4f   "
          "ratio %+.4f  (1/3 = 0.3333)"%(
        mean(iso),mean(ten),mean(ten)/mean(iso) if mean(iso)!=0 else float('nan')))
