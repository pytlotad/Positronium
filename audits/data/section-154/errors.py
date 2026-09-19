# The combination that actually matters is <U> ~ <mu1.mu2> - 3<(mu1.n)(mu2.n)>,
# which the identity of 153c sets to zero.  Reported per trajectory with its
# standard error, because a RATIO of two averages is meaningless when the
# denominator sits near zero (audit 154).
import re,glob,math
def vec(t): return [float(x) for x in t.strip('()').split(',')]
def norm(v): return math.sqrt(sum(x*x for x in v))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
mu=9.2740100657e-24
scale=-4.1150/0.7499            # a mu^2 in GHz per unit structure, from 153
for phenomenon,prepared in ((1,0.7499),(2,-0.2406)):
    rows=[]
    for f in sorted(glob.glob('ph%d_w*_err.txt'%phenomenon)):
        cur=None
        for line in open(f):
            if line.startswith('TRAJ'):
                if cur: rows.append(cur)
                cur=None
            elif line.startswith('SPINPAIR'):
                d=dict(re.findall(r'(\w+)=(\S+)',line))
                cur=(vec(d['mu1']),vec(d['mu2']),vec(d['L']))
        if cur: rows.append(cur)
    combo=[];cosv=[];csq=[]
    for a,b,L in rows:
        ln=norm(L); n=[x/ln for x in L] if ln>0 else [0,0,1]
        iso=dot(a,b)/(mu*mu); ten=dot(a,n)*dot(b,n)/(mu*mu)
        combo.append(iso-3.0*ten)
        cosv.append(dot(a,b)/(norm(a)*norm(b)))
        m=[a[i]+b[i] for i in range(3)]; csq.append((norm(m)/mu)**2)
    k=len(combo); mean=lambda v:sum(v)/len(v)
    sd=lambda v:math.sqrt(sum((x-mean(v))**2 for x in v)/max(1,len(v)-1))
    se=sd(combo)/math.sqrt(k)
    print("phenomenon %d, %d trajectories"%(phenomenon,k))
    print("  <cos> terminal %+.4f  (prepared %+.4f)"%(mean(cosv),prepared))
    print("  <|m|^2/mu^2>   %.4f"%mean(csq))
    print("  <iso> - 3<ten> = %+.4f +- %.4f  ->  <U> = %+.4f +- %.4f GHz "
          "(%.2f sigma from zero)"%(mean(combo),se,mean(combo)*scale,se*abs(scale),
                                    abs(mean(combo))/se if se>0 else 0.0))
