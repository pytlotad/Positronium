import re, math, sys
hbar=1.054571817e-34
rows=[]
pat=re.compile(r"SPINPAIR dt=([0-9.e+-]+) mu1=([-0-9.e+]+),([-0-9.e+]+),([-0-9.e+]+) "
               r"mu2=([-0-9.e+]+),([-0-9.e+]+),([-0-9.e+]+) L=([-0-9.e+]+),([-0-9.e+]+),([-0-9.e+]+) "
               r"ratio1=([-0-9.e+]+) ratio2=([-0-9.e+]+)")
for line in open("trace.txt"):
    m=pat.match(line)
    if m: rows.append([float(x) for x in m.groups()])
print("podkrokow", len(rows))
if not rows: sys.exit(0)
def sub(a,b): return [a[i]-b[i] for i in range(3)]
def add(a,b): return [a[i]+b[i] for i in range(3)]
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def norm(a): return math.sqrt(dot(a,a))
def scale(a,s): return [x*s for x in a]
t=0.0; data=[]
for r in rows:
    dt=r[0]; mu1=r[1:4]; mu2=r[4:7]; L=r[7:10]; g1=r[10]; g2=r[11]
    S1=scale(mu1,1.0/g1); S2=scale(mu2,1.0/g2)
    Stot=add(S1,S2); n=scale(Stot,1.0/max(norm(Stot),1e-300))
    # perpendicular parts are opposite; take S1's, and an azimuth frame from L
    perp1=sub(S1,scale(n,dot(S1,n)))
    e1=sub(L,scale(n,dot(L,n))); e1n=norm(e1)
    if e1n<=0: continue
    e1=scale(e1,1.0/e1n); e2=cross(n,e1)
    psi=math.atan2(dot(perp1,e2),dot(perp1,e1))
    costheta=dot(S1,S2)/(norm(S1)*norm(S2))
    t+=dt
    data.append((t,costheta,psi,norm(Stot)/hbar,norm(S1)/hbar))
# unwrap psi
psi=[d[2] for d in data]
un=[psi[0]]
for k in range(1,len(psi)):
    d=psi[k]-psi[k-1]
    while d> math.pi: d-=2*math.pi
    while d<-math.pi: d+=2*math.pi
    un.append(un[-1]+d)
cos=[d[1] for d in data]; tm=[d[0] for d in data]
print(f"cos theta: start {cos[0]:+.4f} koniec {cos[-1]:+.4f} min {min(cos):+.4f} max {max(cos):+.4f}")
print(f"psi przebiega {un[-1]-un[0]:.3f} rad = {(un[-1]-un[0])/(2*math.pi):.2f} obrotow")
# libration cycles: extrema of cos
ext=[k for k in range(1,len(cos)-1) if (cos[k]-cos[k-1])*(cos[k+1]-cos[k])<0]
print("ekstremow kata:", len(ext))
s=0.5  # |S_i| in hbar
if len(ext)>=3:
    for label,(a,b) in [("pierwszy cykl",(ext[0],ext[2])),("ostatni cykl",(ext[-3],ext[-1]))]:
        period=tm[b]-tm[a]; f=1.0/period if period>0 else float('nan')
        # action J = (1/2pi) integral s cos(theta) dpsi over the cycle, in hbar
        J=0.0
        for k in range(a,b):
            J+=s*0.5*(cos[k]+cos[k+1])*(un[k+1]-un[k])
        J/= (2*math.pi)
        print(f"{label}: okres {period:.4e} s, f_lib {f/1e9:.3f} GHz, "
              f"h f_lib/h = {f/1e9:.3f} GHz wobec 203.3941, dzialanie J = {J:.4f} hbar, "
              f"|d psi| {abs(un[b]-un[a]):.3f} rad")
print(f"|S_pary|/hbar start {data[0][3]:.4f} koniec {data[-1][3]:.4f}; |S_1|/hbar {data[0][4]:.6f}")
