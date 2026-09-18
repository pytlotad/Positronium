# The M1 power of the libration, taken from the trace's own total moment
# (audit 135).  P = mu0 |d2(mu1+mu2)/dt2|^2 / (6 pi c^3), the same formula
# coherentMagneticDipoleOrbitAveragedEmission evaluates analytically.
import re,sys,math
mu0=4.0e-7*math.pi; c=299792458.0; hbar=1.054571817e-34
trace=sys.argv[1]; share=sys.argv[2] if len(sys.argv)>2 else None
t=[];m=[]
now=0.0
for line in open(trace):
    if not line.startswith('SPINPAIR'): continue
    f=dict(re.findall(r'(\w+)=(\S+)',line))
    vec=lambda k:[float(x) for x in f[k].strip('()').split(',')]
    now+=float(f['dt']); t.append(now)
    a,b=vec('mu1'),vec('mu2'); m.append([a[i]+b[i] for i in range(3)])
print("substeps %d, span %.4e s"%(len(t),t[-1]-t[0]))
def secondDerivative(stride):
    out=[]
    for i in range(stride,len(t)-stride):
        h1=t[i]-t[i-stride]; h2=t[i+stride]-t[i]
        if h1<=0 or h2<=0: continue
        d=[2.0*(h2*m[i-stride][k]-(h1+h2)*m[i][k]+h1*m[i+stride][k])
           /(h1*h2*(h1+h2)) for k in range(3)]
        out.append((t[i],sum(x*x for x in d)))
    return out
def power(sq): return mu0*sq/(6.0*math.pi*c**3)
for stride in (1,2):
    d=secondDerivative(stride)
    mean=sum(s for _,s in d)/len(d)
    print("stencil %d: <|m''|^2> %.6e  -> P_M1 %.6e W  (max %.6e W)"%(
        stride,mean,power(mean),power(max(s for _,s in d))))
d1=secondDerivative(1); d2=secondDerivative(2)
mean1=sum(s for _,s in d1)/len(d1); mean2=sum(s for _,s in d2)/len(d2)
print("R135e stencil agreement: %.4f"%(mean2/mean1))
print("mean substep %.4e s; libration period 1.1044e-11 s -> %.1f samples"%(
    (t[-1]-t[0])/len(t),1.1044e-11/((t[-1]-t[0])/len(t))))
P=power(mean1)
omega=2.0*math.pi*90.551e9; J=0.0985*hbar
Elib=J*omega; Tlib=2.0*math.pi/omega
print("E_lib = J omega = %.6e J = %.6e eV"%(Elib,Elib/1.602176634e-19))
print("radiated per libration period %.6e J, i.e. %.3e of E_lib"%(
    P*Tlib,P*Tlib/Elib))
print("damping time E_lib/P = %.6e s = %.3e collapse times (30.672 ps)"%(
    Elib/P,Elib/P/30.672e-12))
if share:
    rows=[]
    for line in open(share):
        if not line.startswith('CREM_M1_SHARE'): continue
        f=dict(re.findall(r'(\w+\|?\w*/?\w*)=(\S+)',line))
        rows.append((float(f['t']),float(f['E1']),float(f['M1'])))
    if rows:
        booked=[r[2] for r in rows]
        print("model books M1 at %d checkpoints: min %.6e max %.6e median %.6e W"
              %(len(booked),min(booked),max(booked),
                sorted(booked)[len(booked)//2]))
        print("E1 median %.6e W; M1/E1 median %.3e"%(
            sorted(r[1] for r in rows)[len(rows)//2],
            sorted(r[2]/r[1] for r in rows if r[1]>0)[len(rows)//2]))
        print("R135a ratio traced/booked (medians): %.4f"%(
            P/max(sorted(booked)[len(booked)//2],1e-300)))
    else:
        print("no CREM_M1_SHARE lines: the model booked no M1 at all")
