# The libration's own M1 power, with the second derivative taken at the
# LIBRATION scale rather than the substep scale (audit 135).  The substep
# stencil is pure integrator jitter; the stencil width is swept until the
# answer stops moving, which is rule R135e.
import re,sys,math
mu0=4.0e-7*math.pi; c=299792458.0; hbar=1.054571817e-34
muB=9.2740100657e-24
trace=sys.argv[1]; tmax=float(sys.argv[2]) if len(sys.argv)>2 else 2.5e-11
t=[];m=[];now=0.0
for line in open(trace):
    if not line.startswith('SPINPAIR'): continue
    f=dict(re.findall(r'(\w+)=(\S+)',line))
    vec=lambda k:[float(x) for x in f[k].strip('()').split(',')]
    now+=float(f['dt'])
    if now>tmax: break
    a,b=vec('mu1'),vec('mu2'); t.append(now); m.append([a[i]+b[i] for i in range(3)])
mag=[math.sqrt(sum(x*x for x in v)) for v in m]
print("window %.4e s, %d substeps, mean spacing %.4e s"%(t[-1],len(t),t[-1]/len(t)))
print("|mu1+mu2|/muB: min %.6f max %.6f mean %.6f, swing %.6f"%(
    min(mag)/muB,max(mag)/muB,sum(mag)/len(mag)/muB,(max(mag)-min(mag))/muB))
def power(sq): return mu0*sq/(6.0*math.pi*c**3)
print("%8s %14s %16s %14s"%("stride","span/T_lib","<|m''|^2>","P_M1 [W]"))
prev=None
for stride in (1,2,5,10,20,40,80,160):
    out=[]
    for i in range(stride,len(t)-stride):
        h1=t[i]-t[i-stride]; h2=t[i+stride]-t[i]
        if h1<=0 or h2<=0: continue
        d=[2.0*(h2*m[i-stride][k]-(h1+h2)*m[i][k]+h1*m[i+stride][k])
           /(h1*h2*(h1+h2)) for k in range(3)]
        out.append(sum(x*x for x in d))
    if not out: continue
    mean=sum(out)/len(out)
    span=2*stride*(t[-1]/len(t))/1.1044e-11
    print("%8d %14.5f %16.6e %14.6e%s"%(stride,span,mean,power(mean),
        "" if prev is None else "   ratio %.4f"%(mean/prev)))
    prev=mean
# The analytic libration estimate: a moment component of amplitude A
# oscillating at omega has |m''| = omega^2 A.
omega=2.0*math.pi*90.551e9
amp=0.5*(max(mag)-min(mag))
print("\nanalytic: omega_lib %.6e rad/s, amplitude %.6e J/T (%.4f muB)"%(
    omega,amp,amp/muB))
print("  |m''| = omega^2 A = %.6e -> P = %.6e W"%(
    omega*omega*amp,power((omega*omega*amp)**2)))
