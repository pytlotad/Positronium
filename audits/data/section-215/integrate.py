import sys,math
S=sys.argv[1]; A=1.05835442e-10; BARRIER=1.826454e-03
def load(f):
    rows=[]
    for ln in open(S+"/"+f):
        t=ln.split()
        if len(t)!=9 or t[0].startswith("#"): continue
        try: rows.append(tuple(float(x) for x in
                              (t[0],t[1],t[2],t[5],t[6])))
        except ValueError: pass
    return sorted(rows)
def mk(xs,ys):
    def g(x):
        i=0
        if x<=xs[0]: i=0
        elif x>=xs[-1]: i=len(xs)-2
        else:
            while xs[i+1]<x: i+=1
        u=(x-xs[i])/(xs[i+1]-xs[i]); return ys[i]+u*(ys[i+1]-ys[i])
    return g
for f,lab,meas in (("f60.txt","delta = 60 deg",+1.13e-03),
                   ("f160.txt","delta = 160 deg",-5.0e-05)):
    rows=load(f); lr=[math.log(r) for r,_,_,_,_ in rows]
    print("\n=== %s ==="%lab)
    print("  %9s %13s %13s"%("r/a_pair","dP/P","(Eo-Ep)/|Ep|"))
    for r,p,o,ep,eo in rows:
        print("  %9.5f %13.4e %13.4e"%(r,(o-p)/p,(eo-ep)/abs(ep)))
    res={}
    for ch,(pi,ei) in enumerate(((1,3),(2,4))):
        lp=mk(lr,[math.log(x[pi]) for x in rows])
        ee=mk(lr,[x[ei] for x in rows])
        lo,hi=max(rows[0][0],BARRIER),rows[-1][0]; N=8000; t=0.0; cum=[]
        for j in range(N-1,-1,-1):
            r0=lo*(hi/lo)**(j/N); r1=lo*(hi/lo)**((j+1)/N)
            dE=ee(math.log(r1))-ee(math.log(r0))
            P=math.exp(lp(math.log(math.sqrt(r0*r1))))
            t+=dE/P; cum.append((r0,t))
        res[ch]=(t,cum)
    tp,cp=res[0]; to,co=res[1]; d=to-tp
    print("  t_para %.6e s   t_ortho %.6e s   (estimator: 3.10e-11 s)"%(tp,to))
    print("  predicted ortho/para - 1  %+.4e"%(to/tp-1.0))
    print("  measured global (214e)    %+.4e"%meas)
    print("  predicted / measured       %.2f"%((to/tp-1.0)/meas))
    print("  %9s %11s %14s"%("r/a_pair","cum t [%]","cum (to-tp) [%]"))
    for tg in (0.7,0.5,0.35,0.25,0.15,0.1,0.05,0.025,0.01,BARRIER):
        x=min(cp,key=lambda c:abs(c[0]-tg)); y=min(co,key=lambda c:abs(c[0]-tg))
        print("  %9.5f %11.4f %14.4f"
              %(x[0],100*x[1]/tp,100*(y[1]-x[1])/d if d else 0.0))
