import sys,math,glob,statistics as st
S=sys.argv[1]
def load(pat,ncol=4):
    d={}
    for f in sorted(glob.glob(S+"/"+pat)):
        for ln in open(f):
            t=ln.split()
            if len(t)==ncol and t[0].isdigit():
                d[int(t[0])]=(float(t[1]),float(t[2]),float(t[3]))
    return d
new=load("s75_*.txt"); old=load("seeds.txt")
dl={}
for ln in open(S+"/delta24.txt"):
    t=ln.split()
    if len(t)==9 and t[0].isdigit(): dl[int(t[0])]=float(t[8])
ks=sorted(new)
print("seeds present at s=0.075: %d of 24"%len(ks))
if len(ks)<24: print("MISSING:",[i for i in range(1,25) if i not in new])
def blk(d,lbl):
    r=[d[k][2] for k in sorted(d)]; p=[d[k][0] for k in sorted(d)]
    o=[d[k][1] for k in sorted(d)]
    n=len(r); m=st.mean(r); s=st.stdev(r); e=s/math.sqrt(n)
    print("%-16s n=%2d mean %+.4e sd %.4e sem %.4e t=%+.2f  pos %d/%d"
          %(lbl,n,m,s,e,m/e,sum(1 for x in r if x>0),n))
    print("%-16s para range/mean %.4e   ortho %.4e"
          %("",(max(p)-min(p))/st.mean(p),(max(o)-min(o))/st.mean(o)))
    return m,e
mn,en=blk(new,"s=0.075"); mo,eo=blk(old,"s=0.30")
print()
lo,hi=mo-2*eo,mo+2*eo
print("OUTCOME C test: 213's 2-sem interval [%.4e, %.4e]"%(lo,hi))
print("  s=0.075 mean %.4e -> %s"%(mn,"OUTSIDE (C holds)" if not lo<=mn<=hi
      else "inside (C fails, 213's bar was honest)"))
def spearman(x,y):
    n=len(x)
    def rk(v):
        o=sorted(range(n),key=lambda i:v[i]); r=[0.0]*n
        i=0
        while i<n:
            j=i
            while j+1<n and v[o[j+1]]==v[o[i]]: j+=1
            for k in range(i,j+1): r[o[k]]=(i+j)/2.0+1
            i=j+1
        return r
    a,b=rk(x),rk(y); ma,mb=st.mean(a),st.mean(b)
    num=sum((a[i]-ma)*(b[i]-mb) for i in range(n))
    den=math.sqrt(sum((a[i]-ma)**2 for i in range(n))
                 *sum((b[i]-mb)**2 for i in range(n)))
    return num/den if den>0 else 0.0
print()
for lbl,sel in (("all 24",lambda d:True),
                ("delta > 90",lambda d:d>90.0),
                ("delta < 90",lambda d:d<90.0)):
    xs=[dl[k] for k in ks if sel(dl[k])]; ys=[new[k][2] for k in ks if sel(dl[k])]
    xo=[dl[k] for k in sorted(old) if sel(dl[k])]
    yo=[old[k][2] for k in sorted(old) if sel(dl[k])]
    print("Spearman rho(ratio, delta)  %-11s s=0.075 %+.3f (n=%2d)"
          "   s=0.30 %+.3f (n=%2d)"%(lbl,spearman(xs,ys),len(xs),
          spearman(xo,yo),len(xo)))
print()
print("%5s %8s %13s %13s %8s"%("seed","delta","s=0.30","s=0.075","ratio"))
for k in sorted(ks,key=lambda k:dl[k]):
    print("%5d %8.2f %13.4e %13.4e %8.3f"
          %(k,dl[k],old[k][2],new[k][2],
            new[k][2]/old[k][2] if old[k][2]!=0 else float('nan')))
print()
print("para column at s=0.075, sorted, 1e-11 s:")
for v,k in sorted((new[k][0],k) for k in ks):
    print("   seed %2d  %.6f"%(k,v*1e11))
