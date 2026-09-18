# |S1+S2| drift, rate asymmetry and mutual-angle motion from a PROJ trace,
# the same quantities the tables of audit 116b and 117d were built from.
import re,sys,glob
print("%-26s %8s %9s %9s %9s %9s %11s %9s %9s"%("run","substeps",
      "S start","S end","drift","range","|w2|/|w1|","cos start","cos end"))
for f in sorted(glob.glob(sys.argv[1])):
    s=[];w1=[];w2=[];cos=[]
    for line in open(f):
        if not line.startswith('PROJ'): continue
        d=dict(re.findall(r'(\S+)=([-\d.e+]+)',line))
        s.append(float(d['Spair_hbar'])); cos.append(float(d['cos12']))
        w1.append(float(d['|w1|'])); w2.append(float(d['|w2|']))
    if not s: print("%-26s no PROJ lines"%f); continue
    ratios=[b/a if a>0 else float('nan') for a,b in zip(w1,w2)]
    ratios.sort()
    print("%-26s %8d %9.6f %9.6f %9.2e %9.2e %11.4f %9.6f %9.6f"%(
        f.replace('_err.txt',''),len(s),s[0],s[-1],abs(s[-1]-s[0]),
        max(s)-min(s),ratios[len(ratios)//2],cos[0],cos[-1]))
