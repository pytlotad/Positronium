# The mutual-moment measure read straight off the section 129 trace: no
# Boltzmann assumption, no temperature, just the time distribution of
# cos(mu1,mu2) that the model's own dynamics produces (audit 130e).
import re,math
rows=[]
for line in open('trace.txt'):
    if not line.startswith('SPINPAIR'): continue
    f=dict(re.findall(r'(\w+)=(\S+)',line))
    vec=lambda k:[float(x) for x in f[k].strip('()').split(',')]
    rows.append((float(f['dt']),vec('mu1'),vec('mu2'),
                 float(f['ratio1']),float(f['ratio2'])))
def cosang(a,b):
    na=math.sqrt(sum(x*x for x in a)); nb=math.sqrt(sum(x*x for x in b))
    return sum(x*y for x,y in zip(a,b))/(na*nb)
tot=sum(r[0] for r in rows)
cm=[cosang(r[1],r[2]) for r in rows]
mean_t=sum(r[0]*c for r,c in zip(rows,cm))/tot
print("substeps",len(rows),"ratio1",rows[0][3],"ratio2",rows[0][4])
print("cos(mu1,mu2): start %+.4f end %+.4f min %+.4f max %+.4f"%(
    cm[0],cm[-1],min(cm),max(cm)))
print("time mean %+.4f, substep mean %+.4f"%(mean_t,sum(cm)/len(cm)))
print("E[w] from this trajectory: time %.4f, substep %.4f"%(
    (1+mean_t)/2,(1+sum(cm)/len(cm))/2))
acc=0.0; tail=[]
for r,c in zip(reversed(rows),reversed(cm)):
    acc+=r[0]; tail.append((r[0],c))
    if acc>0.10*tot: break
wt=sum(d*c for d,c in tail)/sum(d for d,c in tail)
print("last 10%% of the time: <cos_mu> %+.4f -> E[w] %.4f"%(wt,(1+wt)/2))
