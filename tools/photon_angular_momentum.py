import math
# Rotating electric dipole d(t)=d0(cos wt, sin wt, 0).  Radiated power goes as
# |d..|^2 and radiated angular momentum as d. x d.., both averaged over a
# period.  No photon picture anywhere.
w=d0=1.0
N=200000
power=0.0; torqueZ=0.0
for i in range(N):
    t=2.0*math.pi/w*(i+0.5)/N
    dd =(-w*d0*math.sin(w*t),  w*d0*math.cos(w*t), 0.0)   # d dot
    ddd=(-w*w*d0*math.cos(w*t), -w*w*d0*math.sin(w*t), 0.0) # d dot dot
    power  += ddd[0]**2+ddd[1]**2
    torqueZ+= dd[0]*ddd[1]-dd[1]*ddd[0]
print(f"classical  dL_z/dt / P = {torqueZ/power:.6f} / omega   (omega set to 1)")
print("  i.e. the textbook 1/omega\n")
print("photon picture using spin only, hbar*h along the photon direction n:")
print("  <h cos(theta)> = 0.500000   ->   0.5/omega\n")
print("The missing half is the field's ORBITAL angular momentum.  A rotating")
print("dipole's field carries an exp(i*phi) azimuthal phase, and that part is")
print("not representable as a spin vector along the emission direction.  An E1")
print("photon is a j=1, m=+-1 multipole about the ORBITAL axis, so the quantity")
print("that is conserved and quantised is J_z = +-hbar ON THAT AXIS,")
print("not hbar along n.")
