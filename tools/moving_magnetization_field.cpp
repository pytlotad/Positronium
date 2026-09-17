// A uniformly moving magnetic moment's retarded field against the boosted
// rest-frame loop field, B and E, at beta 0, 0.05, 0.3 and 0.1-5 r* (audit
// sections 84b and 103).  e36bf5c failed this at 8% (beta 0.05) and 58%
// (beta 0.3); correct code agrees to ~5e-10.  Validation now carries the same
// comparison as dipole-moving-magnetization.
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/moving_magnetization_field.cpp -o /tmp/movingloop $(root-config --libs)
// Moving magnetic dipole in uniform motion: retarded field vs the boosted rest-frame loop field.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
int main(){
  const double eps=magneticDipoleRadius(), rs=comptonBarrierRadius, mu=firstMagneticMoment;
  for(double beta:{0.0,0.05,0.3}){
    const Vec3 v{beta*c,0,0}; const Vec3 mp{0.2*mu,-0.4*mu,0.8*mu};
    State now{}; now.firstPosition={1e4*rs,0,0}; now.secondPosition={0,0,0}; now.secondVelocity=v; now.secondProperDipole=mp; synchronizeCovariantDipoles(now);
    StateHistory h; const int H=64; const double span=40*rs/c;
    for(int j=H;j>=0;--j){ State p=now; p.time=-span*j/H; p.secondPosition=v*p.time; synchronizeCovariantDipoles(p); h.push_back(p); }
    const double g=1/std::sqrt(1-beta*beta);
    for(double rr:{0.1,1.0,2.0,5.0}){ Vec3 x=Vec3{0.5,0.7,-0.4}*(rr*rs);
      const auto f=retardedMagneticDipoleField(x,0.0,h,now,false);
      Vec3 xr=x+Vec3{1,0,0}*((g-1)*x.x);
      const Vec3 Bp=plummerMagneticDipoleField(xr,mp,eps), Bpoles=plummerDipoleField(xr,mp,eps);
      auto boost=[&](const Vec3& Br,Vec3& El,Vec3& Bl){ Bl=Br*g-v*(g*g/(g+1)*dot(v,Br)/(c*c)); El=cross(v,Br)*(-g); };
      Vec3 El,Bl,Elp,Blp; boost(Bp,El,Bl); boost(Bpoles,Elp,Blp);
      printf("beta %.2f r=%.2f r*: petla |dB|/|B| %.2e |dE|/|E| %.2e   (wzgledem bieguny: |dB|/|B| %.2e)\n",beta,rr,
        (f.magnetic-Bl).norm()/Bl.norm(), beta>0?(f.electric-El).norm()/El.norm():0.0, (f.magnetic-Blp).norm()/Blp.norm());
    }
  }
}
