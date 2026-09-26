// Audit 223: the closed form t = hbar^6 c^3 / (4 mu k^5) is a
// statement about STRUCTURE, so it should survive changing the one
// thing the model lets us change without touching the coupling: the
// reduced mass.  alpha, hbar and c are fixed; mu moves by a factor
// 1836 across the three same-charge-opposite-sign pairs the model
// supports.  Nobody wrote t ~ 1/mu anywhere.
//
// Ratios are reported as well as absolutes, because the step bias of
// s_max = 0.30 is multiplicative (audit 213f) and cancels in a ratio.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    struct Case { const char* name; ParticlePair pair; };
    const Case cases[]={
        {"positronium",{electron,positron}},
        {"true muonium",{muon,antimuon}},
        {"protonium",{proton,antiproton}},
    };
    std::printf("%14s %12s %12s %14s %14s %9s\n",
                "pair","mu/m_e","a_pair [fm]","closed [ps]",
                "measured [ps]","meas/clos");
    double firstClosed=0.0,firstMeasured=0.0;
    for(const Case& cs:cases){
        applyPair(cs.pair);
        gSpinQuantization=true;
        const double mu=firstMass*secondMass/(firstMass+secondMass);
        const double k=pairCoulombStrength;
        const double a=pairBohrRadius(cs.pair);
        const double closed=std::pow(hbar,6.0)*c*c*c/(4.0*mu*std::pow(k,5.0));
        const CremCollapseEstimate e=estimateCremCollapse(1,1,900.0);
        const bool ok=e.calibrationOutcome==SimulationOutcome::ReachedCutoff
                     &&std::isfinite(e.lifetimeSeconds);
        std::printf("%14s %12.4f %12.3f %14.6f %14s %9s\n",
            cs.name,mu/electron.mass,a*1e15,closed*1e12,
            ok?"":"censored","");
        if(ok)
            std::printf("%14s %12s %12s %14s %14.6f %9.4f\n",
                "","","","",e.lifetimeSeconds*1e12,
                e.lifetimeSeconds/closed);
        if(firstClosed==0.0&&ok){ firstClosed=closed;
                                  firstMeasured=e.lifetimeSeconds; }
        else if(ok)
            std::printf("%14s ratio to positronium: closed %10.4f"
                        "  measured %10.4f\n","",
                        firstClosed/closed,firstMeasured/e.lifetimeSeconds);
    }
}
