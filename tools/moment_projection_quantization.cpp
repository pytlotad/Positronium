// Is the projection of a particle's spin magnetic moment on the axis set by
// the PARTNER's net field quantized (audit section 114)?
//
// The axis is the model's own: thomasBmtEffectiveField of the partner's
// instantaneous field at the particle, which is what the Thomas-BMT rotation
// precesses the moment about, and therefore the axis whose projection that
// rotation leaves invariant.  This probe prepares N production states (no
// integration) and histograms cos(mu1, axis) and the spin projection it
// implies in units of hbar.  The trajectory side of the same question is the
// CREM_DEBUG_PROJECTION trace in secular_spin_orbit.hpp.
//
// Usage: moment_projection_quantization <count> [seed base, default 1]
//        [phenomenon, default 1] [quant to impose the MUTUAL angle]
// The last one is the distinction the question turns on: --spin-quantization
// fixes the angle BETWEEN the two moments, not either one's projection on the
// field axis.
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/moment_projection_quantization.cpp -o /tmp/proj $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc,char** argv) {
    const int count=argc>1?atoi(argv[1]):400;
    const unsigned long long base=argc>2?strtoull(argv[2],nullptr,10):1ULL;
    const int phenomenon=argc>3?atoi(argv[3]):1;
    if(argc>4&&std::string(argv[4])=="quant") gSpinQuantization=true;
    std::vector<double> projections,fieldProjections;
    int nearPoles=0;
    for(int index=0;index<count;++index) {
        SimulationOptions options;
        options.frameCount=2;
        options.observationTime=1.0e-24;
        const SimulationResult prepared=simulate(
            splitMix64(base+static_cast<unsigned long long>(index)),
            phenomenon,options);
        if(prepared.frames.empty()) continue;
        // The sampled frame carries the positions and the drawn moments; the
        // velocities are the sharp circular preparation the sampler makes
        // (orbit in the xy plane, separation along x), rebuilt here rather
        // than re-derived from the frame, which does not store them.
        const Frame& frame=prepared.frames.front();
        State state{};
        state.firstPosition=frame.first;
        state.secondPosition=frame.second;
        const Vec3 separationVector=state.firstPosition-state.secondPosition;
        const double radius=separationVector.norm();
        if(!(radius>0.0)) continue;
        const Vec3 radial=separationVector*(1.0/radius);
        const Vec3 tangential=cross(Vec3{0.0,0.0,1.0},radial);
        const double relativeSpeed=
            std::sqrt(pairCoulombStrength/(pairReducedMass*radius));
        const double total=firstMass+secondMass;
        state.firstVelocity=tangential*(relativeSpeed*secondMass/total);
        state.secondVelocity=tangential*(-relativeSpeed*firstMass/total);
        state.firstDipole=frame.firstDipole;
        state.secondDipole=frame.secondDipole;
        state.firstProperDipole=frame.firstDipole;
        state.secondProperDipole=frame.secondDipole;
        synchronizeCovariantDipoles(state);
        const StateHistory history=causalInitialHistory(state);
        const ElectromagneticField partnerField=fieldFromOtherParticleAt(
            state.firstPosition,state.time,state,history,true);
        const Vec3 axis=thomasBmtEffectiveField(state.firstVelocity,
            partnerField,firstSpecies.gFactor);
        const Vec3 moment=state.firstDipole;
        const double axisNorm=axis.norm();
        const double momentNorm=moment.norm();
        if(!(axisNorm>0.0)||!(momentNorm>0.0)) continue;
        const double cosine=dot(moment,axis)/(axisNorm*momentNorm);
        projections.push_back(cosine);
        const double magneticNorm=partnerField.magnetic.norm();
        if(magneticNorm>0.0)
            fieldProjections.push_back(
                dot(moment,partnerField.magnetic)/(magneticNorm*momentNorm));
        if(std::abs(std::abs(cosine)-1.0)<=0.05) ++nearPoles;
    }
    const int total=static_cast<int>(projections.size());
    if(total==0) { std::printf("no prepared states\n"); return 1; }
    std::vector<double> sorted=projections;
    std::sort(sorted.begin(),sorted.end());
    double sum=0.0;
    for(double value:projections) sum+=value;
    std::printf("prepared %d states (base seed %llu, phenomenon %d)\n",
                total,base,phenomenon);
    std::printf("cos(mu1, BMT axis of the partner's field): mean %+.4f "
                "median %+.4f min %+.4f max %+.4f\n",
                sum/total,sorted[static_cast<size_t>(total/2)],
                sorted.front(),sorted.back());
    std::printf("spin projection S1.axis: mean %+.4f hbar, i.e. no pile-up at "
                "+-0.5 unless the next line says so\n",0.5*sum/total);
    std::printf("within 0.05 of +-1 (a quantized axis would be ~100%%): "
                "%d/%d = %.1f%%  (a uniform cos gives 5%%)\n",
                nearPoles,total,100.0*nearPoles/total);
    // Decile occupancy: a continuous measure fills them evenly.
    std::printf("deciles of cos:");
    for(int bin=0;bin<10;++bin) {
        int inBin=0;
        const double low=-1.0+0.2*bin,high=low+0.2;
        for(double value:projections)
            if(value>=low&&(value<high||(bin==9&&value<=high))) ++inBin;
        std::printf(" %d",inBin);
    }
    std::printf("\n");
    if(!fieldProjections.empty()) {
        double fieldSum=0.0;
        for(double value:fieldProjections) fieldSum+=value;
        std::printf("for reference, cos(mu1, B_partner) alone: mean %+.4f "
                    "over %zu states\n",
                    fieldSum/fieldProjections.size(),fieldProjections.size());
    }
}
