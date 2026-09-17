// Is the 2 gamma / 3 gamma split the same division as p-Ps / o-Ps (audit 111)?
//
// The multiplicity is decided by the two-photon weight
// w = |mu1+mu2|^2/(|mu1|+|mu2|)^2 = (1+cos)/2 of the PREPARED moments, drawn
// as a branching ratio (crem_collapse.hpp).  The channel LABEL is the knife
// edge cos >= 0.5 of the same moments.  This probe prepares N states exactly
// as production does -- no integration, the sampled state only -- and reports
// the joint distribution of the label and the drawn multiplicity, plus the
// weight's own distribution.
//
// Usage: annihilation_channel_split <count> [seed base, default 1]
//        [phenomenon, default 1] [quant to impose spin quantization]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/annihilation_channel_split.cpp -o /tmp/split $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc,char** argv) {
    const int count=argc>1?atoi(argv[1]):1000;
    const unsigned long long base=argc>2?strtoull(argv[2],nullptr,10):1ULL;
    const int phenomenon=argc>3?atoi(argv[3]):1;
    // Audit 91 withdrew the imposed mutual angle; with it back on the
    // prepared cosine is exactly +-1 and the weight is exactly 1 or 0.
    if(argc>4&&std::string(argv[4])=="quant") gSpinQuantization=true;
    std::vector<double> cosines,weights;
    int labelPara=0,twoPhoton=0,twoPhotonInPara=0,twoPhotonInOrtho=0;
    int band=0,bandTwoPhoton=0;
    std::uint64_t stream=splitMix64(base^0x9e3779b97f4a7c15ULL);
    for(int index=0;index<count;++index) {
        SimulationOptions options;
        options.frameCount=2;
        options.observationTime=1.0e-24;
        const SimulationResult prepared=simulate(
            splitMix64(base+static_cast<unsigned long long>(index)),
            phenomenon,options);
        if(prepared.frames.empty()) continue;
        const Vec3 first=prepared.frames.front().firstDipole;
        const Vec3 second=prepared.frames.front().secondDipole;
        const double scale=first.norm()+second.norm();
        if(!(scale>0.0)) continue;
        const double cosine=dot(first,second)/(first.norm()*second.norm());
        const double weight=(first+second).squaredNorm()/(scale*scale);
        cosines.push_back(cosine);
        weights.push_back(weight);
        const bool para=cosine>=0.5;
        const bool two=drawUniformUnit(stream)<weight;
        if(para) ++labelPara;
        if(two) ++twoPhoton;
        if(para&&two) ++twoPhotonInPara;
        if(!para&&two) ++twoPhotonInOrtho;
        if(std::abs(cosine)<0.5) { ++band; if(two) ++bandTwoPhoton; }
    }
    const int total=static_cast<int>(weights.size());
    if(total==0) { std::printf("no prepared states\n"); return 1; }
    const auto mean=[](const std::vector<double>& values) {
        double sum=0.0;
        for(double value:values) sum+=value;
        return sum/static_cast<double>(values.size());
    };
    std::vector<double> sorted=cosines;
    std::sort(sorted.begin(),sorted.end());
    const int labelOrtho=total-labelPara;
    std::printf("prepared %d states (phenomenon %d, base seed %llu)\n",
                total,phenomenon,base);
    std::printf("cos: mean %+.4f median %+.4f min %+.4f max %+.4f\n",
                mean(cosines),sorted[static_cast<size_t>(total/2)],
                sorted.front(),sorted.back());
    std::printf("two-photon weight: mean %.4f (= P(2 gamma) over the ensemble)\n",
                mean(weights));
    std::printf("label cos>=0.5 (p-Ps): %d/%d = %.1f%%; below (o-Ps): %d = %.1f%%\n",
                labelPara,total,100.0*labelPara/total,labelOrtho,
                100.0*labelOrtho/total);
    std::printf("drawn 2 gamma: %d/%d = %.1f%%\n",
                twoPhoton,total,100.0*twoPhoton/total);
    if(labelPara>0)
        std::printf("  P(2 gamma | labelled p-Ps) = %d/%d = %.3f\n",
                    twoPhotonInPara,labelPara,
                    static_cast<double>(twoPhotonInPara)/labelPara);
    if(labelOrtho>0)
        std::printf("  P(2 gamma | labelled o-Ps) = %d/%d = %.3f\n",
                    twoPhotonInOrtho,labelOrtho,
                    static_cast<double>(twoPhotonInOrtho)/labelOrtho);
    if(band>0)
        std::printf("  |cos| < 0.5 (neither label is a band edge): %d = %.1f%%, "
                    "of which 2 gamma %.3f\n",band,100.0*band/total,
                    static_cast<double>(bandTwoPhoton)/band);
    std::printf("spin statistics for comparison: P(2 gamma) = 1/4 = 0.250\n");
}
