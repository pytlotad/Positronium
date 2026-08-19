#pragma once

// Optional finite-difference Maxwell backend: Yee grid, AMR, CPML and
// extended charge/current deposition. Included only by the validation build.

// Smooth finite-size source used by the Maxwell-field backend.  The Gaussian
// is spherical in the instantaneous rest frame.  In the laboratory it is
// Lorentz-contracted along the velocity and multiplied by gamma, so its
// spatial integral remains exactly q.  The radius is a numerical UV cutoff,
// not a claim about the physical size of an first.
struct RelativisticChargeCloud {
    double charge = 0.0;
    double restRadius = chargeCloudRestRadius;

    double density(const Vec3& fieldPosition, const Vec3& centre,
                   const Vec3& velocity) const {
        const Vec3 offset = fieldPosition - centre;
        const double speedSquared = velocity.squaredNorm();
        const double lorentzGamma = gamma(velocity);
        double restFrameRadiusSquared = offset.squaredNorm();
        if (speedSquared > 0.0) {
            const double parallelProjection = dot(offset, velocity);
            restFrameRadiusSquared += lorentzGamma*lorentzGamma
                * parallelProjection*parallelProjection / (c*c);
        }
        const double normalization = charge * lorentzGamma
            / (std::pow(2.0*pi, 1.5) * restRadius*restRadius*restRadius);
        return normalization
            * std::exp(-0.5 * restFrameRadiusSquared/(restRadius*restRadius));
    }

    Vec3 currentDensity(const Vec3& fieldPosition, const Vec3& centre,
                        const Vec3& velocity) const {
        return velocity * density(fieldPosition, centre, velocity);
    }

    double supportRadius(const Vec3& velocity) const {
        // At four Gaussian widths the omitted total charge is below 0.12%.
        // AMR deposition renormalizes the retained weights to exactly q.
        (void)velocity;
        return 4.0 * restRadius;
    }
};

struct RankTwoTensor {
    std::array<std::array<double,4>,4> component{};
};

struct AntisymmetricTensor {
    std::array<std::array<double,4>,4> component{};
};

RankTwoTensor lorentzBoostMatrix(const Vec3& velocity) {
    RankTwoTensor boost;
    const double g=gamma(velocity); const Vec3 beta=velocity/c;
    boost.component[0][0]=g;
    for(int i=0;i<3;++i) {
        const double bi=i==0?beta.x:(i==1?beta.y:beta.z);
        boost.component[0][i+1]=g*bi; boost.component[i+1][0]=g*bi;
        for(int j=0;j<3;++j) {
            const double bj=j==0?beta.x:(j==1?beta.y:beta.z);
            boost.component[i+1][j+1]=(i==j?1.0:0.0)
                +(g-1.0)*bi*bj/std::max(beta.squaredNorm(),1.0e-300);
        }
    }
    return boost;
}

FourVector transform(const RankTwoTensor& matrix,const FourVector& vector) {
    const std::array<double,4> input{vector.time,vector.space.x,vector.space.y,vector.space.z};
    std::array<double,4> output{};
    for(int i=0;i<4;++i) for(int j=0;j<4;++j)
        output[i]+=matrix.component[i][j]*input[j];
    return {output[0],{output[1],output[2],output[3]}};
}

RankTwoTensor transform(const RankTwoTensor& boost,const RankTwoTensor& tensor) {
    RankTwoTensor result;
    for(int m=0;m<4;++m) for(int n=0;n<4;++n)
        for(int a=0;a<4;++a) for(int b=0;b<4;++b)
            result.component[m][n]+=boost.component[m][a]*boost.component[n][b]
                                   *tensor.component[a][b];
    return result;
}

AntisymmetricTensor transform(const RankTwoTensor& boost,
                              const AntisymmetricTensor& tensor) {
    AntisymmetricTensor result;
    for(int m=0;m<4;++m) for(int n=0;n<4;++n)
        for(int a=0;a<4;++a) for(int b=0;b<4;++b)
            result.component[m][n]+=boost.component[m][a]*boost.component[n][b]
                                   *tensor.component[a][b];
    return result;
}

// One covariant state for a smooth classical body.  rho0, M^{mu nu}, and
// T_matter^{mu nu} are defined in the instantaneous rest frame and boosted
// together.  pressure is the explicit cohesive (Poincare) stress required by
// this phenomenological rigid Gaussian model; it is not silently included in
// the measured particle mass.
struct CovariantExtendedBody {
    double charge=0.0, physicalMass=0.0, restRadius=chargeCloudRestRadius;
    Vec3 centre, velocity, magneticDipole;

    double electromagneticSelfEnergy() const {
        // Gaussian rho=q exp(-r^2/2 sigma^2)/(2 pi sigma^2)^(3/2).
        return charge*charge/(8.0*std::pow(pi,1.5)*epsilon0*restRadius);
    }
    double electromagneticMass() const { return electromagneticSelfEnergy()/(c*c); }
    double bareMatterMass() const { return physicalMass-electromagneticMass(); }

    double properShape(const Vec3& restOffset) const {
        return std::exp(-0.5*restOffset.squaredNorm()/(restRadius*restRadius))
            /(std::pow(2.0*pi,1.5)*restRadius*restRadius*restRadius);
    }
    double restElectricFieldMagnitude(double radius) const {
        if(radius<1.0e-8*restRadius) return std::abs(charge)*radius
            /(3.0*epsilon0*std::pow(2.0*pi,1.5)*std::pow(restRadius,3));
        const double x=radius/(std::sqrt(2.0)*restRadius);
        const double enclosedFraction=std::erf(x)
            -std::sqrt(2.0/pi)*(radius/restRadius)
             *std::exp(-0.5*radius*radius/(restRadius*restRadius));
        return std::abs(charge)*enclosedFraction
            /(4.0*pi*epsilon0*radius*radius);
    }
    double cohesivePressure(const Vec3& restOffset) const {
        // Static balance: grad(p)=rho E_force with T^{ij}=p delta^{ij}
        // under the (+---) convention used below. Integrate the outward
        // electrostatic force density from r to the Gaussian tail.
        const double lower=restOffset.norm(),upper=8.0*restRadius;
        if(lower>=upper) return 0.0;
        constexpr int intervals=128;
        const double h=(upper-lower)/intervals;
        const auto densityForce=[&](double radius) {
            const double rho=std::abs(charge)*std::exp(
                -0.5*radius*radius/(restRadius*restRadius))
                /(std::pow(2.0*pi,1.5)*std::pow(restRadius,3));
            return rho*restElectricFieldMagnitude(radius);
        };
        double integral=densityForce(lower)+densityForce(upper);
        for(int i=1;i<intervals;++i)
            integral+=(i%2?4.0:2.0)*densityForce(lower+i*h);
        return integral*h/3.0;
    }
    FourVector freeFourCurrent(const Vec3& restOffset) const {
        return transform(lorentzBoostMatrix(velocity),
                         FourVector{c*charge*properShape(restOffset),{}});
    }
    AntisymmetricTensor polarizationMagnetization(const Vec3& restOffset) const {
        AntisymmetricTensor rest;
        const Vec3 m=magneticDipole*properShape(restOffset);
        // Spatial components M^{ij}=epsilon^{ijk} M_k.  Boosting this tensor
        // automatically generates the electric-polarization components M^{0i}.
        rest.component[1][2]=m.z; rest.component[2][1]=-m.z;
        rest.component[2][3]=m.x; rest.component[3][2]=-m.x;
        rest.component[3][1]=m.y; rest.component[1][3]=-m.y;
        return transform(lorentzBoostMatrix(velocity),rest);
    }
    RankTwoTensor matterStressEnergy(const Vec3& restOffset) const {
        RankTwoTensor rest;
        rest.component[0][0]=physicalMass*c*c*properShape(restOffset);
        const double pressure=cohesivePressure(restOffset);
        for(int i=1;i<4;++i) rest.component[i][i]=pressure;
        return transform(lorentzBoostMatrix(velocity),rest);
    }
};

struct MaxwellCell {
    Vec3 electric, magnetic;
    double chargeDensity = 0.0;
    Vec3 currentDensity;
};
struct MaxwellBoundaryFlux {
    double energyRate=0.0;
    Vec3 momentumRate, angularMomentumRate;
};
struct MaxwellVolumeIntegrals {
    double energy=0.0;
    Vec3 momentum, angularMomentum;
};
struct GridDipoleInteraction {
    Vec3 force, torque;
    double fieldPower=0.0;
};
struct ChargeConservingDepositResult {
    double continuityResidual=0.0;
    double longitudinalCurl=0.0;
};

// Finite-volume Yee layout.  Charge lives at cell centres, normal electric
// fields and currents on faces, and magnetic fields on the complementary
// staggered faces.  The trajectory-decomposition deposition below is the
// local charge-conserving construction underlying Esirkepov-type schemes.
class YeeMaxwellBlock {
public:
    YeeMaxwellBlock(int cellsPerAxis,double cellSize,Vec3 centre)
        : n_(cellsPerAxis),dx_(cellSize),origin_(centre-Vec3{
            0.5*cellsPerAxis*cellSize,0.5*cellsPerAxis*cellSize,
            0.5*cellsPerAxis*cellSize}),
          rho_(size()),rhoPrevious_(size()),ex_(size()),ey_(size()),ez_(size()),
          bx_(size()),by_(size()),bz_(size()),jx_(size()),jy_(size()),jz_(size()),
          px_(size()),py_(size()),pz_(size()),pxPrevious_(size()),
          pyPrevious_(size()),pzPrevious_(size()),mx_(size()),my_(size()),mz_(size()) {
        if(n_<6||!(dx_>0.0)) throw std::invalid_argument("invalid Yee block");
        for(auto& memory:electricPmlMemory_) memory.assign(size(),0.0);
        for(auto& memory:magneticPmlMemory_) memory.assign(size(),0.0);
    }
    double courantTimeStep() const { return 0.95*dx_/(c*std::sqrt(3.0)); }
    void clearSources() {
        rhoPrevious_=rho_; std::fill(rho_.begin(),rho_.end(),0.0);
        pxPrevious_=px_; pyPrevious_=py_; pzPrevious_=pz_;
        std::fill(jx_.begin(),jx_.end(),0.0);
        std::fill(jy_.begin(),jy_.end(),0.0);
        std::fill(jz_.begin(),jz_.end(),0.0);
        std::fill(px_.begin(),px_.end(),0.0);
        std::fill(py_.begin(),py_.end(),0.0);
        std::fill(pz_.begin(),pz_.end(),0.0);
        std::fill(mx_.begin(),mx_.end(),0.0);
        std::fill(my_.begin(),my_.end(),0.0);
        std::fill(mz_.begin(),mz_.end(),0.0);
    }
    void depositInitialCic(double charge,const Vec3& position) {
        depositCharge(charge,position,rho_);
    }
    void depositEsirkepovCic(double charge,const Vec3& before,
                             const Vec3& after,double dt) {
        if(!(dt>0.0)) throw std::invalid_argument("invalid Esirkepov step");
        depositCharge(charge,after,rho_);
        Vec3 cursor=before;
        // Symmetric trajectory decomposition removes the leading axis-order
        // bias while every one-dimensional leg conserves charge locally.
        const std::array<std::pair<int,double>,5> legs{{
            {0,0.5*(after.x-before.x)},{1,0.5*(after.y-before.y)},
            {2,after.z-before.z},{1,0.5*(after.y-before.y)},
            {0,0.5*(after.x-before.x)}}};
        for(const auto [axis,delta]:legs) {
            Vec3 target=cursor;
            if(axis==0) target.x+=delta;
            else if(axis==1) target.y+=delta;
            else target.z+=delta;
            depositAxisLeg(charge,cursor,target,axis,dt);
            cursor=target;
        }
    }
    // A smooth relativistic cloud is represented by a fixed, normalized set
    // of material quadrature points.  Depositing every point through the same
    // local trajectory operator retains exact discrete continuity, including
    // the deformation current caused by a changing Lorentz contraction.
    void depositInitialGaussian(double charge,double restRadius,
                                const Vec3& centre,const Vec3& velocity) {
        forEachGaussianMarker(restRadius,centre,velocity,
            [&](double weight,const Vec3& position) {
                depositInitialCic(charge*weight,position);
            });
    }
    void depositGaussianEsirkepov(double charge,double restRadius,
                                  const Vec3& beforeCentre,
                                  const Vec3& beforeVelocity,
                                  const Vec3& afterCentre,
                                  const Vec3& afterVelocity,double dt) {
        std::vector<std::pair<double,Vec3>> before,after;
        forEachGaussianMarker(restRadius,beforeCentre,beforeVelocity,
            [&](double weight,const Vec3& position) {
                before.emplace_back(weight,position);
            });
        forEachGaussianMarker(restRadius,afterCentre,afterVelocity,
            [&](double weight,const Vec3& position) {
                after.emplace_back(weight,position);
            });
        for(std::size_t marker=0;marker<before.size();++marker)
            depositEsirkepovCic(charge*before[marker].first,
                                before[marker].second,after[marker].second,dt);
    }
    // P is stored on electric/current faces. M is stored on complementary
    // edges. Consequently rho_b=-div(P) and J_b=dP/dt+curl(M) obey the same
    // finite-volume continuity equation algebraically, without projection.
    void depositCovariantDipoleYee(const Vec3& centre,const Vec3& velocity,
                                   const Vec3& restDipole,double restRadius) {
        const double g=gamma(velocity);
        const Vec3 beta=velocity/c;
        const Vec3 labM=restDipole*g
            -beta*(dot(beta,restDipole)*g*g/(g+1.0));
        const Vec3 labP=cross(velocity,labM)/(c*c);
        forEachGaussianMarker(restRadius,centre,velocity,
            [&](double weight,const Vec3& position) {
                depositStaggered(px_,labP.x*weight,position,{0.0,0.5,0.5});
                depositStaggered(py_,labP.y*weight,position,{0.5,0.0,0.5});
                depositStaggered(pz_,labP.z*weight,position,{0.5,0.5,0.0});
                depositStaggered(mx_,labM.x*weight,position,{0.5,0.0,0.0});
                depositStaggered(my_,labM.y*weight,position,{0.0,0.5,0.0});
                depositStaggered(mz_,labM.z*weight,position,{0.0,0.0,0.5});
            });
    }
    void finalizeInitialBoundSources() { addBoundSources(0.0,false); }
    void finalizeBoundSources(double dt) {
        if(!(dt>0.0)) throw std::invalid_argument("invalid bound-current step");
        addBoundSources(dt,true);
    }
    double totalCharge() const {
        double total=0.0;
        for(double density:rho_) total+=density*dx_*dx_*dx_;
        return total;
    }
    double maximumContinuityResidual(double dt) const {
        double maximum=0.0;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const double divergence=(jx_[index(i+1,j,k)]-jx_[index(i,j,k)]
                +jy_[index(i,j+1,k)]-jy_[index(i,j,k)]
                +jz_[index(i,j,k+1)]-jz_[index(i,j,k)])/dx_;
            maximum=std::max(maximum,std::abs(
                (rho_[index(i,j,k)]-rhoPrevious_[index(i,j,k)])/dt+divergence));
        }
        return maximum;
    }
    double continuityScale(double dt) const {
        double scale=1.0;
        for(std::size_t p=0;p<size();++p)
            scale=std::max(scale,std::abs(rho_[p]-rhoPrevious_[p])/dt);
        return scale;
    }
    void setPlaneWave(double amplitude,int mode=1) {
        const double waveNumber=2.0*pi*mode/(n_*dx_);
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const double x=origin_.x+(i+0.5)*dx_;
            ey_[index(i,j,k)]=amplitude*std::sin(waveNumber*x);
            bz_[index(i,j,k)]=amplitude*std::sin(waveNumber*(x-0.5*dx_))/c;
        }
    }
    void setPlaneWavePacketX(double amplitude,double width,double centreX=0.0) {
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const double x=origin_.x+(i+0.5)*dx_-centreX;
            const double field=amplitude*std::exp(-0.5*x*x/(width*width));
            ey_[index(i,j,k)]=field;
            // Bz is shifted by half a Yee cell in x.
            const double xb=origin_.x+i*dx_-centreX;
            bz_[index(i,j,k)]=amplitude
                *std::exp(-0.5*xb*xb/(width*width))/c;
        }
    }
    void enableConvolutionalPml(int cells=6,double targetAttenuation=1.0e-8) {
        if(cells<2||2*cells>=n_||!(targetAttenuation>0.0)
           ||!(targetAttenuation<1.0))
            throw std::invalid_argument("invalid Yee CPML configuration");
        absorbingCells_=cells;
        pmlSigmaMaximum_=-4.0*c*std::log(targetAttenuation)/(cells*dx_);
        pmlAlphaMaximum_=0.08*c/dx_;
        for(auto& memory:electricPmlMemory_) memory.assign(size(),0.0);
        for(auto& memory:magneticPmlMemory_) memory.assign(size(),0.0);
    }
    double fieldEnergy(int excludedBoundaryCells=0) const {
        double result=0.0;
        for(int k=excludedBoundaryCells;k<n_-excludedBoundaryCells;++k)
            for(int j=excludedBoundaryCells;j<n_-excludedBoundaryCells;++j)
                for(int i=excludedBoundaryCells;i<n_-excludedBoundaryCells;++i) {
                    const std::size_t p=index(i,j,k);
                    result+=0.5*(epsilon0*(ex_[p]*ex_[p]+ey_[p]*ey_[p]+ez_[p]*ez_[p])
                        +(bx_[p]*bx_[p]+by_[p]*by_[p]+bz_[p]*bz_[p])/mu0)
                        *dx_*dx_*dx_;
                }
        return result;
    }
    std::pair<Vec3,Vec3> interpolateField(const Vec3& position) const {
        return {{interpolateStaggered(ex_,position,{0.0,0.5,0.5}),
                 interpolateStaggered(ey_,position,{0.5,0.0,0.5}),
                 interpolateStaggered(ez_,position,{0.5,0.5,0.0})},
                {interpolateStaggered(bx_,position,{0.5,0.0,0.0}),
                 interpolateStaggered(by_,position,{0.0,0.5,0.0}),
                 interpolateStaggered(bz_,position,{0.0,0.0,0.5})}};
    }
    void setUniformField(const Vec3& electric,const Vec3& magnetic) {
        std::fill(ex_.begin(),ex_.end(),electric.x);
        std::fill(ey_.begin(),ey_.end(),electric.y);
        std::fill(ez_.begin(),ez_.end(),electric.z);
        std::fill(bx_.begin(),bx_.end(),magnetic.x);
        std::fill(by_.begin(),by_.end(),magnetic.y);
        std::fill(bz_.begin(),bz_.end(),magnetic.z);
    }
    void advance(double dt) {
        if(!(dt>0.0)||dt>courantTimeStep())
            throw std::invalid_argument("Yee step violates CFL");
        // Faraday: B^{n+1/2} from edge E^n.
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const std::size_t p=index(i,j,k);
            const double dyEz=pmlDerivative(
                (ez_[index(i,j+1,k)]-ez_[p])/dx_,1,j,dt,
                magneticPmlMemory_[0][p]);
            const double dzEy=pmlDerivative(
                (ey_[index(i,j,k+1)]-ey_[p])/dx_,2,k,dt,
                magneticPmlMemory_[1][p]);
            const double dzEx=pmlDerivative(
                (ex_[index(i,j,k+1)]-ex_[p])/dx_,2,k,dt,
                magneticPmlMemory_[2][p]);
            const double dxEz=pmlDerivative(
                (ez_[index(i+1,j,k)]-ez_[p])/dx_,0,i,dt,
                magneticPmlMemory_[3][p]);
            const double dxEy=pmlDerivative(
                (ey_[index(i+1,j,k)]-ey_[p])/dx_,0,i,dt,
                magneticPmlMemory_[4][p]);
            const double dyEx=pmlDerivative(
                (ex_[index(i,j+1,k)]-ex_[p])/dx_,1,j,dt,
                magneticPmlMemory_[5][p]);
            bx_[p]-=dt*(dyEz-dzEy);
            by_[p]-=dt*(dzEx-dxEz);
            bz_[p]-=dt*(dxEy-dyEx);
        }
        // Ampere-Maxwell: E^{n+1} from face B^{n+1/2} and local face J.
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const std::size_t p=index(i,j,k);
            const double dyBz=pmlDerivative(
                (bz_[p]-bz_[index(i,j-1,k)])/dx_,1,j,dt,
                electricPmlMemory_[0][p]);
            const double dzBy=pmlDerivative(
                (by_[p]-by_[index(i,j,k-1)])/dx_,2,k,dt,
                electricPmlMemory_[1][p]);
            const double dzBx=pmlDerivative(
                (bx_[p]-bx_[index(i,j,k-1)])/dx_,2,k,dt,
                electricPmlMemory_[2][p]);
            const double dxBz=pmlDerivative(
                (bz_[p]-bz_[index(i-1,j,k)])/dx_,0,i,dt,
                electricPmlMemory_[3][p]);
            const double dxBy=pmlDerivative(
                (by_[p]-by_[index(i-1,j,k)])/dx_,0,i,dt,
                electricPmlMemory_[4][p]);
            const double dyBx=pmlDerivative(
                (bx_[p]-bx_[index(i,j-1,k)])/dx_,1,j,dt,
                electricPmlMemory_[5][p]);
            ex_[p]+=dt*(c*c*(dyBz-dzBy)-jx_[p]/epsilon0);
            ey_[p]+=dt*(c*c*(dzBx-dxBz)-jy_[p]/epsilon0);
            ez_[p]+=dt*(c*c*(dxBy-dyBx)-jz_[p]/epsilon0);
        }
    }
    double maximumMagneticDivergence(int excludedBoundaryCells=0) const {
        double maximum=0.0;
        for(int k=excludedBoundaryCells;k<n_-excludedBoundaryCells;++k)
            for(int j=excludedBoundaryCells;j<n_-excludedBoundaryCells;++j)
                for(int i=excludedBoundaryCells;i<n_-excludedBoundaryCells;++i)
            maximum=std::max(maximum,std::abs(
                (bx_[index(i+1,j,k)]-bx_[index(i,j,k)]
                +by_[index(i,j+1,k)]-by_[index(i,j,k)]
                +bz_[index(i,j,k+1)]-bz_[index(i,j,k)])/dx_));
        return maximum;
    }
    double magneticDivergenceScale(int excludedBoundaryCells=0) const {
        double scale=1.0;
        for(int k=excludedBoundaryCells;k<n_-excludedBoundaryCells;++k)
            for(int j=excludedBoundaryCells;j<n_-excludedBoundaryCells;++j)
                for(int i=excludedBoundaryCells;i<n_-excludedBoundaryCells;++i) {
                    const std::size_t p=index(i,j,k);
                    scale=std::max({scale,std::abs(bx_[p])/dx_,
                        std::abs(by_[p])/dx_,std::abs(bz_[p])/dx_});
                }
        return scale;
    }
    // Conservative full-coverage AMR reflux/restriction for refinement ratio
    // two. Cell quantities use volume averages; normal face quantities use
    // the four coincident fine-face areas. These are the finite-volume
    // operators that commute with the Yee divergence.
    std::size_t refluxFromFine(const YeeMaxwellBlock& fine) {
        if(fine.n_!=2*n_||std::abs(fine.dx_*2.0-dx_)>1.0e-12*dx_
           ||(fine.origin_-origin_).norm()>1.0e-12*dx_)
            throw std::invalid_argument("incompatible Yee AMR levels");
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            rho_[index(i,j,k)]=fine.averageCellBlock(fine.rho_,i,j,k);
            rhoPrevious_[index(i,j,k)]=
                fine.averageCellBlock(fine.rhoPrevious_,i,j,k);
            ex_[index(i,j,k)]=fine.averageXFace(fine.ex_,i,j,k);
            ey_[index(i,j,k)]=fine.averageYFace(fine.ey_,i,j,k);
            ez_[index(i,j,k)]=fine.averageZFace(fine.ez_,i,j,k);
            bx_[index(i,j,k)]=fine.averageXFace(fine.bx_,i,j,k);
            by_[index(i,j,k)]=fine.averageYFace(fine.by_,i,j,k);
            bz_[index(i,j,k)]=fine.averageZFace(fine.bz_,i,j,k);
            jx_[index(i,j,k)]=fine.averageXFace(fine.jx_,i,j,k);
            jy_[index(i,j,k)]=fine.averageYFace(fine.jy_,i,j,k);
            jz_[index(i,j,k)]=fine.averageZFace(fine.jz_,i,j,k);
        }
        return size();
    }
    double maximumRestrictionResidual(const YeeMaxwellBlock& fine) const {
        if(fine.n_!=2*n_) return std::numeric_limits<double>::infinity();
        double maximum=0.0;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const std::size_t p=index(i,j,k);
            maximum=std::max({maximum,
                std::abs(rho_[p]-fine.averageCellBlock(fine.rho_,i,j,k)),
                std::abs(ex_[p]-fine.averageXFace(fine.ex_,i,j,k)),
                std::abs(ey_[p]-fine.averageYFace(fine.ey_,i,j,k)),
                std::abs(ez_[p]-fine.averageZFace(fine.ez_,i,j,k)),
                c*std::abs(bx_[p]-fine.averageXFace(fine.bx_,i,j,k)),
                c*std::abs(by_[p]-fine.averageYFace(fine.by_,i,j,k)),
                c*std::abs(bz_[p]-fine.averageZFace(fine.bz_,i,j,k))});
        }
        return maximum;
    }
private:
    int n_; double dx_; Vec3 origin_;
    std::vector<double> rho_,rhoPrevious_,ex_,ey_,ez_,bx_,by_,bz_,jx_,jy_,jz_;
    std::vector<double> px_,py_,pz_,pxPrevious_,pyPrevious_,pzPrevious_;
    std::vector<double> mx_,my_,mz_;
    std::array<std::vector<double>,6> electricPmlMemory_,magneticPmlMemory_;
    int absorbingCells_=0;
    double pmlSigmaMaximum_=0.0,pmlAlphaMaximum_=0.0;
    double pmlKappaMaximum_=6.0;
    std::size_t size() const { return static_cast<std::size_t>(n_)*n_*n_; }
    int wrap(int value) const { value%=n_;return value<0?value+n_:value; }
    std::size_t index(int i,int j,int k) const {
        return (static_cast<std::size_t>(wrap(k))*n_+wrap(j))*n_+wrap(i);
    }
    struct PmlCoordinate { double sigma=0.0,kappa=1.0,alpha=0.0; };
    PmlCoordinate pmlCoordinate(int coordinate) const {
        if(absorbingCells_==0) return {};
        const int edge=std::min(coordinate,n_-1-coordinate);
        if(edge>=absorbingCells_) return {};
        const double depth=static_cast<double>(absorbingCells_-edge)
                          /absorbingCells_;
        const double cubic=depth*depth*depth;
        return {pmlSigmaMaximum_*cubic,
                1.0+(pmlKappaMaximum_-1.0)*cubic,
                pmlAlphaMaximum_*(1.0-depth)};
    }
    double pmlDerivative(double raw,int /*axis*/,int coordinate,double dt,
                         double& memory) {
        const PmlCoordinate profile=pmlCoordinate(coordinate);
        if(profile.sigma==0.0) { memory=0.0; return raw; }
        const double decay=profile.sigma/profile.kappa+profile.alpha;
        const double b=std::exp(-decay*dt);
        const double denominator=profile.sigma+profile.kappa*profile.alpha;
        const double a=denominator>0.0
            ? profile.sigma*(b-1.0)/(profile.kappa*denominator) : 0.0;
        memory=b*memory+a*raw;
        return raw/profile.kappa+memory;
    }
    double averageCellBlock(const std::vector<double>& field,
                            int i,int j,int k) const {
        double result=0.0;
        for(int dk=0;dk<2;++dk) for(int dj=0;dj<2;++dj)
            for(int di=0;di<2;++di)
                result+=field[index(2*i+di,2*j+dj,2*k+dk)]/8.0;
        return result;
    }
    double averageXFace(const std::vector<double>& field,
                        int i,int j,int k) const {
        double result=0.0;
        for(int dk=0;dk<2;++dk) for(int dj=0;dj<2;++dj)
            result+=field[index(2*i,2*j+dj,2*k+dk)]/4.0;
        return result;
    }
    double averageYFace(const std::vector<double>& field,
                        int i,int j,int k) const {
        double result=0.0;
        for(int dk=0;dk<2;++dk) for(int di=0;di<2;++di)
            result+=field[index(2*i+di,2*j,2*k+dk)]/4.0;
        return result;
    }
    double averageZFace(const std::vector<double>& field,
                        int i,int j,int k) const {
        double result=0.0;
        for(int dj=0;dj<2;++dj) for(int di=0;di<2;++di)
            result+=field[index(2*i+di,2*j+dj,2*k)]/4.0;
        return result;
    }
    Vec3 gridCoordinate(const Vec3& position) const {
        return (position-origin_)/dx_-Vec3{0.5,0.5,0.5};
    }
    static Vec3 contractedOffset(const Vec3& restOffset,const Vec3& velocity) {
        const double speed=velocity.norm();
        if(speed==0.0) return restOffset;
        const Vec3 direction=velocity/speed;
        return restOffset+direction*dot(restOffset,direction)
            *(1.0/gamma(velocity)-1.0);
    }
    template<class Function>
    void forEachGaussianMarker(double radius,const Vec3& centre,
                               const Vec3& velocity,Function&& function) {
        if(!(radius>0.0)) throw std::invalid_argument("invalid Gaussian radius");
        constexpr std::array<double,5> node{{-2.0,-1.0,0.0,1.0,2.0}};
        double normalization=0.0;
        for(double x:node) for(double y:node) for(double z:node)
            normalization+=std::exp(-0.5*(x*x+y*y+z*z));
        for(double x:node) for(double y:node) for(double z:node) {
            const double weight=std::exp(-0.5*(x*x+y*y+z*z))/normalization;
            const Vec3 restOffset=Vec3{x,y,z}*radius;
            function(weight,centre+contractedOffset(restOffset,velocity));
        }
    }
    void depositStaggered(std::vector<double>& field,double integratedValue,
                          const Vec3& position,const Vec3& staggering) {
        const Vec3 coordinate=(position-origin_)/dx_-staggering;
        const int i=static_cast<int>(std::floor(coordinate.x));
        const int j=static_cast<int>(std::floor(coordinate.y));
        const int k=static_cast<int>(std::floor(coordinate.z));
        const double fx=coordinate.x-i,fy=coordinate.y-j,fz=coordinate.z-k;
        const double inverseVolume=1.0/(dx_*dx_*dx_);
        for(int dk=0;dk<2;++dk) for(int dj=0;dj<2;++dj)
            for(int di=0;di<2;++di)
                field[index(i+di,j+dj,k+dk)]+=integratedValue*inverseVolume
                    *(di?fx:1-fx)*(dj?fy:1-fy)*(dk?fz:1-fz);
    }
    double interpolateStaggered(const std::vector<double>& field,
                                const Vec3& position,
                                const Vec3& staggering) const {
        const Vec3 coordinate=(position-origin_)/dx_-staggering;
        const int i=static_cast<int>(std::floor(coordinate.x));
        const int j=static_cast<int>(std::floor(coordinate.y));
        const int k=static_cast<int>(std::floor(coordinate.z));
        const double fx=coordinate.x-i,fy=coordinate.y-j,fz=coordinate.z-k;
        double value=0.0;
        for(int dk=0;dk<2;++dk) for(int dj=0;dj<2;++dj)
            for(int di=0;di<2;++di)
                value+=field[index(i+di,j+dj,k+dk)]
                    *(di?fx:1-fx)*(dj?fy:1-fy)*(dk?fz:1-fz);
        return value;
    }
    void addBoundSources(double dt,bool includePolarizationCurrent) {
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const std::size_t p=index(i,j,k);
            rho_[p]-=(px_[index(i+1,j,k)]-px_[p]
                    +py_[index(i,j+1,k)]-py_[p]
                    +pz_[index(i,j,k+1)]-pz_[p])/dx_;
            // Use the same forward incidence orientation as the face-current
            // divergence above. Mixed forward differences commute exactly,
            // hence div(curl(M)) is zero to roundoff on the periodic block.
            jx_[p]+=(mz_[index(i,j+1,k)]-mz_[p]
                    -my_[index(i,j,k+1)]+my_[p])/dx_;
            jy_[p]+=(mx_[index(i,j,k+1)]-mx_[p]
                    -mz_[index(i+1,j,k)]+mz_[p])/dx_;
            jz_[p]+=(my_[index(i+1,j,k)]-my_[p]
                    -mx_[index(i,j+1,k)]+mx_[p])/dx_;
            if(includePolarizationCurrent) {
                jx_[p]+=(px_[p]-pxPrevious_[p])/dt;
                jy_[p]+=(py_[p]-pyPrevious_[p])/dt;
                jz_[p]+=(pz_[p]-pzPrevious_[p])/dt;
            }
        }
    }
    void depositCharge(double charge,const Vec3& position,
                       std::vector<double>& density) {
        const Vec3 coordinate=gridCoordinate(position);
        const int i=static_cast<int>(std::floor(coordinate.x));
        const int j=static_cast<int>(std::floor(coordinate.y));
        const int k=static_cast<int>(std::floor(coordinate.z));
        const double fx=coordinate.x-i,fy=coordinate.y-j,fz=coordinate.z-k;
        const double inverseVolume=1.0/(dx_*dx_*dx_);
        for(int dk=0;dk<2;++dk) for(int dj=0;dj<2;++dj)
            for(int di=0;di<2;++di)
                density[index(i+di,j+dj,k+dk)]+=charge*inverseVolume
                    *(di?fx:1-fx)*(dj?fy:1-fy)*(dk?fz:1-fz);
    }
    void depositAxisLeg(double charge,const Vec3& before,const Vec3& after,
                        int axis,double dt) {
        const Vec3 a=gridCoordinate(before),b=gridCoordinate(after);
        const double displacement=axis==0?b.x-a.x:axis==1?b.y-a.y:b.z-a.z;
        if(displacement==0.0) return;
        const Vec3 midpoint=(a+b)*0.5;
        const int transverseA=axis==0?static_cast<int>(std::floor(midpoint.y))
                              :static_cast<int>(std::floor(midpoint.x));
        const int transverseB=axis==2?static_cast<int>(std::floor(midpoint.y))
                              :static_cast<int>(std::floor(midpoint.z));
        const double fractionA=axis==0?midpoint.y-transverseA
                              :midpoint.x-transverseA;
        const double fractionB=axis==2?midpoint.y-transverseB
                              :midpoint.z-transverseB;
        const double faceCurrent=charge*displacement/(dt*dx_*dx_);
        double start=axis==0?a.x:axis==1?a.y:a.z;
        const double end=axis==0?b.x:axis==1?b.y:b.z;
        const double direction=end>start?1.0:-1.0;
        while(direction*(end-start)>1.0e-14) {
            const int cell=static_cast<int>(std::floor(start));
            const double boundary=direction>0?cell+1.0:cell;
            const double finish=direction>0?std::min(end,boundary)
                                             :std::max(end,boundary);
            const double localFraction=(finish-start)/displacement;
            // CIC charge nodes are the neighbouring cell centres `cell` and
            // `cell+1`; their connecting Yee face has index `cell+1` for
            // either direction of travel.  The current sign carries the
            // direction.
            const int face=cell+1;
            for(int da=0;da<2;++da) for(int db=0;db<2;++db) {
                const double weight=(da?fractionA:1-fractionA)
                                   *(db?fractionB:1-fractionB);
                if(axis==0) jx_[index(face,transverseA+da,transverseB+db)]
                    +=faceCurrent*localFraction*weight;
                else if(axis==1) jy_[index(transverseA+da,face,transverseB+db)]
                    +=faceCurrent*localFraction*weight;
                else jz_[index(transverseA+da,transverseB+db,face)]
                    +=faceCurrent*localFraction*weight;
            }
            start=finish;
        }
    }
};

// Reference two-level Yee hierarchy. The fine level covers the same physical
// volume, making every coarse face an AMR interface for validation. This is a
// deliberately strict reflux test before moving the operator to local patches.
class YeeAmrHierarchy {
public:
    YeeAmrHierarchy(int coarseCells,double coarseSpacing,Vec3 centre)
        : coarse_(coarseCells,coarseSpacing,centre),
          fine_(2*coarseCells,0.5*coarseSpacing,centre) {}
    YeeMaxwellBlock& coarse() { return coarse_; }
    YeeMaxwellBlock& fine() { return fine_; }
    const YeeMaxwellBlock& coarse() const { return coarse_; }
    const YeeMaxwellBlock& fine() const { return fine_; }
    std::size_t reflux() { return coarse_.refluxFromFine(fine_); }
    void enableOuterConvolutionalPml(int cells=6,
                                     double attenuation=1.0e-8) {
        coarse_.enableConvolutionalPml(cells,attenuation);
    }
    void setPlaneWave(double amplitude,int mode=1) {
        coarse_.setPlaneWave(amplitude,mode);
        fine_.setPlaneWave(amplitude,mode);
        reflux();
    }
    void advanceSubcycled(double coarseDt) {
        if(!(coarseDt>0.0)||coarseDt>coarse_.courantTimeStep())
            throw std::invalid_argument("Yee AMR step violates coarse CFL");
        coarse_.advance(coarseDt);
        fine_.advance(0.5*coarseDt);
        fine_.advance(0.5*coarseDt);
        // The area restriction is the accumulated fine flux register for a
        // fully covered parent. Replacing the provisional coarse update is
        // exactly the finite-volume reflux correction on this reference grid.
        reflux();
    }
private:
    YeeMaxwellBlock coarse_,fine_;
};

// Cell-centred finite-difference Maxwell block.  The centred curl operators
// commute, so div(curl(.)) vanishes to roundoff in the periodic interior.
// A Poisson projection restores Gauss's electric constraint after source
// deposition and after every field step.
class MaxwellBlock {
public:
    MaxwellBlock(int cellsPerAxis, double cellSize, Vec3 centre)
        : n_(cellsPerAxis), dx_(cellSize), origin_(centre
          - Vec3{0.5*cellsPerAxis*cellSize, 0.5*cellsPerAxis*cellSize,
                 0.5*cellsPerAxis*cellSize}),
          cells_(static_cast<std::size_t>(n_)*n_*n_),
          scalarWork_(cells_.size()), vectorWork_(cells_.size()),
          magneticWork_(cells_.size()), magnetizationWork_(cells_.size()),
          polarizationWork_(cells_.size()) {
        if (n_ < 6 || !(dx_ > 0.0)) {
            throw std::invalid_argument("Maxwell block requires n>=6 and dx>0");
        }
    }

    int cellsPerAxis() const { return n_; }
    double cellSize() const { return dx_; }
    double extent() const { return n_*dx_; }
    Vec3 centre() const { return origin_+Vec3{0.5*extent(),0.5*extent(),0.5*extent()}; }
    double courantTimeStep() const { return 0.95*dx_/(c*std::sqrt(3.0)); }
    const std::vector<MaxwellCell>& cells() const { return cells_; }
    double absorbedBoundaryEnergy() const { return absorbedBoundaryEnergy_; }
    void enableConvolutionalPml(int cellThickness=4,
                                double targetAttenuation=1.0e-8) {
        if (cellThickness<2 || 2*cellThickness>=n_) {
            throw std::invalid_argument("absorbing layer does not fit the block");
        }
        absorbingCells_=cellThickness;
        // Integral of sigma_max*(x/L)^3/c across the layer equals
        // sigma_max*L/(4c).  Choose sigma_max for the requested attenuation.
        absorbingSigmaMaximum_=-4.0*c*std::log(targetAttenuation)
                              /(absorbingCells_*dx_);
        pmlAlphaMaximum_=0.08*c/dx_;
        for(auto& memory:electricPmlMemory_) memory.assign(cells_.size(),Vec3{});
        for(auto& memory:magneticPmlMemory_) memory.assign(cells_.size(),Vec3{});
    }
    void enableMatchedAbsorbingLayer(int cellThickness=4,
                                     double targetAttenuation=1.0e-8) {
        enableConvolutionalPml(cellThickness,targetAttenuation);
    }
    std::vector<double> chargeDensitySnapshot() const {
        std::vector<double> result; result.reserve(cells_.size());
        for (const MaxwellCell& cell : cells_) result.push_back(cell.chargeDensity);
        return result;
    }

    Vec3 cellPosition(int i, int j, int k) const {
        return origin_ + Vec3{(i + 0.5)*dx_, (j + 0.5)*dx_, (k + 0.5)*dx_};
    }
    Vec3 cellPhase(const Vec3& position) const {
        const Vec3 coordinate=(position-origin_)/dx_-Vec3{0.5,0.5,0.5};
        const auto periodic=[](double value) {
            value-=std::floor(value);
            return value;
        };
        return {periodic(coordinate.x),periodic(coordinate.y),
                periodic(coordinate.z)};
    }
    bool contains(const Vec3& position, double margin=0.0) const {
        const Vec3 upper=origin_+Vec3{extent(),extent(),extent()};
        return position.x>=origin_.x+margin && position.x<=upper.x-margin
            && position.y>=origin_.y+margin && position.y<=upper.y-margin
            && position.z>=origin_.z+margin && position.z<=upper.z-margin;
    }
    bool recenterIfNeeded(const Vec3& requestedCentre,double guardFraction=0.2) {
        const Vec3 displacement=requestedCentre-centre();
        const double threshold=guardFraction*extent();
        if(std::abs(displacement.x)<threshold&&std::abs(displacement.y)<threshold
           &&std::abs(displacement.z)<threshold) return false;
        const Vec3 cellShift{
            std::round(displacement.x/dx_),std::round(displacement.y/dx_),
            std::round(displacement.z/dx_)};
        const Vec3 oldOrigin=origin_;
        const std::vector<MaxwellCell> oldCells=cells_;
        origin_+=cellShift*dx_;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const Vec3 position=cellPosition(i,j,k);
            const Vec3 oldCoordinate=(position-oldOrigin)/dx_-Vec3{0.5,0.5,0.5};
            const int oi=static_cast<int>(std::lround(oldCoordinate.x));
            const int oj=static_cast<int>(std::lround(oldCoordinate.y));
            const int ok=static_cast<int>(std::lround(oldCoordinate.z));
            MaxwellCell replacement;
            if(oi>=0&&oi<n_&&oj>=0&&oj<n_&&ok>=0&&ok<n_)
                replacement=oldCells[(static_cast<std::size_t>(ok)*n_+oj)*n_+oi];
            at(i,j,k)=replacement;
        }
        std::fill(magnetizationWork_.begin(),magnetizationWork_.end(),Vec3{});
        for(auto& memory:electricPmlMemory_)
            std::fill(memory.begin(),memory.end(),Vec3{});
        for(auto& memory:magneticPmlMemory_)
            std::fill(memory.begin(),memory.end(),Vec3{});
        return true;
    }

    void clearSources() {
        for (std::size_t p=0;p<cells_.size();++p) {
            MaxwellCell& cell=cells_[p];
            cell.chargeDensity = 0.0;
            cell.currentDensity = {};
            magnetizationWork_[p]={};
            polarizationWork_[p]={};
        }
    }
    void clearFields() {
        for(MaxwellCell& cell:cells_) {
            cell.electric={};
            cell.magnetic={};
        }
    }
    void addMagneticCurl(const std::vector<Vec3>& vectorPotential) {
        if(vectorPotential.size()!=cells_.size())
            throw std::invalid_argument("invalid magnetic vector potential");
        forEachCell([&](int i,int j,int k,std::size_t p) {
            cells_[p].magnetic+=curlOf(vectorPotential,i,j,k);
        });
    }

    void depositCloud(const RelativisticChargeCloud& cloud, const Vec3& centre,
                      const Vec3& velocity) {
        depositCloudComponents(cloud,centre,velocity,true,true);
    }
    void depositChargeDensity(const RelativisticChargeCloud& cloud,
                              const Vec3& centre) {
        depositCloudComponents(cloud,centre,{},true,false);
    }
    void depositConvectionCurrent(const RelativisticChargeCloud& cloud,
                                  const Vec3& centre,const Vec3& velocity) {
        depositCloudComponents(cloud,centre,velocity,false,true);
    }

private:
    void depositCloudComponents(const RelativisticChargeCloud& cloud,
                                const Vec3& centre,const Vec3& velocity,
                                bool depositCharge,bool depositCurrent) {
        std::vector<std::pair<std::size_t, double>> weights;
        double weightSum = 0.0;
        const double radius = cloud.supportRadius(velocity);
        for (int k = 0; k < n_; ++k) for (int j = 0; j < n_; ++j)
            for (int i = 0; i < n_; ++i) {
                const Vec3 position = cellPosition(i, j, k);
                if ((position - centre).norm() > radius) continue;
                const double density = cloud.density(position, centre, velocity);
                if (density == 0.0) continue;
                weights.emplace_back(index(i, j, k), density);
                weightSum += density*dx_*dx_*dx_;
            }
        if (weights.empty() || weightSum == 0.0) return;
        const double normalization = cloud.charge/weightSum;
        for (const auto& [cellIndex, density] : weights) {
            const double depositedDensity = density*normalization;
            if(depositCharge) cells_[cellIndex].chargeDensity += depositedDensity;
            if(depositCurrent) cells_[cellIndex].currentDensity += velocity*depositedDensity;
        }
    }

public:

    void depositMagneticDipole(const Vec3& centre,const Vec3& dipole,
                               double radius=chargeCloudRestRadius) {
        std::vector<std::pair<std::size_t,double>> weights;
        double integral=0.0;
        const double cutoff=4.0*radius;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const Vec3 offset=cellPosition(i,j,k)-centre;
            if(offset.norm()>cutoff) continue;
            const double weight=std::exp(-0.5*offset.squaredNorm()/(radius*radius));
            weights.emplace_back(index(i,j,k),weight);
            integral+=weight*dx_*dx_*dx_;
        }
        if(integral==0.0) return;
        for(const auto& [p,weight]:weights) magnetizationWork_[p]+=dipole*(weight/integral);
    }

    void depositCovariantDipole(const Vec3& centre,const Vec3& velocity,
                                const Vec3& restDipole,
                                double radius=chargeCloudRestRadius) {
        std::vector<std::pair<std::size_t,double>> weights;
        double integral=0.0; const double cutoff=4.0*radius;
        const double g=gamma(velocity);
        const Vec3 beta=velocity/c;
        const Vec3 labMagnetizationMoment=restDipole*g
            -beta*(dot(beta,restDipole)*g*g/(g+1.0));
        const Vec3 labPolarizationMoment=
            cross(velocity,labMagnetizationMoment)/(c*c);
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const Vec3 offset=cellPosition(i,j,k)-centre;
            if(offset.norm()>cutoff) continue;
            const double weight=std::exp(-0.5*offset.squaredNorm()/(radius*radius));
            weights.emplace_back(index(i,j,k),weight);
            integral+=weight*dx_*dx_*dx_;
        }
        if(integral==0.0) return;
        for(const auto& [p,weight]:weights) {
            magnetizationWork_[p]+=labMagnetizationMoment*(weight/integral);
            polarizationWork_[p]+=labPolarizationMoment*(weight/integral);
        }
    }

    std::vector<Vec3> polarizationSnapshot() const { return polarizationWork_; }

    void finalizeBoundCurrent(const std::vector<Vec3>& previousPolarization,
                              double dt) {
        if(previousPolarization.size()!=cells_.size()||!(dt>0.0))
            throw std::invalid_argument("invalid bound-current history");
        forEachCell([&](int i,int j,int k,std::size_t p) {
            cells_[p].chargeDensity-=divergenceOf(polarizationWork_,i,j,k);
            cells_[p].currentDensity+=curlOf(magnetizationWork_,i,j,k)
                +(polarizationWork_[p]-previousPolarization[p])/dt;
        });
    }
    void finalizeBoundInstantaneous() {
        forEachCell([&](int i,int j,int k,std::size_t p) {
            cells_[p].chargeDensity-=divergenceOf(polarizationWork_,i,j,k);
            cells_[p].currentDensity+=curlOf(magnetizationWork_,i,j,k);
        });
    }

    void finalizeMagnetizationCurrent() {
        forEachCell([&](int i,int j,int k,std::size_t p) {
            cells_[p].currentDensity+=curlOf(magnetizationWork_,i,j,k);
        });
    }

    GridDipoleInteraction dipoleInteraction(const Vec3& centre,const Vec3& dipole,
                                             double radius=chargeCloudRestRadius) const {
        std::vector<Vec3> magnetization(cells_.size());
        double integral=0.0; const double volume=dx_*dx_*dx_;
        const double cutoff=4.0*radius;
        std::vector<std::pair<std::size_t,double>> weights;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const Vec3 offset=cellPosition(i,j,k)-centre;
            if(offset.norm()>cutoff) continue;
            const double weight=std::exp(-0.5*offset.squaredNorm()/(radius*radius));
            weights.emplace_back(index(i,j,k),weight); integral+=weight*volume;
        }
        if(integral==0.0) return {};
        for(const auto& [p,weight]:weights) magnetization[p]=dipole*(weight/integral);
        GridDipoleInteraction result;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const Vec3 current=curlOf(magnetization,i,j,k);
            if(current.squaredNorm()==0.0) continue;
            const Vec3 forceDensity=cross(current,at(i,j,k).magnetic);
            result.force+=forceDensity*volume;
            result.torque+=cross(cellPosition(i,j,k)-centre,forceDensity)*volume;
            result.fieldPower+=dot(current,at(i,j,k).electric)*volume;
        }
        return result;
    }

    void setVacuumPlaneWave(double amplitude, int mode = 1) {
        const double waveNumber = 2.0*pi*mode/extent();
        setVacuumPlaneWaveWavenumber(amplitude,waveNumber);
    }
    void setVacuumPlaneWaveWavenumber(double amplitude,double waveNumber) {
        for (int k = 0; k < n_; ++k) for (int j = 0; j < n_; ++j)
            for (int i = 0; i < n_; ++i) {
                MaxwellCell& cell = at(i,j,k);
                const double phase = waveNumber*cellPosition(i,j,k).x;
                cell.electric = {0.0, amplitude*std::sin(phase), 0.0};
                cell.magnetic = {0.0, 0.0, amplitude*std::sin(phase)/c};
            }
    }
    void setLocalizedWavePacket(double amplitude, double width) {
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            MaxwellCell& cell=at(i,j,k);
            const Vec3 position=cellPosition(i,j,k);
            const double envelope=std::exp(-0.5*position.squaredNorm()/(width*width));
            cell.electric={0.0,amplitude*envelope,0.0};
            cell.magnetic={0.0,0.0,amplitude*envelope/c};
        }
    }
    void setPlaneWavePacketX(double amplitude,double width,double centreX=0.0) {
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            MaxwellCell& cell=at(i,j,k);
            const double x=cellPosition(i,j,k).x-centreX;
            const double field=amplitude*std::exp(-0.5*x*x/(width*width));
            cell.electric={0.0,field,0.0};
            cell.magnetic={0.0,0.0,field/c};
        }
    }
    double planeWaveRelativeError(double amplitude,double waveNumber,
                                  double elapsedTime) const {
        double squaredError=0.0,squaredReference=0.0;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const double reference=amplitude*std::sin(
                waveNumber*(cellPosition(i,j,k).x-c*elapsedTime));
            const MaxwellCell& cell=at(i,j,k);
            squaredError+=(cell.electric.y-reference)*(cell.electric.y-reference)
                         +c*c*(cell.magnetic.z-reference/c)
                              *(cell.magnetic.z-reference/c);
            squaredReference+=2.0*reference*reference;
        }
        return std::sqrt(squaredError/std::max(squaredReference,1.0e-300));
    }
    void setUniformField(const Vec3& electric,const Vec3& magnetic) {
        for(MaxwellCell& cell:cells_) {
            cell.electric=electric;
            cell.magnetic=magnetic;
        }
    }

    std::pair<Vec3,Vec3> interpolateField(const Vec3& position) const {
        const Vec3 coordinate = (position-origin_)/dx_ - Vec3{0.5,0.5,0.5};
        const int i0=static_cast<int>(std::floor(coordinate.x));
        const int j0=static_cast<int>(std::floor(coordinate.y));
        const int k0=static_cast<int>(std::floor(coordinate.z));
        const double fx=coordinate.x-i0, fy=coordinate.y-j0, fz=coordinate.z-k0;
        Vec3 electric, magnetic;
        for(int dk=0;dk<2;++dk) for(int dj=0;dj<2;++dj) for(int di=0;di<2;++di) {
            const double weight=(di?fx:1.0-fx)*(dj?fy:1.0-fy)*(dk?fz:1.0-fz);
            const MaxwellCell& cell=at(i0+di,j0+dj,k0+dk);
            electric += cell.electric*weight;
            magnetic += cell.magnetic*weight;
        }
        return {electric,magnetic};
    }
    std::pair<Vec3,Vec3> numericalSelfField(
        const RelativisticChargeCloud& cloud,const Vec3& position,
        const Vec3& velocity,int projectionIterations=400) const {
        MaxwellBlock isolated(n_,dx_,centre());
        isolated.depositCloud(cloud,position,velocity);
        isolated.projectElectricGaussConstraint(projectionIterations);
        return isolated.interpolateField(position);
    }

    ChargeConservingDepositResult finalizeChargeConservingCurrent(
        const std::vector<double>& previousChargeDensity,double dt,
        int iterations=160) {
        if (previousChargeDensity.size()!=cells_.size() || !(dt>0.0)) {
            throw std::invalid_argument("invalid charge-conserving deposition input");
        }
        // Complete the deposited midpoint convection/bound current by its
        // unique zero-mean longitudinal component.  The transverse current is
        // untouched because the correction is a discrete gradient.
        const std::vector<Vec3> transverseCandidate=[](const auto& cells) {
            std::vector<Vec3> result; result.reserve(cells.size());
            for(const MaxwellCell& cell:cells) result.push_back(cell.currentDensity);
            return result;
        }(cells_);
        std::vector<double> residual(cells_.size());
        forEachCell([&](int i,int j,int k,std::size_t p) {
            residual[p]=(cells_[p].chargeDensity-previousChargeDensity[p])/dt
                       + divergenceCurrent(i,j,k);
        });
        const double mean=std::accumulate(residual.begin(),residual.end(),0.0)
                         /static_cast<double>(residual.size());
        for(double& value:residual) value-=mean;
        solveCompatiblePoisson(residual,iterations);
        forEachCell([&](int i,int j,int k,std::size_t p) {
            cells_[p].currentDensity = cells_[p].currentDensity - Vec3{
                scalar(i+1,j,k)-scalar(i-1,j,k),
                scalar(i,j+1,k)-scalar(i,j-1,k),
                scalar(i,j,k+1)-scalar(i,j,k-1)}*(0.5/dx_);
        });
        ChargeConservingDepositResult result;
        forEachCell([&](int i,int j,int k,std::size_t p) {
            result.continuityResidual=std::max(result.continuityResidual,std::abs(
                (cells_[p].chargeDensity-previousChargeDensity[p])/dt
                +divergenceCurrent(i,j,k)));
        });
        std::vector<Vec3> correction(cells_.size());
        for(std::size_t p=0;p<cells_.size();++p)
            correction[p]=cells_[p].currentDensity-transverseCandidate[p];
        forEachCell([&](int i,int j,int k,std::size_t) {
            result.longitudinalCurl=std::max(
                result.longitudinalCurl,curlOf(correction,i,j,k).norm());
        });
        return result;
    }
    double maximumContinuityResidual(
        const std::vector<double>& previousChargeDensity,double dt) const {
        double maximum=0.0;
        forEachCellConst([&](int i,int j,int k,std::size_t p) {
            maximum=std::max(maximum,std::abs(
                (cells_[p].chargeDensity-previousChargeDensity[p])/dt
                +divergenceCurrent(i,j,k)));
        });
        return maximum;
    }

    void advance(double dt) {
        if (!(dt > 0.0) || dt > courantTimeStep()) {
            throw std::invalid_argument("Maxwell step violates the 3-D CFL limit");
        }
        const double pmlEnergyBefore=pmlLayerEnergy();
        std::vector<Vec3> electricDerivative(cells_.size());
        std::vector<Vec3> magneticDerivative(cells_.size());
        forEachCell([&](int i, int j, int k, std::size_t p) {
            electricDerivative[p] = curlMagnetic(i,j,k)*(c*c)
                                  - cells_[p].currentDensity/epsilon0;
            magneticDerivative[p] = curlElectric(i,j,k)*(-1.0);
        });
        for (std::size_t p=0; p<cells_.size(); ++p) {
            vectorWork_[p] = cells_[p].electric + electricDerivative[p]*(0.5*dt);
            magneticWork_[p] = cells_[p].magnetic + magneticDerivative[p]*(0.5*dt);
        }
        forEachCell([&](int i, int j, int k, std::size_t p) {
            const Vec3 curlB = convolutionalPmlCurl(
                magneticWork_,magneticPmlMemory_,i,j,k,p,dt);
            const Vec3 curlE = convolutionalPmlCurl(
                vectorWork_,electricPmlMemory_,i,j,k,p,dt);
            cells_[p].electric += (curlB*(c*c)
                                  - cells_[p].currentDensity/epsilon0)*dt;
            cells_[p].magnetic += curlE*(-dt);
        });
        // The compatible Poisson solver is periodic.  Applying it across an
        // open CPML boundary would couple opposite faces instantaneously and
        // recreate an incoming wave.  Charge-carrying validation blocks do
        // not use CPML; an eventual open-boundary multigrid projection will
        // replace this branch for source-filled far-field domains.
        if(absorbingCells_==0) projectElectricGaussConstraint(80);
        absorbedBoundaryEnergy_+=std::max(0.0,pmlEnergyBefore-pmlLayerEnergy());
    }

    void projectElectricGaussConstraint(int iterations) {
        // Solve laplacian(phi)=div(E)-rho/eps0, then E <- E-grad(phi).
        std::vector<double> residual(scalarWork_.size(), 0.0);
        forEachCell([&](int i,int j,int k,std::size_t p) {
            residual[p] = divergenceElectric(i,j,k)
                        - cells_[p].chargeDensity/epsilon0;
        });
        // A periodic Poisson equation needs a zero-mean right-hand side.
        const double mean = std::accumulate(residual.begin(), residual.end(), 0.0)
                          / static_cast<double>(residual.size());
        for (double& value : residual) value -= mean;
        solveCompatiblePoisson(residual,iterations);
        forEachCell([&](int i,int j,int k,std::size_t p) {
            cells_[p].electric = cells_[p].electric - Vec3{
                scalar(i+1,j,k)-scalar(i-1,j,k),
                scalar(i,j+1,k)-scalar(i,j-1,k),
                scalar(i,j,k+1)-scalar(i,j,k-1)}*(0.5/dx_);
        });
    }
    void projectMagneticDivergenceConstraint(int iterations) {
        std::vector<double> residual(scalarWork_.size(),0.0);
        forEachCell([&](int i,int j,int k,std::size_t p) {
            residual[p]=divergenceField(i,j,k,false);
        });
        const double mean=std::accumulate(residual.begin(),residual.end(),0.0)
                         /static_cast<double>(residual.size());
        for(double& value:residual) value-=mean;
        solveCompatiblePoisson(residual,iterations);
        forEachCell([&](int i,int j,int k,std::size_t p) {
            cells_[p].magnetic=cells_[p].magnetic-Vec3{
                scalar(i+1,j,k)-scalar(i-1,j,k),
                scalar(i,j+1,k)-scalar(i,j-1,k),
                scalar(i,j,k+1)-scalar(i,j,k-1)}*(0.5/dx_);
        });
    }

    double totalCharge() const {
        double result=0.0;
        for (const MaxwellCell& cell : cells_) result += cell.chargeDensity;
        return result*dx_*dx_*dx_;
    }
    double fieldEnergy() const {
        double result=0.0;
        for (const MaxwellCell& cell : cells_) result +=
            0.5*(epsilon0*cell.electric.squaredNorm()
                + cell.magnetic.squaredNorm()/mu0);
        return result*dx_*dx_*dx_;
    }
    MaxwellVolumeIntegrals volumeIntegrals() const {
        MaxwellVolumeIntegrals result;
        const double volume=dx_*dx_*dx_;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const MaxwellCell& cell=at(i,j,k);
            const Vec3 momentumDensity=cross(cell.electric,cell.magnetic)*epsilon0;
            result.energy+=0.5*(epsilon0*cell.electric.squaredNorm()
                              +cell.magnetic.squaredNorm()/mu0)*volume;
            result.momentum+=momentumDensity*volume;
            result.angularMomentum+=cross(cellPosition(i,j,k),momentumDensity)*volume;
        }
        return result;
    }
    MaxwellBoundaryFlux boundaryFlux() const {
        MaxwellBoundaryFlux result;
        const int lower=absorbingCells_>0 ? absorbingCells_ : 1;
        const int upper=n_-1-lower;
        const double area=dx_*dx_;
        const auto accumulate=[&](int i,int j,int k,const Vec3& normal) {
            const MaxwellCell& cell=at(i,j,k);
            const Vec3 poynting=cross(cell.electric,cell.magnetic)/mu0;
            result.energyRate+=dot(poynting,normal)*area;
            const Vec3 stressOnNormal=
                (cell.electric*dot(cell.electric,normal)
                 -normal*(0.5*cell.electric.squaredNorm()))*epsilon0
              +(cell.magnetic*dot(cell.magnetic,normal)
                 -normal*(0.5*cell.magnetic.squaredNorm()))/mu0;
            const Vec3 outwardMomentum=stressOnNormal*(-area);
            result.momentumRate+=outwardMomentum;
            result.angularMomentumRate+=cross(cellPosition(i,j,k),outwardMomentum);
        };
        for(int a=lower;a<=upper;++a) for(int b=lower;b<=upper;++b) {
            accumulate(lower,a,b,{-1,0,0}); accumulate(upper,a,b,{1,0,0});
            accumulate(a,lower,b,{0,-1,0}); accumulate(a,upper,b,{0,1,0});
            accumulate(a,b,lower,{0,0,-1}); accumulate(a,b,upper,{0,0,1});
        }
        return result;
    }
    double maximumMagneticDivergence() const {
        double result=0.0;
        forEachCellConst([&](int i,int j,int k,std::size_t) {
            result=std::max(result,std::abs(divergenceField(i,j,k,false)));
        });
        return result;
    }
    double maximumMagneticGradientScale() const {
        double result=0.0;
        for(const MaxwellCell& cell:cells_)
            result=std::max(result,cell.magnetic.norm()/dx_);
        return result;
    }
    double maximumElectricGaussResidual() const {
        double result=0.0;
        forEachCellConst([&](int i,int j,int k,std::size_t p) {
            result=std::max(result,std::abs(divergenceElectric(i,j,k)
                       - cells_[p].chargeDensity/epsilon0));
        });
        return result;
    }
    double maximumChargeGaussScale() const {
        double result=0.0;
        for (const MaxwellCell& cell : cells_) {
            result=std::max(result,std::abs(cell.chargeDensity/epsilon0));
        }
        return result;
    }
    double maximumCurrentDivergence() const {
        double result=0.0;
        forEachCellConst([&](int i,int j,int k,std::size_t) {
            result=std::max(result,std::abs(divergenceCurrent(i,j,k)));
        });
        return result;
    }
    double maximumCurrentGradientScale() const {
        double result=0.0;
        for(const MaxwellCell& cell:cells_)
            result=std::max(result,cell.currentDensity.norm()/dx_);
        return result;
    }
    double interiorFieldEnergy(int boundaryCells) const {
        double result=0.0; const double volume=dx_*dx_*dx_;
        for(int k=boundaryCells;k<n_-boundaryCells;++k)
            for(int j=boundaryCells;j<n_-boundaryCells;++j)
                for(int i=boundaryCells;i<n_-boundaryCells;++i) {
                    const MaxwellCell& cell=at(i,j,k);
                    result+=0.5*(epsilon0*cell.electric.squaredNorm()
                               +cell.magnetic.squaredNorm()/mu0)*volume;
                }
        return result;
    }
    std::pair<double,double> xCharacteristicEnergy(int i) const {
        double outgoing=0.0,incoming=0.0;
        const Vec3 normal{1,0,0};
        const double area=dx_*dx_;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) {
            const MaxwellCell& cell=at(i,j,k);
            const Vec3 transverseElectric{0,cell.electric.y,cell.electric.z};
            const Vec3 magneticCharacteristic=cross(normal,cell.magnetic)*c;
            outgoing+=(transverseElectric-magneticCharacteristic).squaredNorm()
                     *(0.25*epsilon0*area);
            incoming+=(transverseElectric+magneticCharacteristic).squaredNorm()
                     *(0.25*epsilon0*area);
        }
        return {outgoing,incoming};
    }
    void replaceFieldAt(const Vec3& position, const Vec3& electric,
                        const Vec3& magnetic) {
        const Vec3 coordinate=(position-origin_)/dx_-Vec3{0.5,0.5,0.5};
        MaxwellCell& cell=at(static_cast<int>(std::lround(coordinate.x)),
                             static_cast<int>(std::lround(coordinate.y)),
                             static_cast<int>(std::lround(coordinate.z)));
        cell.electric=electric;
        cell.magnetic=magnetic;
    }

    void fillBoundaryFromTemporalInterpolation(const MaxwellBlock& before,
                                               const MaxwellBlock& after,
                                               double fraction,
                                               int thickness=2) {
        if(thickness<1||2*thickness>=n_||fraction<0.0||fraction>1.0)
            throw std::invalid_argument("invalid coarse-to-fine boundary fill");
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const int edge=std::min({i,j,k,n_-1-i,n_-1-j,n_-1-k});
            if(edge>=thickness) continue;
            const Vec3 position=cellPosition(i,j,k);
            if(!before.contains(position,before.cellSize())
               ||!after.contains(position,after.cellSize())) continue;
            const auto [electricBefore,magneticBefore]=before.interpolateField(position);
            const auto [electricAfter,magneticAfter]=after.interpolateField(position);
            MaxwellCell& cell=at(i,j,k);
            cell.electric=electricBefore*(1.0-fraction)+electricAfter*fraction;
            cell.magnetic=magneticBefore*(1.0-fraction)+magneticAfter*fraction;
        }
    }

    double maximumBoundaryMismatch(const MaxwellBlock& parent,
                                   int thickness=2) const {
        double maximum=0.0;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const int edge=std::min({i,j,k,n_-1-i,n_-1-j,n_-1-k});
            if(edge>=thickness) continue;
            const Vec3 position=cellPosition(i,j,k);
            if(!parent.contains(position,parent.cellSize())) continue;
            const auto [electric,magnetic]=parent.interpolateField(position);
            maximum=std::max({maximum,(at(i,j,k).electric-electric).norm(),
                              c*(at(i,j,k).magnetic-magnetic).norm()});
        }
        return maximum;
    }

    MaxwellVolumeIntegrals volumeIntegralsOutside(const MaxwellBlock& covered) const {
        MaxwellVolumeIntegrals result;
        const double volume=dx_*dx_*dx_;
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const Vec3 position=cellPosition(i,j,k);
            if(covered.contains(position)) continue;
            const MaxwellCell& cell=at(i,j,k);
            const Vec3 momentumDensity=cross(cell.electric,cell.magnetic)*epsilon0;
            result.energy+=0.5*(epsilon0*cell.electric.squaredNorm()
                              +cell.magnetic.squaredNorm()/mu0)*volume;
            result.momentum+=momentumDensity*volume;
            result.angularMomentum+=cross(position,momentumDensity)*volume;
        }
        return result;
    }
    MaxwellVolumeIntegrals volumeIntegralsInside(const MaxwellBlock& covered) const {
        const MaxwellVolumeIntegrals total=volumeIntegrals();
        const MaxwellVolumeIntegrals outside=volumeIntegralsOutside(covered);
        return {total.energy-outside.energy,total.momentum-outside.momentum,
                total.angularMomentum-outside.angularMomentum};
    }

    std::size_t restrictFieldFrom(const MaxwellBlock& fine) {
        const int ratio=static_cast<int>(std::lround(dx_/fine.cellSize()));
        if(ratio<2||std::abs(dx_/fine.cellSize()-ratio)>1.0e-12)
            throw std::invalid_argument("AMR restriction requires an integer ratio");
        std::size_t restricted=0;
        const double requiredMargin=0.5*dx_;
        const double inverseSamples=1.0/(ratio*ratio*ratio);
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i) {
            const Vec3 position=cellPosition(i,j,k);
            if(!fine.contains(position,requiredMargin)) continue;
            Vec3 electric,magnetic;
            for(int fk=0;fk<ratio;++fk) for(int fj=0;fj<ratio;++fj)
                for(int fi=0;fi<ratio;++fi) {
                    const Vec3 offset{
                        (fi+0.5-0.5*ratio)*fine.cellSize(),
                        (fj+0.5-0.5*ratio)*fine.cellSize(),
                        (fk+0.5-0.5*ratio)*fine.cellSize()};
                    const auto [sampleElectric,sampleMagnetic]=
                        fine.interpolateField(position+offset);
                    electric+=sampleElectric*inverseSamples;
                    magnetic+=sampleMagnetic*inverseSamples;
                }
            replaceFieldAt(position,electric,magnetic);
            ++restricted;
        }
        return restricted;
    }

private:
    int n_; double dx_; Vec3 origin_;
    std::vector<MaxwellCell> cells_;
    std::vector<double> scalarWork_;
    std::vector<Vec3> vectorWork_, magneticWork_, magnetizationWork_,polarizationWork_;
    std::array<std::vector<Vec3>,3> electricPmlMemory_,magneticPmlMemory_;
    int absorbingCells_=0;
    double absorbingSigmaMaximum_=0.0;
    double pmlAlphaMaximum_=0.0;
    double pmlKappaMaximum_=6.0;
    double absorbedBoundaryEnergy_=0.0;
    int wrap(int value) const { value%=n_; return value<0 ? value+n_ : value; }
    std::size_t index(int i,int j,int k) const {
        return (static_cast<std::size_t>(wrap(k))*n_+wrap(j))*n_+wrap(i);
    }
    MaxwellCell& at(int i,int j,int k) { return cells_[index(i,j,k)]; }
    const MaxwellCell& at(int i,int j,int k) const { return cells_[index(i,j,k)]; }
    double scalar(int i,int j,int k) const { return scalarWork_[index(i,j,k)]; }
    template<class F> void forEachCell(F&& f) {
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i)
            f(i,j,k,index(i,j,k));
    }
    template<class F> void forEachCellConst(F&& f) const {
        for(int k=0;k<n_;++k) for(int j=0;j<n_;++j) for(int i=0;i<n_;++i)
            f(i,j,k,index(i,j,k));
    }
    Vec3 curlOf(const std::vector<Vec3>& field,int i,int j,int k) const {
        const auto& xp=field[index(i+1,j,k)]; const auto& xm=field[index(i-1,j,k)];
        const auto& yp=field[index(i,j+1,k)]; const auto& ym=field[index(i,j-1,k)];
        const auto& zp=field[index(i,j,k+1)]; const auto& zm=field[index(i,j,k-1)];
        const double h=0.5/dx_;
        return {(yp.z-ym.z-(zp.y-zm.y))*h,
                (zp.x-zm.x-(xp.z-xm.z))*h,
                (xp.y-xm.y-(yp.x-ym.x))*h};
    }
    double divergenceOf(const std::vector<Vec3>& field,int i,int j,int k) const {
        return ((field[index(i+1,j,k)].x-field[index(i-1,j,k)].x)
              +(field[index(i,j+1,k)].y-field[index(i,j-1,k)].y)
              +(field[index(i,j,k+1)].z-field[index(i,j,k-1)].z))*(0.5/dx_);
    }
    Vec3 curlElectric(int i,int j,int k) const {
        return directCurl(i,j,k,true);
    }
    Vec3 curlMagnetic(int i,int j,int k) const {
        return directCurl(i,j,k,false);
    }
    Vec3 directCurl(int i,int j,int k,bool electric) const {
        const auto value=[&](int a,int b,int d) {
            const MaxwellCell& cell=at(a,b,d); return electric?cell.electric:cell.magnetic;
        };
        const Vec3 xp=value(i+1,j,k),xm=value(i-1,j,k),yp=value(i,j+1,k),
                   ym=value(i,j-1,k),zp=value(i,j,k+1),zm=value(i,j,k-1);
        const double h=0.5/dx_;
        return {(yp.z-ym.z-zp.y+zm.y)*h,(zp.x-zm.x-xp.z+xm.z)*h,
                (xp.y-xm.y-yp.x+ym.x)*h};
    }
    double divergenceField(int i,int j,int k,bool electric) const {
        const auto value=[&](int a,int b,int d) {
            const MaxwellCell& cell=at(a,b,d); return electric?cell.electric:cell.magnetic;
        };
        return ((value(i+1,j,k).x-value(i-1,j,k).x)
              +(value(i,j+1,k).y-value(i,j-1,k).y)
              +(value(i,j,k+1).z-value(i,j,k-1).z))*(0.5/dx_);
    }
    double divergenceElectric(int i,int j,int k) const {
        return divergenceField(i,j,k,true);
    }
    double divergenceCurrent(int i,int j,int k) const {
        return ((at(i+1,j,k).currentDensity.x-at(i-1,j,k).currentDensity.x)
              +(at(i,j+1,k).currentDensity.y-at(i,j-1,k).currentDensity.y)
              +(at(i,j,k+1).currentDensity.z-at(i,j,k-1).currentDensity.z))
              *(0.5/dx_);
    }
    void solveCompatiblePoisson(const std::vector<double>& rightHandSide,
                                int iterations) {
        std::fill(scalarWork_.begin(),scalarWork_.end(),0.0);
        std::vector<double> next(scalarWork_.size(),0.0);
        for(int iteration=0;iteration<iterations;++iteration) {
            forEachCell([&](int i,int j,int k,std::size_t p) {
                next[p]=(scalar(i+2,j,k)+scalar(i-2,j,k)
                        +scalar(i,j+2,k)+scalar(i,j-2,k)
                        +scalar(i,j,k+2)+scalar(i,j,k-2)
                        -4.0*dx_*dx_*rightHandSide[p])/6.0;
            });
            scalarWork_.swap(next);
        }
    }
    struct PmlCoordinate { double sigma=0.0,kappa=1.0,alpha=0.0; };
    PmlCoordinate pmlCoordinate(int coordinate) const {
        if(absorbingCells_==0) return {};
        const int edge=std::min(coordinate,n_-1-coordinate);
        if(edge>=absorbingCells_) return {};
        const double depth=static_cast<double>(absorbingCells_-edge)
                          /absorbingCells_;
        const double cubic=depth*depth*depth;
        return {absorbingSigmaMaximum_*cubic,
                1.0+(pmlKappaMaximum_-1.0)*cubic,
                pmlAlphaMaximum_*(1.0-depth)};
    }
    Vec3 directionalDerivative(const std::vector<Vec3>& field,int axis,
                               int i,int j,int k) const {
        if(axis==0) return (field[index(i+1,j,k)]-field[index(i-1,j,k)])*(0.5/dx_);
        if(axis==1) return (field[index(i,j+1,k)]-field[index(i,j-1,k)])*(0.5/dx_);
        return (field[index(i,j,k+1)]-field[index(i,j,k-1)])*(0.5/dx_);
    }
    Vec3 convolutionalPmlCurl(const std::vector<Vec3>& field,
        std::array<std::vector<Vec3>,3>& memory,int i,int j,int k,
        std::size_t p,double dt) {
        if(absorbingCells_==0) return curlOf(field,i,j,k);
        std::array<Vec3,3> derivative;
        const std::array<int,3> coordinate{i,j,k};
        for(int axis=0;axis<3;++axis) {
            const Vec3 raw=directionalDerivative(field,axis,i,j,k);
            const PmlCoordinate profile=pmlCoordinate(coordinate[axis]);
            if(profile.sigma==0.0) {
                memory[axis][p]={}; derivative[axis]=raw; continue;
            }
            const double decay=profile.sigma/profile.kappa+profile.alpha;
            const double b=std::exp(-decay*dt);
            const double denominator=profile.sigma
                                   +profile.kappa*profile.alpha;
            const double a=denominator>0.0
                ? profile.sigma*(b-1.0)/(profile.kappa*denominator) : 0.0;
            memory[axis][p]=memory[axis][p]*b+raw*a;
            derivative[axis]=raw/profile.kappa+memory[axis][p];
        }
        return {derivative[1].z-derivative[2].y,
                derivative[2].x-derivative[0].z,
                derivative[0].y-derivative[1].x};
    }
    double pmlLayerEnergy() const {
        if(absorbingCells_==0) return 0.0;
        double result=0.0;
        const double volume=dx_*dx_*dx_;
        forEachCellConst([&](int i,int j,int k,std::size_t p) {
            if(std::min({i,j,k,n_-1-i,n_-1-j,n_-1-k})>=absorbingCells_) return;
            result+=0.5*(epsilon0*cells_[p].electric.squaredNorm()
                       +cells_[p].magnetic.squaredNorm()/mu0)*volume;
        });
        return result;
    }
};

class MaxwellAmrHierarchy {
public:
    MaxwellAmrHierarchy(Vec3 centre, int cellsPerAxis = 24,
                        int refinementLevels = 3)
        : refinementRatio_(2) {
        if (refinementLevels < 1) {
            throw std::invalid_argument("AMR hierarchy requires at least one level");
        }
        // The finest level resolves r_c by two cells.  Every coarser level
        // doubles dx and therefore doubles the covered physical extent.
        double spacing = 0.5*chargeCloudRestRadius
                       * std::pow(refinementRatio_, refinementLevels-1);
        for (int level=0; level<refinementLevels; ++level) {
            levels_.emplace_back(cellsPerAxis, spacing, centre);
            spacing /= refinementRatio_;
        }
    }
    std::vector<MaxwellBlock>& levels() { return levels_; }
    const std::vector<MaxwellBlock>& levels() const { return levels_; }
    MaxwellBlock& finest() { return levels_.back(); }
    const MaxwellBlock& finest() const { return levels_.back(); }
    double timeStep() const { return finest().courantTimeStep(); }
    std::pair<Vec3,Vec3> interpolateCompositeField(const Vec3& position) const {
        for(std::size_t level=levels_.size();level-->0;)
            if(levels_[level].contains(position,levels_[level].cellSize()))
                return levels_[level].interpolateField(position);
        return levels_.front().interpolateField(position);
    }
    MaxwellVolumeIntegrals compositeVolumeIntegrals() const {
        MaxwellVolumeIntegrals result;
        for(std::size_t level=0;level<levels_.size();++level) {
            const MaxwellVolumeIntegrals contribution=level+1<levels_.size()
                ? levels_[level].volumeIntegralsOutside(levels_[level+1])
                : levels_[level].volumeIntegrals();
            result.energy+=contribution.energy;
            result.momentum+=contribution.momentum;
            result.angularMomentum+=contribution.angularMomentum;
        }
        return result;
    }
    MaxwellBoundaryFlux farFieldBoundaryFlux() const {
        return levels_.front().boundaryFlux();
    }
    void depositPair(const Vec3& firstPosition, const Vec3& firstVelocity,
                     const Vec3& secondPosition, const Vec3& secondVelocity,
                     const Vec3& firstDipole={},const Vec3& secondDipole={}) {
        const RelativisticChargeCloud first{-eCharge, chargeCloudRestRadius};
        const RelativisticChargeCloud second{eCharge, chargeCloudRestRadius};
        for (MaxwellBlock& level : levels_) {
            level.clearSources();
            level.depositCloud(first, firstPosition, firstVelocity);
            level.depositCloud(second, secondPosition, secondVelocity);
            level.depositMagneticDipole(firstPosition,firstDipole);
            level.depositMagneticDipole(secondPosition,secondDipole);
            level.finalizeMagnetizationCurrent();
        }
    }
    void enableMatchedAbsorbingBoundaries(int cells=4) {
        // Fine-patch boundaries exchange fields with their parent.  Only the
        // common far-field level represents an open physical boundary.
        levels_.front().enableConvolutionalPml(cells);
    }
    std::size_t synchronizeFineToCoarse() {
        // Volume restriction for refinement ratio 2.  Averaging the eight
        // fine-cell centres preserves a constant field and the cell average.
        std::size_t synchronizedCells=0;
        for(std::size_t fineLevel=levels_.size()-1;fineLevel>0;--fineLevel) {
            MaxwellBlock& fine=levels_[fineLevel];
            MaxwellBlock& coarse=levels_[fineLevel-1];
            for(int k=0;k<coarse.cellsPerAxis();++k)
                for(int j=0;j<coarse.cellsPerAxis();++j)
                    for(int i=0;i<coarse.cellsPerAxis();++i) {
                        const Vec3 position=coarse.cellPosition(i,j,k);
                        if(!fine.contains(position,coarse.cellSize())) continue;
                        Vec3 electric,magnetic;
                        for(int dk : {-1,1}) for(int dj : {-1,1})
                            for(int di : {-1,1}) {
                                const Vec3 sample=position+Vec3{
                                    static_cast<double>(di),static_cast<double>(dj),
                                    static_cast<double>(dk)}
                                    *(0.5*fine.cellSize());
                                const auto [sampleElectric,sampleMagnetic]=
                                    fine.interpolateField(sample);
                                electric+=sampleElectric*(1.0/8.0);
                                magnetic+=sampleMagnetic*(1.0/8.0);
                            }
                        coarse.replaceFieldAt(position,electric,magnetic);
                        ++synchronizedCells;
                    }
        }
        return synchronizedCells;
    }
    std::size_t follow(const Vec3& centre) {
        std::size_t moved=0;
        for(MaxwellBlock& level:levels_) if(level.recenterIfNeeded(centre)) ++moved;
        return moved;
    }
    void advanceSubcycled(double coarseTimeStep) {
        if(coarseTimeStep>levels_.front().courantTimeStep())
            throw std::invalid_argument("coarse AMR step violates CFL");
        advanceLevelSubcycled(0,coarseTimeStep);
        synchronizeFineToCoarse();
        // Restriction changes the covered parent cells.  Refill the fine
        // interface once more so both representations end the coarse step at
        // the same physical time and with the same shell values.
        for(std::size_t level=1;level<levels_.size();++level)
            levels_[level].fillBoundaryFromTemporalInterpolation(
                levels_[level-1],levels_[level-1],1.0);
    }
private:
    int refinementRatio_;
    std::vector<MaxwellBlock> levels_;
    void advanceLevelSubcycled(std::size_t level,double dt) {
        const MaxwellBlock before=levels_[level];
        levels_[level].advance(dt);
        if(level+1>=levels_.size()) return;
        for(int substep=0;substep<refinementRatio_;++substep) {
            const double fraction=(substep+0.5)/refinementRatio_;
            levels_[level+1].fillBoundaryFromTemporalInterpolation(
                before,levels_[level],fraction);
            advanceLevelSubcycled(level+1,dt/refinementRatio_);
        }
        levels_[level+1].fillBoundaryFromTemporalInterpolation(
            before,levels_[level],1.0);
    }
};

// Production geometry for a widely separated pair: both particle-centred
// patches exchange radiation through one coarse far-field block.  Keeping the
// patches in one object prevents accidentally evolving two causally isolated
// electromagnetic domains.
class MaxwellPairPatchHierarchy {
public:
    MaxwellPairPatchHierarchy(const Vec3& firstPosition,
                              const Vec3& secondPosition)
        : farField_(32,4.0*chargeCloudRestRadius,
                    (firstPosition+secondPosition)*0.5),
          firstPatch_(24,0.5*chargeCloudRestRadius,firstPosition),
          secondPatch_(24,0.5*chargeCloudRestRadius,secondPosition) {}
    MaxwellBlock& farField() { return farField_; }
    MaxwellBlock& firstPatch() { return firstPatch_; }
    MaxwellBlock& secondPatch() { return secondPatch_; }
    const MaxwellBlock& farField() const { return farField_; }
    const MaxwellBlock& firstPatch() const { return firstPatch_; }
    const MaxwellBlock& secondPatch() const { return secondPatch_; }
    bool coversPair(const Vec3& firstPosition,
                    const Vec3& secondPosition) const {
        const double fineGuard=4.5*chargeCloudRestRadius;
        const double coarseGuard=8.0*chargeCloudRestRadius;
        return firstPatch_.contains(firstPosition,fineGuard)
            &&secondPatch_.contains(secondPosition,fineGuard)
            &&farField_.contains(firstPosition,coarseGuard)
            &&farField_.contains(secondPosition,coarseGuard);
    }
    std::size_t follow(const Vec3& firstPosition,
                       const Vec3& secondPosition) {
        std::size_t moved=0;
        if(firstPatch_.recenterIfNeeded(firstPosition)) ++moved;
        if(secondPatch_.recenterIfNeeded(secondPosition)) ++moved;
        const Vec3 midpoint=(firstPosition+secondPosition)*0.5;
        if(farField_.recenterIfNeeded(midpoint,0.3)) ++moved;
        return moved;
    }
    std::pair<Vec3,Vec3> fieldAt(const Vec3& position) const {
        if(firstPatch_.contains(position,firstPatch_.cellSize()))
            return firstPatch_.interpolateField(position);
        if(secondPatch_.contains(position,secondPatch_.cellSize()))
            return secondPatch_.interpolateField(position);
        return farField_.interpolateField(position);
    }
    void enableFarFieldCpml(int cells=8) {
        farField_.enableConvolutionalPml(cells);
    }
    MaxwellVolumeIntegrals compositeVolumeIntegrals() const {
        const MaxwellVolumeIntegrals far=farField_.volumeIntegrals();
        const MaxwellVolumeIntegrals coveredFirst=
            farField_.volumeIntegralsInside(firstPatch_);
        const MaxwellVolumeIntegrals coveredSecond=
            farField_.volumeIntegralsInside(secondPatch_);
        const MaxwellVolumeIntegrals first=firstPatch_.volumeIntegrals();
        const MaxwellVolumeIntegrals second=secondPatch_.volumeIntegrals();
        return {far.energy-coveredFirst.energy-coveredSecond.energy
                    +first.energy+second.energy,
                far.momentum-coveredFirst.momentum-coveredSecond.momentum
                    +first.momentum+second.momentum,
                far.angularMomentum-coveredFirst.angularMomentum
                    -coveredSecond.angularMomentum+first.angularMomentum
                    +second.angularMomentum};
    }
    struct CouplingStep {
        std::size_t restrictedCells=0;
        double energyDefect=0.0;
        Vec3 momentumDefect,angularMomentumDefect;
    };
    CouplingStep advanceSubcycled(double farTimeStep,bool parallelBranches=true) {
        if(!(farTimeStep>0.0)||farTimeStep>farField_.courantTimeStep())
            throw std::invalid_argument("pair hierarchy violates far-field CFL");
        if((firstPatch_.centre()-secondPatch_.centre()).norm()
            <0.5*(firstPatch_.extent()+secondPatch_.extent()))
            throw std::runtime_error(
                "particle patches overlap; a merged near-field patch is required");
        const MaxwellVolumeIntegrals beforeTotal=compositeVolumeIntegrals();
        const MaxwellBlock farBefore=farField_;
        farField_.advance(farTimeStep);
        const int ratio=static_cast<int>(std::lround(
            farField_.cellSize()/firstPatch_.cellSize()));
        const auto advancePatch=[&](MaxwellBlock& patch) {
            for(int substep=0;substep<ratio;++substep) {
                const double fraction=(substep+0.5)/ratio;
                patch.fillBoundaryFromTemporalInterpolation(
                    farBefore,farField_,fraction);
                patch.advance(farTimeStep/ratio);
            }
            patch.fillBoundaryFromTemporalInterpolation(
                farBefore,farField_,1.0);
        };
        if(parallelBranches) {
            auto firstTask=std::async(std::launch::async,[&] {
                advancePatch(firstPatch_);
            });
            advancePatch(secondPatch_);
            firstTask.get();
        } else {
            advancePatch(firstPatch_);
            advancePatch(secondPatch_);
        }
        CouplingStep result;
        result.restrictedCells+=farField_.restrictFieldFrom(firstPatch_);
        result.restrictedCells+=farField_.restrictFieldFrom(secondPatch_);
        firstPatch_.fillBoundaryFromTemporalInterpolation(
            farField_,farField_,1.0);
        secondPatch_.fillBoundaryFromTemporalInterpolation(
            farField_,farField_,1.0);
        const MaxwellVolumeIntegrals afterTotal=compositeVolumeIntegrals();
        result.energyDefect=afterTotal.energy-beforeTotal.energy;
        result.momentumDefect=afterTotal.momentum-beforeTotal.momentum;
        result.angularMomentumDefect=
            afterTotal.angularMomentum-beforeTotal.angularMomentum;
        return result;
    }
private:
    MaxwellBlock farField_,firstPatch_,secondPatch_;
};

class DynamicSelfFieldCalibration {
public:
    static constexpr int phaseSamples=4;
    DynamicSelfFieldCalibration(const MaxwellBlock& geometry,double charge,
                                int projectionIterations=240)
        : charge_(charge),cellSize_(geometry.cellSize()),centre_(geometry.centre()) {
        const RelativisticChargeCloud cloud{charge,chargeCloudRestRadius};
        samples_.resize(phaseSamples*phaseSamples*phaseSamples);
        for(int k=0;k<phaseSamples;++k) for(int j=0;j<phaseSamples;++j)
            for(int i=0;i<phaseSamples;++i) {
                const Vec3 phase{(i+0.5)/phaseSamples,(j+0.5)/phaseSamples,
                                 (k+0.5)/phaseSamples};
                const Vec3 position=centre_+(phase-Vec3{0.5,0.5,0.5})*cellSize_;
                const auto [electric,magnetic]=geometry.numericalSelfField(
                    cloud,position,{},projectionIterations);
                samples_[index(i,j,k)]={electric,magnetic};
            }
    }
    ElectromagneticField field(const MaxwellBlock& geometry,
                               const Vec3& position,const Vec3& velocity,
                               double requestedCharge) const {
        const Vec3 phase=geometry.cellPhase(position);
        const Vec3 coordinate=phase*phaseSamples-Vec3{0.5,0.5,0.5};
        const int i0=static_cast<int>(std::floor(coordinate.x));
        const int j0=static_cast<int>(std::floor(coordinate.y));
        const int k0=static_cast<int>(std::floor(coordinate.z));
        const double fx=coordinate.x-i0,fy=coordinate.y-j0,fz=coordinate.z-k0;
        Vec3 restElectric,restMagnetic;
        for(int dk=0;dk<2;++dk) for(int dj=0;dj<2;++dj)
            for(int di=0;di<2;++di) {
                const double weight=(di?fx:1.0-fx)*(dj?fy:1.0-fy)
                                   *(dk?fz:1.0-fz);
                const ElectromagneticField& sample=
                    samples_[index(i0+di,j0+dj,k0+dk)];
                restElectric+=sample.electric*weight;
                restMagnetic+=sample.magnetic*weight;
            }
        const double g=1.0/std::sqrt(std::max(1.0e-15,
            1.0-velocity.squaredNorm()/(c*c)));
        const Vec3 beta=velocity/c;
        const double betaDotElectric=beta.x*restElectric.x
            +beta.y*restElectric.y+beta.z*restElectric.z;
        const Vec3 parallel=beta.squaredNorm()>0.0
            ? beta*(betaDotElectric/beta.squaredNorm()) : Vec3{};
        const Vec3 labElectric=parallel+(restElectric-parallel)*g
            -cross(velocity,restMagnetic)*g;
        const Vec3 labMagnetic=restMagnetic*g
            -cross(velocity,restElectric)*(g/(c*c));
        const double chargeScale=requestedCharge/charge_;
        return {labElectric*chargeScale,labMagnetic*chargeScale};
    }
private:
    double charge_,cellSize_; Vec3 centre_;
    std::vector<ElectromagneticField> samples_;
    static int wrap(int value) {
        value%=phaseSamples; return value<0?value+phaseSamples:value;
    }
    static std::size_t index(int i,int j,int k) {
        return (static_cast<std::size_t>(wrap(k))*phaseSamples+wrap(j))
             *phaseSamples+wrap(i);
    }
};
