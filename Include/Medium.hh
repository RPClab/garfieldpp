#ifndef G_MEDIUM_H
#define G_MEDIUM_H

#include <string>
#include <vector>

#include "FundamentalConstants.hh"

namespace Garfield {

/// Abstract base class for media.

class Medium {
 public:
  /// Constructor
  Medium();
  /// Destructor
  virtual ~Medium();

  /// Return the id number of the class instance.
  int GetId() const { return m_id; }
  /// Get the medium name/identifier.
  const std::string& GetName() const { return m_name; }
  /// Is this medium a gas?
  virtual bool IsGas() const { return false; }
  /// Is this medium a semiconductor?
  virtual bool IsSemiconductor() const { return false; }

  /// Set the temperature [K].
  void SetTemperature(const double t);
  /// Get the temperature [K].
  double GetTemperature() const { return m_temperature; }
  // Set the pressure [Torr].
  void SetPressure(const double p);
  // Get the pressure [Torr].
  double GetPressure() const { return m_pressure; }
  /// Set the relative static dielectric constant.
  void SetDielectricConstant(const double eps);
  /// Get the relative static dielectric constant.
  double GetDielectricConstant() const { return m_epsilon; }

  /// Get number of components of the medium.
  unsigned int GetNumberOfComponents() const { return m_nComponents; }
  /// Get the name and fraction of a given component.
  virtual void GetComponent(const unsigned int i, std::string& label,
                            double& f);
  /// Set the effective atomic number.
  virtual void SetAtomicNumber(const double z);
  /// Get the effective atomic number.
  virtual double GetAtomicNumber() const { return m_z; }
  /// Set the effective atomic weight.
  virtual void SetAtomicWeight(const double a);
  /// Get the effective atomic weight.
  virtual double GetAtomicWeight() const { return m_a; }
  /// Set the number density [cm-3].
  virtual void SetNumberDensity(const double n);
  /// Get the number density [cm-3].
  virtual double GetNumberDensity() const { return m_density; }
  /// Set the mass density [g/cm3].
  virtual void SetMassDensity(const double rho);
  /// Get the mass density [g/cm3].
  virtual double GetMassDensity() const;

  /// Switch electron/ion/hole on/off.
  virtual void EnableDrift(const bool on = true) { m_driftable = on; }
  /// Make the medium ionisable.
  virtual void EnablePrimaryIonisation() { m_ionisable = true; }
  /// Make the medium non-ionisable.
  void DisablePrimaryIonisation() { m_ionisable = false; }

  /// Is charge carrier transport enabled in this medium?
  bool IsDriftable() const { return m_driftable; }
  /// Does the medium have electron scattering rates?
  bool IsMicroscopic() const { return m_microscopic; }
  /// Is charge deposition by charged particles/photon enabled in this medium?
  bool IsIonisable() const { return m_ionisable; }

  /// Set the W value (average energy to produce an electron/ion or e/h pair).
  void SetW(const double w) { m_w = w; }
  /// Get the W value.
  double GetW() { return m_w; }
  /// Set the Fano factor.
  void SetFanoFactor(const double f) { m_fano = f; }
  /// Get the Fano factor.
  double GetFanoFactor() { return m_fano; }

  // Transport parameters for electrons
  /// Drift velocity [cm / ns]
  virtual bool ElectronVelocity(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& vx,
                                double& vy, double& vz);
  /// Longitudinal and transverse diffusion coefficients [cm1/2]
  virtual bool ElectronDiffusion(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz, double& dl,
                                 double& dt);
  /// Diffusion tensor: diagonal elements are the diffusion
  /// coefficients [cm] along e, btrans, e x b,
  /// off-diagonal elements are the covariances
  virtual bool ElectronDiffusion(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz,
                                 double cov[3][3]);
  /// Ionisation coefficient [cm-1]
  virtual bool ElectronTownsend(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz,
                                double& alpha);
  /// Attachment coefficient [cm-1]
  virtual bool ElectronAttachment(const double ex, const double ey,
                                  const double ez, const double bx,
                                  const double by, const double bz,
                                  double& eta);
  /// Lorentz angle
  virtual bool ElectronLorentzAngle(const double ex, const double ey,
                                    const double ez, const double bx,
                                    const double by, const double bz,
                                    double& lor);

  // Microscopic electron transport properties

  /// Dispersion relation (energy vs. wave vector)
  virtual double GetElectronEnergy(const double px, const double py,
                                   const double pz, double& vx, double& vy,
                                   double& vz, const int band = 0);
  /// Sample the momentum vector for a given energy 
  /// (only meaningful in semiconductors).
  virtual void GetElectronMomentum(const double e, double& px, double& py,
                                   double& pz, int& band);

  /// Null-collision rate [ns-1]
  virtual double GetElectronNullCollisionRate(const int band = 0);
  /// Collision rate [ns-1] for given electron energy
  virtual double GetElectronCollisionRate(const double e, const int band = 0);
  /// Sample the collision type. Update energy and direction vector.
  virtual bool GetElectronCollision(
      const double e, int& type, int& level, double& e1, double& dx, double& dy,
      double& dz, std::vector<std::pair<int, double> >& secondaries, int& ndxc,
      int& band);
  virtual unsigned int GetNumberOfDeexcitationProducts() const { return 0; }
  virtual bool GetDeexcitationProduct(const unsigned int i, double& t,
                                      double& s, int& type,
                                      double& energy) const;

  // Transport parameters for holes
  virtual bool HoleVelocity(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            double& vx, double& vy, double& vz);
  virtual bool HoleDiffusion(const double ex, const double ey, const double ez,
                             const double bx, const double by, const double bz,
                             double& dl, double& dt);
  virtual bool HoleDiffusion(const double ex, const double ey, const double ez,
                             const double bx, const double by, const double bz,
                             double cov[3][3]);
  virtual bool HoleTownsend(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            double& alpha);
  virtual bool HoleAttachment(const double ex, const double ey, const double ez,
                              const double bx, const double by, const double bz,
                              double& eta);

  // Transport parameters for ions
  virtual bool IonVelocity(const double ex, const double ey, const double ez,
                           const double bx, const double by, const double bz,
                           double& vx, double& vy, double& vz);
  virtual bool IonDiffusion(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            double& dl, double& dt);
  /// Dissociation coefficient
  virtual bool IonDissociation(const double ex, const double ey,
                               const double ez, const double bx,
                               const double by, const double bz, double& diss);

  /// Set the range of fields to be covered by the transport tables.
  void SetFieldGrid(double emin, double emax, const size_t ne, bool logE,
                    double bmin = 0., double bmax = 0., const size_t nb = 1,
                    double amin = HalfPi, double amax = HalfPi, 
                    const size_t na = 1);
  /// Set the fields and E-B angles to be used in the transport tables.
  void SetFieldGrid(const std::vector<double>& efields,
                    const std::vector<double>& bfields,
                    const std::vector<double>& angles);
  /// Get the fields and E-B angles used in the transport tables.
  void GetFieldGrid(std::vector<double>& efields, std::vector<double>& bfields,
                    std::vector<double>& angles);

  /// Get an entry in the table of drift speeds along E.
  bool GetElectronVelocityE(const unsigned int ie, const unsigned int ib,
                            const unsigned int ia, double& v) {
    return GetEntry(ie, ib, ia, "GetElectronVelocityE", m_eVelE, v);
  }
  /// Get an entry in the table of drift speeds along ExB.
  bool GetElectronVelocityExB(const unsigned int ie, const unsigned int ib,
                              const unsigned int ia, double& v) {
    return GetEntry(ie, ib, ia, "GetElectronVelocityExB", m_eVelX, v);
  }
  /// Get an entry in the table of drift speeds along Btrans.
  bool GetElectronVelocityB(const unsigned int ie, const unsigned int ib,
                            const unsigned int ia, double& v) {
    return GetEntry(ie, ib, ia, "GetElectronVelocityB", m_eVelB, v);
  }
  /// Get an entry in the table of longitudinal diffusion coefficients.
  bool GetElectronLongitudinalDiffusion(const unsigned int ie,
                                        const unsigned int ib,
                                        const unsigned int ia, double& dl) {
    return GetEntry(ie, ib, ia, "GetElectronLongitudinalDiffusion", 
                    m_eDifL, dl);
  }
  /// Get an entry in the table of transverse diffusion coefficients.
  bool GetElectronTransverseDiffusion(const unsigned int ie,
                                      const unsigned int ib,
                                      const unsigned int ia, double& dt) {
    return GetEntry(ie, ib, ia, "GetElectronTransverseDiffusion", 
                    m_eDifT, dt);
  }
  /// Get an entry in the table of Townsend coefficients.
  bool GetElectronTownsend(const unsigned int ie, const unsigned int ib,
                           const unsigned int ia, double& alpha) {
    return GetEntry(ie, ib, ia, "GetElectronTownsend", m_eAlp, alpha);
  }
  /// Get an entry in the table of attachment coefficients.
  bool GetElectronAttachment(const unsigned int ie, const unsigned int ib,
                             const unsigned int ia, double& eta) {
    return GetEntry(ie, ib, ia, "GetElectronAttachment", m_eAtt, eta);
  }
 
  /// Get an entry in the table of Lorentz angles.
  bool GetElectronLorentzAngle(const unsigned int ie, const unsigned int ib,
                               const unsigned int ia, double& lor) {
    return GetEntry(ie, ib, ia, "GetElectronLorentzAngle", m_eLor, lor);
  }

  bool GetHoleVelocityE(const unsigned int ie, const unsigned int ib,
                        const unsigned int ia, double& v) {
    return GetEntry(ie, ib, ia, "GetHoleVelocityE", m_hVelE, v);
  }
  bool GetHoleVelocityExB(const unsigned int ie, const unsigned int ib,
                          const unsigned int ia, double& v) {
    return GetEntry(ie, ib, ia, "GetHoleVelocityExB", m_hVelX, v);
  }
  bool GetHoleVelocityB(const unsigned int ie, const unsigned int ib,
                        const unsigned int ia, double& v) {
    return GetEntry(ie, ib, ia, "GetHoleVelocityB", m_hVelB, v);
  }
  bool GetHoleLongitudinalDiffusion(const unsigned int ie,
                                    const unsigned int ib,
                                    const unsigned int ia, double& dl) {
    return GetEntry(ie, ib, ia, "GetHoleLongitudinalDiffusion", m_hDifL, dl);
  }
  bool GetHoleTransverseDiffusion(const unsigned int ie, const unsigned int ib,
                                  const unsigned int ia, double& dt) {
    return GetEntry(ie, ib, ia, "GetHoleTransverseDiffusion", m_hDifT, dt);
  }
  bool GetHoleTownsend(const unsigned int ie, const unsigned int ib,
                       const unsigned int ia, double& alpha) {
    return GetEntry(ie, ib, ia, "GetHoleTownsend", m_hAlp, alpha);
  }
  bool GetHoleAttachment(const unsigned int ie, const unsigned int ib,
                         const unsigned int ia, double& eta) {
    return GetEntry(ie, ib, ia, "GetHoleAttachment", m_hAtt, eta);
  }

  bool GetIonMobility(const unsigned int ie, const unsigned int ib,
                      const unsigned int ia, double& mu) {
    return GetEntry(ie, ib, ia, "GetIonMobility", m_iMob, mu);
  }
  bool GetIonLongitudinalDiffusion(const unsigned int ie, const unsigned int ib,
                                   const unsigned int ia, double& dl) {
    return GetEntry(ie, ib, ia, "GetIonLongitudinalDiffusion", m_iDifL, dl);
  } 
  bool GetIonTransverseDiffusion(const unsigned int ie, const unsigned int ib,
                                 const unsigned int ia, double& dt) {
    return GetEntry(ie, ib, ia, "GetIonTransverseDiffusion", m_iDifT, dt);
  }
  bool GetIonDissociation(const unsigned int ie, const unsigned int ib,
                          const unsigned int ia, double& diss) {
    return GetEntry(ie, ib, ia, "GetIonDissociation", m_iDis, diss);
  }

  /// Reset all tables of transport parameters.
  virtual void ResetTables();

  void ResetElectronVelocity() {
    m_eVelE.clear();
    m_eVelB.clear();
    m_eVelX.clear();
  }
  void ResetElectronDiffusion() {
    m_eDifL.clear();
    m_eDifT.clear();
    m_eDifM.clear();
  }
  void ResetElectronTownsend() { m_eAlp.clear(); }
  void ResetElectronAttachment() { m_eAtt.clear(); }
  void ResetElectronLorentzAngle() { m_eLor.clear(); }

  void ResetHoleVelocity() {
    m_hVelE.clear();
    m_hVelB.clear();
    m_hVelX.clear();
  }
  void ResetHoleDiffusion() {
    m_hDifL.clear();
    m_hDifT.clear();
    m_hDifM.clear();
  }
  void ResetHoleTownsend() { m_hAlp.clear(); }
  void ResetHoleAttachment() { m_hAtt.clear(); }

  void ResetIonMobility() { m_iMob.clear(); }
  void ResetIonDiffusion() {
    m_iDifL.clear();
    m_iDifT.clear();
  }
  void ResetIonDissociation() { m_iDis.clear(); }

  bool SetIonMobility(const unsigned int ie, const unsigned int ib,
                      const unsigned int ia, const double mu);
  bool SetIonMobility(const std::vector<double>& fields,
                      const std::vector<double>& mobilities);

  /// Select the extrapolation method for fields below/above the table range.
  /// Possible options are "constant", "linear", and "exponential".
  void SetExtrapolationMethodVelocity(const std::string& extrLow,
                                      const std::string& extrHigh);
  void SetExtrapolationMethodDiffusion(const std::string& extrLow,
                                       const std::string& extrHigh);
  void SetExtrapolationMethodTownsend(const std::string& extrLow,
                                      const std::string& extrHigh);
  void SetExtrapolationMethodAttachment(const std::string& extrLow,
                                        const std::string& extrHigh);
  void SetExtrapolationMethodIonMobility(const std::string& extrLow,
                                         const std::string& extrHigh);
  void SetExtrapolationMethodIonDissociation(const std::string& extrLow,
                                             const std::string& extrHigh);

  /// Set the degree of polynomial interpolation (usually 2).
  void SetInterpolationMethodVelocity(const unsigned int intrp);
  void SetInterpolationMethodDiffusion(const unsigned int intrp);
  void SetInterpolationMethodTownsend(const unsigned int intrp);
  void SetInterpolationMethodAttachment(const unsigned int intrp);
  void SetInterpolationMethodIonMobility(const unsigned int intrp);
  void SetInterpolationMethodIonDissociation(const unsigned int intrp);

  // Scaling of fields and transport parameters.
  virtual double ScaleElectricField(const double e) const { return e; }
  virtual double UnScaleElectricField(const double e) const { return e; }
  virtual double ScaleVelocity(const double v) const { return v; }
  virtual double ScaleDiffusion(const double d) const { return d; }
  virtual double ScaleDiffusionTensor(const double d) const { return d; }
  virtual double ScaleTownsend(const double alpha) const { return alpha; }
  virtual double ScaleAttachment(const double eta) const { return eta; }
  virtual double ScaleLorentzAngle(const double lor) const { return lor; }
  virtual double ScaleDissociation(const double diss) const { return diss; }

  // Optical properties
  /// Get the energy range [eV] of the available optical data.
  virtual bool GetOpticalDataRange(double& emin, double& emax,
                                   const unsigned int i = 0);
  /// Get the complex dielectric function at a given energy.
  virtual bool GetDielectricFunction(const double e, double& eps1, double& eps2,
                                     const unsigned int i = 0);
  // Get the photoabsorption cross-section [cm2] at a given energy.
  virtual bool GetPhotoAbsorptionCrossSection(const double e, double& sigma,
                                              const unsigned int i = 0);
  virtual double GetPhotonCollisionRate(const double e);
  virtual bool GetPhotonCollision(const double e, int& type, int& level,
                                  double& e1, double& ctheta, int& nsec,
                                  double& esec);

  /// Switch on/off debugging  messages
  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }

 protected:
  std::string m_className = "Medium";

  static int m_idCounter;

  // Id number
  int m_id;
  // Name
  std::string m_name = "";
  // Temperature [K]
  double m_temperature = 293.15;
  // Pressure [Torr]
  double m_pressure = 760.;
  // Static dielectric constant
  double m_epsilon = 1.;
  // Number of components
  unsigned int m_nComponents = 1;
  // (Effective) atomic number Z
  double m_z = 1.;
  // Atomic weight A
  double m_a = 0.;
  // Number density [cm-3]
  double m_density = 0.;

  // Transport flags
  bool m_driftable = false;
  bool m_microscopic = false;
  bool m_ionisable = false;

  // W value
  double m_w = 0.;
  // Fano factor
  double m_fano = 0.;

  // Update flag
  bool m_isChanged = true;

  // Switch on/off debugging messages
  bool m_debug = false;

  // Field grids
  std::vector<double> m_eFields;
  std::vector<double> m_bFields;
  std::vector<double> m_bAngles;

  // Tables of transport parameters
  bool m_tab2d = false;
  // Electrons
  std::vector<std::vector<std::vector<double> > > m_eVelE;
  std::vector<std::vector<std::vector<double> > > m_eVelX;
  std::vector<std::vector<std::vector<double> > > m_eVelB;
  std::vector<std::vector<std::vector<double> > > m_eDifL;
  std::vector<std::vector<std::vector<double> > > m_eDifT;
  std::vector<std::vector<std::vector<double> > > m_eAlp;
  std::vector<std::vector<std::vector<double> > > m_eAtt;
  std::vector<std::vector<std::vector<double> > > m_eLor;

  std::vector<std::vector<std::vector<std::vector<double> > > > m_eDifM;

  // Holes
  std::vector<std::vector<std::vector<double> > > m_hVelE;
  std::vector<std::vector<std::vector<double> > > m_hVelX;
  std::vector<std::vector<std::vector<double> > > m_hVelB;
  std::vector<std::vector<std::vector<double> > > m_hDifL;
  std::vector<std::vector<std::vector<double> > > m_hDifT;
  std::vector<std::vector<std::vector<double> > > m_hAlp;
  std::vector<std::vector<std::vector<double> > > m_hAtt;

  std::vector<std::vector<std::vector<std::vector<double> > > > m_hDifM;

  // Ions
  std::vector<std::vector<std::vector<double> > > m_iMob;
  std::vector<std::vector<std::vector<double> > > m_iDifL;
  std::vector<std::vector<std::vector<double> > > m_iDifT;
  std::vector<std::vector<std::vector<double> > > m_iDis;

  // Thresholds for Townsend, attachment and dissociation coefficients.
  unsigned int m_eThrAlp = 0;
  unsigned int m_eThrAtt = 0;
  unsigned int m_hThrAlp = 0;
  unsigned int m_hThrAtt = 0;
  unsigned int m_iThrDis = 0;

  // Extrapolation methods (TODO: enum).
  std::pair<unsigned int, unsigned int> m_extrVel = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrDif = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrAlp = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrAtt = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrLor = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrMob = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrDis = {0, 1};

  // Interpolation methods
  unsigned int m_intpVel = 2;
  unsigned int m_intpDif = 2;
  unsigned int m_intpAlp = 2;
  unsigned int m_intpAtt = 2;
  unsigned int m_intpLor = 2;
  unsigned int m_intpMob = 2;
  unsigned int m_intpDis = 2;

  bool Velocity(const double ex, const double ey, const double ez,
                const double bx, const double by, const double bz,
                const std::vector<std::vector<std::vector<double> > >& velE,
                const std::vector<std::vector<std::vector<double> > >& velB,
                const std::vector<std::vector<std::vector<double> > >& velX,
                const double q, double& vx, double& vy, double& vz) const;
  bool Diffusion(const double ex, const double ey, const double ez,
                 const double bx, const double by, const double bz,
                 const std::vector<std::vector<std::vector<double> > >& difL,
                 const std::vector<std::vector<std::vector<double> > >& difT,
                 double& dl, double& dt) const;
  bool Diffusion(const double ex, const double ey, const double ez,
    const double bx, const double by, const double bz,
    const std::vector<std::vector<std::vector<std::vector<double> > > >& diff,
    double cov[3][3]) const;
  bool Alpha(const double ex, const double ey, const double ez,
             const double bx, const double by, const double bz,
             const std::vector<std::vector<std::vector<double> > >& tab,
             unsigned int intp, const unsigned int thr, 
             const std::pair<unsigned int, unsigned int>& extr, 
             double& alpha) const; 
  double GetAngle(const double ex, const double ey, const double ez,
                  const double bx, const double by, const double bz,
                  const double e, const double b) const;
  bool Interpolate(const double e, const double b, const double a,
                   const std::vector<std::vector<std::vector<double> > >& table,
                   double& y, const unsigned int intp,
                   const std::pair<unsigned int, unsigned int>& extr) const;

  double Interpolate1D(const double e, const std::vector<double>& table,
                       const std::vector<double>& fields,
                       const unsigned int intpMeth,
                       const std::pair<unsigned int, unsigned int>& extr) const;

  bool GetEntry(const unsigned int i, const unsigned int j, 
                const unsigned int k, const std::string& fcn, 
                const std::vector<std::vector<std::vector<double> > >& tab,
                double& val) const;

  void SetExtrapolationMethod(const std::string& low, const std::string& high,
                              std::pair<unsigned int, unsigned int>& extr,
                              const std::string& fcn);
  bool GetExtrapolationIndex(std::string str, unsigned int& nb) const;
  unsigned int SetThreshold(
      const std::vector<std::vector<std::vector<double> > >& tab) const;

  void Clone(std::vector<std::vector<std::vector<double> > >& tab,
             const std::vector<double>& efields,
             const std::vector<double>& bfields,
             const std::vector<double>& angles, const unsigned int intp,
             const std::pair<unsigned int, unsigned int>& extr,
             const double init, const std::string& label);
  void Clone(
      std::vector<std::vector<std::vector<std::vector<double> > > >& tab,
      const unsigned int n, const std::vector<double>& efields,
      const std::vector<double>& bfields, const std::vector<double>& angles,
      const unsigned int intp,
      const std::pair<unsigned int, unsigned int>& extr, const double init,
      const std::string& label);

  void Init(const size_t nE, const size_t nB, const size_t nA,
            std::vector<std::vector<std::vector<double> > >& tab,
            const double val);
  void Init(
      const size_t nE, const size_t nB, const size_t nA, const size_t nT,
      std::vector<std::vector<std::vector<std::vector<double> > > >& tab,
      const double val);
};
}

#endif
