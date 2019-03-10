#ifndef G_MEDIUM_H
#define G_MEDIUM_H

#include <string>
#include <vector>

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
                    double amin = 0., double amax = 0., const size_t na = 1);
  /// Set the fields and E-B angles to be used in the transport tables.
  void SetFieldGrid(const std::vector<double>& efields,
                    const std::vector<double>& bfields,
                    const std::vector<double>& angles);
  /// Get the fields and E-B angles used in the transport tables.
  void GetFieldGrid(std::vector<double>& efields, std::vector<double>& bfields,
                    std::vector<double>& angles);

  /// Get an entry in the table of drift speeds along E.
  bool GetElectronVelocityE(const unsigned int ie, const unsigned int ib,
                            const unsigned int ia, double& v);
  /// Get an entry in the table of drift speeds along ExB.
  bool GetElectronVelocityExB(const unsigned int ie, const unsigned int ib,
                              const unsigned int ia, double& v);
  /// Get an entry in the table of drift speeds along Btrans.
  bool GetElectronVelocityB(const unsigned int ie, const unsigned int ib,
                            const unsigned int ia, double& v);

  /// Get an entry in the table of longitudinal diffusion coefficients.
  bool GetElectronLongitudinalDiffusion(const unsigned int ie,
                                        const unsigned int ib,
                                        const unsigned int ia, double& dl);
  /// Get an entry in the table of transverse diffusion coefficients.
  bool GetElectronTransverseDiffusion(const unsigned int ie,
                                      const unsigned int ib,
                                      const unsigned int ia, double& dt);
  /// Get an entry in the table of Townsend coefficients.
  bool GetElectronTownsend(const unsigned int ie, const unsigned int ib,
                           const unsigned int ia, double& alpha);
  /// Get an entry in the table of attachment coefficients.
  bool GetElectronAttachment(const unsigned int ie, const unsigned int ib,
                             const unsigned int ia, double& eta);
  /// Get an entry in the table of Lorentz angles.
  bool GetElectronLorentzAngle(const unsigned int ie, const unsigned int ib,
                               const unsigned int ia, double& lor);

  bool GetHoleVelocityE(const unsigned int ie, const unsigned int ib,
                        const unsigned int ia, double& v);
  bool GetHoleVelocityExB(const unsigned int ie, const unsigned int ib,
                          const unsigned int ia, double& v);
  bool GetHoleVelocityB(const unsigned int ie, const unsigned int ib,
                        const unsigned int ia, double& v);
  bool GetHoleLongitudinalDiffusion(const unsigned int ie,
                                    const unsigned int ib,
                                    const unsigned int ia, double& dl);
  bool GetHoleTransverseDiffusion(const unsigned int ie, const unsigned int ib,
                                  const unsigned int ia, double& dt);
  bool GetHoleTownsend(const unsigned int ie, const unsigned int ib,
                       const unsigned int ia, double& alpha);
  bool GetHoleAttachment(const unsigned int ie, const unsigned int ib,
                         const unsigned int ia, double& eta);

  bool GetIonMobility(const unsigned int ie, const unsigned int ib,
                      const unsigned int ia, double& mu);
  bool GetIonLongitudinalDiffusion(const unsigned int ie, const unsigned int ib,
                                   const unsigned int ia, double& dl);
  bool GetIonTransverseDiffusion(const unsigned int ie, const unsigned int ib,
                                 const unsigned int ia, double& dt);
  bool GetIonDissociation(const unsigned int ie, const unsigned int ib,
                          const unsigned int ia, double& diss);

  /// Reset all tables of transport parameters.
  virtual void ResetTables();

  void ResetElectronVelocity() {
    m_eVelocityE.clear();
    m_eVelocityB.clear();
    m_eVelocityExB.clear();
  }
  void ResetElectronDiffusion() {
    m_eDiffLong.clear();
    m_eDiffTrans.clear();
    m_eDiffTens.clear();
  }
  void ResetElectronTownsend() { m_eTownsend.clear(); }
  void ResetElectronAttachment() { m_eAttachment.clear(); }
  void ResetElectronLorentzAngle() { m_eLorentzAngle.clear(); }

  void ResetHoleVelocity() {
    m_hVelocityE.clear();
    m_hVelocityB.clear();
    m_hVelocityExB.clear();
  }
  void ResetHoleDiffusion() {
    m_hDiffLong.clear();
    m_hDiffTrans.clear();
    m_hDiffTens.clear();
  }
  void ResetHoleTownsend() { m_hTownsend.clear(); }
  void ResetHoleAttachment() { m_hAttachment.clear(); }

  void ResetIonMobility() { m_ionMobility.clear(); }
  void ResetIonDiffusion() {
    m_ionDiffLong.clear();
    m_ionDiffTrans.clear();
  }
  void ResetIonDissociation() { m_ionDissociation.clear(); }

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
  bool m_map2d = false;
  // Electrons
  std::vector<std::vector<std::vector<double> > > m_eVelocityE;
  std::vector<std::vector<std::vector<double> > > m_eVelocityExB;
  std::vector<std::vector<std::vector<double> > > m_eVelocityB;
  std::vector<std::vector<std::vector<double> > > m_eDiffLong;
  std::vector<std::vector<std::vector<double> > > m_eDiffTrans;
  std::vector<std::vector<std::vector<double> > > m_eTownsend;
  std::vector<std::vector<std::vector<double> > > m_eAttachment;
  std::vector<std::vector<std::vector<double> > > m_eLorentzAngle;

  std::vector<std::vector<std::vector<std::vector<double> > > > m_eDiffTens;

  // Holes
  std::vector<std::vector<std::vector<double> > > m_hVelocityE;
  std::vector<std::vector<std::vector<double> > > m_hVelocityExB;
  std::vector<std::vector<std::vector<double> > > m_hVelocityB;
  std::vector<std::vector<std::vector<double> > > m_hDiffLong;
  std::vector<std::vector<std::vector<double> > > m_hDiffTrans;
  std::vector<std::vector<std::vector<double> > > m_hTownsend;
  std::vector<std::vector<std::vector<double> > > m_hAttachment;

  std::vector<std::vector<std::vector<std::vector<double> > > > m_hDiffTens;

  // Ions
  std::vector<std::vector<std::vector<double> > > m_ionMobility;
  std::vector<std::vector<std::vector<double> > > m_ionDiffLong;
  std::vector<std::vector<std::vector<double> > > m_ionDiffTrans;
  std::vector<std::vector<std::vector<double> > > m_ionDissociation;

  // Thresholds for Townsend, attachment and dissociation coefficients.
  int thrElectronTownsend = 0;
  int thrElectronAttachment = 0;

  int thrHoleTownsend = 0;
  int thrHoleAttachment = 0;
  int thrIonDissociation = 0;

  // Extrapolation methods (TODO: enum).
  std::pair<unsigned int, unsigned int> m_extrVel = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrDiff = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrTownsend = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrAttachment = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrLorentzAngle = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrMobility = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrDissociation = {0, 1};

  // Interpolation methods
  unsigned int m_intpVel = 2;
  unsigned int m_intpDiff = 2;
  unsigned int m_intpTownsend = 2;
  unsigned int m_intpAttachment = 2;
  unsigned int m_intpLorentzAngle = 2;
  unsigned int m_intpMobility = 2;
  unsigned int m_intpDissociation = 2;

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
  void SetExtrapolationMethod(const std::string& low, const std::string& high,
                              std::pair<unsigned int, unsigned int>& extr,
                              const std::string& fcn);
  bool GetExtrapolationIndex(std::string str, unsigned int& nb) const;

  void CloneTable(std::vector<std::vector<std::vector<double> > >& tab,
                  const std::vector<double>& efields,
                  const std::vector<double>& bfields,
                  const std::vector<double>& angles, const unsigned int intp,
                  const std::pair<unsigned int, unsigned int>& extr,
                  const double init, const std::string& label);
  void CloneTensor(
      std::vector<std::vector<std::vector<std::vector<double> > > >& tab,
      const unsigned int n, const std::vector<double>& efields,
      const std::vector<double>& bfields, const std::vector<double>& angles,
      const unsigned int intp,
      const std::pair<unsigned int, unsigned int>& extr, const double init,
      const std::string& label);

  void InitTable(const size_t nE, const size_t nB, const size_t nA,
                 std::vector<std::vector<std::vector<double> > >& tab,
                 const double val);
  void InitTensor(
      const size_t nE, const size_t nB, const size_t nA, const size_t nT,
      std::vector<std::vector<std::vector<std::vector<double> > > >& tab,
      const double val);
};
}

#endif
