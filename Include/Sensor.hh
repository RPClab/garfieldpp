#ifndef G_SENSOR_H
#define G_SENSOR_H

#include <vector>

#include "ComponentBase.hh"

namespace Garfield {

/// %Sensor

class Sensor {

 public:
  /// Constructor
  Sensor() = default;
  /// Destructor
  ~Sensor() {}

  /// Add a component.
  void AddComponent(ComponentBase* comp);
  /// Get the number of components attached to the sensor.
  unsigned int GetNumberOfComponents() const { return m_components.size(); }
  /// Retrieve the pointer to a given component.
  ComponentBase* GetComponent(const unsigned int i);

  /// Add an electrode.
  void AddElectrode(ComponentBase* comp, const std::string& label);
  /// Get the number of electrodes attached to the sensor.
  unsigned int GetNumberOfElectrodes() const { return m_electrodes.size(); }
  /// Remove all components, electrodes and reset the sensor.
  void Clear();

  /// Get the drift field and potential at (x, y, z).
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& medium,
                     int& status);
  /// Get the drift field at (x, y, z).
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& medium, int& status);

  /// Get the magnetic field at (x, y, z).
  void MagneticField(const double x, const double y, const double z, double& bx,
                     double& by, double& bz, int& status);

  /// Get the weighting field at (x, y, z).
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label);
  /// Get the weighting potential at (x, y, z).
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label);

  /// Get the medium at (x, y, z).
  bool GetMedium(const double x, const double y, const double z,
                 Medium*& medium);

  /// Set the user area to the default.
  bool SetArea();
  /// Set the user area explicitly.
  bool SetArea(const double xmin, const double ymin, const double zmin,
               const double xmax, const double ymax, const double zmax);
  /// Return the current user area.
  bool GetArea(double& xmin, double& ymin, double& zmin, double& xmax,
               double& ymax, double& zmax);
  /// Check if a point is inside the user area.
  bool IsInArea(const double x, const double y, const double z);

  /// Return the voltage range.
  bool GetVoltageRange(double& vmin, double& vmax);

  /// Start a new event, when computing the average signal over multiple events.
  void NewSignal() { ++m_nEvents; }
  /// Reset signals and induced charges of all electrodes.
  void ClearSignal();

  void AddSignal(const double q, const double t, const double dt,
                 const double x, const double y, const double z,
                 const double vx, const double vy, const double vz);
  void AddInducedCharge(const double q, const double x0, const double y0,
                        const double z0, const double x1, const double y1,
                        const double z1);

  /** Set the time window and binning for the signal calculation.
    * \param tstart start time [ns]
    * \param tstep bin width [ns]
    * \param nstep number of bins
    */ 
  void SetTimeWindow(const double tstart, const double tstep, 
                     const unsigned int nsteps);
  /// Retrieve the time window and binning.
  void GetTimeWindow(double& tstart, double& tstep, unsigned int& nsteps) const {
    tstart = m_tStart;
    tstep = m_tStep;
    nsteps = m_nTimeBins;
  }
  /// Retrieve the total signal for a given electrode and time bin. 
  double GetSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the electron signal for a given electrode and time bin.
  double GetElectronSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the ion or hole signal for a given electrode and time bin.
  double GetIonSignal(const std::string& label, const unsigned int bin);
  /// Retrieve the total induced charge for a given electrode, 
  /// calculated using the weighting potentials at the start and end points. 
  double GetInducedCharge(const std::string& label);
  /// Set the function to be used for evaluating the transfer function.
  void SetTransferFunction(double (*f)(double t));
  /// Set the points to be used for interpolating the transfer function. 
  void SetTransferFunction(const std::vector<double>& times,
                           const std::vector<double>& values);
  /// Evaluate the transfer function at a given time.
  double GetTransferFunction(const double t);
  /// Convolute the induced current with the transfer function.
  bool ConvoluteSignal();
  /// Replace the current signal curve by its integral.
  bool IntegrateSignal();
  /// Set the function to be used for evaluating the noise component.
  void SetNoiseFunction(double (*f)(double t));
  /// Add noise to the induced signal.
  void AddNoise(const bool total = true, const bool electron = false, 
                const bool ion = false);
  /** Determine the threshold crossings of the current signal curve. 
    * \param thr threshold value
    * \param label electrode for which to compute the threshold crossings
    * \param n number of threshold crossings
    */ 
  bool ComputeThresholdCrossings(const double thr, const std::string& label,
                                 int& n);
  /// Get the number of threshold crossings 
  /// (after having called ComputeThresholdCrossings).
  unsigned int GetNumberOfThresholdCrossings() const { 
    return m_thresholdCrossings.size(); 
  }
  /** Retrieve the time and type of a given threshold crossing (after having 
    * called ComputeThresholdCrossings.
    * \param i index
    * \param time threshold crossing time [ns]
    * \param level threshold (should correspond to the value requested).
    * \param rise flag whether the crossing is on a rising or falling slope.
    */
  bool GetThresholdCrossing(const unsigned int i, double& time, double& level,
                            bool& rise) const;

  /// Switch debugging messages on/off.
  void EnableDebugging(const bool on = true) { m_debug = on; }

  bool IsWireCrossed(const double x0, const double y0, const double z0,
                     const double x1, const double y1, const double z1,
                     double& xc, double& yc, double& zc);
  bool IsInTrapRadius(const double q0, const double x0, const double y0, 
                      const double z0, double& xw, double& yw, double& rw);

 private:
  std::string m_className = "Sensor";

  // Components
  std::vector<ComponentBase*> m_components;
  ComponentBase* m_lastComponent = nullptr;

  // Electrodes
  struct Electrode {
    ComponentBase* comp;
    std::string label;
    std::vector<double> signal;
    std::vector<double> electronsignal;
    std::vector<double> ionsignal;
    double charge;
  };
  std::vector<Electrode> m_electrodes;

  // Time window for signals
  unsigned int m_nTimeBins = 200;
  double m_tStart = 0.;
  double m_tStep = 10.;
  unsigned int m_nEvents = 0;
  static double m_signalConversion;

  // Transfer function
  bool m_hasTransferFunction = false;
  double (*m_fTransfer)(double t) = nullptr;
  std::vector<double> m_transferFunctionTimes;
  std::vector<double> m_transferFunctionValues;

  // Noise
  bool m_hasNoiseFunction = false;
  double (*m_fNoise)(double t) = nullptr;

  struct ThresholdCrossing {
    double time;
    bool rise;
  };
  std::vector<ThresholdCrossing> m_thresholdCrossings;
  double m_thresholdLevel = 0.;

  // User bounding box
  bool m_hasUserArea = false;
  double m_xMinUser = 0., m_yMinUser = 0., m_zMinUser = 0.;
  double m_xMaxUser = 0., m_yMaxUser = 0., m_zMaxUser = 0.;

  // Switch on/off debugging messages
  bool m_debug = false;

  // Return the current sensor size
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax);

  double InterpolateTransferFunctionTable(const double t) const;
};
}

#endif
