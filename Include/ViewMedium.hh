#ifndef G_VIEW_MEDIUM
#define G_VIEW_MEDIUM

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>

#include "FundamentalConstants.hh"

namespace Garfield {

class Medium;

/// Plot transport coefficients as function of electric and magnetic field.

class ViewMedium {
 public:
  /// Constructor
  ViewMedium();
  /// Destructor
  ~ViewMedium();

  /// Set the canvas to be painted on.
  void SetCanvas(TCanvas* c);
  // Set the medium from which to retrieve the transport coefficients.
  void SetMedium(Medium* m);

  /// Try to choose the x-axis range based on the field grid. 
  void EnableAutoRangeX(const bool on = true) { m_autoRangeX = on; } 
  /// Set the limits of the electric field.
  void SetRangeE(const double emin, const double emax, const bool logscale);
  /// Set the limits of the magnetic field.
  void SetRangeB(const double bmin, const double bmax, const bool logscale);
  /// Set the limits of the angle between electric and magnetic field.
  void SetRangeA(const double amin, const double amax, const bool logscale);
  /// Choose the y-axis range based on the function's minima/maxima.
  void EnableAutoRangeY(const bool on = true) { m_autoRangeY = on; }
  /// Set the range of the function (velocity etc.) to be plotted.
  void SetRangeY(const double ymin, const double ymax, 
                 const bool logscale = false);

  /// Set the electric field to use when plotting as function of B or angle.
  void SetElectricField(const double efield) { m_efield = efield; }
  /// Set the magnetic field to use when plotting as function of E or angle.
  void SetMagneticField(const double bfield) { m_bfield = bfield; }
  /// Set the angle to use when plotting as function of E or B.
  void SetAngle(const double angle) { m_angle = angle; }

  void PlotElectronVelocity(const char xaxis, const bool same = false);
  void PlotHoleVelocity(const char xaxis, const bool same = false);
  void PlotIonVelocity(const char xaxis, const bool same = false);
  void PlotElectronDiffusion(const char xaxis, const bool same = false);
  void PlotHoleDiffusion(const char xaxis, const bool same = false);
  void PlotIonDiffusion(const char xaxis, const bool same = false);
  void PlotElectronTownsend(const char xaxis, const bool same = false);
  void PlotHoleTownsend(const char xaxis, const bool same = false);
  void PlotElectronAttachment(const char xaxis, const bool same = false);
  void PlotHoleAttachment(const char xaxis, const bool same = false);
  void PlotElectronLorentzAngle(const char xaxis, const bool same = false);
  void PlotElectronCrossSections();

  double EvaluateFunction(double* pos, double* par);

  enum Property {
    ElectronVelocityE,
    ElectronTransverseDiffusion,
    ElectronLongitudinalDiffusion,
    ElectronTownsend,
    ElectronAttachment,
    ElectronLorentzAngle,
    HoleVelocityE = 10,
    HoleTransverseDiffusion,
    HoleLongitudinalDiffusion,
    HoleTownsend,
    HoleAttachment,
    IonVelocity = 20,
    IonTransverseDiffusion,
    IonLongitudinalDiffusion,
    ElectronVelocityB,
    ElectronVelocityExB,
    HoleVelocityB,
    HoleVelocityExB
  };

 private:
  std::string m_className = "ViewMedium";

  // Options
  bool m_debug = false;

  // Canvas
  TCanvas* m_canvas = nullptr;
  bool m_hasExternalCanvas = false;

  Medium* m_medium = nullptr;

  // X-axis ranges
  double m_eMin = 100., m_eMax = 100000.;
  double m_bMin = 0., m_bMax = 2.;
  double m_aMin = 0., m_aMax = Pi;
  bool m_logE = true;
  bool m_logB = false;
  bool m_logA = false;
  bool m_autoRangeX = true;
  // Y-axis range
  double m_yMin = 0., m_yMax = 1.;
  bool m_logY = false;
  bool m_autoRangeY = true;

  // E-field to use when plotting as function of B-field or angle.
  double m_efield = 1000.;
  // B-field to use when plotting as function of E-field or angle.
  double m_bfield = 0.;
  // Angle to use when plotting as function of E-field or B-field.
  double m_angle = HalfPi;

  // Functions
  std::vector<TF1> m_functions;
  // Graphs
  std::vector<TGraph*> m_graphs;
  // Labels
  std::vector<std::pair<std::string, int> > m_labels;

  void SetupCanvas();
  void AddFunction(const bool keep, const Property type, const char xaxis);
  void GetAxisRangeX(const std::vector<double>& efields,
                     const std::vector<double>& bfields, 
                     const std::vector<double>& angles, const char xaxis,
                     double& xmin, double& max, bool& logx) const;
  int GetColor(const Property property) const;
  std::string GetLabelX(const char xaxis) const;
  std::string GetLabelY(const Property property) const;
  std::string GetLegend(const Property property) const;
};
}
#endif
