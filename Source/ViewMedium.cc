#include <cmath>
#include <iostream>
#include <string>

#include <TAxis.h>
#include <TLatex.h>

#include "GarfieldConstants.hh"
#include "Medium.hh"
#include "Plotting.hh"
#include "ViewMedium.hh"

namespace {

int FindIndex(const std::vector<double>& fields, const double field,
              const double eps) {

  if (fields.empty()) return -1;
  const int n = fields.size();
  for (int i = 0; i < n; ++i) {
    const double sum = fabs(fields[i]) + fabs(field);
    const double tol = std::max(eps * sum, Garfield::Small);
    if (fabs(fields[i] - field) < tol) return i;
  }
  return -1;
} 

}

namespace Garfield {

ViewMedium::ViewMedium() { plottingEngine.SetDefaultStyle(); }

ViewMedium::~ViewMedium() {
  if (!m_hasExternalCanvas && m_canvas) delete m_canvas;
}

void ViewMedium::SetCanvas(TCanvas* c) {
  if (!c) {
    std::cerr << m_className << "::SetCanvas: Null pointer.\n";
    return;
  }
  if (!m_hasExternalCanvas && m_canvas) {
    delete m_canvas;
    m_canvas = nullptr;
  }
  m_canvas = c;
  m_hasExternalCanvas = true;
}

void ViewMedium::SetMedium(Medium* m) {
  if (!m) {
    std::cerr << m_className << "::SetMedium: Null pointer.\n";
    return;
  }

  m_medium = m;
}

void ViewMedium::SetRangeE(const double emin, const double emax, 
                           const bool logscale) {
  if (emin >= emax || emin < 0.) {
    std::cerr << m_className << "::SetRangeE: Incorrect range.\n";
    return;
  }

  m_eMin = emin;
  m_eMax = emax;
  m_logE = logscale;
}

void ViewMedium::SetRangeB(const double bmin, const double bmax, const bool logscale) {
  if (bmin >= bmax || bmin < 0.) {
    std::cerr << m_className << "::SetRangeB: Incorrect range.\n";
    return;
  }

  m_bMin = bmin;
  m_bMax = bmax;
  m_logB = logscale;
}

void ViewMedium::SetRangeA(const double amin, const double amax, const bool logscale) {
  if (amin >= amax || amin < 0.) {
    std::cerr << m_className << "::SetRangeA: Incorrect range.\n";
    return;
  }

  m_aMin = amin;
  m_aMax = amax;
  m_logA = logscale;
}

void ViewMedium::SetRangeY(const double ymin, const double ymax,
                           const bool logscale) {
  if (ymin >= ymax || ymin < 0.) {
    std::cerr << m_className << "::SetRangeY: Incorrect range.\n";
    return;
  }

  m_yMin = ymin;
  m_yMax = ymax;
  m_logY = logscale;
}

void ViewMedium::PlotElectronVelocity(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, ElectronVelocityE, xaxis);
  AddFunction(true, ElectronVelocityB, xaxis);
  AddFunction(true, ElectronVelocityExB, xaxis);
}

void ViewMedium::PlotHoleVelocity(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, HoleVelocityE, xaxis);
  AddFunction(true, HoleVelocityB, xaxis);
  AddFunction(true, HoleVelocityExB, xaxis);
}

void ViewMedium::PlotIonVelocity(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, IonVelocity, xaxis);
}

void ViewMedium::PlotElectronDiffusion(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, ElectronTransverseDiffusion, xaxis);
  AddFunction(true, ElectronLongitudinalDiffusion, xaxis);
}

void ViewMedium::PlotHoleDiffusion(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, HoleTransverseDiffusion, xaxis);
  AddFunction(true, HoleLongitudinalDiffusion, xaxis);
}

void ViewMedium::PlotIonDiffusion(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, IonTransverseDiffusion, xaxis);
  AddFunction(true, IonLongitudinalDiffusion, xaxis);
}

void ViewMedium::PlotElectronTownsend(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, ElectronTownsend, xaxis);
}

void ViewMedium::PlotHoleTownsend(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, HoleTownsend, xaxis);
}

void ViewMedium::PlotElectronAttachment(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, ElectronAttachment, xaxis);
}

void ViewMedium::PlotHoleAttachment(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, HoleAttachment, xaxis);
}

void ViewMedium::PlotElectronLorentzAngle(const char xaxis, const bool same) {
  SetupCanvas();
  AddFunction(same, ElectronLorentzAngle, xaxis);
}

void ViewMedium::PlotElectronCrossSections() {
  std::cerr << m_className << "::PlotElectronCrossSections: Not implemented.\n";
}

void ViewMedium::SetupCanvas() {
  if (!m_canvas) {
    m_canvas = new TCanvas();
    m_canvas->SetTitle("Medium View");
    if (m_hasExternalCanvas) m_hasExternalCanvas = false;
  }
  m_canvas->cd();
  gPad->SetLeftMargin(0.15);
}

void ViewMedium::AddFunction(const bool keep, const Property type,
                             const char xaxis) {
  // Make sure the medium is set.
  if (!m_medium) {
    std::cerr << m_className << "::AddFunction: Medium is not defined.\n";
    return;
  }

  // Look for an unused function name.
  int idx = 0;
  std::string fname = "fMediumView_0";
  while (gROOT->GetListOfFunctions()->FindObject(fname.c_str())) {
    ++idx;
    fname = "fMediumView_" + std::to_string(idx);
  }
  if (m_debug) {
    std::cout << m_className << "::AddFunction: Adding " << fname << "\n";
  }
  if (!keep) {
    m_functions.clear();
    for (auto graph : m_graphs) delete graph;
    m_graphs.clear();
    m_labels.clear();
  }

  // Get the field grid.
  std::vector<double> efields;
  std::vector<double> bfields;
  std::vector<double> bangles;
  m_medium->GetFieldGrid(efields, bfields, bangles);
  double xmin = 0., xmax = 0.;
  bool logx = false;
  GetAxisRangeX(efields, bfields, bangles, xaxis, xmin, xmax, logx);

  // Create a TF1 and add it to the list of functions.
  TF1 fcn(fname.c_str(), this, &ViewMedium::EvaluateFunction,
          xmin, xmax, 5, "ViewMedium", "EvaluateFunction");
  fcn.SetNpx(1000);
  const std::string xlabel = GetLabelX(xaxis);
  const std::string ylabel = GetLabelY(type);
  const std::string title = ";" + xlabel + ";" + ylabel;
  fcn.SetRange(xmin, xmax);
  if (!m_autoRangeY && (fabs(m_yMax - m_yMin) > 0.)) {
    fcn.SetMinimum(m_yMin);
    fcn.SetMaximum(m_yMax);
  }
  fcn.GetXaxis()->SetTitle(xlabel.c_str());
  fcn.GetXaxis()->SetTitleOffset(1.2);
  fcn.GetYaxis()->SetTitle(ylabel.c_str());
  fcn.SetTitle(title.c_str());
  fcn.SetParameter(0, type);
  fcn.SetParameter(1, xaxis);
  fcn.SetParameter(2, m_efield);
  fcn.SetParameter(3, m_bfield);
  fcn.SetParameter(4, m_angle);
  const auto color = GetColor(type);
  fcn.SetLineColor(color);
  m_functions.push_back(std::move(fcn));
  m_labels.emplace_back(std::make_pair(GetLegend(type), color));
  const unsigned int nE = efields.size();
  const unsigned int nB = bfields.size();
  const unsigned int nA = bangles.size();

  if (!efields.empty() && !bfields.empty() && !bangles.empty()) { 
    // Add a graph with the values at the grid points.
    int nPoints = xaxis == 'e' ? nE : xaxis == 'b' ? nB : nA;
    TGraph* graph = new TGraph(nPoints);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(color);
    bool ok = true;

    constexpr double eps = 1.e-3;
    int ie = FindIndex(efields, m_efield, eps);
    int ib = FindIndex(bfields, m_bfield, eps);
    int ia = FindIndex(bangles, m_angle, eps);
    if ((xaxis == 'e' && (ib < 0 || ia < 0)) ||
        (xaxis == 'b' && (ie < 0 || ia < 0)) ||
        (xaxis == 'a' && (ie < 0 || ib < 0))) {
      ok = false;
    }
    bool nonzero = false;
    for (int j = 0; j < nPoints; ++j) {
      double value = 0.;
      if (!ok) break;
      if (xaxis == 'e') {
        ie = j;
      } else if (xaxis == 'b') {
        ib = j;
      } else if (xaxis == 'a') {
        ia = j;
      } else {
        std::cerr << m_className << "::AddFunction: Unexpected axis.\n";
        break;
      }
      switch (type) {
        case ElectronVelocityE:
          ok = m_medium->GetElectronVelocityE(ie, ib, ia, value);
          value = m_medium->ScaleVelocity(value);
          break;
        case ElectronTransverseDiffusion:
          ok = m_medium->GetElectronTransverseDiffusion(ie, ib, ia, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case ElectronLongitudinalDiffusion:
          ok = m_medium->GetElectronLongitudinalDiffusion(ie, ib, ia, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case ElectronTownsend:
          ok = m_medium->GetElectronTownsend(ie, ib, ia, value);
          value = m_medium->ScaleTownsend(exp(value));
          break;
        case ElectronAttachment:
          ok = m_medium->GetElectronAttachment(ie, ib, ia, value);
          value = m_medium->ScaleAttachment(exp(value));
          break;
        case ElectronLorentzAngle:
          ok = m_medium->GetElectronLorentzAngle(ie, ib, ia, value);
          value = m_medium->ScaleLorentzAngle(value);
          break;
        case HoleVelocityE:
          ok = m_medium->GetHoleVelocityE(ie, ib, ia, value);
          value = m_medium->ScaleVelocity(value);
          break;
        case HoleTransverseDiffusion:
          ok = m_medium->GetHoleTransverseDiffusion(ie, ib, ia, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case HoleLongitudinalDiffusion:
          ok = m_medium->GetHoleLongitudinalDiffusion(ie, ib, ia, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case HoleTownsend:
          ok = m_medium->GetHoleTownsend(ie, ib, ia, value);
          value = m_medium->ScaleTownsend(exp(value));
          break;
        case HoleAttachment:
          ok = m_medium->GetHoleAttachment(ie, ib, ia, value);
          value = m_medium->ScaleAttachment(exp(value));
          break;
        case IonVelocity:
          ok = m_medium->GetIonMobility(ie, ib, ia, value);
          value *= m_medium->UnScaleElectricField(efields[ie]);
          break;
        case IonTransverseDiffusion:
          ok = m_medium->GetIonTransverseDiffusion(ie, ib, ia, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case IonLongitudinalDiffusion:
          ok = m_medium->GetIonLongitudinalDiffusion(ie, ib, ia, value);
          value = m_medium->ScaleDiffusion(value);
          break;
        case ElectronVelocityB:
          ok = m_medium->GetElectronVelocityB(ie, ib, ia, value);
          value = fabs(m_medium->ScaleVelocity(value));
          break;
        case ElectronVelocityExB:
          ok = m_medium->GetElectronVelocityExB(ie, ib, ia, value);
          value = fabs(m_medium->ScaleVelocity(value));
          break;
        case HoleVelocityB:
          ok = m_medium->GetHoleVelocityB(ie, ib, ia, value);
          value = fabs(m_medium->ScaleVelocity(value));
          break;
        case HoleVelocityExB:
          ok = m_medium->GetHoleVelocityExB(ie, ib, ia, value);
          value = fabs(m_medium->ScaleVelocity(value));
          break;
      }
      if (xaxis == 'e')
        graph->SetPoint(j, m_medium->UnScaleElectricField(efields[j]), value);
      else if (xaxis == 'b')
        graph->SetPoint(j, bfields[j], value);
      else if (xaxis == 'a')
        graph->SetPoint(j, bangles[j], value);
      if (fabs(value) > 1.e-10) nonzero = true;
    }
    if (ok && nonzero) {
      m_graphs.push_back(graph);
    } else {
      if (m_debug) {
        std::cerr << m_className << "::AddFunction:\n    Could not retrieve "
                  << "table for " << GetLegend(type) << ".\n";
      }
      delete graph;
    }
  }
  if (m_functions.empty()) return;
  TLatex label;
  if (keep && m_functions.size() > 1) {
    m_functions[0].GetYaxis()->SetTitleOffset(0);
    m_functions[0].Draw("");
    double ymin = m_functions[0].GetMinimum();
    double ymax = m_functions[0].GetMaximum();
    double xLabel = 1.3 * xmin;
    double yLabel = 0.9 * ymax;
    label.SetText(xLabel, yLabel, m_labels[0].first.c_str());
    xLabel += label.GetXsize(); 
    label.SetTextColor(m_labels[0].second);
    label.DrawLatex(xLabel, yLabel, m_labels[0].first.c_str());
    const unsigned int nFunctions = m_functions.size();
    for (unsigned int i = 1; i < nFunctions; ++i) {
      // See if the function has non-zero values.
      const double fmin = m_functions[i].GetMinimum();
      const double fmax = m_functions[i].GetMaximum();
      constexpr double tol = 1.e-10;
      if (fabs(fmin) < tol && fabs(fmax) < tol) continue;
      ymin = std::min(ymin, fmin);
      ymax = std::max(ymax, fmax);
      m_functions[i].Draw("lsame");
      yLabel -= 1.5 * label.GetYsize();
      label.SetTextColor(m_labels[i].second);
      label.DrawLatex(xLabel, yLabel, m_labels[i].first.c_str());
    }
    if (m_autoRangeY) {
      const double dy = ymax - ymin;
      m_functions[0].SetMinimum(std::max(0., ymin - 0.1 * dy));
      m_functions[0].SetMaximum(ymax + 0.1 * dy);
    }
  } else {
    m_functions.back().GetYaxis()->SetTitleOffset(0);
    m_functions.back().Draw("");
  }
  for (auto& graph : m_graphs) graph->Draw("p");
  if (logx) {
    m_canvas->SetLogx(1);
  } else {
    m_canvas->SetLogx(0);
  }
  if (m_logY && !m_autoRangeY) {
    m_canvas->SetLogy(1);
  } else {
    m_canvas->SetLogy(0);
  }
  m_canvas->Update();
}

double ViewMedium::EvaluateFunction(double* pos, double* par) {
  // to be modified to include B and angle

  if (!m_medium) return 0.;
  int type = int(par[0]);
  char xaxis = char(par[1]);
  const double x = pos[0];

  const double ctheta = xaxis == 'a' ? cos(x) : cos(par[4]);
  const double stheta = xaxis == 'a' ? sin(x) : sin(par[4]);
  double ex = xaxis == 'e' ? x : par[2];
  double ey = 0.;
  double bx = xaxis == 'b' ? x : par[3];
  double by = 0.;
  if (type == ElectronVelocityExB || type == HoleVelocityExB ||
             type == ElectronVelocityB || type == HoleVelocityB) {
    ex *= ctheta;
    ey = xaxis == 'e' ? x * stheta : par[2] * stheta; 
  } else {
    bx *= ctheta;
    by = xaxis == 'b' ? x * stheta : par[3] * stheta;
  } 
  // Return value.
  double y = 0.;
  // Auxiliary (dummy) variables.
  double s = 0., t = 0.;
  switch (type) {
    case ElectronVelocityE:
      if (!m_medium->ElectronVelocity(ex, 0, 0, bx, by, 0, y, s, t)) return 0.;
      y = fabs(y);
      break;
    case ElectronTransverseDiffusion:
      if (!m_medium->ElectronDiffusion(ex, 0, 0, bx, by, 0, s, y)) return 0.;
      break;
    case ElectronLongitudinalDiffusion:
      if (!m_medium->ElectronDiffusion(ex, 0, 0, bx, by, 0, y, s)) return 0.;
      break;
    case ElectronTownsend:
      if (!m_medium->ElectronTownsend(ex, 0, 0, bx, by, 0, y)) return 0.;
      break;
    case ElectronAttachment:
      if (!m_medium->ElectronAttachment(ex, 0, 0, bx, by, 0, y)) return 0.;
      break;
    case ElectronLorentzAngle:
      if (!m_medium->ElectronLorentzAngle(ex, 0, 0, bx, by, 0, y)) return 0.;
      break;
    case HoleVelocityE:
      if (!m_medium->HoleVelocity(ex, 0, 0, bx, by, 0, y, s, t)) return 0.;
      break;
    case HoleTransverseDiffusion:
      if (!m_medium->HoleDiffusion(ex, 0, 0, bx, by, 0, s, y)) return 0.;
      break;
    case HoleLongitudinalDiffusion:
      if (!m_medium->HoleDiffusion(ex, 0, 0, bx, by, 0, y, s)) return 0.;
      break;
    case HoleTownsend:
      if (!m_medium->HoleTownsend(ex, 0, 0, bx, by, 0, y)) return 0.;
      break;
    case HoleAttachment:
      if (!m_medium->HoleAttachment(ex, 0, 0, bx, by, 0, y)) return 0.;
      break;
    case IonVelocity:
      if (!m_medium->IonVelocity(ex, 0, 0, bx, by, 0, y, s, t)) return 0.;
      break;
    case IonTransverseDiffusion:
      if (!m_medium->IonDiffusion(ex, 0, 0, bx, by, 0, s, y)) return 0.;
      break;
    case IonLongitudinalDiffusion:
      if (!m_medium->IonDiffusion(ex, 0, 0, bx, by, 0, y, s)) return 0.;
      break;
    case ElectronVelocityB:
      if (!m_medium->ElectronVelocity(ex, ey, 0, bx, 0, 0, y, s, t)) return 0.;
      y = fabs(y);
      break;
    case ElectronVelocityExB:
      if (!m_medium->ElectronVelocity(ex, ey, 0, bx, 0, 0, s, t, y)) return 0.;
      y = fabs(y);
      break;
    case HoleVelocityB:
      if (!m_medium->HoleVelocity(ex, ey, 0, bx, 0, 0, y, s, t)) return 0.;
      y = fabs(y);
      break;
    case HoleVelocityExB:
      if (!m_medium->HoleVelocity(ex, ey, 0, bx, 0, 0, s, t, y)) return 0.;
      y = fabs(y);
      break;
    default:
      std::cerr << m_className << "::EvaluateFunction:\n    "
                << "Unknown type of transport coefficient requested. Bug!\n";
      return 0.;
  }

  return y;
}

void ViewMedium::GetAxisRangeX(const std::vector<double>& efields, 
    const std::vector<double>& bfields, const std::vector<double>& angles,
    const char xaxis, double& xmin, double& xmax, bool& logx) const {

  if (m_autoRangeX) {
    bool ok = false;
    if (xaxis == 'e' && !efields.empty()) {
      xmin = efields.front();
      xmax = efields.back();
      const double dx = xmax - xmin; 
      if (dx > Small && xmax > Small) {
        xmax += 0.1 * dx;
        if (xmin > Small) {
          xmin *= 0.8;
          logx = true;
        } else {
          logx = false;
        }
        ok = true; 
      }
    } else if (xaxis == 'b' && !bfields.empty()) {
      logx = false;
      xmin = bfields.front();
      xmax = bfields.back();
      const double dx = xmax - xmin; 
      if (dx > Small && xmax > Small) {
        if (xmin > Small) xmin *= 0.8;
        xmax += 0.1 * dx;
        ok = true; 
      }
    } else if (xaxis == 'a' && !angles.empty()) {
      logx = false;
      xmin = angles.front();
      xmax = angles.back();
      const double dx = xmax - xmin; 
      if (dx > Small && xmax > Small) {
        ok = true; 
      }
    }
    if (ok) return;
  }
  if (xaxis == 'e') {
    xmin = m_eMin;
    xmax = m_eMax;
    logx = m_logE;
  } else if (xaxis == 'b') {
    xmin = m_bMin;
    xmax = m_bMax;
    logx = m_logB;
  } else if (xaxis == 'a') {
    xmin = m_aMin;
    xmax = m_aMax;
    logx = m_logA;
  }
}

int ViewMedium::GetColor(const Property prop) const {
  if (prop == ElectronLongitudinalDiffusion || prop == ElectronAttachment ||
      prop == ElectronLorentzAngle) {
    return plottingEngine.GetRootColorLine1();
  } else if (prop == HoleLongitudinalDiffusion || prop == HoleAttachment ||
             prop == IonLongitudinalDiffusion) {
    return plottingEngine.GetRootColorLine2();
  } else if (prop < HoleVelocityE) {
    return plottingEngine.GetRootColorElectron();
  } else if (prop == ElectronVelocityB || prop == HoleVelocityB) {
    return kGreen + 2;
  } else if (prop == ElectronVelocityExB || prop == HoleVelocityExB) {
    return kRed + 1;
  } else if (prop < IonVelocity) {
    return plottingEngine.GetRootColorHole();
  }
  return plottingEngine.GetRootColorIon();
}

std::string ViewMedium::GetLabelX(const char xaxis) const {

  if (xaxis == 'e') {
    return "electric field [V/cm]";
  } else if (xaxis == 'b') {
    return "magnetic field [T]";
  } else if (xaxis == 'a') {
    return "angle between #bf{E} and #bf{B} [rad]";
  }
  return "";
}

std::string ViewMedium::GetLabelY(const Property prop) const {

  switch (prop) {
    case ElectronVelocityE:
    case ElectronVelocityB:
    case ElectronVelocityExB:
    case HoleVelocityE:
    case HoleVelocityB:
    case HoleVelocityExB:
    case IonVelocity:
      return "drift velocity [cm/ns]";
      break;
    case ElectronLongitudinalDiffusion:
    case ElectronTransverseDiffusion:
    case HoleLongitudinalDiffusion: 
    case HoleTransverseDiffusion:
    case IonLongitudinalDiffusion: 
    case IonTransverseDiffusion:
      return "diffusion coefficient [#sqrt{cm}]";
      break;
    case ElectronTownsend:
    case HoleTownsend:
      return "#alpha [1/cm]";
      break;
    case ElectronAttachment:
    case HoleAttachment:
      return "#eta [1/cm]";
      break;
    case ElectronLorentzAngle:
      return "Angle between #bf{v} and #bf{E} [rad]";
      break;
    default:
      return "unknown";
  }
  return "";
}

std::string ViewMedium::GetLegend(const Property prop) const {
  
  std::string label = "";
  switch (prop) {
    case ElectronVelocityE:
    case HoleVelocityE:
      label += "velocity #parallel #bf{E}";
      break;
    case ElectronVelocityB:
    case HoleVelocityB:
      label += "velocity #parallel #bf{B}";
      break;
    case ElectronVelocityExB:
    case HoleVelocityExB:
      label += "velocity #parallel #bf{E}#times #bf{B}";
      break;
    case IonVelocity:
      label += "ion velocity";
      break;
    case ElectronLongitudinalDiffusion:
    case HoleLongitudinalDiffusion:
    case IonLongitudinalDiffusion: 
      label += "longitudinal diffusion";
      break; 
    case ElectronTransverseDiffusion:
    case HoleTransverseDiffusion:
    case IonTransverseDiffusion:
      label += "transverse diffusion";
      break; 
    case ElectronTownsend:
    case HoleTownsend:
      label += "Townsend";
      break;
    case ElectronAttachment:
    case HoleAttachment:
      label += "Attachment";
      break;
    case ElectronLorentzAngle:
      label += "Lorenz angle";
      break;
  }
  return label;
}

}
