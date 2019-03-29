#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>

#include <TMath.h>

#include "FundamentalConstants.hh"
#include "GarfieldConstants.hh"
#include "MagboltzInterface.hh"
#include "MediumMagboltz.hh"
#include "OpticalData.hh"
#include "Random.hh"

namespace {

void PrintErrorMixer(const std::string& fcn) {
  std::cerr << fcn << ": Error calculating the collision rates table.\n";
}

std::string GetDescription(const unsigned int index,
                           char scrpt[][Garfield::Magboltz::nCharDescr]) {
  return std::string(scrpt[index],
                     scrpt[index] + Garfield::Magboltz::nCharDescr);
}
}

namespace Garfield {

const int MediumMagboltz::DxcTypeRad = 0;
const int MediumMagboltz::DxcTypeCollIon = 1;
const int MediumMagboltz::DxcTypeCollNonIon = -1;

MediumMagboltz::MediumMagboltz()
    : MediumGas(),
      m_eFinal(40.),
      m_eStep(m_eFinal / Magboltz::nEnergySteps),
      m_eHigh(1.e4),
      m_eHighLog(log(m_eHigh)),
      m_lnStep(1.),
      m_eFinalGamma(20.),
      m_eStepGamma(m_eFinalGamma / nEnergyStepsGamma) {
  m_className = "MediumMagboltz";

  // Set physical constants in Magboltz common blocks.
  Magboltz::cnsts_.echarg = ElementaryCharge * 1.e-15;
  Magboltz::cnsts_.emass = ElectronMassGramme;
  Magboltz::cnsts_.amu = AtomicMassUnit;
  Magboltz::cnsts_.pir2 = BohrRadius * BohrRadius * Pi;
  Magboltz::inpt_.ary = RydbergEnergy;

  // Set parameters in Magboltz common blocks.
  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = Magboltz::nEnergySteps;
  // Select the scattering model.
  Magboltz::inpt_.nAniso = 2;
  // Max. energy [eV]
  Magboltz::inpt_.efinal = m_eFinal;
  // Energy step size [eV]
  Magboltz::inpt_.estep = m_eStep;
  // Temperature and pressure
  Magboltz::inpt_.akt = BoltzmannConstant * m_temperature;
  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;
  // Disable Penning transfer.
  Magboltz::inpt_.ipen = 0;

  m_description.assign(Magboltz::nMaxLevels,
                       std::string(Magboltz::nCharDescr, ' '));

  m_cfTot.assign(Magboltz::nEnergySteps, 0.);
  m_cfTotLog.assign(nEnergyStepsLog, 0.);
  m_cf.assign(Magboltz::nEnergySteps,
              std::vector<double>(Magboltz::nMaxLevels, 0.));
  m_cfLog.assign(nEnergyStepsLog,
                 std::vector<double>(Magboltz::nMaxLevels, 0.));

  m_isChanged = true;

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = true;

  m_scaleExc.fill(1.);
}

bool MediumMagboltz::SetMaxElectronEnergy(const double e) {
  if (e <= Small) {
    std::cerr << m_className << "::SetMaxElectronEnergy: Invalid energy.\n";
    return false;
  }
  m_eFinal = e;

  // Determine the energy interval size.
  m_eStep = std::min(m_eFinal, m_eHigh) / Magboltz::nEnergySteps;

  // Set max. energy and step size also in Magboltz common block.
  Magboltz::inpt_.efinal = m_eFinal;
  Magboltz::inpt_.estep = m_eStep;

  // Force recalculation of the scattering rates table.
  m_isChanged = true;

  return true;
}

bool MediumMagboltz::SetMaxPhotonEnergy(const double e) {
  if (e <= Small) {
    std::cerr << m_className << "::SetMaxPhotonEnergy: Invalid energy.\n";
    return false;
  }
  m_eFinalGamma = e;

  // Determine the energy interval size.
  m_eStepGamma = m_eFinalGamma / nEnergyStepsGamma;

  // Force recalculation of the scattering rates table.
  m_isChanged = true;

  return true;
}

void MediumMagboltz::SetSplittingFunctionOpalBeaty() {
  m_useOpalBeaty = true;
  m_useGreenSawada = false;
}

void MediumMagboltz::SetSplittingFunctionGreenSawada() {
  m_useOpalBeaty = false;
  m_useGreenSawada = true;
  if (m_isChanged) return;

  bool allset = true;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (!m_hasGreenSawada[i]) {
      if (allset) {
        std::cout << m_className << "::SetSplittingFunctionGreenSawada:\n";
        allset = false;
      }
      std::cout << "    Fit parameters for " << m_gas[i] << " not available.\n"
                << "    Using Opal-Beaty formula instead.\n";
    }
  }
}

void MediumMagboltz::SetSplittingFunctionFlat() {
  m_useOpalBeaty = false;
  m_useGreenSawada = false;
}

void MediumMagboltz::EnableDeexcitation() {
  if (m_usePenning) {
    std::cout << m_className << "::EnableDeexcitation:\n"
              << "    Penning transfer will be switched off.\n";
  }
  // if (m_useRadTrap) {
  //   std::cout << "    Radiation trapping is switched on.\n";
  // } else {
  //   std::cout << "    Radiation trapping is switched off.\n";
  // }
  m_usePenning = false;
  m_useDeexcitation = true;
  m_isChanged = true;
  m_dxcProducts.clear();
}

void MediumMagboltz::EnableRadiationTrapping() {
  m_useRadTrap = true;
  if (!m_useDeexcitation) {
    std::cout << m_className << "::EnableRadiationTrapping:\n    "
              << "Radiation trapping is enabled but de-excitation is not.\n";
  } else {
    m_isChanged = true;
  }
}

bool MediumMagboltz::EnablePenningTransfer(const double r,
                                           const double lambda) {
   
  if (!MediumGas::EnablePenningTransfer(r, lambda)) return false;
 
  m_rPenning.fill(0.);
  m_lambdaPenning.fill(0.);

  // Make sure that the collision rate table is updated.
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::EnablePenningTransfer");
      return false;
    }
    m_isChanged = false;
  }
  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation) {
      ++nLevelsFound;
    }
    m_rPenning[i] = m_rPenningGlobal;
    m_lambdaPenning[i] = m_lambdaPenningGlobal;
  }

  if (nLevelsFound > 0) {
    std::cout << m_className << "::EnablePenningTransfer:\n    "
              << "Updated Penning transfer parameters for " << nLevelsFound
              << " excitation cross-sections.\n";
    if (nLevelsFound != m_excLevels.size() && !m_excLevels.empty()) {
      std::cerr << m_className << "::EnablePenningTransfer:\n    Warning: "
                << "mismatch between number of excitation cross-sections ("
                << nLevelsFound << ")\n    and number of excitation rates in "
                << "the gas table (" << m_excLevels.size() << ").\n    "
                << "The gas table was probably calculated using a different "
                << "version of Magboltz.\n";
    }
  } else {
    std::cerr << m_className << "::EnablePenningTransfer:\n    "
              << "No excitation cross-sections in the present energy range.\n";
  }

  if (m_useDeexcitation) {
    std::cout << m_className << "::EnablePenningTransfer:\n    "
              << "Deexcitation handling will be switched off.\n";
  }
  m_usePenning = true;
  return true;
}

bool MediumMagboltz::EnablePenningTransfer(const double r, const double lambda,
                                           std::string gasname) {

  if (!MediumGas::EnablePenningTransfer(r, lambda, gasname)) return false;

  // Look (again) for this gas in the present mixture.
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      iGas = i;
      break;
    }
  }

  // Make sure that the collision rate table is updated.
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::EnablePenningTransfer");
      return false;
    }
    m_isChanged = false;
  }
  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (int(m_csType[i] / nCsTypes) != iGas) continue;
    if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation) {
      ++nLevelsFound;
    }
    m_rPenning[i] = m_rPenningGas[iGas];
    m_lambdaPenning[i] = m_lambdaPenningGas[iGas];
  }

  if (nLevelsFound > 0) {
    std::cout << m_className << "::EnablePenningTransfer:\n"
              << "    Penning transfer parameters for " << nLevelsFound
              << " excitation levels set to:\n"
              << "      r      = " << m_rPenningGas[iGas] << "\n"
              << "      lambda = " << m_lambdaPenningGas[iGas] << " cm\n";
  } else {
    std::cerr << m_className << "::EnablePenningTransfer:\n"
              << "    Specified gas (" << gasname
              << ") has no excitation levels in the present energy range.\n";
  }

  m_usePenning = true;
  return true;
}

void MediumMagboltz::DisablePenningTransfer() {

  MediumGas::DisablePenningTransfer();
  m_rPenning.fill(0.);
  m_lambdaPenning.fill(0.);

  m_usePenning = false;
}

bool MediumMagboltz::DisablePenningTransfer(std::string gasname) {

  if (!MediumGas::DisablePenningTransfer(gasname)) return false;
  // Get the "standard" name of this gas.
  gasname = GetGasName(gasname);
  if (gasname.empty()) return false;

  // Look (again) for this gas in the present gas mixture.
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      iGas = i;
      break;
    }
  }

  if (iGas < 0) return false;

  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (int(m_csType[i] / nCsTypes) == iGas) {
      m_rPenning[i] = 0.;
      m_lambdaPenning[i] = 0.;
    } else {
      if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation &&
          m_rPenning[i] > Small) {
        ++nLevelsFound;
      }
    }
  }

  if (nLevelsFound == 0) {
    // There are no more excitation levels with r > 0.
    std::cout << m_className << "::DisablePenningTransfer:\n"
              << "    Penning transfer switched off for all excitations.\n";
    m_usePenning = false;
  }
  return true;
}

void MediumMagboltz::SetExcitationScaling(const double r, std::string gasname) {
  if (r <= 0.) {
    std::cerr << m_className << "::SetExcitationScaling: Incorrect value.\n";
    return;
  }

  // Get the "standard" name of this gas.
  gasname = GetGasName(gasname);
  if (gasname.empty()) {
    std::cerr << m_className << "::SetExcitationScaling: Unknown gas name.\n";
    return;
  }

  // Look for this gas in the present gas mixture.
  bool found = false;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      m_scaleExc[i] = r;
      found = true;
      break;
    }
  }

  if (!found) {
    std::cerr << m_className << "::SetExcitationScaling:\n"
              << "    Specified gas (" << gasname
              << ") is not part of the present gas mixture.\n";
    return;
  }

  // Make sure that the collision rate table is updated.
  m_isChanged = true;
}

bool MediumMagboltz::Initialise(const bool verbose) {
  if (!m_isChanged) {
    if (m_debug) {
      std::cerr << m_className << "::Initialise: Nothing changed.\n";
    }
    return true;
  }
  if (!Mixer(verbose)) {
    PrintErrorMixer(m_className + "::Initialise");
    return false;
  }
  m_isChanged = false;
  return true;
}

void MediumMagboltz::PrintGas() {
  MediumGas::PrintGas();

  if (m_isChanged) {
    if (!Initialise()) return;
  }

  std::cout << "    Electron cross-sections:\n";
  int igas = -1; 
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    // Collision type
    int type = m_csType[i] % nCsTypes;
    if (igas != int(m_csType[i] / nCsTypes)) {
      igas = int(m_csType[i] / nCsTypes);
      std::cout << "      " << m_gas[igas] << "\n";
    }
    // Description (from Magboltz)
    // Threshold energy
    double e = m_rgas[igas] * m_energyLoss[i];
    std::cout << "        Level " << i << ": " << m_description[i] << "\n";
    std::cout << "          Type " << type;
    if (type == ElectronCollisionTypeElastic) {
      std::cout << " (elastic)\n";
    } else if (type == ElectronCollisionTypeIonisation) {
      std::cout << " (ionisation). Ionisation threshold: " << e << " eV.\n";
    } else if (type == ElectronCollisionTypeAttachment) {
      std::cout << " (attachment)\n";
    } else if (type == ElectronCollisionTypeInelastic) {
      std::cout << " (inelastic). Energy loss: " << e << " eV.\n";
    } else if (type == ElectronCollisionTypeExcitation) {
      std::cout << " (excitation). Excitation energy: " << e << " eV.\n";
    } else if (type == ElectronCollisionTypeSuperelastic) {
      std::cout << " (super-elastic). Energy gain: " << -e << " eV.\n";
    } else if (type == ElectronCollisionTypeVirtual) {
      std::cout << " (virtual)\n";
    } else {
      std::cout << " (unknown)\n";
    }
    if (type == ElectronCollisionTypeExcitation && m_usePenning &&
        e > m_minIonPot) {
      std::cout << "          Penning transfer coefficient: " 
                << m_rPenning[i] << "\n";
    } else if (type == ElectronCollisionTypeExcitation && m_useDeexcitation) {
      const int idxc = m_iDeexcitation[i];
      if (idxc < 0 || idxc >= (int)m_deexcitations.size()) {
        std::cout << "          Deexcitation cascade not implemented.\n";
        continue;
      }
      const auto& dxc = m_deexcitations[idxc];
      if (dxc.osc > 0.) {
        std::cout << "          Oscillator strength: " << dxc.osc << "\n";
      }
      std::cout << "          Decay channels:\n";
      const int nChannels = dxc.type.size();
      for (int j = 0; j < nChannels; ++j) {
        if (dxc.type[j] == DxcTypeRad) {
          std::cout << "          Radiative decay to ";
          if (dxc.final[j] < 0) {
            std::cout << "ground state: ";
          } else {
            std::cout << m_deexcitations[dxc.final[j]].label << ": ";
          }
        } else if (dxc.type[j] == DxcTypeCollIon) {
          if (dxc.final[j] < 0) {
            std::cout << "          Penning ionisation: ";
          } else {
            std::cout << "          Associative ionisation: ";
          }
        } else if (dxc.type[j] == DxcTypeCollNonIon) {
          if (dxc.final[j] >= 0) {
            std::cout << "          Collision-induced transition to "
                      << m_deexcitations[dxc.final[j]].label << ": ";
          } else {
            std::cout << "          Loss: ";
          }
        }
        const double br = j == 0 ? dxc.p[j] : dxc.p[j] - dxc.p[j - 1];
        std::cout << std::setprecision(5) << br * 100. << "%\n";
      }
    }
  }
}

double MediumMagboltz::GetElectronNullCollisionRate(const int band) {
  // If necessary, update the collision rates table.
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::GetElectronNullCollisionRate");
      return 0.;
    }
    m_isChanged = false;
  }

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronNullCollisionRate: Band > 0.\n";
  }

  return m_cfNull;
}

double MediumMagboltz::GetElectronCollisionRate(const double e,
                                                const int band) {
  // Check if the electron energy is within the currently set range.
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate: Invalid energy.\n";
    return m_cfTot[0];
  }
  if (e > m_eFinal && m_useAutoAdjust) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n    Rate at " << e
              << " eV is not included in the current table.\n    "
              << "Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  // If necessary, update the collision rates table.
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::GetElectronCollisionRate");
      return 0.;
    }
    m_isChanged = false;
  }

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronCollisionRate: Band > 0.\n";
  }

  // Get the energy interval.
  if (e <= m_eHigh) {
    // Linear binning
    constexpr int iemax = Magboltz::nEnergySteps - 1;
    const int iE = std::min(std::max(int(e / m_eStep), 0), iemax);
    return m_cfTot[iE];
  }

  // Logarithmic binning
  const double eLog = log(e);
  int iE = int((eLog - m_eHighLog) / m_lnStep);
  // Calculate the collision rate by log-log interpolation.
  const double fmax = m_cfTotLog[iE];
  const double fmin = iE == 0 ? log(m_cfTot.back()) : m_cfTotLog[iE - 1];
  const double emin = m_eHighLog + iE * m_lnStep;
  const double f = fmin + (eLog - emin) * (fmax - fmin) / m_lnStep;
  return exp(f);
}

double MediumMagboltz::GetElectronCollisionRate(const double e,
                                                const unsigned int level,
                                                const int band) {
  // Check if the electron energy is within the currently set range.
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate: Invalid energy.\n";
    return 0.;
  }

  // Check if the level exists.
  if (level >= m_nTerms) {
    std::cerr << m_className << "::GetElectronCollisionRate: Invalid level.\n";
    return 0.;
  }

  // Get the total scattering rate.
  double rate = GetElectronCollisionRate(e, band);
  // Get the energy interval.
  if (e <= m_eHigh) {
    // Linear binning
    constexpr int iemax = Magboltz::nEnergySteps - 1;
    const int iE = std::min(std::max(int(e / m_eStep), 0), iemax);
    if (level == 0) {
      rate *= m_cf[iE][0];
    } else {
      rate *= m_cf[iE][level] - m_cf[iE][level - 1];
    }
  } else {
    // Logarithmic binning
    const int iE = int((log(e) - m_eHighLog) / m_lnStep);
    if (level == 0) {
      rate *= m_cfLog[iE][0];
    } else {
      rate *= m_cfLog[iE][level] - m_cfLog[iE][level - 1];
    }
  }
  return rate;
}

bool MediumMagboltz::GetElectronCollision(
    const double e, int& type, int& level, double& e1, double& dx, double& dy,
    double& dz, std::vector<std::pair<int, double> >& secondaries, int& ndxc,
    int& band) {
  ndxc = 0;
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollision: Invalid energy.\n";
    return false;
  }
  // Check if the electron energy is within the currently set range.
  if (e > m_eFinal && m_useAutoAdjust) {
    std::cerr << m_className << "::GetElectronCollision:\n    Provided energy ("
              << e << " eV) exceeds current energy range.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  // If necessary, update the collision rates table.
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::GetElectronCollision");
      return false;
    }
    m_isChanged = false;
  }

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronCollision: Band > 0.\n";
  }

  double angCut = 1.;
  double angPar = 0.5;

  if (e <= m_eHigh) {
    // Linear binning
    // Get the energy interval.
    constexpr int iemax = Magboltz::nEnergySteps - 1;
    const int iE = std::min(std::max(int(e / m_eStep), 0), iemax);

    // Sample the scattering process.
    const double r = RndmUniform();
    if (r <= m_cf[iE][0]) {
      level = 0;
    } else if (r >= m_cf[iE][m_nTerms - 1]) {
      level = m_nTerms - 1;
    } else {
      const auto begin = m_cf[iE].cbegin();
      level = std::lower_bound(begin, begin + m_nTerms, r) - begin;
    }
    // Get the angular distribution parameters.
    angCut = m_scatCut[iE][level];
    angPar = m_scatPar[iE][level];
  } else {
    // Logarithmic binning
    // Get the energy interval.
    const int iE = std::min(std::max(int(log(e / m_eHigh) / m_lnStep), 0),
                            nEnergyStepsLog - 1);
    // Sample the scattering process.
    const double r = RndmUniform();
    if (r <= m_cfLog[iE][0]) {
      level = 0;
    } else if (r >= m_cfLog[iE][m_nTerms - 1]) {
      level = m_nTerms - 1;
    } else {
      const auto begin = m_cfLog[iE].cbegin();
      level = std::lower_bound(begin, begin + m_nTerms, r) - begin;
    }
    // Get the angular distribution parameters.
    angCut = m_scatCutLog[iE][level];
    angPar = m_scatParLog[iE][level];
  }

  // Extract the collision type.
  type = m_csType[level] % nCsTypes;
  const int igas = int(m_csType[level] / nCsTypes);
  // Increase the collision counters.
  ++m_nCollisions[type];
  ++m_nCollisionsDetailed[level];

  // Get the energy loss for this process.
  double loss = m_energyLoss[level];

  if (type == ElectronCollisionTypeVirtual) return true;

  if (type == ElectronCollisionTypeIonisation) {
    // Sample the secondary electron energy according to
    // the Opal-Beaty-Peterson parameterisation.
    double esec = 0.;
    if (m_useOpalBeaty) {
      // Get the splitting parameter.
      const double w = m_wOpalBeaty[level];
      esec = w * tan(RndmUniform() * atan(0.5 * (e - loss) / w));
      // Rescaling (SST)
      // esec = w * pow(esec / w, 0.9524);
    } else if (m_useGreenSawada) {
      const double gs = m_parGreenSawada[igas][0];
      const double gb = m_parGreenSawada[igas][1];
      const double w = gs * e / (e + gb);
      const double ts = m_parGreenSawada[igas][2];
      const double ta = m_parGreenSawada[igas][3];
      const double tb = m_parGreenSawada[igas][4];
      const double esec0 = ts - ta / (e + tb);
      const double r = RndmUniform();
      esec = esec0 +
             w * tan((r - 1.) * atan(esec0 / w) +
                     r * atan((0.5 * (e - loss) - esec0) / w));
    } else {
      esec = RndmUniform() * (e - loss);
    }
    if (esec <= 0) esec = Small;
    loss += esec;
    // Add the secondary electron.
    secondaries.emplace_back(std::make_pair(IonProdTypeElectron, esec));
    // Add the ion.
    secondaries.emplace_back(std::make_pair(IonProdTypeIon, 0.));
  } else if (type == ElectronCollisionTypeExcitation) {
    // if (m_gas[igas] == "CH4" && loss * m_rgas[igas] < 13.35 && e > 12.65) {
    //   if (RndmUniform() < 0.5) {
    //     loss = 8.55 + RndmUniform() * (13.3 - 8.55);
    //     loss /= m_rgas[igas];
    //   } else {
    //     loss = std::max(Small, RndmGaussian(loss * m_rgas[igas], 1.));
    //     loss /= m_rgas[igas];
    //   }
    // }
    // Follow the de-excitation cascade (if switched on).
    if (m_useDeexcitation && m_iDeexcitation[level] >= 0) {
      int fLevel = 0;
      ComputeDeexcitationInternal(m_iDeexcitation[level], fLevel);
      ndxc = m_dxcProducts.size();
    } else if (m_usePenning) {
      m_dxcProducts.clear();
      // Simplified treatment of Penning ionisation.
      // If the energy threshold of this level exceeds the
      // ionisation potential of one of the gases,
      // create a new electron (with probability rPenning).
      if (m_energyLoss[level] * m_rgas[igas] > m_minIonPot &&
          RndmUniform() < m_rPenning[level]) {
        // The energy of the secondary electron is assumed to be given by
        // the difference of excitation and ionisation threshold.
        double esec = m_energyLoss[level] * m_rgas[igas] - m_minIonPot;
        if (esec <= 0) esec = Small;
        // Add the secondary electron to the list.
        dxcProd newDxcProd;
        newDxcProd.t = 0.;
        newDxcProd.s = 0.;
        if (m_lambdaPenning[level] > Small) {
          // Uniform distribution within a sphere of radius lambda
          newDxcProd.s =
              m_lambdaPenning[level] * pow(RndmUniformPos(), 1. / 3.);
        }
        newDxcProd.energy = esec;
        newDxcProd.type = DxcProdTypeElectron;
        m_dxcProducts.push_back(std::move(newDxcProd));
        ndxc = 1;
        ++m_nPenning;
      }
    }
  }

  // Make sure the energy loss is smaller than the energy.
  if (e < loss) loss = e - 0.0001;

  // Determine the scattering angle.
  double ctheta0 = 1. - 2. * RndmUniform();
  if (m_useAnisotropic) {
    switch (m_scatModel[level]) {
      case 0:
        break;
      case 1:
        ctheta0 = 1. - RndmUniform() * angCut;
        if (RndmUniform() > angPar) ctheta0 = -ctheta0;
        break;
      case 2:
        ctheta0 = (ctheta0 + angPar) / (1. + angPar * ctheta0);
        break;
      default:
        std::cerr << m_className << "::GetElectronCollision:\n"
                  << "    Unknown scattering model.\n"
                  << "    Using isotropic distribution.\n";
        break;
    }
  }

  const double s1 = m_rgas[igas];
  const double s2 = (s1 * s1) / (s1 - 1.);
  const double theta0 = acos(ctheta0);
  const double arg = std::max(1. - s1 * loss / e, Small);
  const double d = 1. - ctheta0 * sqrt(arg);

  // Update the energy.
  e1 = std::max(e * (1. - loss / (s1 * e) - 2. * d / s2), Small);
  double q = std::min(sqrt((e / e1) * arg) / s1, 1.);
  const double theta = asin(q * sin(theta0));
  double ctheta = cos(theta);
  if (ctheta0 < 0.) {
    const double u = (s1 - 1.) * (s1 - 1.) / arg;
    if (ctheta0 * ctheta0 > u) ctheta = -ctheta;
  }
  const double stheta = sin(theta);
  // Calculate the direction after the collision.
  dz = std::min(dz, 1.);
  const double argZ = sqrt(dx * dx + dy * dy);

  // Azimuth is chosen at random.
  const double phi = TwoPi * RndmUniform();
  const double cphi = cos(phi);
  const double sphi = sin(phi);
  if (argZ == 0.) {
    dz = ctheta;
    dx = cphi * stheta;
    dy = sphi * stheta;
  } else {
    const double a = stheta / argZ;
    const double dz1 = dz * ctheta + argZ * stheta * sphi;
    const double dy1 = dy * ctheta + a * (dx * cphi - dy * dz * sphi);
    const double dx1 = dx * ctheta - a * (dy * cphi + dx * dz * sphi);
    dz = dz1;
    dy = dy1;
    dx = dx1;
  }

  return true;
}

bool MediumMagboltz::GetDeexcitationProduct(const unsigned int i, double& t,
                                            double& s, int& type,
                                            double& energy) const {
  if (i >= m_dxcProducts.size() || !(m_useDeexcitation || m_usePenning)) {
    return false;
  }
  t = m_dxcProducts[i].t;
  s = m_dxcProducts[i].s;
  type = m_dxcProducts[i].type;
  energy = m_dxcProducts[i].energy;
  return true;
}

double MediumMagboltz::GetPhotonCollisionRate(const double e) {
  if (e <= 0.) {
    std::cerr << m_className << "::GetPhotonCollisionRate: Invalid  energy.\n";
    return m_cfTotGamma[0];
  }
  if (e > m_eFinalGamma && m_useAutoAdjust) {
    std::cerr << m_className << "::GetPhotonCollisionRate:\n    Rate at " << e
              << " eV is not included in the current table.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxPhotonEnergy(1.05 * e);
  }

  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::GetPhotonCollisionRate");
      return 0.;
    }
    m_isChanged = false;
  }

  const int iE =
      std::min(std::max(int(e / m_eStepGamma), 0), nEnergyStepsGamma - 1);

  double cfSum = m_cfTotGamma[iE];
  if (m_useDeexcitation && m_useRadTrap && !m_deexcitations.empty()) {
    // Loop over the excitations.
    for (const auto& dxc : m_deexcitations) {
      if (dxc.cf > 0. && fabs(e - dxc.energy) <= dxc.width) {
        cfSum += dxc.cf *
                 TMath::Voigt(e - dxc.energy, dxc.sDoppler, 2 * dxc.gPressure);
      }
    }
  }

  return cfSum;
}

bool MediumMagboltz::GetPhotonCollision(const double e, int& type, int& level,
                                        double& e1, double& ctheta, int& nsec,
                                        double& esec) {
  if (e <= 0.) {
    std::cerr << m_className << "::GetPhotonCollision: Invalid energy.\n";
    return false;
  }
  if (e > m_eFinalGamma && m_useAutoAdjust) {
    std::cerr << m_className << "::GetPhotonCollision:\n    Provided energy ("
              << e << " eV) exceeds current energy range.\n"
              << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxPhotonEnergy(1.05 * e);
  }

  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::GetPhotonCollision");
      return false;
    }
    m_isChanged = false;
  }

  // Energy interval
  const int iE =
      std::min(std::max(int(e / m_eStepGamma), 0), nEnergyStepsGamma - 1);

  double r = m_cfTotGamma[iE];
  if (m_useDeexcitation && m_useRadTrap && !m_deexcitations.empty()) {
    int nLines = 0;
    std::vector<double> pLine(0);
    std::vector<int> iLine(0);
    // Loop over the excitations.
    const unsigned int nDeexcitations = m_deexcitations.size();
    for (unsigned int i = 0; i < nDeexcitations; ++i) {
      const auto& dxc = m_deexcitations[i];
      if (dxc.cf > 0. && fabs(e - dxc.energy) <= dxc.width) {
        r += dxc.cf *
             TMath::Voigt(e - dxc.energy, dxc.sDoppler, 2 * dxc.gPressure);
        pLine.push_back(r);
        iLine.push_back(i);
        ++nLines;
      }
    }
    r *= RndmUniform();
    if (nLines > 0 && r >= m_cfTotGamma[iE]) {
      // Photon is absorbed by a discrete line.
      for (int i = 0; i < nLines; ++i) {
        if (r <= pLine[i]) {
          ++m_nPhotonCollisions[PhotonCollisionTypeExcitation];
          int fLevel = 0;
          ComputeDeexcitationInternal(iLine[i], fLevel);
          type = PhotonCollisionTypeExcitation;
          nsec = m_dxcProducts.size();
          return true;
        }
      }
      std::cerr << m_className << "::GetPhotonCollision:\n";
      std::cerr << "    Random sampling of deexcitation line failed.\n";
      std::cerr << "    Program bug!\n";
      return false;
    }
  } else {
    r *= RndmUniform();
  }

  if (r <= m_cfGamma[iE][0]) {
    level = 0;
  } else if (r >= m_cfGamma[iE][m_nPhotonTerms - 1]) {
    level = m_nPhotonTerms - 1;
  } else {
    const auto begin = m_cfGamma[iE].cbegin();
    level = std::lower_bound(begin, begin + m_nPhotonTerms, r) - begin;
  }

  nsec = 0;
  esec = e1 = 0.;
  type = csTypeGamma[level];
  // Collision type
  type = type % nCsTypesGamma;
  int ngas = int(csTypeGamma[level] / nCsTypesGamma);
  ++m_nPhotonCollisions[type];
  // Ionising collision
  if (type == 1) {
    esec = std::max(e - m_ionPot[ngas], Small);
    nsec = 1;
  }

  // Determine the scattering angle
  ctheta = 2 * RndmUniform() - 1.;

  return true;
}

void MediumMagboltz::ResetCollisionCounters() {
  m_nCollisions.fill(0);
  m_nCollisionsDetailed.assign(m_nTerms, 0);
  m_nPenning = 0;
  m_nPhotonCollisions.fill(0);
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions() const {
  return std::accumulate(std::begin(m_nCollisions), std::end(m_nCollisions), 0);
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions(
    unsigned int& nElastic, unsigned int& nIonisation,
    unsigned int& nAttachment, unsigned int& nInelastic,
    unsigned int& nExcitation, unsigned int& nSuperelastic) const {
  nElastic = m_nCollisions[ElectronCollisionTypeElastic];
  nIonisation = m_nCollisions[ElectronCollisionTypeIonisation];
  nAttachment = m_nCollisions[ElectronCollisionTypeAttachment];
  nInelastic = m_nCollisions[ElectronCollisionTypeInelastic];
  nExcitation = m_nCollisions[ElectronCollisionTypeExcitation];
  nSuperelastic = m_nCollisions[ElectronCollisionTypeSuperelastic];
  return nElastic + nIonisation + nAttachment + nInelastic + nExcitation +
         nSuperelastic;
}

unsigned int MediumMagboltz::GetNumberOfLevels() {
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::GetNumberOfLevels");
      return 0;
    }
    m_isChanged = false;
  }

  return m_nTerms;
}

bool MediumMagboltz::GetLevel(const unsigned int i, int& ngas, int& type,
                              std::string& descr, double& e) {
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::GetLevel");
      return false;
    }
    m_isChanged = false;
  }

  if (i >= m_nTerms) {
    std::cerr << m_className << "::GetLevel: Index out of range.\n";
    return false;
  }

  // Collision type
  type = m_csType[i] % nCsTypes;
  ngas = int(m_csType[i] / nCsTypes);
  // Description (from Magboltz)
  descr = m_description[i];
  // Threshold energy
  e = m_rgas[ngas] * m_energyLoss[i];
  if (m_debug) {
    std::cout << m_className << "::GetLevel:\n"
              << "    Level " << i << ": " << descr << "\n"
              << "    Type " << type << "\n"
              << "    Threshold energy: " << e << " eV\n";
    if (type == ElectronCollisionTypeExcitation && m_usePenning &&
        e > m_minIonPot) {
      std::cout << "    Penning transfer coefficient: " << m_rPenning[i]
                << "\n";
    } else if (type == ElectronCollisionTypeExcitation && m_useDeexcitation) {
      const int idxc = m_iDeexcitation[i];
      if (idxc < 0 || idxc >= (int)m_deexcitations.size()) {
        std::cout << "    Deexcitation cascade not implemented.\n";
        return true;
      }
      const auto& dxc = m_deexcitations[idxc];
      if (dxc.osc > 0.) {
        std::cout << "    Oscillator strength: " << dxc.osc << "\n";
      }
      std::cout << "    Decay channels:\n";
      const int nChannels = dxc.type.size();
      for (int j = 0; j < nChannels; ++j) {
        if (dxc.type[j] == DxcTypeRad) {
          std::cout << "      Radiative decay to ";
          if (dxc.final[j] < 0) {
            std::cout << "ground state: ";
          } else {
            std::cout << m_deexcitations[dxc.final[j]].label << ": ";
          }
        } else if (dxc.type[j] == DxcTypeCollIon) {
          if (dxc.final[j] < 0) {
            std::cout << "      Penning ionisation: ";
          } else {
            std::cout << "      Associative ionisation: ";
          }
        } else if (dxc.type[j] == DxcTypeCollNonIon) {
          if (dxc.final[j] >= 0) {
            std::cout << "      Collision-induced transition to "
                      << m_deexcitations[dxc.final[j]].label << ": ";
          } else {
            std::cout << "      Loss: ";
          }
        }
        const double br = j == 0 ? dxc.p[j] : dxc.p[j] - dxc.p[j - 1];
        std::cout << std::setprecision(5) << br * 100. << "%\n";
      }
    }
  }

  return true;
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions(
    const unsigned int level) const {
  if (level >= m_nTerms) {
    std::cerr << m_className << "::GetNumberOfElectronCollisions: "
              << "Level " << level << " does not exist.\n";
    return 0;
  }
  return m_nCollisionsDetailed[level];
}

unsigned int MediumMagboltz::GetNumberOfPhotonCollisions() const {
  return std::accumulate(std::begin(m_nPhotonCollisions),
                         std::end(m_nPhotonCollisions), 0);
}

unsigned int MediumMagboltz::GetNumberOfPhotonCollisions(
    unsigned int& nElastic, unsigned int& nIonising,
    unsigned int& nInelastic) const {
  nElastic = m_nPhotonCollisions[0];
  nIonising = m_nPhotonCollisions[1];
  nInelastic = m_nPhotonCollisions[2];
  return nElastic + nIonising + nInelastic;
}

int MediumMagboltz::GetGasNumberMagboltz(const std::string& input) const {

  if (input.empty()) return 0;

  if (input == "CF4") {
    return 1;
  } else if (input == "Ar") {
    return 2;
  } else if (input == "He" || input == "He-4") {
    // Helium 4
    return 3;
  } else if (input == "He-3") {
    // Helium 3
    return 4;
  } else if (input == "Ne") {
    return 5;
  } else if (input == "Kr") {
    return 6;
  } else if (input == "Xe") {
    return 7;
  } else if (input == "CH4") {
    // Methane
    return 8;
  } else if (input == "C2H6") {
    // Ethane
    return 9;
  } else if (input == "C3H8") {
    // Propane
    return 10;
  } else if (input == "iC4H10") {
    // Isobutane
    return 11;
  } else if (input == "CO2") {
    return 12;
  } else if (input == "neoC5H12") {
    // Neopentane
    return 13;
  } else if (input == "H2O") {
    return 14;
  } else if (input == "O2") {
    return 15;
  } else if (input == "N2") {
    return 16;
  } else if (input == "NO") {
    // Nitric oxide (NO)
    return 17;
  } else if (input == "N2O") {
    // Nitrous oxide (N2O)
    return 18;
  } else if (input == "C2H4") {
    // Ethene (C2H4)
    return 19;
  } else if (input == "C2H2") {
    // Acetylene (C2H2)
    return 20;
  } else if (input == "H2") {
    // Hydrogen
    return 21;
  } else if (input == "D2") {
    // Deuterium
    return 22;
  } else if (input == "CO") {
    // Carbon monoxide (CO)
    return 23;
  } else if (input == "Methylal") {
    // Methylal (dimethoxymethane, CH3-O-CH2-O-CH3, "hot" version)
    return 24;
  } else if (input == "DME") {
    return 25;
  } else if (input == "Reid-Step") {
    return 26;
  } else if (input == "Maxwell-Model") {
    return 27;
  } else if (input == "Reid-Ramp") {
    return 28;
  } else if (input == "C2F6") {
    return 29;
  } else if (input == "SF6") {
    return 30;
  } else if (input == "NH3") {
    return 31;
  } else if (input == "C3H6") {
    // Propene
    return 32;
  } else if (input == "cC3H6") {
    // Cyclopropane
    return 33;
  } else if (input == "CH3OH") {
    // Methanol
    return 34;
  } else if (input == "C2H5OH") {
    // Ethanol
    return 35;
  } else if (input == "C3H7OH") {
    // Propanol
    return 36;
  } else if (input == "Cs") {
    return 37;
  } else if (input == "F2") {
    // Fluorine
    return 38;
  } else if (input == "CS2") {
    return 39;
  } else if (input == "COS") {
    return 40;
  } else if (input == "CD4") {
    // Deuterated methane
    return 41;
  } else if (input == "BF3") {
    return 42;
  } else if (input == "C2HF5" || input == "C2H2F4") {
    return 43;
  } else if (input == "TMA") {
    return 44;
  } else if (input == "paraH2") {
    // Para hydrogen
    return 45;
  } else if (input == "nC3H7OH") {
    // n-propanol
    return 46;
  } else if (input == "CHF3") {
    return 50;
  } else if (input == "CF3Br") {
    return 51;
  } else if (input == "C3F8") {
    return 52;
  } else if (input == "O3") {
    // Ozone
    return 53;
  } else if (input == "Hg") {
    // Mercury
    return 54;
  } else if (input == "H2S") {
    return 55;
  } else if (input == "nC4H10") {
    // n-Butane
    return 56;
  } else if (input == "nC5H12") {
    // n-Pentane
    return 57;
  } else if (input == "N2 (Phelps)") {
    return 58;
  } else if (input == "GeH4") {
    // Germane, GeH4
    return 59;
  } else if (input == "SiH4") {
    // Silane, SiH4
    return 60;
  }

  std::cerr << m_className << "::GetGasNumberMagboltz:\n"
            << "    Gas " << input << " is not defined.\n";
  return 0;
}

bool MediumMagboltz::Mixer(const bool verbose) {
  // Set constants and parameters in Magboltz common blocks.
  Magboltz::cnsts_.echarg = ElementaryCharge * 1.e-15;
  Magboltz::cnsts_.emass = ElectronMassGramme;
  Magboltz::cnsts_.amu = AtomicMassUnit;
  Magboltz::cnsts_.pir2 = BohrRadius * BohrRadius * Pi;
  Magboltz::inpt_.ary = RydbergEnergy;

  Magboltz::inpt_.akt = BoltzmannConstant * m_temperature;
  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;

  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = Magboltz::nEnergySteps;
  Magboltz::inpt_.nAniso = m_useAnisotropic ? 2 : 0;

  for (unsigned int i = 0; i < Magboltz::nEnergySteps; ++i) {
    const double en = (i + 0.5) * m_eStep;
    Magboltz::mix2_.eg[i] = en;
    Magboltz::mix2_.eroot[i] = sqrt(en);
    Magboltz::dens_.den[i] = 0.;
  }
  constexpr int iemax = Magboltz::nEnergySteps - 1;

  // Calculate the atomic density (ideal gas law).
  const double dens = GetNumberDensity();
  // Prefactor for calculation of scattering rate from cross-section.
  const double prefactor = dens * SpeedOfLight * sqrt(2. / ElectronMass);

  m_rgas.fill(1.);

  m_ionPot.fill(-1.);
  m_minIonPot = -1.;

  m_parGreenSawada.fill({1., 0., 0., 0., 0.});
  m_hasGreenSawada.fill(false);

  m_wOpalBeaty.fill(1.);
  m_energyLoss.fill(0.);
  m_csType.fill(0);
  m_scatModel.fill(0);

  m_rPenning.fill(0.);
  m_lambdaPenning.fill(0.);

  m_deexcitations.clear();
  m_iDeexcitation.fill(-1);

  // Reset the collision rates.
  m_cfTot.assign(Magboltz::nEnergySteps, 0.);
  m_cfTotLog.assign(nEnergyStepsLog, 0.);

  m_cf.assign(Magboltz::nEnergySteps,
              std::vector<double>(Magboltz::nMaxLevels, 0.));
  m_cfLog.assign(nEnergyStepsLog,
                 std::vector<double>(Magboltz::nMaxLevels, 0.));

  m_scatPar.assign(Magboltz::nEnergySteps,
                   std::vector<double>(Magboltz::nMaxLevels, 0.5));
  m_scatCut.assign(Magboltz::nEnergySteps,
                   std::vector<double>(Magboltz::nMaxLevels, 1.));

  m_scatParLog.assign(nEnergyStepsLog,
                      std::vector<double>(Magboltz::nMaxLevels, 0.5));
  m_scatCutLog.assign(nEnergyStepsLog,
                      std::vector<double>(Magboltz::nMaxLevels, 1.));

  // Cross-sections
  // 0: total, 1: elastic,
  // 2: ionisation, 3: attachment,
  // 4, 5: unused
  static double q[Magboltz::nEnergySteps][6];
  // Inelastic cross-sections
  static double qIn[Magboltz::nEnergySteps][Magboltz::nMaxInelasticTerms];
  // Ionisation cross-sections
  static double qIon[Magboltz::nEnergySteps][Magboltz::nMaxIonisationTerms];
  // Attachment cross-sections
  static double qAtt[Magboltz::nEnergySteps][Magboltz::nMaxAttachmentTerms];
  // "Null-collision" cross-sections
  static double qNull[Magboltz::nEnergySteps][Magboltz::nMaxNullTerms];
  // Parameters for scattering angular distribution
  static double pEqEl[Magboltz::nEnergySteps][6];
  // Parameters for angular distribution in inelastic collisions
  static double pEqIn[Magboltz::nEnergySteps][Magboltz::nMaxInelasticTerms];
  // Parameters for angular distribution in ionising collisions
  static double pEqIon[Magboltz::nEnergySteps][Magboltz::nMaxIonisationTerms];
  // Penning transfer parameters
  static double penFra[Magboltz::nMaxInelasticTerms][3];
  // Description of cross-section terms
  static char scrpt[Magboltz::nMaxLevelsPerComponent][Magboltz::nCharDescr];
  // Description of "null-collision" cross-section terms
  static char scrptn[Magboltz::nMaxNullTerms][Magboltz::nCharDescr];

  // Check the gas composition and establish the gas numbers.
  int gasNumber[m_nMaxGases];
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const int ng = GetGasNumberMagboltz(m_gas[i]);
    if (ng <= 0) {
      std::cerr << m_className << "::Mixer:\n    Gas " << m_gas[i]
                << " does not have a gas number in Magboltz.\n";
      return false;
    }
    gasNumber[i] = ng;
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n    " << Magboltz::nEnergySteps
              << " linear energy steps between 0 and "
              << std::min(m_eFinal, m_eHigh) << " eV.\n";
    if (m_eFinal > m_eHigh) {
      std::cout << "    " << nEnergyStepsLog << " logarithmic steps between "
                << m_eHigh << " and " << m_eFinal << " eV\n";
    }
  }
  m_nTerms = 0;

  std::ofstream outfile;
  if (m_useCsOutput) {
    outfile.open("cs.txt", std::ios::out);
    outfile << "# energy [eV] vs. cross-section [cm2]\n";
  }

  // Loop over the gases in the mixture.
  for (unsigned int iGas = 0; iGas < m_nComponents; ++iGas) {
    Magboltz::inpt_.efinal = std::min(m_eFinal, m_eHigh);
    Magboltz::inpt_.estep = m_eStep;
    Magboltz::mix2_.eg[iemax] = (iemax + 0.5) * m_eStep;
    Magboltz::mix2_.eroot[iemax] = sqrt((iemax + 0.5) * m_eStep);
    char name[] = "                         ";
    // Number of inelastic cross-sections
    long long nIn = 0;
    // Number of ionisation cross-sections
    long long nIon = 0;
    // Number of attachment cross-sections
    long long nAtt = 1;
    // Number of "null-collision" cross-sections
    long long nNull = 0;
    // Virial coefficient (not used)
    double virial = 0.;
    // Thresholds/characteristic energies.
    std::array<double, 6> e;
    // Energy losses for inelastic cross-sections.
    std::array<double, Magboltz::nMaxInelasticTerms> eIn;
    // Ionisation thresholds.
    std::array<double, Magboltz::nMaxIonisationTerms> eIon;
    // Scattering algorithms
    std::array<long long, Magboltz::nMaxInelasticTerms> kIn;
    std::array<long long, 6> kEl;
    // Opal-Beaty parameter
    std::array<double, Magboltz::nMaxIonisationTerms> eoby;
    // Scaling factor for "null-collision" terms
    std::array<double, Magboltz::nMaxNullTerms> scln;
    // Parameters for simulation of Auger and fluorescence processes.
    std::array<long long, Magboltz::nMaxIonisationTerms> nc0;
    std::array<long long, Magboltz::nMaxIonisationTerms> ng1;
    std::array<long long, Magboltz::nMaxIonisationTerms> ng2;
    std::array<double, Magboltz::nMaxIonisationTerms> ec0;
    std::array<double, Magboltz::nMaxIonisationTerms> wklm;
    std::array<double, Magboltz::nMaxIonisationTerms> efl;
    std::array<double, Magboltz::nMaxIonisationTerms> eg1;
    std::array<double, Magboltz::nMaxIonisationTerms> eg2;

    // Retrieve the cross-section data for this gas from Magboltz.
    long long ngs = gasNumber[iGas];
    Magboltz::gasmix_(
        &ngs, q[0], qIn[0], &nIn, e.data(), eIn.data(), name, &virial,
        eoby.data(), pEqEl[0], pEqIn[0], penFra[0], kEl.data(), kIn.data(),
        qIon[0], pEqIon[0], eIon.data(), &nIon, qAtt[0], &nAtt, qNull[0],
        &nNull, scln.data(), nc0.data(), ec0.data(), wklm.data(), efl.data(),
        ng1.data(), eg1.data(), ng2.data(), eg2.data(), scrpt, scrptn);
    if (m_debug || verbose) {
      const double m = (2. / e[1]) * ElectronMass / AtomicMassUnitElectronVolt;
      std::cout << "    " << name << "\n"
                << "      mass: " << m << " amu\n";
      if (nIon > 1) {
        std::cout << "      ionisation threshold: " << eIon[0] << " eV\n";
      } else {
        std::cout << "      ionisation threshold: " << e[2] << " eV\n";
      }
      if (e[3] > 0. && e[4] > 0.) {
        std::cout << "      cross-sections at minimum ionising energy:\n"
                  << "        excitation: " << e[3] * 1.e18 << " Mbarn\n"
                  << "        ionisation: " << e[4] * 1.e18 << " Mbarn\n";
      }
    }
    int np0 = m_nTerms;
    // Make sure there is still sufficient space.
    if (np0 + nIn + nIon + nAtt + nNull >= Magboltz::nMaxLevels) {
      std::cerr << m_className << "::Mixer:\n"
                << "    Max. number of levels (" << Magboltz::nMaxLevels
                << ") exceeded.\n";
      return false;
    }
    const double van = m_fraction[iGas] * prefactor;
    int np = np0;
    if (m_useCsOutput) {
      outfile << "# cross-sections for " << name << "\n";
      outfile << "# cross-section types:\n";
      outfile << "# elastic\n";
    }
    // Elastic scattering
    ++m_nTerms;
    m_scatModel[np] = kEl[1];
    const double r = 1. + 0.5 * e[1];
    m_rgas[iGas] = r;
    m_energyLoss[np] = 0.;
    m_description[np] = GetDescription(1, scrpt);
    m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeElastic;
    bool withIon = false;
    // Ionisation
    if (nIon > 1) {
      for (int j = 0; j < nIon; ++j) {
        if (m_eFinal < eIon[j]) continue;
        withIon = true;
        ++m_nTerms;
        ++np;
        m_scatModel[np] = kEl[2];
        m_energyLoss[np] = eIon[j] / r;
        m_wOpalBeaty[np] = eoby[j];
        m_description[np] = GetDescription(2 + j, scrpt);
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeIonisation;
        if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
      }
      m_parGreenSawada[iGas][0] = eoby[0];
      m_parGreenSawada[iGas][4] = 2 * eIon[0];
      m_ionPot[iGas] = eIon[0];
    } else {
      if (m_eFinal >= e[2]) {
        withIon = true;
        ++m_nTerms;
        ++np;
        m_scatModel[np] = kEl[2];
        m_energyLoss[np] = e[2] / r;
        m_wOpalBeaty[np] = eoby[0];
        m_parGreenSawada[iGas][0] = eoby[0];
        m_parGreenSawada[iGas][4] = 2 * e[2];
        m_ionPot[iGas] = e[2];
        m_description[np] = GetDescription(2, scrpt);
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeIonisation;
        if (m_useCsOutput) outfile << "# ionisation (gross)\n";
      }
    }
    // Attachment
    for (int j = 0; j < nAtt; ++j) {
      ++m_nTerms;
      ++np;
      m_scatModel[np] = 0;
      m_energyLoss[np] = 0.;
      m_description[np] = GetDescription(2 + nIon + j, scrpt);
      m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeAttachment;
      if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
    }
    // Inelastic terms
    int nExc = 0, nSuperEl = 0;
    for (int j = 0; j < nIn; ++j) {
      ++np;
      m_scatModel[np] = kIn[j];
      m_energyLoss[np] = eIn[j] / r;
      m_description[np] = GetDescription(4 + nIon + nAtt + j, scrpt);
      if ((m_description[np][1] == 'E' && m_description[np][2] == 'X') ||
          (m_description[np][0] == 'E' && m_description[np][1] == 'X') ||
          (m_gas[iGas] == "N2" && eIn[j] > 6.)) {
        // Excitation
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeExcitation;
        ++nExc;
      } else if (eIn[j] < 0.) {
        // Super-elastic collision
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeSuperelastic;
        ++nSuperEl;
      } else {
        // Inelastic collision
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeInelastic;
      }
      if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
    }
    m_nTerms += nIn;
    if (nNull > 0) {
      for (int j = 0; j < nNull; ++j) {
        ++m_nTerms;
        ++np;
        m_scatModel[np] = 0;
        m_energyLoss[np] = 0.;
        m_description[np] = GetDescription(j, scrptn);
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeVirtual;
        if (m_useCsOutput) outfile << "# " << m_description[np] << "\n";
      }
    }
    // Loop over the energy table.
    for (unsigned int iE = 0; iE < Magboltz::nEnergySteps; ++iE) {
      np = np0;
      if (m_useCsOutput) {
        outfile << (iE + 0.5) * m_eStep << "  " << q[iE][1] << "  ";
      }
      // Elastic scattering
      m_cf[iE][np] = q[iE][1] * van;
      SetScatteringParameters(m_scatModel[np], pEqEl[iE][1], m_scatCut[iE][np],
                              m_scatPar[iE][np]);
      // Ionisation
      if (withIon) {
        if (nIon > 1) {
          for (int j = 0; j < nIon; ++j) {
            if (m_eFinal < eIon[j]) continue;
            ++np;
            m_cf[iE][np] = qIon[iE][j] * van;
            SetScatteringParameters(m_scatModel[np], pEqIon[iE][j],
                                    m_scatCut[iE][np], m_scatPar[iE][np]);
            if (m_useCsOutput) outfile << qIon[iE][j] << "  ";
          }
        } else {
          ++np;
          m_cf[iE][np] = q[iE][2] * van;
          SetScatteringParameters(m_scatModel[np], pEqEl[iE][2],
                                  m_scatCut[iE][np], m_scatPar[iE][np]);
          if (m_useCsOutput) outfile << q[iE][2] << "  ";
        }
      }
      // Attachment
      for (int j = 0; j < nAtt; ++j) {
        ++np;
        m_cf[iE][np] = qAtt[iE][j] * van;
        // m_cf[iE][np] = q[iE][3] * van;
        m_scatPar[iE][np] = 0.5;
        if (m_useCsOutput) outfile << qAtt[iE][j] << "  ";
      }
      // Inelastic terms
      for (int j = 0; j < nIn; ++j) {
        ++np;
        if (m_useCsOutput) outfile << qIn[iE][j] << "  ";
        m_cf[iE][np] = qIn[iE][j] * van;
        // Scale the excitation cross-sections (for error estimates).
        m_cf[iE][np] *= m_scaleExc[iGas];
        if (m_cf[iE][np] < 0.) {
          std::cerr << m_className << "::Mixer:\n"
                    << "    Negative inelastic cross-section at "
                    << (iE + 0.5) * m_eStep << " eV. Set to zero.\n";
          m_cf[iE][np] = 0.;
        }
        SetScatteringParameters(m_scatModel[np], pEqIn[iE][j],
                                m_scatCut[iE][np], m_scatPar[iE][np]);
      }
      if ((m_debug || verbose) && nIn > 0 && iE == iemax) {
        std::cout << "      " << nIn << " inelastic terms (" << nExc
                  << " excitations, " << nSuperEl << " superelastic, "
                  << nIn - nExc - nSuperEl << " other)\n";
      }
      if (nNull > 0) {
        for (int j = 0; j < nNull; ++j) {
          ++np;
          m_cf[iE][np] = qNull[iE][j] * van * scln[j];
          ;
          if (m_useCsOutput) outfile << qNull[iE][j] << "  ";
        }
      }
      if (m_useCsOutput) outfile << "\n";
    }
    if (m_eFinal <= m_eHigh) continue;
    // Fill the high-energy part (logarithmic binning).
    // Calculate the growth factor.
    const double rLog = pow(m_eFinal / m_eHigh, 1. / nEnergyStepsLog);
    m_lnStep = log(rLog);
    // Set the upper limit of the first bin.
    double emax = m_eHigh * rLog;

    for (int iE = 0; iE < nEnergyStepsLog; ++iE) {
      Magboltz::inpt_.estep = emax / (Magboltz::nEnergySteps - 0.5);
      Magboltz::inpt_.efinal = emax + 0.5 * Magboltz::inpt_.estep;
      Magboltz::mix2_.eg[iemax] = emax;
      Magboltz::mix2_.eroot[iemax] = sqrt(emax);
      Magboltz::gasmix_(
          &ngs, q[0], qIn[0], &nIn, e.data(), eIn.data(), name, &virial,
          eoby.data(), pEqEl[0], pEqIn[0], penFra[0], kEl.data(), kIn.data(),
          qIon[0], pEqIon[0], eIon.data(), &nIon, qAtt[0], &nAtt, qNull[0],
          &nNull, scln.data(), nc0.data(), ec0.data(), wklm.data(), efl.data(),
          ng1.data(), eg1.data(), ng2.data(), eg2.data(), scrpt, scrptn);
      np = np0;
      if (m_useCsOutput) outfile << emax << "  " << q[iemax][1] << "  ";
      // Elastic scattering
      m_cfLog[iE][np] = q[iemax][1] * van;
      SetScatteringParameters(m_scatModel[np], pEqEl[iemax][1],
                              m_scatCutLog[iE][np], m_scatParLog[iE][np]);
      // Ionisation
      if (withIon) {
        if (nIon > 1) {
          for (int j = 0; j < nIon; ++j) {
            if (m_eFinal < eIon[j]) continue;
            ++np;
            m_cfLog[iE][np] = qIon[iemax][j] * van;
            SetScatteringParameters(m_scatModel[np], pEqIon[iemax][j],
                                    m_scatCutLog[iE][np], m_scatParLog[iE][np]);
            if (m_useCsOutput) outfile << qIon[iemax][j] << "  ";
          }
        } else {
          ++np;
          // Gross cross-section
          m_cfLog[iE][np] = q[iemax][2] * van;
          // Counting cross-section
          // m_cfLog[iE][np] = q[iemax][4] * van;
          SetScatteringParameters(m_scatModel[np], pEqEl[iemax][2],
                                  m_scatCutLog[iE][np], m_scatParLog[iE][np]);
          if (m_useCsOutput) outfile << q[iemax][2] << "  ";
        }
      }
      // Attachment
      for (int j = 0; j < nAtt; ++j) {
        ++np;
        m_cfLog[iE][np] = qAtt[iemax][j] * van;
        // m_cfLog[iE][np] = q[iemax][3] * van;
        if (m_useCsOutput) outfile << qAtt[iemax][j] << "  ";
      }
      // Inelastic terms
      for (int j = 0; j < nIn; ++j) {
        ++np;
        if (m_useCsOutput) outfile << qIn[iemax][j] << "  ";
        m_cfLog[iE][np] = qIn[iemax][j] * van;
        // Scale the excitation cross-sections (for error estimates).
        m_cfLog[iE][np] *= m_scaleExc[iGas];
        if (m_cfLog[iE][np] < 0.) {
          std::cerr << m_className << "::Mixer:\n"
                    << "    Negative inelastic cross-section at " << emax
                    << " eV. Set to zero.\n";
          m_cfLog[iE][np] = 0.;
        }
        SetScatteringParameters(m_scatModel[np], pEqIn[iemax][j],
                                m_scatCutLog[iE][np], m_scatParLog[iE][np]);
      }
      if (nNull > 0) {
        for (int j = 0; j < nNull; ++j) {
          ++np;
          m_cfLog[iE][np] = qNull[iemax][j] * van * scln[j];
          if (m_useCsOutput) outfile << qNull[iemax][j] << "  ";
        }
      }
      if (m_useCsOutput) outfile << "\n";
      // Increase the energy.
      emax *= rLog;
    }
  }
  if (m_useCsOutput) outfile.close();

  // Find the smallest ionisation threshold.
  auto it = std::min_element(std::begin(m_ionPot),
                             std::begin(m_ionPot) + m_nComponents);
  m_minIonPot = *it;
  std::string minIonPotGas = m_gas[std::distance(std::begin(m_ionPot), it)];

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n"
              << "    Lowest ionisation threshold in the mixture: "
              << m_minIonPot << " eV (" << minIonPotGas << ")\n";
  }

  for (unsigned int iE = 0; iE < Magboltz::nEnergySteps; ++iE) {
    // Calculate the total collision frequency.
    for (unsigned int k = 0; k < m_nTerms; ++k) {
      if (m_cf[iE][k] < 0.) {
        std::cerr << m_className << "::Mixer:\n"
                  << "    Negative collision rate at " << (iE + 0.5) * m_eStep
                  << " eV, cross-section " << k << ". Set to zero.\n";
        std::cout << m_description[k] << "\n";
        m_cf[iE][k] = 0.;
      }
      m_cfTot[iE] += m_cf[iE][k];
    }
    // Normalise the collision probabilities.
    if (m_cfTot[iE] > 0.) {
      for (unsigned int k = 0; k < m_nTerms; ++k) m_cf[iE][k] /= m_cfTot[iE];
    }
    for (unsigned int k = 1; k < m_nTerms; ++k) {
      m_cf[iE][k] += m_cf[iE][k - 1];
    }
    const double ekin = m_eStep * (iE + 0.5);
    m_cfTot[iE] *= sqrt(ekin);
    // Use relativistic expression at high energies.
    if (ekin > 1.e3) {
      const double re = ekin / ElectronMass;
      m_cfTot[iE] *= sqrt(1. + 0.5 * re) / (1. + re);
    }
  }

  if (m_eFinal > m_eHigh) {
    const double rLog = pow(m_eFinal / m_eHigh, 1. / nEnergyStepsLog);
    for (int iE = 0; iE < nEnergyStepsLog; ++iE) {
      // Calculate the total collision frequency.
      for (unsigned int k = 0; k < m_nTerms; ++k) {
        if (m_cfLog[iE][k] < 0.) m_cfLog[iE][k] = 0.;
        m_cfTotLog[iE] += m_cfLog[iE][k];
      }
      // Normalise the collision probabilities.
      if (m_cfTotLog[iE] > 0.) {
        for (unsigned int k = 0; k < m_nTerms; ++k) {
          m_cfLog[iE][k] /= m_cfTotLog[iE];
        }
      }
      for (unsigned int k = 1; k < m_nTerms; ++k) {
        m_cfLog[iE][k] += m_cfLog[iE][k - 1];
      }
      const double ekin = m_eHigh * pow(rLog, iE + 1);
      const double re = ekin / ElectronMass;
      m_cfTotLog[iE] *= sqrt(ekin) * sqrt(1. + re) / (1. + re);
      // Store the logarithm (for log-log interpolation)
      m_cfTotLog[iE] = log(m_cfTotLog[iE]);
    }
  }

  // Determine the null collision frequency.
  m_cfNull = 0.;
  for (unsigned int j = 0; j < Magboltz::nEnergySteps; ++j) {
    if (m_cfTot[j] > m_cfNull) m_cfNull = m_cfTot[j];
  }
  if (m_eFinal > m_eHigh) {
    for (int j = 0; j < nEnergyStepsLog; ++j) {
      const double r = exp(m_cfTotLog[j]);
      if (r > m_cfNull) m_cfNull = r;
    }
  }

  // Reset the collision counters.
  m_nCollisionsDetailed.assign(m_nTerms, 0);
  m_nCollisions.fill(0);

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n"
              << "    Energy [eV]    Collision Rate [ns-1]\n";
    const double emax = std::min(m_eHigh, m_eFinal);
    for (int i = 0; i < 8; ++i) {
      const double en = (2 * i + 1) * emax / 16;
      const double cf = m_cfTot[(i + 1) * Magboltz::nEnergySteps / 16];
      std::printf("    %10.2f    %18.2f\n", en, cf);
    }
  }

  // Set up the de-excitation channels.
  if (m_useDeexcitation) {
    ComputeDeexcitationTable(verbose);
    for (const auto& dxc : m_deexcitations) {
      if (dxc.p.size() == dxc.final.size() && dxc.p.size() == dxc.type.size())
        continue;
      std::cerr << m_className << "::Mixer:\n"
                << "    Mismatch in deexcitation channel count. Program bug!\n"
                << "    Deexcitation handling is switched off.\n";
      m_useDeexcitation = false;
      break;
    }
  }

  // Fill the photon collision rates table.
  if (!ComputePhotonCollisionTable(verbose)) {
    std::cerr << m_className << "::Mixer:\n";
    std::cerr << "    Photon collision rates could not be calculated.\n";
    if (m_useDeexcitation) {
      std::cerr << "    Deexcitation handling is switched off.\n";
      m_useDeexcitation = false;
    }
  }

  // Reset the Penning transfer parameters.
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    m_rPenning[i] = m_rPenningGlobal;
    int iGas = int(m_csType[i] / nCsTypes);
    if (m_rPenningGas[iGas] > Small) {
      m_rPenning[i] = m_rPenningGas[iGas];
      m_lambdaPenning[i] = m_lambdaPenningGas[iGas];
    }
  }

  // Set the Green-Sawada splitting function parameters.
  SetupGreenSawada();

  return true;
}

void MediumMagboltz::SetupGreenSawada() {
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const double ta = 1000.;
    const double tb = m_parGreenSawada[i][4];
    m_hasGreenSawada[i] = true;
    if (m_gas[i] == "He" || m_gas[i] == "He-3") {
      m_parGreenSawada[i] = {15.5, 24.5, -2.25, ta, tb};
    } else if (m_gas[i] == "Ne") {
      m_parGreenSawada[i] = {24.3, 21.6, -6.49, ta, tb};
    } else if (m_gas[i] == "Ar") {
      m_parGreenSawada[i] = {6.92, 7.85, 6.87, ta, tb};
    } else if (m_gas[i] == "Kr") {
      m_parGreenSawada[i] = {7.95, 13.5, 3.90, ta, tb};
    } else if (m_gas[i] == "Xe") {
      m_parGreenSawada[i] = {7.93, 11.5, 3.81, ta, tb};
    } else if (m_gas[i] == "H2" || m_gas[i] == "D2") {
      m_parGreenSawada[i] = {7.07, 7.7, 1.87, ta, tb};
    } else if (m_gas[i] == "N2") {
      m_parGreenSawada[i] = {13.8, 15.6, 4.71, ta, tb};
    } else if (m_gas[i] == "O2") {
      m_parGreenSawada[i] = {18.5, 12.1, 1.86, ta, tb};
    } else if (m_gas[i] == "CH4") {
      m_parGreenSawada[i] = {7.06, 12.5, 3.45, ta, tb};
    } else if (m_gas[i] == "H2O") {
      m_parGreenSawada[i] = {12.8, 12.6, 1.28, ta, tb};
    } else if (m_gas[i] == "CO") {
      m_parGreenSawada[i] = {13.3, 14.0, 2.03, ta, tb};
    } else if (m_gas[i] == "C2H2") {
      m_parGreenSawada[i] = {9.28, 5.8, 1.37, ta, tb};
    } else if (m_gas[i] == "NO") {
      m_parGreenSawada[i] = {10.4, 9.5, -4.30, ta, tb};
    } else if (m_gas[i] == "CO2") {
      m_parGreenSawada[i] = {12.3, 13.8, -2.46, ta, tb};
    } else {
      m_parGreenSawada[i][3] = 0.;
      m_hasGreenSawada[i] = false;
      if (m_useGreenSawada) {
        std::cout << m_className << "::SetupGreenSawada:\n"
                  << "    Fit parameters for " << m_gas[i]
                  << " not available.\n"
                  << "    Opal-Beaty formula is used instead.\n";
      }
    }
  }
}

void MediumMagboltz::SetScatteringParameters(const int model,
                                             const double parIn, double& cut,
                                             double& parOut) const {
  cut = 1.;
  parOut = 0.5;
  if (model <= 0) return;

  if (model >= 2) {
    parOut = parIn;
    return;
  }

  // Set cuts on angular distribution and
  // renormalise forward scattering probability.

  if (parIn <= 1.) {
    parOut = parIn;
    return;
  }

  constexpr double rads = 2. / Pi;
  const double cns = parIn - 0.5;
  const double thetac = asin(2. * sqrt(cns - cns * cns));
  const double fac = (1. - cos(thetac)) / pow(sin(thetac), 2.);
  parOut = cns * fac + 0.5;
  cut = thetac * rads;
}

void MediumMagboltz::ComputeDeexcitationTable(const bool verbose) {
  m_iDeexcitation.fill(-1);
  m_deexcitations.clear();

  // Optical data (for quencher photoabsorption cs and ionisation yield)
  OpticalData optData;

  // Indices of "de-excitable" gases (only Ar for the time being).
  int iAr = -1;

  // Map Magboltz level names to internal ones.
  std::map<std::string, std::string> levelNamesAr = {
      {"1S5    ", "Ar_1S5"},     {"1S4    ", "Ar_1S4"},
      {"1S3    ", "Ar_1S3"},     {"1S2    ", "Ar_1S2"},
      {"2P10   ", "Ar_2P10"},    {"2P9    ", "Ar_2P9"},
      {"2P8    ", "Ar_2P8"},     {"2P7    ", "Ar_2P7"},
      {"2P6    ", "Ar_2P6"},     {"2P5    ", "Ar_2P5"},
      {"2P4    ", "Ar_2P4"},     {"2P3    ", "Ar_2P3"},
      {"2P2    ", "Ar_2P2"},     {"2P1    ", "Ar_2P1"},
      {"3D6    ", "Ar_3D6"},     {"3D5    ", "Ar_3D5"},
      {"3D3    ", "Ar_3D3"},     {"3D4!   ", "Ar_3D4!"},
      {"3D4    ", "Ar_3D4"},     {"3D1!!  ", "Ar_3D1!!"},
      {"2S5    ", "Ar_2S5"},     {"2S4    ", "Ar_2S4"},
      {"3D1!   ", "Ar_3D1!"},    {"3D2    ", "Ar_3D2"},
      {"3S1!!!!", "Ar_3S1!!!!"}, {"3S1!!  ", "Ar_3S1!!"},
      {"3S1!!! ", "Ar_3S1!!!"},  {"2S3    ", "Ar_2S3"},
      {"2S2    ", "Ar_2S2"},     {"3S1!   ", "Ar_3S1!"},
      {"4D5    ", "Ar_4D5"},     {"3S4    ", "Ar_3S4"},
      {"4D2    ", "Ar_4D2"},     {"4S1!   ", "Ar_4S1!"},
      {"3S2    ", "Ar_3S2"},     {"5D5    ", "Ar_5D5"},
      {"4S4    ", "Ar_4S4"},     {"5D2    ", "Ar_5D2"},
      {"6D5    ", "Ar_6D5"},     {"5S1!   ", "Ar_5S1!"},
      {"4S2    ", "Ar_4S2"},     {"5S4    ", "Ar_5S4"},
      {"6D2    ", "Ar_6D2"},     {"HIGH   ", "Ar_Higher"}};

  std::map<std::string, int> mapLevels;
  // Make a mapping of all excitation levels.
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    // Check if the level is an excitation.
    if (m_csType[i] % nCsTypes != ElectronCollisionTypeExcitation) continue;
    // Extract the index of the gas.
    const int ngas = int(m_csType[i] / nCsTypes);
    if (m_gas[ngas] == "Ar") {
      // Argon
      if (iAr < 0) iAr = ngas;
      // Get the level description (as specified in Magboltz).
      std::string level = "       ";
      for (int j = 0; j < 7; ++j) level[j] = m_description[i][5 + j];
      if (levelNamesAr.find(level) != levelNamesAr.end()) {
        mapLevels[levelNamesAr[level]] = i;
      } else {
        std::cerr << m_className << "::ComputeDeexcitationTable:\n"
                  << "    Unknown Ar excitation level: " << level << "\n";
      }
    }
  }

  // Count the excitation levels.
  unsigned int nDeexcitations = 0;
  std::map<std::string, int> lvl;
  for (auto it = mapLevels.cbegin(), end = mapLevels.cend(); it != end; ++it) {
    std::string level = (*it).first;
    lvl[level] = nDeexcitations;
    m_iDeexcitation[(*it).second] = nDeexcitations;
    ++nDeexcitations;
  }

  // Conversion factor from oscillator strength to transition rate.
  constexpr double f2A =
      2. * SpeedOfLight * FineStructureConstant / (3. * ElectronMass * HbarC);

  // Radiative de-excitation channels
  // Transition rates (unless indicated otherwise) are taken from:
  //     NIST Atomic Spectra Database
  // Transition rates for lines missing in the NIST database:
  //     O. Zatsarinny and K. Bartschat, J. Phys. B 39 (2006), 2145-2158
  // Oscillator strengths not included in the NIST database:
  //     J. Berkowitz, Atomic and Molecular Photoabsorption (2002)
  //     C.-M. Lee and K. T. Lu, Phys. Rev. A 8 (1973), 1241-1257
  for (auto it = mapLevels.cbegin(), end = mapLevels.cend(); it != end; ++it) {
    std::string level = (*it).first;
    Deexcitation dxc;
    dxc.gas = int(m_csType[(*it).second] / nCsTypes);
    dxc.level = (*it).second;
    dxc.label = level;
    // Excitation energy
    dxc.energy = m_energyLoss[(*it).second] * m_rgas[dxc.gas];
    // Oscillator strength
    dxc.osc = dxc.cf = 0.;
    dxc.sDoppler = dxc.gPressure = dxc.width = 0.;
    const std::vector<int> levelsAr4s = {lvl["Ar_1S5"], lvl["Ar_1S4"],
                                         lvl["Ar_1S3"], lvl["Ar_1S2"]};
    if (level == "Ar_1S5" || level == "Ar_1S3") {
      // Metastables
    } else if (level == "Ar_1S4") {
      dxc.osc = 0.0609;  // NIST
      // Berkowitz: f = 0.058
      dxc.p = {0.119};
      dxc.final = {-1};
    } else if (level == "Ar_1S2") {
      dxc.osc = 0.25;  // NIST
      // Berkowitz: 0.2214
      dxc.p = {0.51};
      dxc.final = {-1};
    } else if (level == "Ar_2P10") {
      dxc.p = {0.0189, 5.43e-3, 9.8e-4, 1.9e-4};
      dxc.final = levelsAr4s;
    } else if (level == "Ar_2P9") {
      dxc.p = {0.0331};
      dxc.final = {lvl["Ar_1S5"]};
    } else if (level == "Ar_2P8") {
      dxc.p = {9.28e-3, 0.0215, 1.47e-3};
      dxc.final = {lvl["Ar_1S5"], lvl["Ar_1S4"], lvl["Ar_1S2"]};
    } else if (level == "Ar_2P7") {
      dxc.p = {5.18e-3, 0.025, 2.43e-3, 1.06e-3};
      dxc.final = levelsAr4s;
    } else if (level == "Ar_2P6") {
      dxc.p = {0.0245, 4.9e-3, 5.03e-3};
      dxc.final = {lvl["Ar_1S5"], lvl["Ar_1S4"], lvl["Ar_1S2"]};
    } else if (level == "Ar_2P5") {
      dxc.p = {0.0402};
      dxc.final = {lvl["Ar_1S4"]};
    } else if (level == "Ar_2P4") {
      dxc.p = {6.25e-4, 2.2e-5, 0.0186, 0.0139};
      dxc.final = levelsAr4s;
    } else if (level == "Ar_2P3") {
      dxc.p = {3.8e-3, 8.47e-3, 0.0223};
      dxc.final = {lvl["Ar_1S5"], lvl["Ar_1S4"], lvl["Ar_1S2"]};
    } else if (level == "Ar_2P2") {
      dxc.p = {6.39e-3, 1.83e-3, 0.0117, 0.0153};
      dxc.final = levelsAr4s;
    } else if (level == "Ar_2P1") {
      dxc.p = {2.36e-4, 0.0445};
      dxc.final = {lvl["Ar_1S4"], lvl["Ar_1S2"]};
    } else if (level == "Ar_3D6") {
      // Additional line (2P7) from Bartschat
      dxc.p = {8.1e-3, 7.73e-4, 1.2e-4, 3.6e-4};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P7"], lvl["Ar_2P4"], lvl["Ar_2P2"]};
    } else if (level == "Ar_3D5") {
      dxc.osc = 0.0011;  // Berkowitz
      // Additional lines (2P7, 2P6, 2P5, 2P1) from Bartschat
      // Transition probability to ground state calculated from osc. strength
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {7.4e-3, 3.9e-5, 3.09e-4, 1.37e-3, 5.75e-4,
               3.2e-5, 1.4e-4, 1.7e-4,  2.49e-6, p0};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P8"],
                   lvl["Ar_2P7"],  lvl["Ar_2P6"],
                   lvl["Ar_2P5"],  lvl["Ar_2P4"],
                   lvl["Ar_2P3"],  lvl["Ar_2P2"],
                   lvl["Ar_2P1"],  -1};
    } else if (level == "Ar_3D3") {
      // Additional lines (2P9, 2P4) from Bartschat
      dxc.p = {4.9e-3, 9.82e-5, 1.2e-4, 2.6e-4,
               2.5e-3, 9.41e-5, 3.9e-4, 1.1e-4};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P7"],
                   lvl["Ar_2P6"],  lvl["Ar_2P4"], lvl["Ar_2P3"], lvl["Ar_2P2"]};
    } else if (level == "Ar_3D4!") {
      // Transition probability for 2P9 transition from Bartschat
      dxc.p = {0.01593};
      dxc.final = {lvl["Ar_2P9"]};
    } else if (level == "Ar_3D4") {
      // Additional lines (2P9, 2P3) from Bartschat
      dxc.p = {2.29e-3, 0.011, 8.8e-5, 2.53e-6};
      dxc.final = {lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P6"], lvl["Ar_2P3"]};
    } else if (level == "Ar_3D1!!") {
      // Additional lines (2P10, 2P6, 2P4 - 2P2) from Bartschat
      dxc.p = {5.85e-6, 1.2e-4,  5.7e-3,  7.3e-3,
               2.e-4,   1.54e-6, 2.08e-5, 6.75e-7};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P7"],
                   lvl["Ar_2P6"],  lvl["Ar_2P4"], lvl["Ar_2P3"], lvl["Ar_2P2"]};
    } else if (level == "Ar_2S5") {
      dxc.p = {4.9e-3, 0.011, 1.1e-3, 4.6e-4, 3.3e-3, 5.9e-5, 1.2e-4, 3.1e-4};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P7"],
                   lvl["Ar_2P6"],  lvl["Ar_2P4"], lvl["Ar_2P3"], lvl["Ar_2P2"]};
    } else if (level == "Ar_2S4") {
      dxc.osc = 0.027;  // NIST
      // Berkowitz: f = 0.026;
      dxc.p = {0.077,  2.44e-3, 8.9e-3, 4.6e-3, 2.7e-3,
               1.3e-3, 4.5e-4,  2.9e-5, 3.e-5,  1.6e-4};
      dxc.final = {-1,
                   lvl["Ar_2P10"],
                   lvl["Ar_2P8"],
                   lvl["Ar_2P7"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P4"],
                   lvl["Ar_2P3"],
                   lvl["Ar_2P2"],
                   lvl["Ar_2P1"]};
    } else if (level == "Ar_3D1!") {
      // Additional line (2P6) from Bartschat
      dxc.p = {3.1e-3, 2.e-3, 0.015, 9.8e-6};
      dxc.final = {lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P6"], lvl["Ar_2P3"]};
    } else if (level == "Ar_3D2") {
      dxc.osc = 0.0932;  // NIST
      // Berkowitz: f = 0.09
      // Additional lines (2P10, 2P6, 2P4-2P1) from Bartschat
      dxc.p = {0.27,   1.35e-5, 9.52e-4, 0.011,   4.01e-5,
               4.3e-3, 8.96e-4, 4.45e-5, 5.87e-5, 8.77e-4};
      dxc.final = {-1,
                   lvl["Ar_2P10"],
                   lvl["Ar_2P8"],
                   lvl["Ar_2P7"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P4"],
                   lvl["Ar_2P3"],
                   lvl["Ar_2P2"],
                   lvl["Ar_2P1"]};
    } else if (level == "Ar_3S1!!!!") {
      // Additional lines (2P10, 2P9, 2P7, 2P6, 2P2) from Bartschat
      dxc.p = {7.51e-6, 4.3e-5, 8.3e-4, 5.01e-5,
               2.09e-4, 0.013,  2.2e-3, 3.35e-6};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P7"],
                   lvl["Ar_2P6"],  lvl["Ar_2P4"], lvl["Ar_2P3"], lvl["Ar_2P2"]};
    } else if (level == "Ar_3S1!!") {
      // Additional lines (2P10 - 2P8, 2P4, 2P3)
      dxc.p = {1.89e-4, 1.52e-4, 7.21e-4, 3.69e-4,
               3.76e-3, 1.72e-4, 5.8e-4,  6.2e-3};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P7"],
                   lvl["Ar_2P6"],  lvl["Ar_2P4"], lvl["Ar_2P3"], lvl["Ar_2P2"]};
    } else if (level == "Ar_3S1!!!") {
      // Additional lines (2P9, 2P8, 2P6) from Bartschat
      dxc.p = {7.36e-4, 4.2e-5, 9.3e-5, 0.015};
      dxc.final = {lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P6"], lvl["Ar_2P3"]};
    } else if (level == "Ar_2S3") {
      dxc.p = {3.26e-3, 2.22e-3, 0.01, 5.1e-3};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P7"], lvl["Ar_2P4"], lvl["Ar_2P2"]};
    } else if (level == "Ar_2S2") {
      dxc.osc = 0.0119;  // NIST
      // Berkowitz: f = 0.012;
      dxc.p = {0.035,  1.76e-3, 2.1e-4, 2.8e-4, 1.39e-3,
               3.8e-4, 2.0e-3,  8.9e-3, 3.4e-3, 1.9e-3};
      dxc.final = {-1,
                   lvl["Ar_2P10"],
                   lvl["Ar_2P8"],
                   lvl["Ar_2P7"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P4"],
                   lvl["Ar_2P3"],
                   lvl["Ar_2P2"],
                   lvl["Ar_2P1"]};
    } else if (level == "Ar_3S1!") {
      dxc.osc = 0.106;  // NIST
      // Berkowitz: f = 0.106
      // Additional lines (2P10, 2P8, 2P7, 2P3) from Bartschat
      dxc.p = {0.313,  2.05e-5, 8.33e-5, 3.9e-4, 3.96e-4,
               4.2e-4, 4.5e-3,  4.84e-5, 7.1e-3, 5.2e-3};
      dxc.final = {-1,
                   lvl["Ar_2P10"],
                   lvl["Ar_2P8"],
                   lvl["Ar_2P7"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P4"],
                   lvl["Ar_2P3"],
                   lvl["Ar_2P2"],
                   lvl["Ar_2P1"]};
    } else if (level == "Ar_4D5") {
      dxc.osc = 0.0019;  // Berkowitz
      // Transition probability to ground state calculated from osc. strength
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {2.78e-3, 2.8e-4, 8.6e-4, 9.2e-4, 4.6e-4, 1.6e-4, p0};
      dxc.final = {lvl["Ar_2P10"],
                   lvl["Ar_2P8"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P3"],
                   lvl["Ar_2P2"],
                   -1};
    } else if (level == "Ar_3S4") {
      dxc.osc = 0.0144;  // Berkowitz
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {4.21e-4, 2.e-3,  1.7e-3, 7.2e-4, 3.5e-4,
               1.2e-4,  4.2e-6, 3.3e-5, 9.7e-5, p0};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P8"],
                   lvl["Ar_2P7"],  lvl["Ar_2P6"],
                   lvl["Ar_2P5"],  lvl["Ar_2P4"],
                   lvl["Ar_2P3"],  lvl["Ar_2P2"],
                   lvl["Ar_2P1"],  -1};
    } else if (level == "Ar_4D2") {
      dxc.osc = 0.048;  // Berkowitz
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {1.7e-4, p0};
      dxc.final = {lvl["Ar_2P7"], -1};
    } else if (level == "Ar_4S1!") {
      dxc.osc = 0.0209;  // Berkowitz
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {1.05e-3, 3.1e-5, 2.5e-5, 4.0e-4, 5.8e-5, 1.2e-4, p0};
      dxc.final = {lvl["Ar_2P10"],
                   lvl["Ar_2P8"],
                   lvl["Ar_2P7"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P3"],
                   -1};
    } else if (level == "Ar_3S2") {
      dxc.osc = 0.0221;  // Berkowitz
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {2.85e-4, 5.1e-5,  5.3e-5, 1.6e-4,  1.5e-4,
               6.0e-4,  2.48e-3, 9.6e-4, 3.59e-4, p0};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P8"],
                   lvl["Ar_2P7"],  lvl["Ar_2P6"],
                   lvl["Ar_2P5"],  lvl["Ar_2P4"],
                   lvl["Ar_2P3"],  lvl["Ar_2P2"],
                   lvl["Ar_2P1"],  -1};
    } else if (level == "Ar_5D5") {
      dxc.osc = 0.0041;  // Berkowitz
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {2.2e-3, 1.1e-4, 7.6e-5, 4.2e-4, 2.4e-4,
               2.1e-4, 2.4e-4, 1.2e-4, p0};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P8"], lvl["Ar_2P7"],
                   lvl["Ar_2P6"],  lvl["Ar_2P5"], lvl["Ar_2P4"],
                   lvl["Ar_2P3"],  lvl["Ar_2P2"], -1};
    } else if (level == "Ar_4S4") {
      dxc.osc = 0.0139;  // Berkowitz
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {1.9e-4, 1.1e-3, 5.2e-4, 5.1e-4, 9.4e-5, 5.4e-5, p0};
      dxc.final = {lvl["Ar_2P10"],
                   lvl["Ar_2P8"],
                   lvl["Ar_2P7"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P4"],
                   -1};
    } else if (level == "Ar_5D2") {
      dxc.osc = 0.0426;  // Berkowitz
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {5.9e-5, 9.0e-6, 1.5e-4, 3.1e-5, p0};
      dxc.final = {lvl["Ar_2P8"], lvl["Ar_2P7"], lvl["Ar_2P5"], lvl["Ar_2P2"],
                   -1};
    } else if (level == "Ar_6D5") {
      dxc.osc = 0.00075;  // Lee and Lu
      // Berkowitz estimates f = 0.0062 for the sum of
      // all "weak" nd levels with n = 6 and higher.
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {1.9e-3, 4.2e-4, 3.e-4, 5.1e-5, 6.6e-5, 1.21e-4, p0};
      dxc.final = {lvl["Ar_2P10"],
                   lvl["Ar_2P6"],
                   lvl["Ar_2P5"],
                   lvl["Ar_2P4"],
                   lvl["Ar_2P3"],
                   lvl["Ar_2P1"],
                   -1};
    } else if (level == "Ar_5S1!") {
      dxc.osc = 0.00051;  // Lee and Lu
      // Berkowitz estimates f = 0.0562 for the sum
      // of all nd' levels with n = 5 and higher.
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {7.7e-5, p0};
      dxc.final = {lvl["Ar_2P5"], -1};
    } else if (level == "Ar_4S2") {
      dxc.osc = 0.00074;  // Lee and Lu
      // Berkowitz estimates f = 0.0069 for the sum over all
      // ns' levels with n = 7 and higher.
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {4.5e-4, 2.e-4, 2.1e-4, 1.2e-4, 1.8e-4, 9.e-4, 3.3e-4, p0};
      dxc.final = {lvl["Ar_2P10"], lvl["Ar_2P8"], lvl["Ar_2P7"], lvl["Ar_2P5"],
                   lvl["Ar_2P4"],  lvl["Ar_2P3"], lvl["Ar_2P2"], -1};
    } else if (level == "Ar_5S4") {
      // dxc.osc = 0.0130; // Lee and Lu
      // Berkowitz estimates f = 0.0211 for the sum of all
      // ns levels with n = 8 and higher.
      dxc.osc = 0.0211;
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {3.6e-4, 1.2e-4, 1.5e-4, 1.4e-4, 7.5e-5, p0};
      dxc.final = {lvl["Ar_2P8"], lvl["Ar_2P6"], lvl["Ar_2P4"],
                   lvl["Ar_2P3"], lvl["Ar_2P2"], -1};
    } else if (level == "Ar_6D2") {
      // dxc.osc = 0.0290; // Lee and Lu
      // Berkowitz estimates f = 0.0574 for the sum of all
      // "strong" nd levels with n = 6 and higher.
      dxc.osc = 0.0574;
      // Additional line: 2P7
      const double p0 = f2A * dxc.energy * dxc.energy * dxc.osc;
      dxc.p = {3.33e-3, p0};
      dxc.final = {lvl["Ar_2P7"], -1};
    } else if (level == "Ar_Higher") {
      dxc.osc = 0.;
      // This (artificial) level represents the sum of higher J = 1 states.
      // The deeexcitation cascade is simulated by allocating it
      // with equal probability to one of the five nearest levels below.
      dxc.type.assign(5, DxcTypeCollNonIon);
      dxc.p = {100., 100., 100., 100., 100.};
      dxc.final = {lvl["Ar_6D5"], lvl["Ar_5S1!"], lvl["Ar_4S2"], lvl["Ar_5S4"],
                   lvl["Ar_6D2"]};
    } else {
      std::cerr << m_className << "::ComputeDeexcitationTable:\n"
                << "    Missing de-excitation data for level " << level
                << ". Program bug!\n";
      return;
    }
    if (level != "Ar_Higher") dxc.type.assign(dxc.p.size(), DxcTypeRad);
    m_deexcitations.push_back(std::move(dxc));
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::ComputeDeexcitationTable:\n";
    std::cout << "    Found " << m_deexcitations.size() << " levels "
              << "with available radiative de-excitation data.\n";
  }

  // Collisional de-excitation channels
  if (iAr >= 0) {
    // Add the Ar dimer ground state.
    Deexcitation dimer;
    dimer.label = "Ar_Dimer";
    dimer.level = -1;
    dimer.gas = iAr;
    dimer.energy = 14.71;
    dimer.osc = dimer.cf = 0.;
    dimer.sDoppler = dimer.gPressure = dimer.width = 0.;
    lvl["Ar_Dimer"] = m_deexcitations.size();
    m_deexcitations.push_back(std::move(dimer));
    ++nDeexcitations;
    // Add an Ar excimer level.
    Deexcitation excimer;
    excimer.label = "Ar_Excimer";
    excimer.level = -1;
    excimer.gas = iAr;
    excimer.energy = 14.71;
    excimer.osc = excimer.cf = 0.;
    excimer.sDoppler = excimer.gPressure = excimer.width = 0.;
    lvl["Ar_Excimer"] = m_deexcitations.size();
    m_deexcitations.push_back(std::move(excimer));
    ++nDeexcitations;
    const double nAr = GetNumberDensity() * m_fraction[iAr];
    // Flags for two-body and three-body collision rate constants.
    // Three-body collisions lead to excimer formation.
    // Two-body collisions give rise to collisional mixing.
    constexpr bool useTachibanaData = false;
    constexpr bool useCollMixing = true;
    for (auto& dxc : m_deexcitations) {
      const std::string level = dxc.label;
      if (level == "Ar_1S5") {
        // K. Tachibana, Phys. Rev. A 34 (1986), 1007-1015
        // Kolts and Setser, J. Chem. Phys. 68 (1978), 4848-4859
        constexpr double k3b = useTachibanaData ? 1.4e-41 : 1.1e-41;
        dxc.p.push_back(k3b * nAr * nAr);
        dxc.final.push_back(lvl["Ar_Excimer"]);
        if (useCollMixing) {
          constexpr double k2b = useTachibanaData ? 2.3e-24 : 2.1e-24;
          dxc.p.push_back(k2b * nAr);
          dxc.final.push_back(lvl["Ar_1S4"]);
          dxc.type.push_back(DxcTypeCollNonIon);
        }
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_1S3") {
        // K. Tachibana, Phys. Rev. A 34 (1986), 1007-1015
        // Kolts and Setser, J. Chem. Phys. 68 (1978), 4848-4859
        constexpr double k3b = useTachibanaData ? 1.5e-41 : 0.83e-41;
        dxc.p.push_back(k3b * nAr * nAr);
        dxc.final.push_back(lvl["Ar_Excimer"]);
        if (useCollMixing) {
          constexpr double k2b = useTachibanaData ? 4.3e-24 : 5.3e-24;
          dxc.p.push_back(k2b * nAr);
          dxc.final.push_back(lvl["Ar_1S4"]);
        }
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      }
      const std::vector<int> levels4s = {lvl["Ar_1S5"], lvl["Ar_1S4"],
                                         lvl["Ar_1S3"], lvl["Ar_1S2"]};
      if (level == "Ar_2P1") {
        // Transfer to 4s states
        // Inoue, Setser, and Sadeghi, J. Chem. Phys. 75 (1982), 977-983
        // constexpr double k4s = 2.9e-20;
        // Sadeghi et al. J. Chem. Phys. 115 (2001), 3144-3154
        constexpr double k4s = 1.6e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P2") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k23 = 0.5e-21;
        dxc.p.push_back(k23 * nAr);
        dxc.final.push_back(lvl["Ar_2P3"]);
        // Transfer to 4s states
        // Inoue, Setser, and Sadeghi, J. Chem. Phys. 75 (1982), 977-983
        // constexpr double k4s = 3.8e-20;
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        constexpr double k4s = 5.3e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P3") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k34 = 27.5e-21;
        constexpr double k35 = 0.3e-21;
        constexpr double k36 = 44.0e-21;
        constexpr double k37 = 1.4e-21;
        constexpr double k38 = 1.9e-21;
        constexpr double k39 = 0.8e-21;
        dxc.p.insert(dxc.p.end(), {k34 * nAr, k35 * nAr, k36 * nAr, k37 * nAr,
                                   k38 * nAr, k39 * nAr});
        dxc.final.insert(dxc.final.end(),
                         {lvl["Ar_2P4"], lvl["Ar_2P5"], lvl["Ar_2P6"],
                          lvl["Ar_2P7"], lvl["Ar_2P8"], lvl["Ar_2P9"]});
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        constexpr double k4s = 4.7e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P4") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k43 = 23.0e-21;
        constexpr double k45 = 0.7e-21;
        constexpr double k46 = 4.8e-21;
        constexpr double k47 = 3.2e-21;
        constexpr double k48 = 1.4e-21;
        constexpr double k49 = 3.3e-21;
        dxc.p.insert(dxc.p.end(), {k43 * nAr, k45 * nAr, k46 * nAr, k47 * nAr,
                                   k48 * nAr, k49 * nAr});
        dxc.final.insert(dxc.final.end(),
                         {lvl["Ar_2P3"], lvl["Ar_2P5"], lvl["Ar_2P6"],
                          lvl["Ar_2P7"], lvl["Ar_2P8"], lvl["Ar_2P9"]});
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        constexpr double k4s = 3.9e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P5") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k54 = 1.7e-21;
        constexpr double k56 = 11.3e-21;
        constexpr double k58 = 9.5e-21;
        dxc.p.insert(dxc.p.end(), {k54 * nAr, k56 * nAr, k58 * nAr});
        dxc.final.insert(dxc.final.end(),
                         {lvl["Ar_2P4"], lvl["Ar_2P6"], lvl["Ar_2P8"]});
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P6") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k67 = 4.1e-21;
        constexpr double k68 = 6.0e-21;
        constexpr double k69 = 1.0e-21;
        dxc.p.insert(dxc.p.end(), {k67 * nAr, k68 * nAr, k69 * nAr});
        dxc.final.insert(dxc.final.end(),
                         {lvl["Ar_2P7"], lvl["Ar_2P8"], lvl["Ar_2P9"]});
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P7") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k76 = 2.5e-21;
        constexpr double k78 = 14.3e-21;
        constexpr double k79 = 23.3e-21;
        dxc.p.insert(dxc.p.end(), {k76 * nAr, k78 * nAr, k79 * nAr});
        dxc.final.insert(dxc.final.end(),
                         {lvl["Ar_2P6"], lvl["Ar_2P8"], lvl["Ar_2P9"]});
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        constexpr double k4s = 5.5e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P8") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k86 = 0.3e-21;
        constexpr double k87 = 0.8e-21;
        constexpr double k89 = 18.2e-21;
        constexpr double k810 = 1.0e-21;
        dxc.p.insert(dxc.p.end(),
                     {k86 * nAr, k87 * nAr, k89 * nAr, k810 * nAr});
        dxc.final.insert(dxc.final.end(), {lvl["Ar_2P6"], lvl["Ar_2P7"],
                                           lvl["Ar_2P9"], lvl["Ar_2P10"]});
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        constexpr double k4s = 3.e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P9") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        constexpr double k98 = 6.8e-21;
        constexpr double k910 = 5.1e-21;
        dxc.p.insert(dxc.p.end(), {k98 * nAr, k910 * nAr});
        dxc.final.insert(dxc.final.end(), {lvl["Ar_2P8"], lvl["Ar_2P10"]});
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        constexpr double k4s = 3.5e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_2P10") {
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        constexpr double k4s = 2.0e-20;
        dxc.p.resize(dxc.p.size() + levels4s.size(), 0.25 * k4s * nAr);
        dxc.final.insert(dxc.final.end(), levels4s.begin(), levels4s.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      }
      const std::vector<int> levels4p = {
          lvl["Ar_2P10"], lvl["Ar_2P9"], lvl["Ar_2P8"], lvl["Ar_2P7"],
          lvl["Ar_2P6"],  lvl["Ar_2P5"], lvl["Ar_2P4"], lvl["Ar_2P3"],
          lvl["Ar_2P2"],  lvl["Ar_2P1"]};
      if (level == "Ar_3D6" || level == "Ar_3D5" || level == "Ar_3D3" ||
          level == "Ar_3D4!" || level == "Ar_3D4" || level == "Ar_3D1!!" ||
          level == "Ar_3D1!" || level == "Ar_3D2" || level == "Ar_3S1!!!!" ||
          level == "Ar_3S1!!" || level == "Ar_3S1!!!" || level == "Ar_3S1!" ||
          level == "Ar_2S5" || level == "Ar_2S4" || level == "Ar_2S3" ||
          level == "Ar_2S2") {
        // 3d and 5s levels
        // Transfer to 4p levels
        // Parameter to be tuned (order of magnitude guess).
        constexpr double k4p = 1.e-20;
        dxc.p.resize(dxc.p.size() + levels4p.size(), 0.1 * k4p * nAr);
        dxc.final.insert(dxc.final.end(), levels4p.begin(), levels4p.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      } else if (level == "Ar_4D5" || level == "Ar_3S4" || level == "Ar_4D2" ||
                 level == "Ar_4S1!" || level == "Ar_3S2" || level == "Ar_5D5" ||
                 level == "Ar_4S4" || level == "Ar_5D2" || level == "Ar_6D5" ||
                 level == "Ar_5S1!" || level == "Ar_4S2" || level == "Ar_5S4" ||
                 level == "Ar_6D2") {
        // Transfer to 4p levels
        // Parameter to be tuned (order of magnitude guess).
        constexpr double k4p = 1.e-20;
        dxc.p.resize(dxc.p.size() + levels4p.size(), 0.1 * k4p * nAr);
        dxc.final.insert(dxc.final.end(), levels4p.begin(), levels4p.end());
        dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
        // Hornbeck-Molnar ionisation
        // P. Becker and F. Lampe, J. Chem. Phys. 42 (1965), 3857-3863
        // A. Bogaerts and R. Gijbels, Phys. Rev. A 52 (1995), 3743-3751
        // This value seems high, to be checked!
        constexpr double kHM = 2.e-18;
        constexpr bool useHornbeckMolnar = true;
        if (useHornbeckMolnar) {
          dxc.p.push_back(kHM * nAr);
          dxc.final.push_back(lvl["Ar_Dimer"]);
          dxc.type.push_back(DxcTypeCollIon);
        }
      }
    }
  }

  // Collisional deexcitation by quenching gases.
  int iCO2 = -1;
  int iCH4 = -1;
  int iC2H6 = -1;
  int iIso = -1;
  int iC2H2 = -1;
  int iCF4 = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == "CO2")
      iCO2 = i;
    else if (m_gas[i] == "CH4")
      iCH4 = i;
    else if (m_gas[i] == "C2H6")
      iC2H6 = i;
    else if (m_gas[i] == "C2H2")
      iC2H2 = i;
    else if (m_gas[i] == "CF4")
      iCF4 = i;
    else if (m_gas[i] == "iC4H10")
      iIso = i;
  }

  // Collision radii for hard-sphere approximation.
  constexpr double rAr3d = 436.e-10;
  constexpr double rAr5s = 635.e-10;

  if (iAr >= 0 && iCO2 >= 0) {
    // Partial density of CO2
    const double nQ = GetNumberDensity() * m_fraction[iCO2];
    // Collision radius
    constexpr double rCO2 = 165.e-10;
    for (auto& dxc : m_deexcitations) {
      std::string level = dxc.label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      optData.GetPhotoabsorptionCrossSection("CO2", dxc.energy, pacs, eta);
      const double pPenningWK = pow(eta, 0.4);
      if (level == "Ar_1S5") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        constexpr double kQ = 5.3e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        constexpr double kQ = 5.0e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_1S3") {
        constexpr double kQ = 5.9e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_1S2") {
        constexpr double kQ = 7.4e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        constexpr double kQ = 6.4e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_2P6") {
        // Rate constant from Sadeghi et al.
        constexpr double kQ = 6.1e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_2P5") {
        // Rate constant from Sadeghi et al.
        constexpr double kQ = 6.6e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_2P1") {
        // Rate constant from Sadeghi et al.
        constexpr double kQ = 6.2e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of 4p rate constants from Sadeghi et al.
        constexpr double kQ = 6.33e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (dxc.osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double kQ = RateConstantWK(dxc.energy, dxc.osc, pacs, iAr, iCO2);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        const double kQ = RateConstantHardSphere(rAr3d, rCO2, iAr, iCO2);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        const double kQ = RateConstantHardSphere(rAr5s, rCO2, iAr, iCO2);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      }
      dxc.final.resize(dxc.p.size(), -1);
    }
  }
  if (iAr >= 0 && iCH4 >= 0) {
    // Partial density of methane
    const double nQ = GetNumberDensity() * m_fraction[iCH4];
    // Collision radius
    constexpr double rCH4 = 190.e-10;
    for (auto& dxc : m_deexcitations) {
      std::string level = dxc.label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      optData.GetPhotoabsorptionCrossSection("CH4", dxc.energy, pacs, eta);
      const double pPenningWK = pow(eta, 0.4);
      if (level == "Ar_1S5") {
        // Rate constant from Chen and Setser, J. Phys. Chem. 95 (1991)
        constexpr double kQ = 4.55e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        constexpr double kQ = 4.5e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_1S3") {
        // Rate constant from Chen and Setser
        constexpr double kQ = 5.30e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_1S2") {
        // Rate constant from Velazco et al.
        constexpr double kQ = 5.7e-19;
        dxc.p.push_back(kQ * nQ);
        dxc.type.push_back(DxcTypeCollNonIon);
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        constexpr double kQ = 7.4e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P6") {
        constexpr double kQ = 3.4e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P5") {
        constexpr double kQ = 6.0e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P1") {
        constexpr double kQ = 9.3e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of rate constants given by Sadeghi et al.
        constexpr double kQ = 6.53e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (dxc.osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double kQ = RateConstantWK(dxc.energy, dxc.osc, pacs, iAr, iCH4);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        const double kQ = RateConstantHardSphere(rAr3d, rCH4, iAr, iCH4);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        const double kQ = RateConstantHardSphere(rAr5s, rCH4, iAr, iCH4);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      }
      dxc.final.resize(dxc.p.size(), -1);
    }
  }
  if (iAr >= 0 && iC2H6 >= 0) {
    // Partial density of ethane
    const double nQ = GetNumberDensity() * m_fraction[iC2H6];
    // Collision radius
    constexpr double rC2H6 = 195.e-10;
    for (auto& dxc : m_deexcitations) {
      std::string level = dxc.label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      optData.GetPhotoabsorptionCrossSection("C2H6", dxc.energy, pacs, eta);
      const double pPenningWK = pow(eta, 0.4);
      if (level == "Ar_1S5") {
        // Rate constant from Chen and Setser, J. Phys. Chem. 95 (1991)
        constexpr double kQ = 5.29e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        constexpr double kQ = 6.2e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_1S3") {
        // Rate constant from Chen and Setser
        constexpr double kQ = 6.53e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_1S2") {
        // Rate constant from Velazco et al.
        constexpr double kQ = 10.7e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        constexpr double kQ = 9.2e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P6") {
        constexpr double kQ = 4.8e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P5") {
        constexpr double kQ = 9.9e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P1") {
        constexpr double kQ = 11.0e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of rate constants given by Sadeghi et al.
        constexpr double kQ = 8.7e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (dxc.osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double kQ = RateConstantWK(dxc.energy, dxc.osc, pacs, iAr, iC2H6);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        const double kQ = RateConstantHardSphere(rAr3d, rC2H6, iAr, iC2H6);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        const double kQ = RateConstantHardSphere(rAr5s, rC2H6, iAr, iC2H6);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      }
      dxc.final.resize(dxc.p.size(), -1);
    }
  }
  if (iAr >= 0 && iIso >= 0) {
    // Partial density of isobutane
    const double nQ = GetNumberDensity() * m_fraction[iIso];
    // Collision radius
    constexpr double rIso = 250.e-10;
    // For the 4p levels, the rate constants are estimated by scaling
    // the values for ethane.
    // Ar radius [pm]
    constexpr double r4p = 340.;
    // Molecular radii are 195 pm for ethane, 250 pm for isobutane.
    constexpr double fr = (r4p + 250.) / (r4p + 195.);
    // Masses [amu]
    constexpr double mAr = 39.9;
    constexpr double mEth = 30.1;
    constexpr double mIso = 58.1;
    // Scaling factor.
    const double f4p =
        fr * fr * sqrt((mEth / mIso) * (mAr + mIso) / (mAr + mEth));
    for (auto& dxc : m_deexcitations) {
      std::string level = dxc.label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      // Use n-butane as approximation for isobutane.
      optData.GetPhotoabsorptionCrossSection("nC4H10", dxc.energy, pacs, eta);
      const double pPenningWK = pow(eta, 0.4);
      if (level == "Ar_1S5") {
        // Rate constant from
        // Piper et al., J. Chem. Phys. 59 (1973), 3323-3340
        constexpr double kQ = 7.1e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_1S4") {
        // Rate constant from Piper et al.
        constexpr double kQ = 6.1e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_1S3") {
        // Rate constant for n-butane from
        // Velazco et al., J. Chem. Phys. 69 (1978)
        constexpr double kQ = 8.5e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_1S2") {
        // Rate constant from Piper et al.
        constexpr double kQ = 11.0e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P8") {
        constexpr double kEth = 9.2e-19;  // ethane
        AddPenningDeexcitation(dxc, f4p * kEth * nQ, pPenningWK);
      } else if (level == "Ar_2P6") {
        constexpr double kEth = 4.8e-19;  // ethane
        AddPenningDeexcitation(dxc, f4p * kEth * nQ, pPenningWK);
      } else if (level == "Ar_2P5") {
        const double kEth = 9.9e-19;  // ethane
        AddPenningDeexcitation(dxc, f4p * kEth * nQ, pPenningWK);
      } else if (level == "Ar_2P1") {
        const double kEth = 11.0e-19;  // ethane
        AddPenningDeexcitation(dxc, f4p * kEth * nQ, pPenningWK);
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        constexpr double kEth = 5.5e-19;  // ethane
        AddPenningDeexcitation(dxc, f4p * kEth * nQ, pPenningWK);
      } else if (dxc.osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double kQ = RateConstantWK(dxc.energy, dxc.osc, pacs, iAr, iIso);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        const double kQ = RateConstantHardSphere(rAr3d, rIso, iAr, iIso);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        const double kQ = RateConstantHardSphere(rAr5s, rIso, iAr, iIso);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      }
      dxc.final.resize(dxc.p.size(), -1);
    }
  }
  if (iAr >= 0 && iC2H2 >= 0) {
    // Partial density of acetylene
    const double nQ = GetNumberDensity() * m_fraction[iC2H2];
    // Collision radius
    constexpr double rC2H2 = 165.e-10;
    for (auto& dxc : m_deexcitations) {
      std::string level = dxc.label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      optData.GetPhotoabsorptionCrossSection("C2H2", dxc.energy, pacs, eta);
      const double pPenningWK = pow(eta, 0.4);
      if (level == "Ar_1S5") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        constexpr double kQ = 5.6e-19;
        // Branching ratio for ionization according to
        // Jones et al., J. Phys. Chem. 89 (1985)
        // p = 0.61, p = 0.74 (agrees roughly with WK estimate)
        AddPenningDeexcitation(dxc, kQ * nQ, 0.61);
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al.
        constexpr double kQ = 4.6e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_1S3") {
        constexpr double kQ = 5.6e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, 0.61);
      } else if (level == "Ar_1S2") {
        // Rate constant from Velazco et al.
        constexpr double kQ = 8.7e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        constexpr double kQ = 5.0e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, 0.3);
      } else if (level == "Ar_2P6") {
        constexpr double kQ = 5.7e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, 0.3);
      } else if (level == "Ar_2P5") {
        constexpr double kQ = 6.0e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, 0.3);
      } else if (level == "Ar_2P1") {
        constexpr double kQ = 5.3e-19;  // Sadeghi
        AddPenningDeexcitation(dxc, kQ * nQ, 0.3);
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of rate constants given by Sadeghi et al.
        constexpr double kQ = 5.5e-19;
        AddPenningDeexcitation(dxc, kQ * nQ, 0.3);
      } else if (dxc.osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double kQ = RateConstantWK(dxc.energy, dxc.osc, pacs, iAr, iC2H2);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        const double kQ = RateConstantHardSphere(rAr3d, rC2H2, iAr, iC2H2);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        const double kQ = RateConstantHardSphere(rAr5s, rC2H2, iAr, iC2H2);
        AddPenningDeexcitation(dxc, kQ * nQ, pPenningWK);
      }
      dxc.final.resize(dxc.p.size(), -1);
    }
  }
  if (iAr >= 0 && iCF4 >= 0) {
    // Partial density of CF4
    const double nQ = GetNumberDensity() * m_fraction[iCF4];
    // Collision radius
    constexpr double rCF4 = 235.e-10;
    for (auto& dxc : m_deexcitations) {
      std::string level = dxc.label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      optData.GetPhotoabsorptionCrossSection("CF4", dxc.energy, pacs, eta);
      if (level == "Ar_1S5") {
        // Rate constant from Chen and Setser
        constexpr double kQ = 0.33e-19;
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_1S3") {
        // Rate constant from Chen and Setser
        constexpr double kQ = 0.26e-19;
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al.
        constexpr double kQ = 1.7e-19;
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_2P6") {
        constexpr double kQ = 1.7e-19;  // Sadeghi
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_2P5") {
        constexpr double kQ = 1.6e-19;  // Sadeghi
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_2P1") {
        constexpr double kQ = 2.2e-19;  // Sadeghi
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of 4p rate constants from Sadeghi et al.
        constexpr double kQ = 1.8e-19;
        dxc.p.push_back(kQ * nQ);
      } else if (dxc.osc > 0.) {
        // Resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double kQ = RateConstantWK(dxc.energy, dxc.osc, pacs, iAr, iCF4);
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        const double kQ = RateConstantHardSphere(rAr3d, rCF4, iAr, iCF4);
        dxc.p.push_back(kQ * nQ);
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        const double kQ = RateConstantHardSphere(rAr5s, rCF4, iAr, iCF4);
        dxc.p.push_back(kQ * nQ);
      }
      dxc.type.resize(dxc.p.size(), DxcTypeCollNonIon);
      dxc.final.resize(dxc.p.size(), -1);
    }
  }

  if ((m_debug || verbose) && nDeexcitations > 0) {
    std::cout << m_className << "::ComputeDeexcitationTable:\n"
              << "      Level  Energy [eV]                    Lifetimes [ns]\n"
              << "                            Total    Radiative       "
              << "     Collisional\n"
              << "                               "
              << "                Ionisation  Transfer      Loss\n";
  }

  for (auto& dxc : m_deexcitations) {
    // Calculate the total decay rate of each level.
    dxc.rate = 0.;
    double fRad = 0.;
    double fCollIon = 0., fCollTransfer = 0., fCollLoss = 0.;
    const unsigned int nChannels = dxc.type.size();
    for (unsigned int j = 0; j < nChannels; ++j) {
      dxc.rate += dxc.p[j];
      if (dxc.type[j] == DxcTypeRad) {
        fRad += dxc.p[j];
      } else if (dxc.type[j] == DxcTypeCollIon) {
        fCollIon += dxc.p[j];
      } else if (dxc.type[j] == DxcTypeCollNonIon) {
        if (dxc.final[j] < 0) {
          fCollLoss += dxc.p[j];
        } else {
          fCollTransfer += dxc.p[j];
        }
      } else {
        std::cerr << m_className << "::ComputeDeexcitationTable:\n    "
                  << "Unknown type of deexcitation channel (level " << dxc.label
                  << "). Program bug!\n";
      }
    }
    if (dxc.rate > 0.) {
      // Print the radiative and collisional decay rates.
      if (m_debug || verbose) {
        std::cout << std::setw(12) << dxc.label << "  " << std::fixed
                  << std::setprecision(3) << std::setw(7) << dxc.energy << "  "
                  << std::setw(10) << 1. / dxc.rate << "  ";
        if (fRad > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fRad << " ";
        } else {
          std::cout << "---------- ";
        }
        if (fCollIon > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fCollIon << " ";
        } else {
          std::cout << "---------- ";
        }
        if (fCollTransfer > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fCollTransfer << " ";
        } else {
          std::cout << "---------- ";
        }
        if (fCollLoss > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fCollLoss << "\n";
        } else {
          std::cout << "---------- \n";
        }
      }
      // Normalise the decay branching ratios.
      for (unsigned int j = 0; j < nChannels; ++j) {
        dxc.p[j] /= dxc.rate;
        if (j > 0) dxc.p[j] += dxc.p[j - 1];
      }
    }
  }
}

double MediumMagboltz::RateConstantWK(const double energy, const double osc,
                                      const double pacs, const int igas1,
                                      const int igas2) const {
  // Calculate rate constant from Watanabe-Katsuura formula.
  const double m1 = ElectronMassGramme / (m_rgas[igas1] - 1.);
  const double m2 = ElectronMassGramme / (m_rgas[igas2] - 1.);
  // Compute the reduced mass.
  double mR = (m1 * m2 / (m1 + m2)) / AtomicMassUnit;
  const double uA = (RydbergEnergy / energy) * osc;
  const double uQ = (2 * RydbergEnergy / energy) * pacs /
                    (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
  return 2.591e-19 * pow(uA * uQ, 0.4) * pow(m_temperature / mR, 0.3);
}

double MediumMagboltz::RateConstantHardSphere(const double r1, const double r2,
                                              const int igas1,
                                              const int igas2) const {
  // Hard sphere cross-section
  const double r = r1 + r2;
  const double sigma = r * r * Pi;
  // Reduced mass
  const double m1 = ElectronMass / (m_rgas[igas1] - 1.);
  const double m2 = ElectronMass / (m_rgas[igas2] - 1.);
  const double mR = m1 * m2 / (m1 + m2);
  // Relative velocity
  const double vel =
      SpeedOfLight * sqrt(8. * BoltzmannConstant * m_temperature / (Pi * mR));
  return sigma * vel;
}

void MediumMagboltz::ComputeDeexcitation(int iLevel, int& fLevel) {
  if (!m_useDeexcitation) {
    std::cerr << m_className << "::ComputeDeexcitation: Not enabled.\n";
    return;
  }

  // Make sure that the tables are updated.
  if (m_isChanged) {
    if (!Mixer()) {
      PrintErrorMixer(m_className + "::ComputeDeexcitation");
      return;
    }
    m_isChanged = false;
  }

  if (iLevel < 0 || iLevel >= (int)m_nTerms) {
    std::cerr << m_className << "::ComputeDeexcitation: Index out of range.\n";
    return;
  }

  iLevel = m_iDeexcitation[iLevel];
  if (iLevel < 0 || iLevel >= (int)m_deexcitations.size()) {
    std::cerr << m_className << "::ComputeDeexcitation:\n"
              << "    Level is not deexcitable.\n";
    return;
  }

  ComputeDeexcitationInternal(iLevel, fLevel);
  if (fLevel >= 0 && fLevel < (int)m_deexcitations.size()) {
    fLevel = m_deexcitations[fLevel].level;
  }
}

void MediumMagboltz::ComputeDeexcitationInternal(int iLevel, int& fLevel) {
  m_dxcProducts.clear();

  double t = 0.;
  fLevel = iLevel;
  while (iLevel >= 0 && iLevel < (int)m_deexcitations.size()) {
    const auto& dxc = m_deexcitations[iLevel];
    const int nChannels = dxc.p.size();
    if (dxc.rate <= 0. || nChannels <= 0) {
      // This level is a dead end.
      fLevel = iLevel;
      return;
    }
    // Determine the de-excitation time.
    t += -log(RndmUniformPos()) / dxc.rate;
    // Select the transition.
    fLevel = -1;
    int type = DxcTypeRad;
    const double r = RndmUniform();
    for (int j = 0; j < nChannels; ++j) {
      if (r <= dxc.p[j]) {
        fLevel = dxc.final[j];
        type = dxc.type[j];
        break;
      }
    }
    if (type == DxcTypeRad) {
      // Radiative decay
      dxcProd photon;
      photon.s = 0.;
      photon.t = t;
      photon.type = DxcProdTypePhoton;
      photon.energy = dxc.energy;
      if (fLevel >= 0) {
        // Decay to a lower lying excited state.
        photon.energy -= m_deexcitations[fLevel].energy;
        if (photon.energy < Small) photon.energy = Small;
        m_dxcProducts.push_back(std::move(photon));
        // Proceed with the next level in the cascade.
        iLevel = fLevel;
      } else {
        // Decay to ground state.
        double delta = RndmVoigt(0., dxc.sDoppler, dxc.gPressure);
        while (photon.energy + delta < Small || fabs(delta) >= dxc.width) {
          delta = RndmVoigt(0., dxc.sDoppler, dxc.gPressure);
        }
        photon.energy += delta;
        m_dxcProducts.push_back(std::move(photon));
        // Deexcitation cascade is over.
        fLevel = iLevel;
        return;
      }
    } else if (type == DxcTypeCollIon) {
      // Ionisation electron
      dxcProd electron;
      electron.s = 0.;
      electron.t = t;
      electron.type = DxcProdTypeElectron;
      electron.energy = dxc.energy;
      if (fLevel >= 0) {
        // Associative ionisation
        electron.energy -= m_deexcitations[fLevel].energy;
        if (electron.energy < Small) electron.energy = Small;
        ++m_nPenning;
        m_dxcProducts.push_back(std::move(electron));
        // Proceed with the next level in the cascade.
        iLevel = fLevel;
      } else {
        // Penning ionisation
        electron.energy -= m_minIonPot;
        if (electron.energy < Small) electron.energy = Small;
        ++m_nPenning;
        m_dxcProducts.push_back(std::move(electron));
        // Deexcitation cascade is over.
        fLevel = iLevel;
        return;
      }
    } else if (type == DxcTypeCollNonIon) {
      // Proceed with the next level in the cascade.
      iLevel = fLevel;
    } else {
      std::cerr << m_className << "::ComputeDeexcitationInternal:\n"
                << "    Unknown deexcitation type (" << type << "). Bug!\n";
      // Abort the calculation.
      fLevel = iLevel;
      return;
    }
  }
}

bool MediumMagboltz::ComputePhotonCollisionTable(const bool verbose) {
  OpticalData data;
  double cs;
  double eta;

  // Atomic density
  const double dens = GetNumberDensity();

  // Reset the collision rate arrays.
  m_cfTotGamma.assign(nEnergyStepsGamma, 0.);
  m_cfGamma.assign(nEnergyStepsGamma, std::vector<double>());
  csTypeGamma.clear();

  m_nPhotonTerms = 0;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const double prefactor = dens * SpeedOfLight * m_fraction[i];
    // Check if optical data for this gas is available.
    std::string gasname = m_gas[i];
    if (gasname == "iC4H10") {
      gasname = "nC4H10";
      if (m_debug || verbose) {
        std::cout << m_className << "::ComputePhotonCollisionTable:\n"
                  << "    Photoabsorption cross-section for "
                  << "iC4H10 not available.\n"
                  << "    Using n-butane cross-section instead.\n";
      }
    }
    if (!data.IsAvailable(gasname)) return false;
    csTypeGamma.push_back(i * nCsTypesGamma + PhotonCollisionTypeIonisation);
    csTypeGamma.push_back(i * nCsTypesGamma + PhotonCollisionTypeInelastic);
    m_nPhotonTerms += 2;
    for (int j = 0; j < nEnergyStepsGamma; ++j) {
      // Retrieve total photoabsorption cross-section and ionisation yield.
      data.GetPhotoabsorptionCrossSection(gasname, (j + 0.5) * m_eStepGamma, cs,
                                          eta);
      m_cfTotGamma[j] += cs * prefactor;
      // Ionisation
      m_cfGamma[j].push_back(cs * prefactor * eta);
      // Inelastic absorption
      m_cfGamma[j].push_back(cs * prefactor * (1. - eta));
    }
  }

  // If requested, write the cross-sections to file.
  if (m_useCsOutput) {
    std::ofstream csfile;
    csfile.open("csgamma.txt", std::ios::out);
    for (int j = 0; j < nEnergyStepsGamma; ++j) {
      csfile << (j + 0.5) * m_eStepGamma << "  ";
      for (unsigned int i = 0; i < m_nPhotonTerms; ++i)
        csfile << m_cfGamma[j][i] << "  ";
      csfile << "\n";
    }
    csfile.close();
  }

  // Calculate the cumulative rates.
  for (int j = 0; j < nEnergyStepsGamma; ++j) {
    for (unsigned int i = 0; i < m_nPhotonTerms; ++i) {
      if (i > 0) m_cfGamma[j][i] += m_cfGamma[j][i - 1];
    }
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::ComputePhotonCollisionTable:\n";
    std::cout << "    Energy [eV]      Mean free path [um]\n";
    for (int i = 0; i < 10; ++i) {
      const int j = (2 * i + 1) * nEnergyStepsGamma / 20;
      const double en = (2 * i + 1) * m_eFinalGamma / 20;
      const double imfp = m_cfTotGamma[j] / SpeedOfLight;
      if (imfp > 0.) {
        printf("    %10.2f    %18.4f\n", en, 1.e4 / imfp);
      } else {
        printf("    %10.2f          ------------\n", en);
      }
    }
  }

  if (!m_useDeexcitation) return true;

  // Conversion factor from oscillator strength to cross-section
  constexpr double f2cs =
      FineStructureConstant * 2 * Pi2 * HbarC * HbarC / ElectronMass;
  // Discrete absorption lines
  int nResonanceLines = 0;
  for (auto& dxc : m_deexcitations) {
    if (dxc.osc < Small) continue;
    const double prefactor = dens * SpeedOfLight * m_fraction[dxc.gas];
    dxc.cf = prefactor * f2cs * dxc.osc;
    // Compute the line width due to Doppler broadening.
    const double mgas = ElectronMass / (m_rgas[dxc.gas] - 1.);
    const double wDoppler = sqrt(BoltzmannConstant * m_temperature / mgas);
    dxc.sDoppler = wDoppler * dxc.energy;
    // Compute the half width at half maximum due to resonance broadening.
    //   A. W. Ali and H. R. Griem, Phys. Rev. 140, 1044
    //   A. W. Ali and H. R. Griem, Phys. Rev. 144, 366
    const double kResBroad = 1.92 * Pi * sqrt(1. / 3.);
    dxc.gPressure = kResBroad * FineStructureConstant * pow(HbarC, 3) *
                    dxc.osc * dens * m_fraction[dxc.gas] /
                    (ElectronMass * dxc.energy);
    // Make an estimate for the width within which a photon can be
    // absorbed by the line
    constexpr double nWidths = 1000.;
    // Calculate the FWHM of the Voigt distribution according to the
    // approximation formula given in
    // Olivero and Longbothum, J. Quant. Spectr. Rad. Trans. 17, 233-236
    const double fwhmGauss = dxc.sDoppler * sqrt(2. * log(2.));
    const double fwhmLorentz = dxc.gPressure;
    const double fwhmVoigt =
        0.5 * (1.0692 * fwhmLorentz + sqrt(0.86639 * fwhmLorentz * fwhmLorentz +
                                           4 * fwhmGauss * fwhmGauss));
    dxc.width = nWidths * fwhmVoigt;
    ++nResonanceLines;
  }

  if (nResonanceLines <= 0) {
    std::cerr << m_className << "::ComputePhotonCollisionTable:\n"
              << "    No resonance lines found.\n";
    return true;
  }

  if (!(m_debug || verbose)) return true;
  std::cout << m_className << "::ComputePhotonCollisionTable:\n    "
            << "Discrete absorption lines:\n   Energy [eV]   "
            << "Line width (FWHM) [eV]    Mean free path [um]\n        "
            << "              Doppler    Pressure         (peak)\n";
  for (const auto& dxc : m_deexcitations) {
    if (dxc.osc < Small) continue;
    const double wp = 2 * dxc.gPressure;
    const double wd = 2 * sqrt(2 * log(2.)) * dxc.sDoppler;
    const double imfpP =
        (dxc.cf / SpeedOfLight) * TMath::Voigt(0., dxc.sDoppler, wp);
    if (imfpP > 0.) {
      printf("  %6.3f +/- %6.1e  %6.2e  %6.3e  %14.4f\n", dxc.energy,
             dxc.width, wd, wp, 1.e4 / imfpP);
    } else {
      printf("  %6.3f +/- %6.1e  %6.2e  %6.3e  -------------\n", dxc.energy,
             dxc.width, wd, wp);
    }
  }
  return true;
}

void MediumMagboltz::RunMagboltz(
    const double e, const double bmag, const double btheta, const int ncoll,
    bool verbose, double& vx, double& vy, double& vz, double& dl, double& dt,
    double& alpha, double& eta, double& lor, double& vxerr, double& vyerr,
    double& vzerr, double& dlerr, double& dterr, double& alphaerr,
    double& etaerr, double& lorerr, double& alphatof,
    std::array<double, 6>& difftens) {
  // Initialize the values.
  vx = vy = vz = 0.;
  dl = dt = 0.;
  alpha = eta = alphatof = 0.;
  lor = 0.;
  vxerr = vyerr = vzerr = 0.;
  dlerr = dterr = 0.;
  alphaerr = etaerr = 0.;
  lorerr = 0.;

  // Set input parameters in Magboltz common blocks.
  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = 4000;
  Magboltz::inpt_.nAniso = 2;

  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;
  Magboltz::inpt_.ipen = 0;
  Magboltz::setp_.nmax = ncoll;

  Magboltz::thrm_.ithrm = m_useGasMotion ? 1 : 0;

  Magboltz::setp_.efield = e;
  // Convert from Tesla to kGauss.
  Magboltz::bfld_.bmag = bmag * 10.;
  // Convert from radians to degree.
  Magboltz::bfld_.btheta = btheta * RadToDegree;

  // Set the gas composition in Magboltz.
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const int ng = GetGasNumberMagboltz(m_gas[i]);
    if (ng <= 0) {
      std::cerr << m_className << "::RunMagboltz:\n    Gas " << m_gas[i]
                << " does not have a gas number in Magboltz.\n";
      return;
    }
    Magboltz::gasn_.ngasn[i] = ng;
    Magboltz::ratio_.frac[i] = 100 * m_fraction[i];
  }

  // Run Magboltz.
  Magboltz::magboltz_();

  // Velocities. Convert to cm / ns.
  vx = Magboltz::vel_.wx * 1.e-9;
  vxerr = Magboltz::velerr_.dwx;
  vy = Magboltz::vel_.wy * 1.e-9;
  vyerr = Magboltz::velerr_.dwy;
  vz = Magboltz::vel_.wz * 1.e-9;
  vzerr = Magboltz::velerr_.dwz;

  // Calculate the Lorentz angle.
  const double vt = sqrt(vx * vx + vy * vy);
  const double v2 = (vx * vx + vy * vy + vz * vz);
  lor = atan2(vt, vz);
  if (vt > 0. && v2 > 0. && fabs(lor) > 0.) {
    const double dvx = vx * vxerr;
    const double dvy = vy * vyerr;
    const double dvz = vz * vzerr;
    const double a = vz / vt;
    lorerr = sqrt(a * a * (vx * vx * dvx * dvx + vy * vy * dvy * dvy) +
                  vt * vt * dvz * dvz) /
             v2;
    lorerr /= lor;
  }

  // Diffusion coefficients.
  dt = sqrt(0.2 * 0.5 * (Magboltz::diflab_.difxx + Magboltz::diflab_.difyy) /
            vz) *
       1.e-4;
  dterr = Magboltz::diferl_.dfter;
  dl = sqrt(0.2 * Magboltz::diflab_.difzz / vz) * 1.e-4;
  dlerr = Magboltz::diferl_.dfler;
  // Diffusion tensor.
  difftens[0] = 0.2e-4 * Magboltz::diflab_.difzz / vz;
  difftens[1] = 0.2e-4 * Magboltz::diflab_.difxx / vz;
  difftens[2] = 0.2e-4 * Magboltz::diflab_.difyy / vz;
  difftens[3] = 0.2e-4 * Magboltz::diflab_.difxz / vz;
  difftens[4] = 0.2e-4 * Magboltz::diflab_.difyz / vz;
  difftens[5] = 0.2e-4 * Magboltz::diflab_.difxy / vz;
  // Townsend and attachment coefficients.
  alpha = Magboltz::ctowns_.alpha;
  alphaerr = Magboltz::ctwner_.alper;
  eta = Magboltz::ctowns_.att;
  etaerr = Magboltz::ctwner_.atter;

  // Calculate effective Townsend SST coefficient from TOF results. 
  if (fabs(Magboltz::tofout_.tofdl) > 0.) {
    const double wrzn = 1.e5 * Magboltz::tofout_.tofwr;
    const double fc1 = 0.5 * wrzn / Magboltz::tofout_.tofdl;
    const double fc2 = (Magboltz::tofout_.ralpha - 
                        Magboltz::tofout_.rattof) * 1.e12 / 
                       Magboltz::tofout_.tofdl;
    alphatof = fc1 - sqrt(fc1 * fc1 - fc2);
  }
  // Print the results.
  if (!(m_debug || verbose)) return;
  std::cout << m_className << "::RunMagboltz: Results:\n";
  printf("    Drift velocity along E:   %12.8f cm/ns +/- %5.2f%%\n", vz, vzerr);
  printf("    Drift velocity along Bt:  %12.8f cm/ns +/- %5.2f%%\n", vx, vxerr);
  printf("    Drift velocity along ExB: %12.8f cm/ns +/- %5.2f%%\n", vy, vyerr);
  printf("    Lorentz angle:            %12.3f degree\n", lor * RadToDegree);
  printf("    Longitudinal diffusion:   %12.8f cm1/2 +/- %5.2f%%\n", dl, dlerr);
  printf("    Transverse diffusion:     %12.8f cm1/2 +/- %5.2f%%\n", dt, dterr);
  printf("    Townsend coefficient:     %12.4f cm-1  +/- %5.2f%%\n", alpha,
         alphaerr);
  printf("    Attachment coefficient:   %12.4f cm-1  +/- %5.2f%%\n", eta,
         etaerr);
  if (alphatof > 0.) {
    printf("    TOF effective Townsend:   %12.4f cm-1 (alpha - eta)\n",
           alphatof);
  }
}

void MediumMagboltz::GenerateGasTable(const int numColl, const bool verbose) {
  // Set the reference pressure and temperature.
  m_pressureTable = m_pressure;
  m_temperatureTable = m_temperature;

  // Initialize the parameter arrays.
  const unsigned int nEfields = m_eFields.size();
  const unsigned int nBfields = m_bFields.size();
  const unsigned int nAngles = m_bAngles.size();
  Init(nEfields, nBfields, nAngles, m_eVelE, 0.);
  Init(nEfields, nBfields, nAngles, m_eVelB, 0.);
  Init(nEfields, nBfields, nAngles, m_eVelX, 0.);
  Init(nEfields, nBfields, nAngles, m_eDifL, 0.);
  Init(nEfields, nBfields, nAngles, m_eDifT, 0.);
  Init(nEfields, nBfields, nAngles, m_eLor, 0.);
  Init(nEfields, nBfields, nAngles, m_eAlp, -30.);
  Init(nEfields, nBfields, nAngles, m_eAlp0, -30.);
  Init(nEfields, nBfields, nAngles, m_eAtt, -30.);
  Init(nEfields, nBfields, nAngles, 6, m_eDifM, 0.);

  m_excRates.clear();
  m_ionRates.clear();
  m_excLevels.clear();
  m_ionLevels.clear();
  std::vector<unsigned int> excLevelIndex;
  std::vector<unsigned int> ionLevelIndex;

  double vx = 0., vy = 0., vz = 0.;
  double difl = 0., dift = 0.;
  double alpha = 0., eta = 0.;
  double lor = 0.;
  double vxerr = 0., vyerr = 0., vzerr = 0.;
  double diflerr = 0., difterr = 0.;
  double alphaerr = 0., etaerr = 0.;
  double alphatof = 0.;
  double lorerr = 0.;
  std::array<double, 6> difftens;

  // Run through the grid of E- and B-fields and angles.
  for (unsigned int i = 0; i < nEfields; ++i) {
    const double e = m_eFields[i];
    for (unsigned int j = 0; j < nAngles; ++j) {
      const double a = m_bAngles[j];
      for (unsigned int k = 0; k < nBfields; ++k) {
        const double b = m_bFields[k];
        std::cout << m_className << "::GenerateGasTable: E = " << e
                  << " V/cm, B = " << b << " T, angle: " << a << " rad\n";
        RunMagboltz(e, b, a, numColl, verbose, vx, vy, vz, difl, dift, alpha,
                    eta, lor, vxerr, vyerr, vzerr, diflerr, difterr, alphaerr,
                    etaerr, lorerr, alphatof, difftens);
        m_eVelE[j][k][i] = vz;
        m_eVelX[j][k][i] = vy;
        m_eVelB[j][k][i] = vx;
        m_eDifL[j][k][i] = difl;
        m_eDifT[j][k][i] = dift;
        m_eLor[j][k][i] = lor;
        m_eAlp[j][k][i] = alpha > 0. ? log(alpha) : -30.;
        m_eAlp0[j][k][i] = alpha > 0. ? log(alpha) : -30.;
        m_eAtt[j][k][i] = eta > 0. ? log(eta) : -30.;
        for (unsigned int l = 0; l < 6; ++l) {
          m_eDifM[l][j][k][i] = difftens[l];
        }
        // If not done yet, retrieve the excitation and ionisation levels.
        if (m_excLevels.empty() && m_ionLevels.empty()) {
          for (long long il = 0; il < Magboltz::nMaxLevels; ++il) {
            if (Magboltz::large_.iarry[il] <= 0) break;
            // Skip levels that are not ionisations or inelastic collisions.
            const int cstype = (Magboltz::large_.iarry[il] - 1) % 5;
            if (cstype != 1 && cstype != 3) continue;
            const int igas = int((Magboltz::large_.iarry[il] - 1) / 5);
            std::string descr = GetDescription(il, Magboltz::scrip_.dscrpt);
            if (cstype == 3) {
              // Skip levels that are not excitations.
              if (!(descr[1] == 'E' && descr[2] == 'X') &&
                  !(descr[0] == 'E' && descr[1] == 'X'))
                continue;
            }
            descr = m_gas[igas] + descr;
            if (cstype == 3) {
              ExcLevel exc;
              exc.label = descr;
              exc.energy = Magboltz::large_.ein[il];
              exc.prob = 0.;
              exc.rms = 0.;
              exc.dt = 0.;
              m_excLevels.push_back(std::move(exc));
              excLevelIndex.push_back(il);
            } else {
              IonLevel ion;
              ion.label = descr;
              ion.energy = Magboltz::large_.ein[il];
              m_ionLevels.push_back(std::move(ion));
              ionLevelIndex.push_back(il);
            }
          }
          std::cout << m_className << "::GenerateGasTable: Found "
                    << m_excLevels.size() << " excitations and "
                    << m_ionLevels.size() << " ionisations.\n";
          for (const auto& exc : m_excLevels) {
            std::cout << "    " << exc.label << ", energy = " << exc.energy
                      << " eV.\n";
          }
          for (const auto& ion : m_ionLevels) {
            std::cout << "    " << ion.label << ", energy = " << ion.energy
                      << " eV.\n";
          }
          Init(nEfields, nBfields, nAngles, m_excLevels.size(),
               m_excRates, 0.);
          Init(nEfields, nBfields, nAngles, m_ionLevels.size(),
               m_ionRates, 0.);
        }
        // Retrieve the excitation and ionisation rates.
        const unsigned int nExc = m_excLevels.size();
        for (unsigned int ie = 0; ie < nExc; ++ie) {
          const unsigned int level = excLevelIndex[ie];
          m_excRates[ie][j][k][i] = Magboltz::outpt_.icoln[level];
        }
        const unsigned int nIon = m_ionLevels.size();
        for (unsigned int ii = 0; ii < nIon; ++ii) {
          const unsigned int level = ionLevelIndex[ii];
          m_ionRates[ii][j][k][i] = Magboltz::outpt_.icoln[level];
        }
      }
    }
  }
  // Set the threshold indices.
  SetThreshold(m_eAlp);
  SetThreshold(m_eAtt);
}
}
