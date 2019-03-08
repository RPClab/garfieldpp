#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <ctime>

#include "Utilities.hh"
#include "MediumGas.hh"
#include "OpticalData.hh"
#include "FundamentalConstants.hh"
#include "GarfieldConstants.hh"

namespace {

std::string FormatFloat(const double x, 
    const unsigned int width = 15, const unsigned int precision = 8) {

  char buffer[256];
  std::snprintf(buffer, width + 1, "%*.*E", width, precision, x);
  return std::string(buffer);
}

std::string FormatInt(const int n, const unsigned int width) {

  char buffer[256];
  std::snprintf(buffer, width + 1, "%*d", width, n); 
  return std::string(buffer);
}

void PrintArray(const std::vector<double>& values, std::ofstream& outfile, 
    int& col, const int ncols) {

  for (const auto value : values) {
    outfile << FormatFloat(value);
    ++col;
    if (col % ncols == 0) outfile << "\n";
  }
}

void PrintExtrapolation(const std::pair<unsigned int, unsigned int>& extr) {

  std::cout << "        Low field extrapolation: ";
  if (extr.first == 0)
    std::cout << " constant\n";
  else if (extr.first == 1)
    std::cout << " linear\n";
  else if (extr.first == 2)
    std::cout << " exponential\n";
  else
    std::cout << " unknown\n";
  std::cout << "        High field extrapolation: ";
  if (extr.second == 0)
    std::cout << " constant\n";
  else if (extr.second == 1)
    std::cout << " linear\n";
  else if (extr.second == 2)
    std::cout << " exponential\n";
  else
    std::cout << " unknown\n";
}

}

namespace Garfield {

MediumGas::MediumGas() : Medium(),
      m_pressureTable(m_pressure),
      m_temperatureTable(m_temperature) {

  m_className = "MediumGas";

  m_gas.fill("");
  m_fraction.fill(0.);
  m_atWeight.fill(0.);
  m_atNum.fill(0.);
  // Default gas mixture: pure argon
  m_gas[0] = "Ar";
  m_fraction[0] = 1.;
  m_name = m_gas[0];
  GetGasInfo(m_gas[0], m_atWeight[0], m_atNum[0]);

  m_rPenningGas.fill(0.);
  m_lambdaPenningGas.fill(0.);
 
  m_isChanged = true;

  EnableDrift();
  EnablePrimaryIonisation();

}

bool MediumGas::SetComposition(const std::string& gas1, const double f1,
                               const std::string& gas2, const double f2,
                               const std::string& gas3, const double f3,
                               const std::string& gas4, const double f4,
                               const std::string& gas5, const double f5,
                               const std::string& gas6, const double f6) {

  std::array<std::string, 6> gases = {gas1, gas2, gas3, gas4, gas5, gas6};
  std::array<double, 6> fractions = {f1, f2, f3, f4, f5, f6};

  // Make a backup copy of the gas composition.
  const std::array<std::string, m_nMaxGases> gasOld = m_gas;
  const unsigned int nComponentsOld = m_nComponents;

  m_nComponents = 0;
  m_gas.fill("");
  m_fraction.fill(0.);
  m_atWeight.fill(0.);
  m_atNum.fill(0.);
  for (unsigned int i = 0; i < 6; ++i) {
    // Find the gas name corresponding to the input string.
    std::string gasname = "";
    if (fractions[i] > 0. && GetGasName(gases[i], gasname)) {
      m_gas[m_nComponents] = gasname;
      m_fraction[m_nComponents] = fractions[i];
      ++m_nComponents;
    }
  }

  // Check if at least one valid ingredient was specified.
  if (m_nComponents == 0) {
    std::cerr << m_className << "::SetComposition:\n"
              << "    Error setting the composition. No valid components.\n";
    return false;
  }

  // Establish the name.
  m_name = "";
  double sum = 0.;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (i > 0) m_name += "/";
    m_name += m_gas[i];
    sum += m_fraction[i];
  }
  // Normalise the fractions to one.
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    m_fraction[i] /= sum;
  }

  // Set the atomic weight and number.
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    GetGasInfo(m_gas[i], m_atWeight[i], m_atNum[i]);
  }

  // Print the composition.
  std::cout << m_className << "::SetComposition:\n    " << m_name;
  if (m_nComponents > 1) {
    std::cout << " (" << m_fraction[0] * 100;
    for (unsigned int i = 1; i < m_nComponents; ++i) {
      std::cout << "/" << m_fraction[i] * 100;
    }
    std::cout << ")";
  }
  std::cout << "\n";

  // Force a recalculation of the collision rates.
  m_isChanged = true;

  // Copy the previous Penning transfer parameters.
  std::array<double, m_nMaxGases> rPenningGasOld;
  std::array<double, m_nMaxGases> lambdaPenningGasOld;
  rPenningGasOld.swap(m_rPenningGas);
  lambdaPenningGasOld.swap(m_lambdaPenningGas);
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    for (unsigned int j = 0; j < nComponentsOld; ++j) {
      if (m_gas[i] != gasOld[j]) continue;
      if (rPenningGasOld[j] < Small) continue;
      m_rPenningGas[i] = rPenningGasOld[j];
      m_lambdaPenningGas[i] = lambdaPenningGasOld[i];
      std::cout << m_className << "::SetComposition:\n"
                << "    Using Penning transfer parameters for " 
                << m_gas[i] << " from previous mixture.\n"
                << "      r      = " << m_rPenningGas[i] << "\n"
                << "      lambda = " << m_lambdaPenningGas[i] << " cm\n";
    }
  }
  return true;
}

void MediumGas::GetComposition(std::string& gas1, double& f1, std::string& gas2,
                               double& f2, std::string& gas3, double& f3,
                               std::string& gas4, double& f4, std::string& gas5,
                               double& f5, std::string& gas6, double& f6) {

  gas1 = m_gas[0];
  gas2 = m_gas[1];
  gas3 = m_gas[2];
  gas4 = m_gas[3];
  gas5 = m_gas[4];
  gas6 = m_gas[5];
  f1 = m_fraction[0];
  f2 = m_fraction[1];
  f3 = m_fraction[2];
  f4 = m_fraction[3];
  f5 = m_fraction[4];
  f6 = m_fraction[5];
}

void MediumGas::GetComponent(const unsigned int i, 
                             std::string& label, double& f) {

  if (i >= m_nComponents) {
    std::cerr << m_className << "::GetComponent: Index out of range.\n";
    label = "";
    f = 0.;
    return;
  }

  label = m_gas[i];
  f = m_fraction[i];
}

void MediumGas::SetAtomicNumber(const double z) {

  std::cerr << m_className << "::SetAtomicNumber:\n"
            << "    Effective Z cannot be changed directly to " << z << ".\n"
            << "    Use SetComposition to define the gas mixture.\n";
}

void MediumGas::SetAtomicWeight(const double a) {

  std::cerr << m_className << "::SetAtomicWeight:\n"
            << "    Effective A cannot be changed directly to " << a << ".\n"
            << "    Use SetComposition to define the gas mixture.\n";
}

void MediumGas::SetNumberDensity(const double n) {

  std::cerr << m_className << "::SetNumberDensity:\n"
            << "    Density cannot directly be changed to " << n << ".\n"
            << "    Use SetTemperature and SetPressure.\n";
}

void MediumGas::SetMassDensity(const double rho) {

  std::cerr << m_className << "::SetMassDensity:\n"
            << "    Density cannot directly be changed to " << rho << ".\n"
            << "    Use SetTemperature, SetPressure and SetComposition.\n";
}

double MediumGas::GetAtomicWeight() const {

  // Effective A, weighted by the fractions of the components.
  double a = 0.;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    a += m_atWeight[i] * m_fraction[i];
  }
  return a;
}

double MediumGas::GetNumberDensity() const {

  // Ideal gas law.
  return LoschmidtNumber * (m_pressure / AtmosphericPressure) *
         (ZeroCelsius / m_temperature);
}

double MediumGas::GetMassDensity() const {

  return GetNumberDensity() * GetAtomicWeight() * AtomicMassUnit;
}

double MediumGas::GetAtomicNumber() const {

  // Effective Z, weighted by the fractions of the components.
  double z = 0.;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    z += m_atNum[i] * m_fraction[i];
  }
  return z;
}

bool MediumGas::LoadGasFile(const std::string& filename) {

  std::ifstream gasfile;
  // Open the file.
  gasfile.open(filename.c_str());
  // Make sure the file could be opened.
  if (!gasfile.is_open()) {
    std::cerr << m_className << "::LoadGasFile:\n"
              << "    Cannot open file " << filename << ".\n";
    return false;
  }
  std::cout << m_className << "::LoadGasFile: Reading " << filename << ".\n";

  // GASOK bits
  std::string gasBits = "";

  // Gas composition
  constexpr int nMagboltzGases = 60;
  std::array<double, nMagboltzGases> mixture;
  mixture.fill(0.);

  int nE = 1;
  int nB = 1;
  int nA = 1;

  int version = 12;

  m_excLevels.clear();
  m_ionLevels.clear();

  // Start reading the data.
  if (m_debug) std::cout << m_className << "::LoadGasFile: Header...\n";
  bool atTables = false;
  while (!atTables) {
    char line[256];
    gasfile.getline(line, 256);
    if (strncmp(line, " The gas tables follow:", 8) == 0 ||
        strncmp(line, "The gas tables follow:", 7) == 0) {
      atTables = true;
      break;
    }
    char* token = strtok(line, " :,%");
    while (token) {
      if (strcmp(token, "Version") == 0) {
        token = strtok(NULL, " :,%");
        version = atoi(token);
        // Check the version number.
        if (version != 10 && version != 11 && version != 12) {
          std::cerr << m_className << "::LoadGasFile:\n"
                    << "    The file has version number " << version << ".\n"
                    << "    Files written in this format cannot be read.\n";
          gasfile.close();
          return false;
        } else {
          std::cout << m_className << "::LoadGasFile:\n"
                    << "    Version: " << version << "\n";
        }
      } else if (strcmp(token, "GASOK") == 0) {
        // Get the GASOK bits indicating if a parameter
        // is present in the table (T) or not (F).
        token = strtok(NULL, " :,%\t");
        token = strtok(NULL, " :,%\t");
        gasBits += token;
        if (m_debug) std::cout << "    GASOK bits: " << gasBits << "\n";
      } else if (strcmp(token, "Identifier") == 0) {
        // Get the identification string.
        std::string identifier = "";
        token = strtok(NULL, "\n");
        if (token != NULL) identifier += token;
        if (m_debug) std::cout << "    Identifier: " << identifier << "\n";
      } else if (strcmp(token, "Dimension") == 0) {
        token = strtok(NULL, " :,%\t");
        if (strcmp(token, "F") == 0) {
          m_map2d = false;
        } else {
          m_map2d = true;
        }
        token = strtok(NULL, " :,%\t");
        nE = atoi(token);
        // Check the number of E points.
        if (nE <= 0) {
          std::cerr << m_className << "::LoadGasFile:\n"
                    << "    Number of E fields out of range.\n";
          gasfile.close();
          return false;
        }
        token = strtok(NULL, " :,%\t");
        nA = atoi(token);
        // Check the number of angles.
        if (m_map2d && nA <= 0) {
          std::cerr << m_className << "::LoadGasFile:\n"
                    << "    Number of E-B angles out of range.\n";
          gasfile.close();
          return false;
        }

        token = strtok(NULL, " :,%\t");
        nB = atoi(token);
        // Check the number of B points.
        if (m_map2d && nB <= 0) {
          std::cerr << m_className << "::LoadGasFile:\n"
                    << "    Number of B fields out of range.\n";
          gasfile.close();
          return false;
        }

        m_eFields.resize(nE);
        m_bFields.resize(nB);
        m_bAngles.resize(nA);

        // Fill in the excitation/ionisation structs
        // Excitation
        token = strtok(NULL, " :,%\t");
        const int nexc = atoi(token);
        // Ionization
        token = strtok(NULL, " :,%\t");
        const int nion = atoi(token);
        if (m_debug) {
          std::cout << "    " << nexc << " excitations, " 
                    << nion << " ionisations.\n";
        }
      } else if (strcmp(token, "E") == 0) {
        token = strtok(NULL, " :,%");
        if (strcmp(token, "fields") == 0) {
          for (int i = 0; i < nE; ++i) gasfile >> m_eFields[i];
        }
      } else if (strcmp(token, "E-B") == 0) {
        token = strtok(NULL, " :,%");
        if (strcmp(token, "angles") == 0) {
          for (int i = 0; i < nA; ++i) gasfile >> m_bAngles[i];
        }
      } else if (strcmp(token, "B") == 0) {
        token = strtok(NULL, " :,%");
        if (strcmp(token, "fields") == 0) {
          double bstore = 0.;
          for (int i = 0; i < nB; i++) {
            gasfile >> bstore;
            m_bFields[i] = bstore / 100.;
          }
        }
      } else if (strcmp(token, "Mixture") == 0) {
        for (int i = 0; i < nMagboltzGases; ++i) {
          gasfile >> mixture[i];
        }
      } else if (strcmp(token, "Excitation") == 0) {
        // Skip number.
        token = strtok(NULL, " :,%");
        // Get label.
        token = strtok(NULL, " :,%");
        ExcLevel exc;
        exc.label = token;
        // Get energy.
        token = strtok(NULL, " :,%");
        exc.energy = atof(token);
        // Get Penning probability.
        token = strtok(NULL, " :,%");
        exc.prob = atof(token);
        exc.rms = 0.;
        exc.dt = 0.;
        if (version >= 11) {
          // Get Penning rms distance.
          token = strtok(NULL, " :,%");
          if (token) {
            exc.rms = atof(token);
            // Get decay time.
            token = strtok(NULL, " :,%");
            if (token) exc.dt = atof(token);
          }
        }
        m_excLevels.push_back(std::move(exc));
      } else if (strcmp(token, "Ionisation") == 0) {
        // Skip number.
        token = strtok(NULL, " :,%");
        // Get label.
        token = strtok(NULL, " :,%");
        IonLevel ion;
        ion.label += token;
        // Get energy.
        token = strtok(NULL, " :,%");
        ion.energy = atof(token);
        m_ionLevels.push_back(std::move(ion));
      }
      token = strtok(NULL, " :,%");
    }
  }

  // Decode the GASOK bits.
  // GASOK(I)   : .TRUE. if present
  // (1)  electron drift velocity || E
  // (2)  ion mobility,
  // (3)  longitudinal diffusion || E
  // (4)  Townsend coefficient,
  // (5)  cluster size distribution.
  // (6)  attachment coefficient,
  // (7)  Lorentz angle,
  // (8)  transverse diffusion || ExB and Bt
  // (9)  electron drift velocity || Bt
  // (10) electron drift velocity || ExB
  // (11) diffusion tensor
  // (12) ion dissociation
  // (13) allocated for SRIM data (not used)
  // (14) allocated for HEED data (not used)
  // (15) excitation rates
  // (16) ionisation rates

  if (m_debug) {
    std::cout << m_className << "::LoadGasFile:\n    "
              << nE << " electric fields, " 
              << nB << " magnetic fields, " 
              << nA << " angles.\n    ";
  }
  m_eVelocityE.clear();
  if (gasBits[0] == 'T') InitTable(nE, nB, nA, m_eVelocityE, 0.);
  m_ionMobility.clear();
  if (gasBits[1] == 'T') InitTable(nE, nB, nA, m_ionMobility, 0.);
  m_eDiffLong.clear();
  if (gasBits[2] == 'T') InitTable(nE, nB, nA, m_eDiffLong, 0.);
  m_eTownsend.clear();
  m_eTownsendNoPenning.clear();
  if (gasBits[3] == 'T') {
    InitTable(nE, nB, nA, m_eTownsend, -30.);
    InitTable(nE, nB, nA, m_eTownsendNoPenning, -30.);
  }
  // gasBits[4]: cluster size distribution; skipped
  m_eAttachment.clear();
  if (gasBits[5] == 'T') InitTable(nE, nB, nA, m_eAttachment, -30.);
  m_eLorentzAngle.clear();
  if (gasBits[6] == 'T') InitTable(nE, nB, nA, m_eLorentzAngle, -30.);
  m_eDiffTrans.clear();
  if (gasBits[7] == 'T') InitTable(nE, nB, nA, m_eDiffTrans, 0.);
  m_eVelocityB.clear();
  if (gasBits[8] == 'T') InitTable(nE, nB, nA, m_eVelocityB, 0.);
  m_eVelocityExB.clear();
  if (gasBits[9] == 'T') InitTable(nE, nB, nA, m_eVelocityExB, 0.);
  m_eDiffTens.clear();
  if (gasBits[10] == 'T') InitTensor(nE, nB, nA, 6, m_eDiffTens, 0.);
  m_ionDissociation.clear();
  if (gasBits[11] == 'T') InitTable(nE, nB, nA, m_ionDissociation, -30.);
  // gasBits[12]: SRIM; skipped
  // gasBits[13]: HEED; skipped
  m_excRates.clear();
  if (gasBits[14] == 'T') {
    InitTensor(nE, nB, nA, m_excLevels.size(), m_excRates, 0.);
  }
  m_ionRates.clear();
  if (gasBits[15] == 'T') {
    InitTensor(nE, nB, nA, m_ionLevels.size(), m_ionRates, 0.);
  }

  // Check the gas mixture.
  std::vector<std::string> gasnames;
  std::vector<double> percentages;
  bool gasMixOk = true;
  unsigned int gasCount = 0;
  for (int i = 0; i < nMagboltzGases; ++i) {
    if (mixture[i] < Small) continue;
    std::string gasname = "";
    if (!GetGasName(i + 1, version, gasname)) {
      std::cerr << m_className << "::LoadGasFile:\n"
                << "    Unknown gas (gas number " << i + 1 << ").\n";
      gasMixOk = false;
      break;
    }
    gasnames.push_back(gasname);
    percentages.push_back(mixture[i]);
    ++gasCount;
  }
  if (gasCount > m_nMaxGases) {
    std::cerr << m_className << "::LoadGasFile:\n"
              << "    Gas mixture has " << gasCount << " components.\n"
              << "    Number of gases is limited to " << m_nMaxGases << ".\n";
    gasMixOk = false;
  } else if (gasCount == 0) {
    std::cerr << m_className << "::LoadGasFile:\n"
              << "    Gas mixture is not defined (zero components).\n";
    gasMixOk = false;
  }
  double sum = 0.;
  for (unsigned int i = 0; i < gasCount; ++i) sum += percentages[i];
  if (gasMixOk && sum != 100.) {
    std::cout << m_className << "::LoadGasFile:\n"
              << "    Renormalizing the percentages.\n";
    for (unsigned int i = 0; i < gasCount; ++i) percentages[i] *= 100. / sum;
  }

  // Force re-initialisation of collision rates etc.
  m_isChanged = true;

  if (gasMixOk) {
    m_name = "";
    m_nComponents = gasCount;
    for (unsigned int i = 0; i < m_nComponents; ++i) {
      if (i > 0) m_name += "/";
      m_name += gasnames[i];
      m_gas[i] = gasnames[i];
      m_fraction[i] = percentages[i] / 100.;
      GetGasInfo(m_gas[i], m_atWeight[i], m_atNum[i]);
    }
    std::cout << m_className << "::LoadGasFile:\n"
              << "    Gas composition set to " << m_name;
    if (m_nComponents > 1) {
      std::cout << " (" << m_fraction[0] * 100;
      for (unsigned int i = 1; i < m_nComponents; ++i) {
        std::cout << "/" << m_fraction[i] * 100;
      }
      std::cout << ")";
    }
    std::cout << "\n";
  } else {
    std::cerr << m_className << "::LoadGasFile:\n"
              << "    Gas composition could not be established.\n";
  }

  if (m_debug) std::cout << m_className << "::LoadGasFile: Gas tables...\n";

  if (m_map2d) {
    if (m_debug) {
      std::cout << m_className << "::LoadGasFile: Loading 3D table.\n";
    }
    for (int i = 0; i < nE; i++) {
      for (int j = 0; j < nA; j++) {
        for (int k = 0; k < nB; k++) {
          // Drift velocity along E, Bt and ExB
          double ve = 0., vb = 0., vexb = 0.;
          gasfile >> ve >> vb >> vexb;
          // Convert from cm / us to cm / ns
          ve *= 1.e-3;
          vb *= 1.e-3;
          vexb *= 1.e-3;
          if (!m_eVelocityE.empty()) m_eVelocityE[j][k][i] = ve;
          if (!m_eVelocityB.empty()) m_eVelocityB[j][k][i] = vb;
          if (!m_eVelocityExB.empty()) m_eVelocityExB[j][k][i] = vexb;
          // Longitudinal and transverse diffusion coefficients
          double dl = 0., dt = 0.;
          gasfile >> dl >> dt;
          if (!m_eDiffLong.empty()) m_eDiffLong[j][k][i] = dl;
          if (!m_eDiffTrans.empty()) m_eDiffTrans[j][k][i] = dt;
          // Townsend and attachment coefficients
          double alpha = 0., alpha0 = 0., eta = 0.;
          gasfile >> alpha >> alpha0 >> eta;
          if (!m_eTownsend.empty()) {
            m_eTownsend[j][k][i] = alpha;
            m_eTownsendNoPenning[j][k][i] = alpha0;
          }
          if (!m_eAttachment.empty()) m_eAttachment[j][k][i] = eta;
          // Ion mobility
          double mu = 0.;
          gasfile >> mu;
          // Convert from cm2 / (V us) to cm2 / (V ns)
          mu *= 1.e-3;
          if (!m_ionMobility.empty()) m_ionMobility[j][k][i] = mu;
          // Lorentz angle
          double lor = 0.;
          gasfile >> lor;
          if (!m_eLorentzAngle.empty()) m_eLorentzAngle[j][k][i] = lor;
          // Ion dissociation
          double diss = 0.;
          gasfile >> diss;
          if (!m_ionDissociation.empty()) m_ionDissociation[j][k][i] = diss;
          // Diffusion tensor
          for (int l = 0; l < 6; l++) {
            double diff = 0.;
            gasfile >> diff;
            if (!m_eDiffTens.empty()) m_eDiffTens[l][j][k][i] = diff;
          }
          // Excitation rates
          const unsigned int nexc = m_excLevels.size();
          for (unsigned int l = 0; l < nexc; ++l) {
            double rate = 0.;
            gasfile >> rate;
            if (!m_excRates.empty()) m_excRates[l][j][k][i] = rate;
          }
          // Ionization rates
          const unsigned int nion = m_ionLevels.size();
          for (unsigned int l = 0; l < nion; ++l) {
            double rate = 0.;
            gasfile >> rate;
            if (!m_ionRates.empty()) m_ionRates[l][j][k][i] = rate;
          }
        }
      }
    }
  } else {
    if (m_debug) {
      std::cout << m_className << "::LoadGasFile: Reading 1D table.\n";
    }
    for (int i = 0; i < nE; i++) {
      double waste = 0.;
      // Drift velocity along E, Bt, ExB
      double ve = 0., vb = 0., vexb = 0.;
      gasfile >> ve >> waste >> vb >> waste >> vexb >> waste;
      if (!m_eVelocityE.empty()) m_eVelocityE[0][0][i] = 1.e-3 * ve;
      if (!m_eVelocityB.empty()) m_eVelocityB[0][0][i] = 1.e-3 * vb;
      if (!m_eVelocityExB.empty()) m_eVelocityExB[0][0][i] = 1.e-3 * vexb;
      // Longitudinal and transferse diffusion coefficients
      double dl = 0., dt = 0.;
      gasfile >> dl >> waste >> dt >> waste;
      if (!m_eDiffLong.empty()) m_eDiffLong[0][0][i] = dl;
      if (!m_eDiffTrans.empty()) m_eDiffTrans[0][0][i] = dt;
      // Townsend and attachment coefficients
      double alpha = 0., alpha0 = 0., eta = 0.;
      gasfile >> alpha >> waste >> alpha0 >> eta >> waste;
      if (!m_eTownsend.empty()) {
        m_eTownsend[0][0][i] = alpha;
        m_eTownsendNoPenning[0][0][i] = alpha0;
      }
      if (!m_eAttachment.empty()) {
        m_eAttachment[0][0][i] = eta;
      }
      // Ion mobility
      double mu = 0.;
      gasfile >> mu >> waste;
      if (!m_ionMobility.empty()) m_ionMobility[0][0][i] = 1.e-3 * mu;
      // Lorentz angle
      double lor = 0.;
      gasfile >> lor >> waste;
      if (!m_eLorentzAngle.empty()) m_eLorentzAngle[0][0][i] = lor;
      // Ion dissociation
      double diss = 0.;
      gasfile >> diss >> waste;
      if (!m_ionDissociation.empty()) m_ionDissociation[0][0][i] = diss;
      // Diffusion tensor
      for (int j = 0; j < 6; j++) {
        double diff = 0.;
        gasfile >> diff >> waste;
        if (!m_eDiffTens.empty()) m_eDiffTens[j][0][0][i] = diff;
      }
      // Excitation rates
      const unsigned int nexc = m_excLevels.size();
      for (unsigned int j = 0; j < nexc; ++j) {
        double rate = 0.;
        gasfile >> rate >> waste;
        if (!m_excRates.empty()) m_excRates[j][0][0][i] = rate;
      }
      // Ionization rates
      const unsigned int nion = m_ionLevels.size();
      for (unsigned int j = 0; j < nion; ++j) {
        double rate = 0.;
        gasfile >> rate >> waste;
        if (!m_ionRates.empty()) m_ionRates[j][0][0][i] = rate;
      }
    }
  }

  // Extrapolation methods
  std::array<unsigned int, 13> hExtrap = {{0}};
  std::array<unsigned int, 13> lExtrap = {{1}};
  // Interpolation methods
  std::array<unsigned int, 13> interpMeth = {{2}};
  // Ion diffusion coefficients.
  double ionDiffL = 0.;
  double ionDiffT = 0.;
  // Moving on to the file footer
  gasfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  if (m_debug) std::cout << m_className << "::LoadGasFile: Footer...\n";
  bool done = false;
  while (!done) {
    char line[256];
    gasfile.getline(line, 256);
    char* token = strtok(line, " :,%=\t");
    while (token) {
      std::cout << "TOKEN: " << token << "\n";
      if (strcmp(token, "H") == 0) {
        token = strtok(NULL, " :,%=\t");
        for (int i = 0; i < 13; i++) {
          token = strtok(NULL, " :,%=\t");
          if (token != NULL) hExtrap[i] = atoi(token);
        }
      } else if (strcmp(token, "L") == 0) {
        token = strtok(NULL, " :,%=\t");
        for (int i = 0; i < 13; i++) {
          token = strtok(NULL, " :,%=\t");
          if (token != NULL) lExtrap[i] = atoi(token);
        }
      } else if (strcmp(token, "Thresholds") == 0) {
        token = strtok(NULL, " :,%=\t");
        if (token != NULL) thrElectronTownsend = atoi(token);
        token = strtok(NULL, " :,%=\t");
        if (token != NULL) thrElectronAttachment = atoi(token);
        token = strtok(NULL, " :,%=\t");
        if (token != NULL) thrIonDissociation = atoi(token);
      } else if (strcmp(token, "Interp") == 0) {
        for (int i = 0; i < 13; i++) {
          token = strtok(NULL, " :,%=\t");
          if (token != NULL) interpMeth[i] = atoi(token);
        }
      } else if (strcmp(token, "A") == 0) {
        token = strtok(NULL, " :,%=\t");
        // Parameter for energy loss distribution, currently not used
        // double a;
        // if (token != NULL) a = atof(token);
      } else if (strcmp(token, "Z") == 0) {
        // Parameter for energy loss distribution, currently not used
        token = strtok(NULL, " :,%=\t");
        // double z;
        // if (token != NULL) z = atof(token);
      } else if (strcmp(token, "EMPROB") == 0) {
        // Parameter for energy loss distribution, currently not used
        token = strtok(NULL, " :,%=\t");
        // double emprob;
        // if (token != NULL) emprob = atof(token);
      } else if (strcmp(token, "EPAIR") == 0) {
        // Parameter for energy loss distribution, currently not used
        token = strtok(NULL, " :,%=\t");
        // double epair;
        // if (token != NULL) epair = atof(token);
      } else if (strcmp(token, "Ion") == 0) {
        // Ion diffusion coefficients
        token = strtok(NULL, " :,%=\t");
        token = strtok(NULL, " :,%=\t");
        if (token != NULL) ionDiffL = atof(token);
        token = strtok(NULL, " :,%=\t");
        if (token != NULL) ionDiffT = atof(token);
      } else if (strcmp(token, "CMEAN") == 0) {
        // Cluster parameter, currently not used
        token = strtok(NULL, " :,%=\t");
        // double clsPerCm;
        // if (token != NULL) clsPerCm = atof(token);
      } else if (strcmp(token, "RHO") == 0) {
        // Parameter for energy loss distribution, currently not used
        token = strtok(NULL, " :,%=\t");
        // double rho;
        // if (token != NULL) rho = atof(token);
      } else if (strcmp(token, "PGAS") == 0) {
        token = strtok(NULL, " :,%=\t");
        double pTorr = 760.;
        if (token != NULL) pTorr = atof(token);
        if (pTorr > 0.) m_pressure = pTorr;
      } else if (strcmp(token, "TGAS") == 0) {
        token = strtok(NULL, " :,%=\t");
        double tKelvin = 293.15;
        if (token != NULL) tKelvin = atof(token);
        if (tKelvin > 0.) m_temperature = tKelvin;
        done = true;
        break;
      } else {
        done = true;
        break;
      }
      token = strtok(NULL, " :,%=\t");
    }
  }

  gasfile.close();

  // Set the reference pressure and temperature.
  m_pressureTable = m_pressure;
  m_temperatureTable = m_temperature;

  // Multiply the E/p values by the pressure.
  for (auto& field : m_eFields) field *= m_pressureTable;

  // Scale the parameters.
  const double sqrp = sqrt(m_pressureTable);
  const double logp = log(m_pressureTable);
  for (int i = nE; i--;) {
    for (int j = nA; j--;) {
      for (int k = nB; k--;) {
        if (!m_eDiffLong.empty()) m_eDiffLong[j][k][i] /= sqrp;
        if (!m_eDiffTrans.empty()) m_eDiffTrans[j][k][i] /= sqrp;
        if (!m_eDiffTens.empty()) {
          for (int l = 6; l--;) m_eDiffTens[l][j][k][i] /= m_pressureTable;
        }
        if (!m_eTownsend.empty()) m_eTownsend[j][k][i] += logp;
        if (!m_eAttachment.empty()) m_eAttachment[j][k][i] += logp;
        if (!m_ionDissociation.empty()) m_ionDissociation[j][k][i] += logp;
      }
    }
  }

  // Decode the extrapolation and interpolation tables.
  m_extrVel = {lExtrap[0], hExtrap[0]};
  m_intpVel = interpMeth[0];
  // Indices 1 and 2 correspond to velocities along Bt and ExB.
  m_extrDiff = {lExtrap[3], hExtrap[3]};
  m_intpDiff = interpMeth[3];
  m_extrTownsend = {lExtrap[4], hExtrap[4]};
  m_intpTownsend = interpMeth[4];
  m_extrAttachment = {lExtrap[5], hExtrap[5]};
  m_intpAttachment = interpMeth[5];
  m_extrMobility = {lExtrap[6], hExtrap[6]};
  m_intpMobility = interpMeth[6];
  m_extrLorentzAngle = {lExtrap[7], hExtrap[7]};
  m_intpLorentzAngle = interpMeth[7];
  // Index 8: transv. diff.
  m_extrDissociation = {lExtrap[9], hExtrap[9]};
  m_intpDissociation = interpMeth[9];
  // Index 10: diff. tensor
  m_extrExcRates = {lExtrap[11], hExtrap[11]};
  m_intpExcRates = interpMeth[11];
  m_extrIonRates = {lExtrap[12], hExtrap[12]};
  m_intpIonRates = interpMeth[12];

  // Ion diffusion
  m_ionDiffLong.clear();
  if (ionDiffL > 0.) InitTable(nE, nB, nA, m_ionDiffLong, ionDiffL);
  m_ionDiffTrans.clear();
  if (ionDiffT > 0.) InitTable(nE, nB, nA, m_ionDiffTrans, ionDiffT);

  if (m_debug) std::cout << m_className << "::LoadGasFile: Done.\n";

  return true;
}

bool MediumGas::WriteGasFile(const std::string& filename) {

  // Set the gas mixture.
  constexpr int nMagboltzGases = 60;
  std::vector<double> mixture(nMagboltzGases, 0.);
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    int ng = 0;
    if (!GetGasNumberGasFile(m_gas[i], ng) || ng == 0) {
      std::cerr << m_className << "::WriteGasFile:\n"
                << "    Error retrieving gas number for " << m_gas[i] << ".\n";
    } else {
      mixture[ng - 1] = m_fraction[i] * 100.;
    }
  }

  if (m_debug) {
    std::cout << m_className << "::WriteGasFile:\n"
              << "    Writing gas tables to file " << filename << "\n";
  }

  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out);
  if (!outfile.is_open()) {
    std::cerr << m_className << "::WriteGasFile:\n"
              << "    Cannot open file " << filename << ".\n";
    outfile.close();
    return false;
  }

  // Assemble the GASOK bits.
  std::string gasBits = "FFFFFFFFFFFFFFFFFFFF";
  if (!m_eVelocityE.empty()) gasBits[0] = 'T';
  if (!m_ionMobility.empty()) gasBits[1] = 'T';
  if (!m_eDiffLong.empty()) gasBits[2] = 'T';
  if (!m_eTownsend.empty()) gasBits[3] = 'T';
  // Cluster size distribution; skipped
  if (!m_eAttachment.empty()) gasBits[5] = 'T';
  if (!m_eLorentzAngle.empty()) gasBits[6] = 'T';
  if (!m_eDiffTrans.empty()) gasBits[7] = 'T';
  if (!m_eVelocityB.empty()) gasBits[8] = 'T';
  if (!m_eVelocityExB.empty()) gasBits[9] = 'T';
  if (!m_eDiffTens.empty()) gasBits[10] = 'T';
  if (!m_ionDissociation.empty()) gasBits[11] = 'T';
  // SRIM, HEED; skipped
  if (!m_excRates.empty()) gasBits[14] = 'T';
  if (!m_ionRates.empty()) gasBits[15] = 'T';

  // Get the current time.
  time_t rawtime = time(0);
  tm timeinfo = *localtime(&rawtime);
  char datebuf[80] = {0};
  char timebuf[80] = {0};
  // Format date and time.
  strftime(datebuf, sizeof(datebuf) - 1, "%d/%m/%y", &timeinfo);
  strftime(timebuf, sizeof(timebuf) - 1, "%H.%M.%S", &timeinfo);
  // Set the member name.
  std::string member = "< none >";
  // Write the header.
  outfile << "*----.----1----.----2----.----3----.----4----.----"
          << "5----.----6----.----7----.----8----.----9----.---"
          << "10----.---11----.---12----.---13--\n";
  outfile << "% Created " << datebuf << " at " << timebuf << " ";
  outfile << member << " GAS      ";
  // Add remark.
  std::string buffer;
  outfile << "\"none" << std::string(25, ' ') << "\"\n";
  const int version = 12;
  outfile << " Version   : " << version << "\n";
  outfile << " GASOK bits: " << gasBits << "\n";
  std::stringstream ids;
  ids.str("");
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    ids << m_gas[i] << " " << 100. * m_fraction[i] << "%, ";
  }
  ids << "T=" << m_temperatureTable << " K, "
      << "p=" << m_pressureTable / AtmosphericPressure << " atm";
  outfile << " Identifier: " << std::setw(80) << std::left << ids.str() << "\n";
  outfile << " Clusters  : " << std::string(80, ' ') << "\n";
  outfile << " Dimension : ";
  if (m_map2d) {
    outfile << "T ";
  } else {
    outfile << "F ";
  }

  const unsigned int nE = m_eFields.size();
  const unsigned int nB = m_bFields.size();
  const unsigned int nA = m_bAngles.size();
  outfile << FormatInt(nE, 9) << " " << FormatInt(nA, 9) << " "
          << FormatInt(nB, 9) << " " 
          << FormatInt(m_excLevels.size(), 9) << " " 
          << FormatInt(m_ionLevels.size(), 9) << "\n"; 
  // Store reduced electric fields (E/p).
  outfile << " E fields   \n";
  std::vector<double> efields = m_eFields;
  for (auto& field: efields) field /= m_pressure;
  int cnt = 0;
  // List 5 values, then new line.
  PrintArray(efields, outfile, cnt, 5); 
  if (nE % 5 != 0) outfile << "\n";
  // Store angles.
  outfile << " E-B angles \n";
  cnt = 0;
  PrintArray(m_bAngles, outfile, cnt, 5);
  if (nA % 5 != 0) outfile << "\n";
  // Store B fields (convert to hGauss).
  outfile << " B fields   \n";
  std::vector<double> bfields = m_bFields;
  for (auto& field: bfields) field *= 100.;
  cnt = 0;
  PrintArray(bfields, outfile, cnt, 5);
  if (nB % 5 != 0) outfile << "\n";

  // Store the gas composition.
  outfile << " Mixture:   \n";
  cnt = 0;
  PrintArray(mixture, outfile, cnt, 5);
  if (nMagboltzGases % 5 != 0) outfile << "\n";

  cnt = 0;
  for (const auto& exc: m_excLevels) {
    ++cnt;
    outfile << " Excitation " << FormatInt(cnt, 5) << ": " << std::setw(45)
            << exc.label << "  " << FormatFloat(exc.energy)
            << FormatFloat(exc.prob) << FormatFloat(exc.rms)
            << FormatFloat(exc.dt) << "\n";
  }
  cnt = 0;
  for (const auto& ion : m_ionLevels) {
    ++cnt;
    outfile << " Ionisation " << FormatInt(cnt, 5) << ": " << std::setw(45)
            << ion.label << "  " << FormatFloat(ion.energy) << "\n";
  }

  const double sqrp = sqrt(m_pressureTable);
  const double logp = log(m_pressureTable);

  outfile << " The gas tables follow:\n";
  cnt = 0;
  for (unsigned int i = 0; i < nE; ++i) {
    for (unsigned int j = 0; j < nA; ++j) {
      for (unsigned int k = 0; k < nB; ++k) {
        // Get the velocities.
        double ve = m_eVelocityE.empty() ? 0. : m_eVelocityE[j][k][i];
        double vb = m_eVelocityB.empty() ? 0. : m_eVelocityB[j][k][i];
        double vexb = m_eVelocityExB.empty() ? 0. : m_eVelocityExB[j][k][i];
        // Convert from cm / ns to cm / us.
        ve *= 1.e3;
        vb *= 1.e3;
        vexb *= 1.e3;
        // Make a list of the values to be written, start with the velocities.
        std::vector<double> val;
        if (m_map2d) {
          val = {ve, vb, vexb};
        } else {
          // Add dummy spline values in case of a 1D table.
          val = {ve, 0., vb, 0., vexb, 0.};
        }
        // Get the diffusion coefficients.
        double dl = m_eDiffLong.empty() ? 0. : m_eDiffLong[j][k][i] * sqrp;
        double dt = m_eDiffTrans.empty() ? 0. : m_eDiffTrans[j][k][i] * sqrp;
        // Get the Townsend and attachment coefficients.
        double alpha = m_eTownsend.empty() ? -30. : m_eTownsend[j][k][i] - logp;
        double alpha0 = m_eTownsendNoPenning.empty() ? -30. : m_eTownsendNoPenning[j][k][i] - logp;
        double eta = m_eAttachment.empty() ? -30. : m_eAttachment[j][k][i] - logp;
        // Add them to the list.
        if (m_map2d) {
          val.insert(val.end(), {dl, dt, alpha, alpha0, eta});
        } else {
          val.insert(val.end(), {dl, 0., dt, 0., alpha, 0., alpha0 ,eta, 0.});
        } 
        // Get the ion mobility and convert from cm2 / (V ns) to cm2 / (V us).
        double mu = m_ionMobility.empty() ? 0. :  1.e3 * m_ionMobility[j][k][i];
        // Get the Lorentz angle.
        double lor = m_eLorentzAngle.empty() ? 0 : m_eLorentzAngle[j][k][i];
        // Get the dissociation coefficient.
        double diss = m_ionDissociation.empty() ? -30. : m_ionDissociation[j][k][i] - logp;
        // Add them to the list.
        if (m_map2d) {
          val.insert(val.end(), {mu, lor, diss});
        } else {
          val.insert(val.end(), {mu, 0., lor, 0., diss, 0.});
        } 
        // Get the components of the diffusion tensor.
        for (int l = 0; l < 6; ++l) {
          if (!m_eDiffTens.empty()) { 
            const double cov = m_eDiffTens[l][j][k][i] * m_pressureTable;
            val.push_back(cov);
          } else {
            val.push_back(0.);
          }
          if (!m_map2d) val.push_back(0.);
        }
        // Get the excitation and ionisation rates.
        for (const auto& rexc : m_excRates) {
          val.push_back(rexc[j][k][i]);
          if (!m_map2d) val.push_back(0.);
        }
        for (const auto& rion : m_ionRates) {
          val.push_back(rion[j][k][i]);
          if (!m_map2d) val.push_back(0.);
        }
        PrintArray(val, outfile, cnt, 8);
      }
    }
    if (cnt % 8 != 0) outfile << "\n";
    cnt = 0;
  }

  // Extrapolation methods
  int hExtrap[13], lExtrap[13];
  int interpMeth[13];

  lExtrap[0] = lExtrap[1] = lExtrap[2] = m_extrVel.first;
  hExtrap[0] = hExtrap[1] = hExtrap[2] = m_extrVel.second;
  interpMeth[0] = interpMeth[1] = interpMeth[2] = m_intpVel;
  lExtrap[3] = lExtrap[8] = lExtrap[10] = m_extrDiff.first;
  hExtrap[3] = hExtrap[8] = hExtrap[10] = m_extrDiff.second;
  interpMeth[3] = interpMeth[8] = interpMeth[10] = m_intpDiff;
  lExtrap[4] = m_extrTownsend.first;
  hExtrap[4] = m_extrTownsend.second;
  interpMeth[4] = m_intpTownsend;
  lExtrap[5] = m_extrAttachment.first;
  hExtrap[5] = m_extrAttachment.second;
  interpMeth[5] = m_intpAttachment;
  lExtrap[6] = m_extrMobility.first;
  hExtrap[6] = m_extrMobility.second;
  interpMeth[6] = m_intpMobility;
  // Lorentz angle
  lExtrap[7] = m_extrLorentzAngle.first;
  hExtrap[7] = m_extrLorentzAngle.second;
  interpMeth[7] = m_intpLorentzAngle;
  lExtrap[9] = m_extrDissociation.first;
  hExtrap[9] = m_extrDissociation.second;
  interpMeth[9] = m_intpDissociation;
  lExtrap[11] = m_extrExcRates.first;
  hExtrap[11] = m_extrExcRates.second;
  interpMeth[11] = m_intpExcRates;
  lExtrap[12] = m_extrIonRates.first;
  hExtrap[12] = m_extrIonRates.second;
  interpMeth[12] = m_intpIonRates;

  outfile << " H Extr: ";
  for (int i = 0; i < 13; i++) outfile << FormatInt(hExtrap[i], 5);
  outfile << "\n";
  outfile << " L Extr: ";
  for (int i = 0; i < 13; i++) outfile << FormatInt(lExtrap[i], 5);
  outfile << "\n";
  outfile << " Thresholds: " << FormatInt(thrElectronTownsend, 10)
          << FormatInt(thrElectronAttachment, 10)  
          << FormatInt(thrIonDissociation, 10) << "\n";
  outfile << " Interp: ";
  for (int i = 0; i < 13; i++) outfile << FormatInt(interpMeth[i], 5);;
  outfile << "\n";
  outfile << " A     =" << FormatFloat(0.) << ","
          << " Z     =" << FormatFloat(0.) << ","
          << " EMPROB=" << FormatFloat(0.) << ","
          << " EPAIR =" << FormatFloat(0.) << "\n";
  const double dli = m_ionDiffLong.empty() ? 0. : m_ionDiffLong[0][0][0];
  const double dti = m_ionDiffTrans.empty() ? 0. : m_ionDiffTrans[0][0][0];
  outfile << " Ion diffusion: " << FormatFloat(dli) << FormatFloat(dti) << "\n";
  outfile << " CMEAN =" << FormatFloat(0.) << ","
          << " RHO   =" << FormatFloat(0.) << ","
          << " PGAS  =" << FormatFloat(m_pressureTable) << ","
          << " TGAS  =" << FormatFloat(m_temperatureTable) << "\n";
  outfile << " CLSTYP    : NOT SET   \n"
          << " FCNCLS    : " << std::string(80, ' ') << "\n"
          << " NCLS      : " << FormatInt(0, 10) << "\n"
          << " Average   : " << FormatFloat(0., 25, 18) << "\n"
          << "  Heed initialisation done: F\n"
          << "  SRIM initialisation done: F\n";
  outfile.close();

  return true;
}

void MediumGas::PrintGas() {

  // Print a summary.
  std::cout << m_className << "::PrintGas:\n"
            << "    Gas composition: " << m_name;
  if (m_nComponents > 1) {
    std::cout << " (" << m_fraction[0] * 100;
    for (unsigned int i = 1; i < m_nComponents; ++i) {
      std::cout << "/" << m_fraction[i] * 100;
    }
    std::cout << ")";
  }
  std::cout << "\n";
  std::cout << "    Pressure:    " << m_pressure << " Torr\n";
  std::cout << "    Temperature: " << m_temperature << " K\n";
  std::cout << "    Gas file:\n";
  std::cout << "      Pressure:    " << m_pressureTable << " Torr\n";
  std::cout << "      Temperature: " << m_temperatureTable << " K\n";
  if (m_eFields.size() > 1) {
    std::cout << "    Electric field range:  " << m_eFields[0] << " - "
              << m_eFields.back() << " V/cm in " << m_eFields.size() - 1
              << " steps.\n";
  } else if (m_eFields.size() == 1) {
    std::cout << "    Electric field:        " << m_eFields[0] << " V/cm\n";
  } else {
    std::cout << "    Electric field range: not set\n";
  }
  if (m_bFields.size() > 1) {
    std::cout << "    Magnetic field range:  " << m_bFields[0] << " - "
              << m_bFields.back() << " T in " << m_bFields.size() - 1
              << " steps.\n";
  } else if (m_bFields.size() == 1) {
    std::cout << "    Magnetic field:        " << m_bFields[0] << "\n";
  } else {
    std::cout << "    Magnetic field range: not set\n";
  }
  if (m_bAngles.size() > 1) {
    std::cout << "    Angular range:         " << m_bAngles[0] << " - "
              << m_bAngles.back() << " in " << m_bAngles.size() - 1 
              << " steps.\n";
  } else if (m_bAngles.size() == 1) {
    std::cout << "    Angle between E and B: " << m_bAngles[0] << "\n";
  } else {
    std::cout << "    Angular range: not set\n";
  }

  std::cout << "    Available electron transport data:\n";
  if (!m_eVelocityE.empty()) {
    std::cout << "      Velocity along E\n";
  }
  if (!m_eVelocityB.empty()) {
    std::cout << "      Velocity along Bt\n";
  }
  if (!m_eVelocityExB.empty()) {
    std::cout << "      Velocity along ExB\n";
  }
  if (!m_eVelocityE.empty() || !m_eVelocityB.empty() || !m_eVelocityExB.empty()) {
    PrintExtrapolation(m_extrVel);
    std::cout << "        Interpolation order: " << m_intpVel << "\n";
  }
  if (!m_eDiffLong.empty()) {
    std::cout << "      Longitudinal diffusion coefficient\n";
  }
  if (!m_eDiffTrans.empty()) {
    std::cout << "      Transverse diffusion coefficient\n";
  }
  if (!m_eDiffTens.empty()) {
    std::cout << "      Diffusion tensor\n";
  }
  if (!m_eDiffLong.empty() || !m_eDiffTrans.empty() || !m_eDiffTens.empty()) {
    PrintExtrapolation(m_extrDiff);
    std::cout << "        Interpolation order: " << m_intpDiff << "\n";
  }
  if (!m_eTownsend.empty()) {
    std::cout << "      Townsend coefficient\n";
    PrintExtrapolation(m_extrTownsend);
    std::cout << "        Interpolation order: " << m_intpTownsend << "\n";
  }
  if (!m_eAttachment.empty()) {
    std::cout << "      Attachment coefficient\n";
    PrintExtrapolation(m_extrAttachment);
    std::cout << "        Interpolation order: " << m_intpAttachment << "\n";
  }
  if (!m_eLorentzAngle.empty()) {
    std::cout << "      Lorentz Angle\n";
    PrintExtrapolation(m_extrLorentzAngle);
    std::cout << "        Interpolation order: " << m_intpLorentzAngle << "\n";
  }
  if (!m_excRates.empty()) {
    std::cout << "      Excitation rates\n";
    PrintExtrapolation(m_extrExcRates);
    std::cout << "        Interpolation order: " << m_intpExcRates << "\n";
  }
  if (!m_ionRates.empty()) {
    std::cout << "      Ionisation rates\n";
    PrintExtrapolation(m_extrIonRates);
    std::cout << "        Interpolation order: " << m_intpIonRates << "\n";
  }
  if (m_eVelocityE.empty() && m_eVelocityB.empty() && m_eVelocityExB.empty() &&
      m_eDiffLong.empty() && m_eDiffTrans.empty() && m_eDiffTens.empty() && 
      m_eTownsend.empty() && m_eAttachment.empty() && 
      m_excRates.empty() && m_ionRates.empty() && m_eLorentzAngle.empty()) {
    std::cout << "      none\n";
  }

  std::cout << "    Available ion transport data:\n";
  if (!m_ionMobility.empty()) {
    std::cout << "      Mobility\n";
    PrintExtrapolation(m_extrMobility);
    std::cout << "        Interpolation order: " << m_intpMobility << "\n";
  }
  if (!m_ionDiffLong.empty()) {
    std::cout << "      Longitudinal diffusion coefficient\n";
  }
  if (!m_ionDiffTrans.empty()) {
    std::cout << "      Transverse diffusion coefficient\n";
  }
  if (!m_ionDiffLong.empty() || !m_ionDiffTrans.empty()) {
    PrintExtrapolation(m_extrDiff);
    std::cout << "        Interpolation order: " << m_intpDiff << "\n";
  }
  if (!m_ionDissociation.empty()) {
    std::cout << "      Dissociation coefficient\n";
    PrintExtrapolation(m_extrDissociation);
    std::cout << "        Interpolation order: " << m_intpDissociation << "\n";
  }
  if (m_ionMobility.empty() && m_ionDiffLong.empty() && m_ionDiffTrans.empty() &&
      m_ionDissociation.empty()) {
    std::cout << "      none\n";
  }
}

bool MediumGas::LoadIonMobility(const std::string& filename) {

  // Open the file.
  std::ifstream infile;
  infile.open(filename.c_str(), std::ios::in);
  // Make sure the file could actually be opened.
  if (!infile) {
    std::cerr << m_className << "::LoadIonMobility:\n"
              << "    Error opening file " << filename << ".\n";;
    return false;
  }

  double field = -1., mu = -1.;
  double lastField = field;
  std::vector<double> efields;
  std::vector<double> mobilities;

  // Read the file line by line.
  char line[100];

  int i = 0;
  while (!infile.eof()) {
    ++i;
    // Read the next line.
    infile.getline(line, 100);

    char* token = strtok(line, " ,\t");
    if (!token) {
      break;
    } else if (strcmp(token, "#") == 0 || strcmp(token, "*") == 0 ||
               strcmp(token, "//") == 0) {
      continue;
    } else {
      field = atof(token);
      token = strtok(NULL, " ,\t");
      if (!token) {
        std::cerr << m_className << "::LoadIonMobility:\n"
                  << "    Found E/N but no mobility before the end-of-line.\n"
                  << "    Skipping line " << i << ".\n";
        continue;
      }
      mu = atof(token);
    }
    token = strtok(NULL, " ,\t");
    if (token && strcmp(token, "//") != 0) {
      std::cerr << m_className << "::LoadIonMobility:\n"
                << "    Unexpected non-comment characters after the mobility.\n"
                << "    Skipping line " << i << ".\n";
      continue;
    }
    if (m_debug) {
      std::cout << "    E/N = " << field << " Td: mu = " << mu << " cm2/(Vs)\n";
    }
    // Check if the data has been read correctly.
    if (infile.fail() && !infile.eof()) {
      std::cerr << m_className << "::LoadIonMobility:\n";
      std::cerr << "    Error reading file\n";
      std::cerr << "    " << filename << " (line " << i << ").\n";
      return false;
    }
    // Make sure the values make sense.
    // Negative field values are not allowed.
    if (field < 0.) {
      std::cerr << m_className << "::LoadIonMobility:\n";
      std::cerr << "    Negative electric field (line " << i << ").\n";
      return false;
    }
    // The table has to be in ascending order.
    if (field <= lastField) {
      std::cerr << m_className << "::LoadIonMobility:\n";
      std::cerr << "    Table is not in ascending order (line " << i << ").\n";
      return false;
    }
    // Add the values to the list.
    efields.push_back(field);
    mobilities.push_back(mu);
    lastField = field;
  }

  const int ne = efields.size();
  if (ne <= 0) {
    std::cerr << m_className << "::LoadIonMobilities:\n";
    std::cerr << "    No valid data found.\n";
    return false;
  }

  // The E/N values in the file are supposed to be in Td (10^-17 V cm2).
  const double scaleField = 1.e-17 * GetNumberDensity();
  // The reduced mobilities in the file are supposed to be in cm2/(V s).
  const double scaleMobility =
      1.e-9 * (AtmosphericPressure / m_pressure) * (m_temperature / ZeroCelsius);
  for (int j = ne; j--;) {
    // Scale the fields and mobilities.
    efields[j] *= scaleField;
    mobilities[j] *= scaleMobility;
  }

  std::cout << m_className << "::LoadIonMobility:\n";
  std::cout << "    Read " << ne << " values from file " << filename << "\n";

  return SetIonMobility(efields, mobilities);
}

void MediumGas::SetExtrapolationMethodExcitationRates(
    const std::string& low, const std::string& high) {

  SetExtrapolationMethod(low, high, m_extrExcRates, "ExcitationRates");
}

void MediumGas::SetExtrapolationMethodIonisationRates(
    const std::string& low, const std::string& high) {

  SetExtrapolationMethod(low, high, m_extrIonRates, "IonisationRates");
}

void MediumGas::SetInterpolationMethodExcitationRates(const int intrp) {

  if (intrp > 0) m_intpExcRates = intrp;
}

void MediumGas::SetInterpolationMethodIonisationRates(const int intrp) {

  if (intrp > 0) m_intpIonRates = intrp;
}

bool MediumGas::GetGasInfo(const std::string& gasname, double& a,
                           double& z) const {

  if (gasname == "CF4") {
    a = 12.0107 + 4 * 18.9984032;
    z = 6 + 4 * 9;
    return true;
  } else if (gasname == "Ar") {
    a = 39.948;
    z = 18;
  } else if (gasname == "He") {
    a = 4.002602;
    z = 2;
  } else if (gasname == "He-3") {
    a = 3.01602931914;
    z = 2;
  } else if (gasname == "Ne") {
    a = 20.1797;
    z = 10;
  } else if (gasname == "Kr") {
    a = 37.798;
    z = 36;
  } else if (gasname == "Xe") {
    a = 131.293;
    z = 54;
  } else if (gasname == "CH4") {
    a = 12.0107 + 4 * 1.00794;
    z = 6 + 4;
  } else if (gasname == "C2H6") {
    a = 2 * 12.0107 + 6 * 1.00794;
    z = 2 * 6 + 6;
  } else if (gasname == "C3H8") {
    a = 3 * 12.0107 + 8 * 1.00794;
    z = 3 * 6 + 8;
  } else if (gasname == "iC4H10") {
    a = 4 * 12.0107 + 10 * 1.00794;
    z = 4 * 6 + 10;
  } else if (gasname == "CO2") {
    a = 12.0107 + 2 * 15.9994;
    z = 6 + 2 * 8;
  } else if (gasname == "neoC5H12") {
    a = 5 * 12.0107 + 12 * 1.00794;
    z = 5 * 6 + 12;
  } else if (gasname == "H2O") {
    a = 2 * 1.00794 + 15.9994;
    z = 2 + 8;
  } else if (gasname == "O2") {
    a = 2 * 15.9994;
    z = 2 * 8;
  } else if (gasname == "N2") {
    a = 2 * 14.0067;
    z = 2 * 7;
  } else if (gasname == "NO") {
    a = 14.0067 + 15.9994;
    z = 7 + 8;
  } else if (gasname == "N2O") {
    a = 2 * 14.0067 + 15.9994;
    z = 2 * 7 + 8;
  } else if (gasname == "C2H4") {
    a = 2 * 12.0107 + 4 * 1.00794;
    z = 2 * 6 + 4;
  } else if (gasname == "C2H2") {
    a = 2 * 12.0107 + 2 * 1.00794;
    z = 2 * 6 + 2;
  } else if (gasname == "H2") {
    a = 2 * 1.00794;
    z = 2;
  } else if (gasname == "D2") {
    a = 2 * 2.01410177785;
    z = 2;
  } else if (gasname == "CO") {
    a = 12.0107 + 15.9994;
    z = 6 + 8;
  } else if (gasname == "Methylal") {
    a = 3 * 12.0107 + 8 * 1.00794 + 2 * 15.9994;
    z = 3 * 6 + 8 + 2 * 8;
  } else if (gasname == "DME") {
    a = 4 * 12.0107 + 10 * 1.00794 + 2 * 15.9994;
    z = 4 * 6 + 10 + 2 * 8;
  } else if (gasname == "Reid-Step" || gasname == "Mawell-Model" ||
             gasname == "Reid-Ramp") {
    a = 1.;
    z = 1.;
  } else if (gasname == "C2F6") {
    a = 2 * 12.0107 + 6 * 18.9984032;
    z = 2 * 6 + 6 * 9;
  } else if (gasname == "SF6") {
    a = 32.065 + 6 * 18.9984032;
    z = 16 + 6 * 9;
  } else if (gasname == "NH3") {
    a = 14.0067 + 3 * 1.00794;
    z = 7 + 3;
  } else if (gasname == "C3H6") {
    a = 3 * 12.0107 + 6 * 1.00794;
    z = 3 * 6 + 6;
  } else if (gasname == "cC3H6") {
    a = 3 * 12.0107 + 6 * 1.00794;
    z = 3 * 6 + 6;
  } else if (gasname == "CH3OH") {
    a = 12.0107 + 4 * 1.00794 + 15.9994;
    z = 6 + 4 + 8;
  } else if (gasname == "C2H5OH") {
    a = 2 * 12.0107 + 6 * 1.00794 + 15.9994;
    z = 2 * 6 + 6 + 8;
  } else if (gasname == "C3H7OH") {
    a = 3 * 12.0107 + 8 * 1.00794 + 15.9994;
    z = 3 * 6 + 8 * 8;
  } else if (gasname == "Cs") {
    a = 132.9054519;
    z = 55;
  } else if (gasname == "F2") {
    a = 2 * 18.9984032;
    z = 2 * 9;
  } else if (gasname == "CS2") {
    a = 12.0107 + 2 * 32.065;
    z = 6 + 2 * 16;
  } else if (gasname == "COS") {
    a = 12.0107 + 15.9994 + 32.065;
    z = 6 + 8 + 16;
  } else if (gasname == "CD4") {
    a = 12.0107 + 4 * 2.01410177785;
    z = 6 + 4;
  } else if (gasname == "BF3") {
    a = 10.811 + 3 * 18.9984032;
    z = 5 + 3 * 9;
  } else if (gasname == "C2H2F4") {
    a = 2 * 12.0107 + 2 * 1.00794 + 4 * 18.9984032;
    z = 2 * 6 + 2 + 4 * 9;
  } else if (gasname == "CHF3") {
    a = 12.0107 + 1.00794 + 3 * 18.9984032;
    z = 6 + 1 + 3 * 9;
  } else if (gasname == "CF3Br") {
    a = 12.0107 + 3 * 18.9984032 + 79.904;
    z = 6 + 3 * 9 + 35;
  } else if (gasname == "C3F8") {
    a = 3 * 12.0107 + 8 * 18.9984032;
    z = 3 * 6 + 8 * 9;
  } else if (gasname == "O3") {
    a = 3 * 15.9994;
    z = 3 * 8;
  } else if (gasname == "Hg") {
    a = 2 * 200.59;
    z = 80;
  } else if (gasname == "H2S") {
    a = 2 * 1.00794 + 32.065;
    z = 2 + 16;
  } else if (gasname == "nC4H10") {
    a = 4 * 12.0107 + 10 * 1.00794;
    z = 4 * 6 + 10;
  } else if (gasname == "nC5H12") {
    a = 5 * 12.0107 + 12 * 1.00794;
    z = 5 * 6 + 12;
  } else if (gasname == "N2") {
    a = 2 * 14.0067;
    z = 2 * 7;
  } else if (gasname == "GeH4") {
    a = 72.64 + 4 * 1.00794;
    z = 32 + 4;
  } else if (gasname == "SiH4") {
    a = 28.0855 + 4 * 1.00794;
    z = 14 + 4;
  } else {
    a = 0.;
    z = 0.;
    return false;
  }

  return true;
}

bool MediumGas::GetGasName(const int gasnumber, const int version,
                           std::string& gasname) {

  switch (gasnumber) {
    case 1:
      gasname = "CF4";
      break;
    case 2:
      gasname = "Ar";
      break;
    case 3:
      gasname = "He";
      break;
    case 4:
      gasname = "He-3";
      break;
    case 5:
      gasname = "Ne";
      break;
    case 6:
      gasname = "Kr";
      break;
    case 7:
      gasname = "Xe";
      break;
    case 8:
      gasname = "CH4";
      break;
    case 9:
      gasname = "C2H6";
      break;
    case 10:
      gasname = "C3H8";
      break;
    case 11:
      gasname = "iC4H10";
      break;
    case 12:
      gasname = "CO2";
      break;
    case 13:
      gasname = "neoC5H12";
      break;
    case 14:
      gasname = "H2O";
      break;
    case 15:
      gasname = "O2";
      break;
    case 16:
      gasname = "N2";
      break;
    case 17:
      gasname = "NO";
      break;
    case 18:
      gasname = "N2O";
      break;
    case 19:
      gasname = "C2H4";
      break;
    case 20:
      gasname = "C2H2";
      break;
    case 21:
      gasname = "H2";
      break;
    case 22:
      gasname = "D2";
      break;
    case 23:
      gasname = "CO";
      break;
    case 24:
      gasname = "Methylal";
      break;
    case 25:
      gasname = "DME";
      break;
    case 26:
      gasname = "Reid-Step";
      break;
    case 27:
      gasname = "Maxwell-Model";
      break;
    case 28:
      gasname = "Reid-Ramp";
      break;
    case 29:
      gasname = "C2F6";
      break;
    case 30:
      gasname = "SF6";
      break;
    case 31:
      gasname = "NH3";
      break;
    case 32:
      gasname = "C3H6";
      break;
    case 33:
      gasname = "cC3H6";
      break;
    case 34:
      gasname = "CH3OH";
      break;
    case 35:
      gasname = "C2H5OH";
      break;
    case 36:
      gasname = "C3H7OH";
      break;
    case 37:
      gasname = "Cs";
      break;
    case 38:
      gasname = "F2";
      break;
    case 39:
      gasname = "CS2";
      break;
    case 40:
      gasname = "COS";
      break;
    case 41:
      gasname = "CD4";
      break;
    case 42:
      gasname = "BF3";
      break;
    case 43:
      gasname = "C2H2F4";
      break;
    case 44:
      if (version <= 11) {
        gasname = "He-3";
      } else {
        gasname = "TMA";
      }
      break;
    case 45:
      gasname = "He";
      break;
    case 46:
      gasname = "Ne";
      break;
    case 47:
      gasname = "Ar";
      break;
    case 48:
      gasname = "Kr";
      break;
    case 49:
      gasname = "Xe";
      break;
    case 50:
      gasname = "CHF3";
      break;
    case 51:
      gasname = "CF3Br";
      break;
    case 52:
      gasname = "C3F8";
      break;
    case 53:
      gasname = "O3";
      break;
    case 54:
      gasname = "Hg";
      break;
    case 55:
      gasname = "H2S";
      break;
    case 56:
      gasname = "nC4H10";
      break;
    case 57:
      gasname = "nC5H12";
      break;
    case 58:
      gasname = "N2";
      break;
    case 59:
      gasname = "GeH4";
      break;
    case 60:
      gasname = "SiH4";
      break;
    default:
      gasname = "";
      return false;
      break;
  }
  return true;
}

bool MediumGas::GetGasName(std::string input, std::string& gasname) const {

  // Convert to upper-case
  for (unsigned int i = 0; i < input.length(); ++i) {
    input[i] = toupper(input[i]);
  }

  gasname = "";

  if (input == "") return false;

  // CF4
  if (input == "CF4" || input == "FREON" || input == "FREON-14" ||
      input == "TETRAFLUOROMETHANE") {
    gasname = "CF4";
    return true;
  }
  // Argon
  if (input == "AR" || input == "ARGON") {
    gasname = "Ar";
    return true;
  }
  // Helium 4
  if (input == "HE" || input == "HELIUM" || input == "HE-4" ||
      input == "HE 4" || input == "HE4" || input == "4-HE" || input == "4 HE" ||
      input == "4HE" || input == "HELIUM-4" || input == "HELIUM 4" ||
      input == "HELIUM4") {
    gasname = "He";
    return true;
  }
  // Helium 3
  if (input == "HE-3" || input == "HE3" || input == "HELIUM-3" ||
      input == "HELIUM 3" || input == "HELIUM3") {
    gasname = "He-3";
    return true;
  }
  // Neon
  if (input == "NE" || input == "NEON") {
    gasname = "Ne";
    return true;
  }
  // Krypton
  if (input == "KR" || input == "KRYPTON") {
    gasname = "Kr";
    return true;
  }
  // Xenon
  if (input == "XE" || input == "XENON") {
    gasname = "Xe";
    return true;
  }
  // Methane
  if (input == "CH4" || input == "METHANE") {
    gasname = "CH4";
    return true;
  }
  // Ethane
  if (input == "C2H6" || input == "ETHANE") {
    gasname = "C2H6";
    return true;
  }
  // Propane
  if (input == "C3H8" || input == "PROPANE") {
    gasname = "C3H8";
    return true;
  }
  // Isobutane
  if (input == "C4H10" || input == "ISOBUTANE" || input == "ISO" ||
      input == "IC4H10" || input == "ISO-C4H10" || input == "ISOC4H10") {
    gasname = "iC4H10";
    return true;
  }
  // Carbon dioxide (CO2)
  if (input == "CO2" || input == "CARBON-DIOXIDE" ||
      input == "CARBON DIOXIDE" || input == "CARBONDIOXIDE") {
    gasname = "CO2";
    return true;
  }
  // Neopentane
  if (input == "NEOPENTANE" || input == "NEO-PENTANE" || input == "NEO-C5H12" ||
      input == "NEOC5H12" || input == "DIMETHYLPROPANE" || input == "C5H12") {
    gasname = "neoC5H12";
    return true;
  }
  // Water
  if (input == "H2O" || input == "WATER" || input == "WATER-VAPOUR" ||
      input == "WATER VAPOUR") {
    gasname = "H2O";
    return true;
  }
  // Oxygen
  if (input == "O2" || input == "OXYGEN") {
    gasname = "O2";
    return true;
  }
  // Nitrogen
  if (input == "NI" || input == "NITRO" || input == "N2" ||
      input == "NITROGEN") {
    gasname = "N2";
    return true;
  }
  // Nitric oxide (NO)
  if (input == "NO" || input == "NITRIC-OXIDE" || input == "NITRIC OXIDE" ||
      input == "NITROGEN-MONOXIDE" || input == "NITROGEN MONOXIDE") {
    gasname = "NO";
    return true;
  }
  // Nitrous oxide (N2O)
  if (input == "N2O" || input == "NITROUS-OXIDE" || input == "NITROUS OXIDE" ||
      input == "DINITROGEN-MONOXIDE" || input == "LAUGHING-GAS") {
    gasname = "N2O";
    return true;
  }
  // Ethene (C2H4)
  if (input == "C2H4" || input == "ETHENE" || input == "ETHYLENE") {
    gasname = "C2H4";
    return true;
  }
  // Acetylene (C2H2)
  if (input == "C2H2" || input == "ACETYL" || input == "ACETYLENE" ||
      input == "ETHYNE") {
    gasname = "C2H2";
    return true;
  }
  // Hydrogen
  if (input == "H2" || input == "HYDROGEN") {
    gasname = "H2";
    return true;
  }
  // Deuterium
  if (input == "D2" || input == "DEUTERIUM") {
    gasname = "D2";
    return true;
  }
  // Carbon monoxide (CO)
  if (input == "CO" || input == "CARBON-MONOXIDE" ||
      input == "CARBON MONOXIDE") {
    gasname = "CO";
    return true;
  }
  // Methylal (dimethoxymethane, CH3-O-CH2-O-CH3, "hot" version)
  if (input == "METHYLAL" || input == "METHYLAL-HOT" || input == "DMM" ||
      input == "DIMETHOXYMETHANE" || input == "FORMAL" || input == "C3H8O2") {
    gasname = "Methylal";
    return true;
  }
  // DME
  if (input == "DME" || input == "DIMETHYL-ETHER" || input == "DIMETHYLETHER" ||
      input == "DIMETHYL ETHER" || input == "METHYL ETHER" ||
      input == "METHYL-ETHER" || input == "METHYLETHER" ||
      input == "WOOD-ETHER" || input == "WOODETHER" || input == "WOOD ETHER" ||
      input == "DIMETHYL OXIDE" || input == "DIMETHYL-OXIDE" ||
      input == "DEMEON" || input == "METHOXYMETHANE" || input == "C4H10O2") {
    gasname = "DME";
    return true;
  }
  // Reid step
  if (input == "REID-STEP") {
    gasname = "Reid-Step";
    return true;
  }
  // Maxwell model
  if (input == "MAXWELL-MODEL") {
    gasname = "Maxwell-Model";
    return true;
  }
  // Reid ramp
  if (input == "REID-RAMP") {
    gasname = "Reid-Ramp";
    return true;
  }
  // C2F6
  if (input == "C2F6" || input == "FREON-116" || input == "ZYRON-116" ||
      input == "ZYRON-116-N5" || input == "HEXAFLUOROETHANE") {
    gasname = "C2F6";
    return true;
  }
  // SF6
  if (input == "SF6" || input == "SULPHUR-HEXAFLUORIDE" ||
      input == "SULFUR-HEXAFLUORIDE" || input == "SULPHUR HEXAFLUORIDE" ||
      input == "SULFUR HEXAFLUORIDE") {
    gasname = "SF6";
    return true;
  }
  // NH3
  if (input == "NH3" || input == "AMMONIA") {
    gasname = "NH3";
    return true;
  }
  // Propene
  if (input == "C3H6" || input == "PROPENE" || input == "PROPYLENE") {
    gasname = "C3H6";
    return true;
  }
  // Cyclopropane
  if (input == "C-PROPANE" || input == "CYCLO-PROPANE" ||
      input == "CYCLO PROPANE" || input == "CYCLOPROPANE" ||
      input == "C-C3H6" || input == "CC3H6" || input == "CYCLO-C3H6") {
    gasname = "cC3H6";
    return true;
  }
  // Methanol
  if (input == "METHANOL" || input == "METHYL-ALCOHOL" ||
      input == "METHYL ALCOHOL" || input == "WOOD ALCOHOL" ||
      input == "WOOD-ALCOHOL" || input == "CH3OH") {
    gasname = "CH3OH";
    return true;
  }
  // Ethanol
  if (input == "ETHANOL" || input == "ETHYL-ALCOHOL" ||
      input == "ETHYL ALCOHOL" || input == "GRAIN ALCOHOL" ||
      input == "GRAIN-ALCOHOL" || input == "C2H5OH") {
    gasname = "C2H5OH";
    return true;
  }
  // Propanol
  if (input == "PROPANOL" || input == "2-PROPANOL" || input == "ISOPROPYL" ||
      input == "ISO-PROPANOL" || input == "ISOPROPANOL" ||
      input == "ISOPROPYL ALCOHOL" || input == "ISOPROPYL-ALCOHOL" ||
      input == "C3H7OH") {
    gasname = "C3H7OH";
    return true;
  }
  // Cesium / Caesium.
  if (input == "CS" || input == "CESIUM" || input == "CAESIUM") {
    gasname = "Cs";
    return true;
  }
  // Fluorine
  if (input == "F2" || input == "FLUOR" || input == "FLUORINE") {
    gasname = "F2";
    return true;
  }
  // CS2
  if (input == "CS2" || input == "CARBON-DISULPHIDE" ||
      input == "CARBON-DISULFIDE" || input == "CARBON DISULPHIDE" ||
      input == "CARBON DISULFIDE") {
    gasname = "CS2";
    return true;
  }
  // COS
  if (input == "COS" || input == "CARBONYL-SULPHIDE" ||
      input == "CARBONYL-SULFIDE" || input == "CARBONYL SULFIDE") {
    gasname = "COS";
    return true;
  }
  // Deuterated methane
  if (input == "DEUT-METHANE" || input == "DEUTERIUM-METHANE" ||
      input == "DEUTERATED-METHANE" || input == "DEUTERATED METHANE" ||
      input == "DEUTERIUM METHANE" || input == "CD4") {
    gasname = "CD4";
    return true;
  }
  // BF3
  if (input == "BF3" || input == "BORON-TRIFLUORIDE" ||
      input == "BORON TRIFLUORIDE") {
    gasname = "BF3";
    return true;
  }
  // C2H2F4 (and C2HF5).
  if (input == "C2HF5" || input == "C2H2F4" || input == "C2F5H" ||
      input == "C2F4H2" || input == "FREON 134" || input == "FREON 134A" ||
      input == "FREON-134" || input == "FREON-134-A" || input == "FREON 125" ||
      input == "ZYRON 125" || input == "FREON-125" || input == "ZYRON-125" ||
      input == "TETRAFLUOROETHANE" || input == "PENTAFLUOROETHANE") {
    gasname = "C2H2F4";
    return true;
  }
  // TMA
  if (input == "TMA" || input == "TRIMETHYLAMINE" || input == "N(CH3)3" ||
      input == "N-(CH3)3") {
    gasname = "TMA";
    return true;
  }
  // CHF3
  if (input == "CHF3" || input == "FREON-23" || input == "TRIFLUOROMETHANE" ||
      input == "FLUOROFORM") {
    gasname = "CHF3";
    return true;
  }
  // CF3Br
  if (input == "CF3BR" || input == "TRIFLUOROBROMOMETHANE" ||
      input == "BROMOTRIFLUOROMETHANE" || input == "HALON-1301" ||
      input == "HALON 1301" || input == "FREON-13B1" || input == "FREON 13BI") {
    gasname = "CF3Br";
    return true;
  }
  // C3F8
  if (input == "C3F8" || input == "OCTAFLUOROPROPANE" || input == "R218" ||
      input == "R-218" || input == "FREON 218" || input == "FREON-218" ||
      input == "PERFLUOROPROPANE" || input == "RC 218" || input == "PFC 218" ||
      input == "RC-218" || input == "PFC-218" || input == "FLUTEC PP30" ||
      input == "GENETRON 218") {
    gasname = "C3F8";
    return true;
  }
  // Ozone
  if (input == "OZONE" || input == "O3") {
    gasname = "O3";
    return true;
  }
  // Mercury
  if (input == "MERCURY" || input == "HG" || input == "HG2") {
    gasname = "Hg";
    return true;
  }
  // H2S
  if (input == "H2S" || input == "HYDROGEN SULPHIDE" || input == "SEWER GAS" ||
      input == "HYDROGEN-SULPHIDE" || input == "SEWER-GAS" ||
      input == "HYDROGEN SULFIDE" || input == "HEPATIC ACID" ||
      input == "HYDROGEN-SULFIDE" || input == "HEPATIC-ACID" ||
      input == "SULFUR HYDRIDE" || input == "DIHYDROGEN MONOSULFIDE" ||
      input == "SULFUR-HYDRIDE" || input == "DIHYDROGEN-MONOSULFIDE" ||
      input == "DIHYDROGEN MONOSULPHIDE" || input == "SULPHUR HYDRIDE" ||
      input == "DIHYDROGEN-MONOSULPHIDE" || input == "SULPHUR-HYDRIDE" ||
      input == "STINK DAMP" || input == "SULFURATED HYDROGEN" ||
      input == "STINK-DAMP" || input == "SULFURATED-HYDROGEN") {
    gasname = "H2S";
    return true;
  }
  // n-Butane
  if (input == "N-BUTANE" || input == "N-C4H10" || input == "NBUTANE" ||
      input == "NC4H10") {
    gasname = "nC4H10";
    return true;
  }
  // n-Pentane
  if (input == "N-PENTANE" || input == "N-C5H12" || input == "NPENTANE" ||
      input == "NC5H12") {
    gasname = "nC5H12";
    return true;
  }
  // Nitrogen
  if (input == "NI-PHELPS" || input == "NI PHELPS" ||
      input == "NITROGEN-PHELPS" || input == "NITROGEN PHELPHS" ||
      input == "N2-PHELPS" || input == "N2 PHELPS" || input == "N2 (PHELPS)") {
    gasname = "N2 (Phelps)";
    return true;
  }
  // Germane, GeH4
  if (input == "GERMANE" || input == "GERM" || input == "GERMANIUM-HYDRIDE" ||
      input == "GERMANIUM HYDRIDE" || input == "GERMANIUM TETRAHYDRIDE" ||
      input == "GERMANIUM-TETRAHYDRIDE" || input == "GERMANOMETHANE" ||
      input == "MONOGERMANE" || input == "GEH4") {
    gasname = "GeH4";
    return true;
  }
  // Silane, SiH4
  if (input == "SILANE" || input == "SIL" || input == "SILICON-HYDRIDE" ||
      input == "SILICON HYDRIDE" || input == "SILICON-TETRAHYDRIDE" ||
      input == "SILICANE" || input == "MONOSILANE" || input == "SIH4") {
    gasname = "SiH4";
    return true;
  }

  std::cerr << m_className << "::GetGasName:\n";
  std::cerr << "    Gas " << input << " is not defined.\n";
  return false;
}

bool MediumGas::GetGasNumberGasFile(const std::string& input,
                                    int& number) const {

  if (input == "") {
    number = 0;
    return false;
  }

  // CF4
  if (input == "CF4") {
    number = 1;
    return true;
  }
  // Argon
  if (input == "Ar") {
    number = 2;
    return true;
  }
  // Helium 4
  if (input == "He" || input == "He-4") {
    number = 3;
    return true;
  }
  // Helium 3
  if (input == "He-3") {
    number = 4;
    return true;
  }
  // Neon
  if (input == "Ne") {
    number = 5;
    return true;
  }
  // Krypton
  if (input == "Kr") {
    number = 6;
    return true;
  }
  // Xenon
  if (input == "Xe") {
    number = 7;
    return true;
  }
  // Methane
  if (input == "CH4") {
    number = 8;
    return true;
  }
  // Ethane
  if (input == "C2H6") {
    number = 9;
    return true;
  }
  // Propane
  if (input == "C3H8") {
    number = 10;
    return true;
  }
  // Isobutane
  if (input == "iC4H10") {
    number = 11;
    return true;
  }
  // Carbon dioxide (CO2)
  if (input == "CO2") {
    number = 12;
    return true;
  }
  // Neopentane
  if (input == "neoC5H12") {
    number = 13;
    return true;
  }
  // Water
  if (input == "H2O") {
    number = 14;
    return true;
  }
  // Oxygen
  if (input == "O2") {
    number = 15;
    return true;
  }
  // Nitrogen
  if (input == "N2") {
    number = 16;
    return true;
  }
  // Nitric oxide (NO)
  if (input == "NO") {
    number = 17;
    return true;
  }
  // Nitrous oxide (N2O)
  if (input == "N2O") {
    number = 18;
    return true;
  }
  // Ethene (C2H4)
  if (input == "C2H4") {
    number = 19;
    return true;
  }
  // Acetylene (C2H2)
  if (input == "C2H2") {
    number = 20;
    return true;
  }
  // Hydrogen
  if (input == "H2") {
    number = 21;
    return true;
  }
  // Deuterium
  if (input == "D2") {
    number = 22;
    return true;
  }
  // Carbon monoxide (CO)
  if (input == "CO") {
    number = 23;
    return true;
  }
  // Methylal (dimethoxymethane, CH3-O-CH2-O-CH3, "hot" version)
  if (input == "Methylal") {
    number = 24;
    return true;
  }
  // DME
  if (input == "DME") {
    number = 25;
    return true;
  }
  // Reid step
  if (input == "Reid-Step") {
    number = 26;
    return true;
  }
  // Maxwell model
  if (input == "Maxwell-Model") {
    number = 27;
    return true;
  }
  // Reid ramp
  if (input == "Reid-Ramp") {
    number = 28;
    return true;
  }
  // C2F6
  if (input == "C2F6") {
    number = 29;
    return true;
  }
  // SF6
  if (input == "SF6") {
    number = 30;
    return true;
  }
  // NH3
  if (input == "NH3") {
    number = 31;
    return true;
  }
  // Propene
  if (input == "C3H6") {
    number = 32;
    return true;
  }
  // Cyclopropane
  if (input == "cC3H6") {
    number = 33;
    return true;
  }
  // Methanol
  if (input == "CH3OH") {
    number = 34;
    return true;
  }
  // Ethanol
  if (input == "C2H5OH") {
    number = 35;
    return true;
  }
  // Propanol
  if (input == "C3H7OH") {
    number = 36;
    return true;
  }
  // Cesium / Caesium.
  if (input == "Cs") {
    number = 37;
    return true;
  }
  // Fluorine
  if (input == "F2") {
    number = 38;
    return true;
  }
  // CS2
  if (input == "CS2") {
    number = 39;
    return true;
  }
  // COS
  if (input == "COS") {
    number = 40;
    return true;
  }
  // Deuterated methane
  if (input == "CD4") {
    number = 41;
    return true;
  }
  // BF3
  if (input == "BF3") {
    number = 42;
    return true;
  }
  // C2HF5 and C2H2F4.
  if (input == "C2HF5" || input == "C2H2F4") {
    number = 43;
    return true;
  }
  // TMA
  if (input == "TMA") {
    number = 44;
    return true;
  }
  // CHF3
  if (input == "CHF3") {
    number = 50;
    return true;
  }
  // CF3Br
  if (input == "CF3Br") {
    number = 51;
    return true;
  }
  // C3F8
  if (input == "C3F8") {
    number = 52;
    return true;
  }
  // Ozone
  if (input == "O3") {
    number = 53;
    return true;
  }
  // Mercury
  if (input == "Hg") {
    number = 54;
    return true;
  }
  // H2S
  if (input == "H2S") {
    number = 55;
    return true;
  }
  // n-butane
  if (input == "nC4H10") {
    number = 56;
    return true;
  }
  // n-pentane
  if (input == "nC5H12") {
    number = 57;
    return true;
  }
  // Nitrogen
  if (input == "N2 (Phelps)") {
    number = 58;
    return true;
  }
  // Germane, GeH4
  if (input == "GeH4") {
    number = 59;
    return true;
  }
  // Silane, SiH4
  if (input == "SiH4") {
    number = 60;
    return true;
  }

  std::cerr << m_className << "::GetGasNumberGasFile:\n";
  std::cerr << "    Gas " << input << " is not defined.\n";
  return false;
}

bool MediumGas::GetPhotoabsorptionCrossSection(const double& e, double& sigma,
                                               const unsigned int& i) {

  if (i >= m_nMaxGases) {
    std::cerr << m_className << "::GetPhotoabsorptionCrossSection:\n";
    std::cerr << "    Index (" << i << ") out of range.\n";
    return false;
  }

  OpticalData optData;
  if (!optData.IsAvailable(m_gas[i])) return false;
  double eta = 0.;
  return optData.GetPhotoabsorptionCrossSection(m_gas[i], e, sigma, eta);
}
}
