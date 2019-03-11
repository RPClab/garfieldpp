#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "FundamentalConstants.hh"
#include "GarfieldConstants.hh"
#include "MediumGas.hh"
#include "OpticalData.hh"
#include "Utilities.hh"

namespace {

std::string FmtFloat(const double x, const unsigned int width = 15,
                        const unsigned int precision = 8) {
  char buffer[256];
  std::snprintf(buffer, width + 1, "%*.*E", width, precision, x);
  return std::string(buffer);
}

std::string FmtInt(const int n, const unsigned int width) {
  char buffer[256];
  std::snprintf(buffer, width + 1, "%*d", width, n);
  return std::string(buffer);
}

void PrintArray(const std::vector<double>& values, std::ofstream& outfile,
                int& col, const int ncols) {
  for (const auto value : values) {
    outfile << FmtFloat(value);
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

MediumGas::MediumGas()
    : Medium(), m_pressureTable(m_pressure), m_temperatureTable(m_temperature) {
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

  // Reset all transport parameters.
  ResetTables();
  // Force a recalculation of the collision rates.
  m_isChanged = true;

  m_nComponents = 0;
  m_gas.fill("");
  m_fraction.fill(0.);
  m_atWeight.fill(0.);
  m_atNum.fill(0.);
  for (unsigned int i = 0; i < 6; ++i) {
    if (fractions[i] < Small) continue;
    // Find the gas name corresponding to the input string.
    const std::string gasname = GetGasName(gases[i]);
    if (!gasname.empty()) {
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
                << "    Using Penning transfer parameters for " << m_gas[i]
                << " from previous mixture.\n"
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

void MediumGas::GetComponent(const unsigned int i, std::string& label,
                             double& f) {
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

  ResetTables();

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
          std::cout << m_className << "::LoadGasFile: Version " 
                    << version << "\n";
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
          m_tab2d = false;
        } else {
          m_tab2d = true;
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
        if (m_tab2d && nA <= 0) {
          std::cerr << m_className << "::LoadGasFile:\n"
                    << "    Number of E-B angles out of range.\n";
          gasfile.close();
          return false;
        }

        token = strtok(NULL, " :,%\t");
        nB = atoi(token);
        // Check the number of B points.
        if (m_tab2d && nB <= 0) {
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
          std::cout << "    " << nexc << " excitations, " << nion
                    << " ionisations.\n";
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
    std::cout << m_className << "::LoadGasFile:\n    " << nE
              << " electric fields, " << nB << " magnetic fields, " << nA
              << " angles.\n    ";
  }
  if (gasBits[0] == 'T') InitTable(nE, nB, nA, m_eVelE, 0.);
  if (gasBits[1] == 'T') InitTable(nE, nB, nA, m_iMob, 0.);
  if (gasBits[2] == 'T') InitTable(nE, nB, nA, m_eDifL, 0.);
  if (gasBits[3] == 'T') {
    InitTable(nE, nB, nA, m_eAlp, -30.);
    InitTable(nE, nB, nA, m_eAlpNoPenning, -30.);
  }
  // gasBits[4]: cluster size distribution; skipped
  if (gasBits[5] == 'T') InitTable(nE, nB, nA, m_eAtt, -30.);
  if (gasBits[6] == 'T') InitTable(nE, nB, nA, m_eLor, -30.);
  if (gasBits[7] == 'T') InitTable(nE, nB, nA, m_eDifT, 0.);
  if (gasBits[8] == 'T') InitTable(nE, nB, nA, m_eVelB, 0.);
  if (gasBits[9] == 'T') InitTable(nE, nB, nA, m_eVelX, 0.);
  if (gasBits[10] == 'T') InitTensor(nE, nB, nA, 6, m_eDifM, 0.);
  if (gasBits[11] == 'T') InitTable(nE, nB, nA, m_iDis, -30.);
  // gasBits[12]: SRIM; skipped
  // gasBits[13]: HEED; skipped
  if (gasBits[14] == 'T') {
    InitTensor(nE, nB, nA, m_excLevels.size(), m_excRates, 0.);
  }
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
    const std::string gasname = GetGasName(i + 1, version);
    if (gasname.empty()) {
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

  if (m_tab2d) {
    if (m_debug) {
      std::cout << m_className << "::LoadGasFile: Loading 3D table.\n";
    }
    for (int i = 0; i < nE; i++) {
      for (int j = 0; j < nA; j++) {
        for (int k = 0; k < nB; k++) {
          // Drift velocity along E, Bt and ExB
          double ve = 0., vb = 0., vx = 0.;
          gasfile >> ve >> vb >> vx;
          // Convert from cm / us to cm / ns and add to the table.
          if (!m_eVelE.empty()) m_eVelE[j][k][i] = 1.e-3 * ve;
          if (!m_eVelB.empty()) m_eVelB[j][k][i] = 1.e-3 * vb;
          if (!m_eVelX.empty()) m_eVelX[j][k][i] = 1.e-3 * vx;
          // Longitudinal and transverse diffusion coefficients
          double dl = 0., dt = 0.;
          gasfile >> dl >> dt;
          if (!m_eDifL.empty()) m_eDifL[j][k][i] = dl;
          if (!m_eDifT.empty()) m_eDifT[j][k][i] = dt;
          // Townsend and attachment coefficients
          double alpha = 0., alpha0 = 0., eta = 0.;
          gasfile >> alpha >> alpha0 >> eta;
          if (!m_eAlp.empty()) {
            m_eAlp[j][k][i] = alpha;
            m_eAlpNoPenning[j][k][i] = alpha0;
          }
          if (!m_eAtt.empty()) m_eAtt[j][k][i] = eta;
          // Ion mobility
          double mu = 0.;
          gasfile >> mu;
          // Convert from cm2 / (V us) to cm2 / (V ns)
          if (!m_iMob.empty()) m_iMob[j][k][i] = 1.e-3 * mu;
          // Lorentz angle
          double lor = 0.;
          gasfile >> lor;
          if (!m_eLor.empty()) m_eLor[j][k][i] = lor;
          // Ion dissociation
          double dis = 0.;
          gasfile >> dis;
          if (!m_iDis.empty()) m_iDis[j][k][i] = dis;
          // Diffusion tensor
          for (int l = 0; l < 6; l++) {
            double diff = 0.;
            gasfile >> diff;
            if (!m_eDifM.empty()) m_eDifM[l][j][k][i] = diff;
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
      double ve = 0., vb = 0., vx = 0.;
      gasfile >> ve >> waste >> vb >> waste >> vx >> waste;
      if (!m_eVelE.empty()) m_eVelE[0][0][i] = 1.e-3 * ve;
      if (!m_eVelB.empty()) m_eVelB[0][0][i] = 1.e-3 * vb;
      if (!m_eVelX.empty()) m_eVelX[0][0][i] = 1.e-3 * vx;
      // Longitudinal and transferse diffusion coefficients
      double dl = 0., dt = 0.;
      gasfile >> dl >> waste >> dt >> waste;
      if (!m_eDifL.empty()) m_eDifL[0][0][i] = dl;
      if (!m_eDifT.empty()) m_eDifT[0][0][i] = dt;
      // Townsend and attachment coefficients
      double alpha = 0., alpha0 = 0., eta = 0.;
      gasfile >> alpha >> waste >> alpha0 >> eta >> waste;
      if (!m_eAlp.empty()) {
        m_eAlp[0][0][i] = alpha;
        m_eAlpNoPenning[0][0][i] = alpha0;
      }
      if (!m_eAtt.empty()) {
        m_eAtt[0][0][i] = eta;
      }
      // Ion mobility
      double mu = 0.;
      gasfile >> mu >> waste;
      if (!m_iMob.empty()) m_iMob[0][0][i] = 1.e-3 * mu;
      // Lorentz angle
      double lor = 0.;
      gasfile >> lor >> waste;
      if (!m_eLor.empty()) m_eLor[0][0][i] = lor;
      // Ion dissociation
      double diss = 0.;
      gasfile >> diss >> waste;
      if (!m_iDis.empty()) m_iDis[0][0][i] = diss;
      // Diffusion tensor
      for (int j = 0; j < 6; j++) {
        double diff = 0.;
        gasfile >> diff >> waste;
        if (!m_eDifM.empty()) m_eDifM[j][0][0][i] = diff;
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
        if (token != NULL) m_eThrAlp = atoi(token);
        token = strtok(NULL, " :,%=\t");
        if (token != NULL) m_eThrAtt = atoi(token);
        token = strtok(NULL, " :,%=\t");
        if (token != NULL) m_iThrDis = atoi(token);
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
        if (!m_eDifL.empty()) m_eDifL[j][k][i] /= sqrp;
        if (!m_eDifT.empty()) m_eDifT[j][k][i] /= sqrp;
        if (!m_eDifM.empty()) {
          for (int l = 6; l--;) m_eDifM[l][j][k][i] /= m_pressureTable;
        }
        if (!m_eAlp.empty()) m_eAlp[j][k][i] += logp;
        if (!m_eAtt.empty()) m_eAtt[j][k][i] += logp;
        if (!m_iDis.empty()) m_iDis[j][k][i] += logp;
        /*
        for (auto& exc : m_excRates) {
          exc[j][k][i] /= m_pressureTable;
        }
        for (auto& ion : m_ionRates) {
          ion[j][k][i] /= m_pressureTable;
        }
        */
      }
    }
  }

  // Decode the extrapolation and interpolation tables.
  m_extrVel = {lExtrap[0], hExtrap[0]};
  m_intpVel = interpMeth[0];
  // Indices 1 and 2 correspond to velocities along Bt and ExB.
  m_extrDif = {lExtrap[3], hExtrap[3]};
  m_intpDif = interpMeth[3];
  m_extrAlp = {lExtrap[4], hExtrap[4]};
  m_intpAlp = interpMeth[4];
  m_extrAtt = {lExtrap[5], hExtrap[5]};
  m_intpAtt = interpMeth[5];
  m_extrMob = {lExtrap[6], hExtrap[6]};
  m_intpMob = interpMeth[6];
  m_extrLor = {lExtrap[7], hExtrap[7]};
  m_intpLor = interpMeth[7];
  // Index 8: transv. diff.
  m_extrDis = {lExtrap[9], hExtrap[9]};
  m_intpDis = interpMeth[9];
  // Index 10: diff. tensor
  m_extrExc = {lExtrap[11], hExtrap[11]};
  m_intpExc = interpMeth[11];
  m_extrIon = {lExtrap[12], hExtrap[12]};
  m_intpIon = interpMeth[12];

  // Ion diffusion
  if (ionDiffL > 0.) InitTable(nE, nB, nA, m_iDifL, ionDiffL);
  if (ionDiffT > 0.) InitTable(nE, nB, nA, m_iDifT, ionDiffT);

  if (m_debug) std::cout << m_className << "::LoadGasFile: Done.\n";

  return true;
}

bool MediumGas::WriteGasFile(const std::string& filename) {
  // Set the gas mixture.
  constexpr int nMagboltzGases = 60;
  std::vector<double> mixture(nMagboltzGases, 0.);
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const int ng = GetGasNumberGasFile(m_gas[i]);
    if (ng <= 0) {
      std::cerr << m_className << "::WriteGasFile:\n"
                << "    Error retrieving gas number for " << m_gas[i] << ".\n";
      continue;
    }
    mixture[ng - 1] = m_fraction[i] * 100.;
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
  if (!m_eVelE.empty()) gasBits[0]  = 'T';
  if (!m_iMob.empty())  gasBits[1]  = 'T';
  if (!m_eDifL.empty()) gasBits[2]  = 'T';
  if (!m_eAlp.empty())  gasBits[3]  = 'T';
  // Cluster size distribution; skipped
  if (!m_eAtt.empty())  gasBits[5]  = 'T';
  if (!m_eLor.empty())  gasBits[6]  = 'T';
  if (!m_eDifT.empty()) gasBits[7]  = 'T';
  if (!m_eVelB.empty()) gasBits[8]  = 'T';
  if (!m_eVelX.empty()) gasBits[9]  = 'T';
  if (!m_eDifM.empty()) gasBits[10] = 'T';
  if (!m_iDis.empty())  gasBits[11] = 'T';
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
  if (m_tab2d) {
    outfile << "T ";
  } else {
    outfile << "F ";
  }

  const unsigned int nE = m_eFields.size();
  const unsigned int nB = m_bFields.size();
  const unsigned int nA = m_bAngles.size();
  outfile << FmtInt(nE, 9) << " " << FmtInt(nA, 9) << " "
          << FmtInt(nB, 9) << " " << FmtInt(m_excLevels.size(), 9) << " "
          << FmtInt(m_ionLevels.size(), 9) << "\n";
  // Store reduced electric fields (E/p).
  outfile << " E fields   \n";
  std::vector<double> efields = m_eFields;
  for (auto& field : efields) field /= m_pressure;
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
  for (auto& field : bfields) field *= 100.;
  cnt = 0;
  PrintArray(bfields, outfile, cnt, 5);
  if (nB % 5 != 0) outfile << "\n";

  // Store the gas composition.
  outfile << " Mixture:   \n";
  cnt = 0;
  PrintArray(mixture, outfile, cnt, 5);
  if (nMagboltzGases % 5 != 0) outfile << "\n";

  cnt = 0;
  for (const auto& exc : m_excLevels) {
    ++cnt;
    outfile << " Excitation " << FmtInt(cnt, 5) << ": " << std::setw(45)
            << exc.label << "  " << FmtFloat(exc.energy) << FmtFloat(exc.prob) 
            << FmtFloat(exc.rms) << FmtFloat(exc.dt) << "\n";
  }
  cnt = 0;
  for (const auto& ion : m_ionLevels) {
    ++cnt;
    outfile << " Ionisation " << FmtInt(cnt, 5) << ": " << std::setw(45)
            << ion.label << "  " << FmtFloat(ion.energy) << "\n";
  }

  const double sqrp = sqrt(m_pressureTable);
  const double logp = log(m_pressureTable);

  outfile << " The gas tables follow:\n";
  cnt = 0;
  for (unsigned int i = 0; i < nE; ++i) {
    for (unsigned int j = 0; j < nA; ++j) {
      for (unsigned int k = 0; k < nB; ++k) {
        // Get the velocities.
        double ve = m_eVelE.empty() ? 0. : m_eVelE[j][k][i];
        double vb = m_eVelB.empty() ? 0. : m_eVelB[j][k][i];
        double vx = m_eVelX.empty() ? 0. : m_eVelX[j][k][i];
        // Convert from cm / ns to cm / us.
        ve *= 1.e3;
        vb *= 1.e3;
        vx *= 1.e3;
        // Make a list of the values to be written, start with the velocities.
        std::vector<double> val;
        if (m_tab2d) {
          val = {ve, vb, vx};
        } else {
          // Add dummy spline values in case of a 1D table.
          val = {ve, 0., vb, 0., vx, 0.};
        }
        // Get the diffusion coefficients.
        double dl = m_eDifL.empty() ? 0. : m_eDifL[j][k][i] * sqrp;
        double dt = m_eDifT.empty() ? 0. : m_eDifT[j][k][i] * sqrp;
        // Get the Townsend and attachment coefficients.
        double alpha = m_eAlp.empty() ? -30. : m_eAlp[j][k][i] - logp;
        double alpha0 = m_eAlpNoPenning.empty() ? -30. : m_eAlpNoPenning[j][k][i] - logp;
        double eta = m_eAtt.empty() ? -30. : m_eAtt[j][k][i] - logp;
        // Add them to the list.
        if (m_tab2d) {
          val.insert(val.end(), {dl, dt, alpha, alpha0, eta});
        } else {
          val.insert(val.end(), {dl, 0., dt, 0., alpha, 0., alpha0, eta, 0.});
        }
        // Get the ion mobility and convert from cm2 / (V ns) to cm2 / (V us).
        double mu = m_iMob.empty() ? 0. : 1.e3 * m_iMob[j][k][i];
        // Get the Lorentz angle.
        double lor = m_eLor.empty() ? 0 : m_eLor[j][k][i];
        // Get the dissociation coefficient.
        double diss = m_iDis.empty() ? -30. : m_iDis[j][k][i] - logp;
        // Add them to the list.
        if (m_tab2d) {
          val.insert(val.end(), {mu, lor, diss});
        } else {
          val.insert(val.end(), {mu, 0., lor, 0., diss, 0.});
        }
        // Get the components of the diffusion tensor.
        for (int l = 0; l < 6; ++l) {
          if (!m_eDifM.empty()) {
            const double cov = m_eDifM[l][j][k][i] * m_pressureTable;
            val.push_back(cov);
          } else {
            val.push_back(0.);
          }
          if (!m_tab2d) val.push_back(0.);
        }
        // Get the excitation and ionisation rates.
        for (const auto& rexc : m_excRates) {
          val.push_back(rexc[j][k][i]);
          if (!m_tab2d) val.push_back(0.);
        }
        for (const auto& rion : m_ionRates) {
          val.push_back(rion[j][k][i]);
          if (!m_tab2d) val.push_back(0.);
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
  lExtrap[3] = lExtrap[8] = lExtrap[10] = m_extrDif.first;
  hExtrap[3] = hExtrap[8] = hExtrap[10] = m_extrDif.second;
  interpMeth[3] = interpMeth[8] = interpMeth[10] = m_intpDif;
  lExtrap[4] = m_extrAlp.first;
  hExtrap[4] = m_extrAlp.second;
  interpMeth[4] = m_intpAlp;
  lExtrap[5] = m_extrAtt.first;
  hExtrap[5] = m_extrAtt.second;
  interpMeth[5] = m_intpAtt;
  lExtrap[6] = m_extrMob.first;
  hExtrap[6] = m_extrMob.second;
  interpMeth[6] = m_intpMob;
  // Lorentz angle
  lExtrap[7] = m_extrLor.first;
  hExtrap[7] = m_extrLor.second;
  interpMeth[7] = m_intpLor;
  lExtrap[9] = m_extrDis.first;
  hExtrap[9] = m_extrDis.second;
  interpMeth[9] = m_intpDis;
  lExtrap[11] = m_extrExc.first;
  hExtrap[11] = m_extrExc.second;
  interpMeth[11] = m_intpExc;
  lExtrap[12] = m_extrIon.first;
  hExtrap[12] = m_extrIon.second;
  interpMeth[12] = m_intpIon;

  outfile << " H Extr: ";
  for (int i = 0; i < 13; i++) outfile << FmtInt(hExtrap[i], 5);
  outfile << "\n";
  outfile << " L Extr: ";
  for (int i = 0; i < 13; i++) outfile << FmtInt(lExtrap[i], 5);
  outfile << "\n";
  outfile << " Thresholds: " << FmtInt(m_eThrAlp, 10)
          << FmtInt(m_eThrAtt, 10) << FmtInt(m_iThrDis, 10) << "\n";
  outfile << " Interp: ";
  for (int i = 0; i < 13; i++) outfile << FmtInt(interpMeth[i], 5);
  outfile << "\n";
  outfile << " A     =" << FmtFloat(0.) << ", Z     =" << FmtFloat(0.) << ","
          << " EMPROB=" << FmtFloat(0.) << ", EPAIR =" << FmtFloat(0.) << "\n";
  const double dli = m_iDifL.empty() ? 0. : m_iDifL[0][0][0];
  const double dti = m_iDifT.empty() ? 0. : m_iDifT[0][0][0];
  outfile << " Ion diffusion: " << FmtFloat(dli) << FmtFloat(dti) << "\n";
  outfile << " CMEAN =" << FmtFloat(0.) << ", RHO   =" << FmtFloat(0.) << ","
          << " PGAS  =" << FmtFloat(m_pressureTable) << ","
          << " TGAS  =" << FmtFloat(m_temperatureTable) << "\n";
  outfile << " CLSTYP    : NOT SET   \n"
          << " FCNCLS    : " << std::string(80, ' ') << "\n"
          << " NCLS      : " << FmtInt(0, 10) << "\n"
          << " Average   : " << FmtFloat(0., 25, 18) << "\n"
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
  std::cout << "    Pressure:    " << m_pressure << " Torr\n"
            << "    Temperature: " << m_temperature << " K\n"
            << "    Gas file:\n"
            << "      Pressure:    " << m_pressureTable << " Torr\n"
            << "      Temperature: " << m_temperatureTable << " K\n";
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
    std::cout << "    Magnetic field:        " << m_bFields[0] << " T\n";
  } else {
    std::cout << "    Magnetic field range: not set\n";
  }
  if (m_bAngles.size() > 1) {
    std::cout << "    Angular range:         " << m_bAngles[0] << " - "
              << m_bAngles.back() << " rad in " << m_bAngles.size() - 1
              << " steps.\n";
  } else if (m_bAngles.size() == 1) {
    std::cout << "    Angle between E and B: " << m_bAngles[0] << " rad\n";
  } else {
    std::cout << "    Angular range: not set\n";
  }

  std::cout << "    Available electron transport data:\n";
  if (!m_eVelE.empty()) {
    std::cout << "      Velocity along E\n";
  }
  if (!m_eVelB.empty()) {
    std::cout << "      Velocity along Bt\n";
  }
  if (!m_eVelX.empty()) {
    std::cout << "      Velocity along ExB\n";
  }
  if (!(m_eVelE.empty() && m_eVelB.empty() && 
        m_eVelX.empty())) {
    PrintExtrapolation(m_extrVel);
    std::cout << "        Interpolation order: " << m_intpVel << "\n";
  }
  if (!m_eDifL.empty()) {
    std::cout << "      Longitudinal diffusion coefficient\n";
  }
  if (!m_eDifT.empty()) {
    std::cout << "      Transverse diffusion coefficient\n";
  }
  if (!m_eDifM.empty()) {
    std::cout << "      Diffusion tensor\n";
  }
  if (!m_eDifL.empty() || !m_eDifT.empty() || !m_eDifM.empty()) {
    PrintExtrapolation(m_extrDif);
    std::cout << "        Interpolation order: " << m_intpDif << "\n";
  }
  if (!m_eAlp.empty()) {
    std::cout << "      Townsend coefficient\n";
    PrintExtrapolation(m_extrAlp);
    std::cout << "        Interpolation order: " << m_intpAlp << "\n";
  }
  if (!m_eAtt.empty()) {
    std::cout << "      Attachment coefficient\n";
    PrintExtrapolation(m_extrAtt);
    std::cout << "        Interpolation order: " << m_intpAtt << "\n";
  }
  if (!m_eLor.empty()) {
    std::cout << "      Lorentz Angle\n";
    PrintExtrapolation(m_extrLor);
    std::cout << "        Interpolation order: " << m_intpLor << "\n";
  }
  if (!m_excRates.empty()) {
    std::cout << "      Excitation rates\n";
    for (const auto& exc : m_excLevels) {
      std::cout << "        " << exc.label << "\n";
      std::cout << "          Energy = " << exc.energy << " eV";
      if (exc.prob > 0.) {
        std::cout << ", Penning transfer probability = " << exc.prob;
      } 
      std::cout << "\n";
    }
    PrintExtrapolation(m_extrExc);
    std::cout << "        Interpolation order: " << m_intpExc << "\n";
  }
  if (!m_ionRates.empty()) {
    std::cout << "      Ionisation rates\n";
    for (const auto& ion : m_ionLevels) {
      std::cout << "        " << ion.label << "\n";
      std::cout << "          Threshold = " << ion.energy << " eV\n";
    }
    PrintExtrapolation(m_extrIon);
    std::cout << "        Interpolation order: " << m_intpIon << "\n";
  }
  if (m_eVelE.empty() && m_eVelB.empty() && m_eVelX.empty() &&
      m_eDifL.empty() && m_eDifT.empty() && m_eDifM.empty() &&
      m_eAlp.empty() && m_eAtt.empty() && m_excRates.empty() &&
      m_ionRates.empty() && m_eLor.empty()) {
    std::cout << "      none\n";
  }

  std::cout << "    Available ion transport data:\n";
  if (!m_iMob.empty()) {
    std::cout << "      Mobility\n";
    PrintExtrapolation(m_extrMob);
    std::cout << "        Interpolation order: " << m_intpMob << "\n";
  }
  if (!m_iDifL.empty()) {
    std::cout << "      Longitudinal diffusion coefficient\n";
  }
  if (!m_iDifT.empty()) {
    std::cout << "      Transverse diffusion coefficient\n";
  }
  if (!m_iDifL.empty() || !m_iDifT.empty()) {
    PrintExtrapolation(m_extrDif);
    std::cout << "        Interpolation order: " << m_intpDif << "\n";
  }
  if (!m_iDis.empty()) {
    std::cout << "      Dissociation coefficient\n";
    PrintExtrapolation(m_extrDis);
    std::cout << "        Interpolation order: " << m_intpDis << "\n";
  }
  if (m_iMob.empty() && m_iDifL.empty() && m_iDifT.empty() && m_iDis.empty()) {
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
              << "    Error opening file " << filename << ".\n";
    ;
    return false;
  }

  double field = -1., mu = -1.;
  double lastField = field;
  std::vector<double> efields;
  std::vector<double> mobilities;

  // Read the file line by line.
  int i = 0;
  while (!infile.eof()) {
    ++i;
    // Read the next line.
    char line[100];
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
  const double scaleMobility = 1.e-9 * (AtmosphericPressure / m_pressure) *
                               (m_temperature / ZeroCelsius);
  for (int j = ne; j--;) {
    // Scale the fields and mobilities.
    efields[j] *= scaleField;
    mobilities[j] *= scaleMobility;
  }

  std::cout << m_className << "::LoadIonMobility:\n";
  std::cout << "    Read " << ne << " values from file " << filename << "\n";

  return SetIonMobility(efields, mobilities);
}

void MediumGas::ResetTables() {

  Medium::ResetTables();
  m_eAlpNoPenning.clear();
  m_excLevels.clear();
  m_ionLevels.clear();
  m_excRates.clear();
  m_ionRates.clear();
}

bool MediumGas::EnablePenningTransfer(const double r,
                                      const double lambda) {

  if (r < 0. || r > 1.) {
    std::cerr << m_className << "::EnablePenningTransfer:\n"
              << "    Transfer probability must be in the range [0, 1].\n";
    return false;
  }

  m_rPenningGlobal = r;
  m_lambdaPenningGlobal = lambda > Small ? lambda : 0.;

  std::cout << m_className << "::EnablePenningTransfer:\n"
            << "    Global Penning transfer parameters set to:\n"
            << "    r      = " << m_rPenningGlobal << "\n"
            << "    lambda = " << m_lambdaPenningGlobal << " cm\n";

  // Find the min. ionisation energy.
  if (m_ionLevels.empty()) {
    std::cerr << m_className << "::EnablePenningTransfer:\n    Warning: present"
              << " gas table has no ionisation rates.\n    Ignore this message "
              << "if you are using microscopic tracking only.\n";
    return true;
  }
  double minIonPot = -1.;
  for (const auto& ion : m_ionLevels) {
    if (minIonPot < 0.) {
      minIonPot = ion.energy;
    } else {
      minIonPot = std::min(minIonPot, ion.energy);
    }
  }

  // Update the transfer probabilities of the excitation levels in the table.
  unsigned int nLevelsFound = 0;
  for (auto& exc : m_excLevels) {
    if (exc.energy < minIonPot) continue;
    exc.prob = m_rPenningGlobal;
    exc.rms = m_lambdaPenningGlobal;
    ++nLevelsFound; 
  }
  if (nLevelsFound > 0) {
    std::cout << m_className << "::EnablePenningTransfer:\n"
              << "    Updated transfer probabilities for " << nLevelsFound 
              << " excitation rates.\n";
    AdjustTownsendCoefficient();
  } else {
    std::cerr << m_className << "::EnablePenningTransfer:\n    Warning: present"
              << " gas table has no eligible excitation rates.\n    Ignore this"
              << " message if you are using microscopic tracking only.\n";
  }
  return true;
}

bool MediumGas::EnablePenningTransfer(const double r, const double lambda,
                                      std::string gasname) {

  if (r < 0. || r > 1.) {
    std::cerr << m_className << "::EnablePenningTransfer:\n"
              << "    Transfer probability must be in the range [0, 1].\n";
    return false;
  }

  // Get the "standard" name of this gas.
  gasname = GetGasName(gasname);
  if (gasname.empty()) {
    std::cerr << m_className << "::EnablePenningTransfer: Unknown gas name.\n";
    return false;
  }

  // Look for this gas in the present gas mixture.
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      m_rPenningGas[i] = r;
      m_lambdaPenningGas[i] = lambda > Small ? lambda : 0.;
      iGas = i;
      break;
    }
  }

  if (iGas < 0) {
    std::cerr << m_className << "::EnablePenningTransfer:\n"
              << "    Requested gas (" << gasname
              << ") is not part of the present gas mixture.\n";
    return false;
  }

  // Find the min. ionisation energy.
  if (m_ionLevels.empty()) {
    std::cerr << m_className << "::EnablePenningTransfer:\n    Warning: present"
              << " gas table has no ionisation rates.\n    Ignore this message"
              << " if you are using microscopic tracking only.\n";
    return true;
  }
  double minIonPot = -1.;
  for (const auto& ion : m_ionLevels) {
    if (minIonPot < 0.) {
      minIonPot = ion.energy;
    } else {
      minIonPot = std::min(minIonPot, ion.energy);
    }
  }
  // Update the transfer probabilities of the excitation levels in the table.
  unsigned int nLevelsFound = 0;
  for (auto& exc : m_excLevels) {
    if (exc.energy < minIonPot) continue;
    // Try to extract the gas name from the label. 
    // TODO: test if this works for all gases and excitation levels.
    const auto pos = exc.label.find('-');
    if (pos == std::string::npos) continue;
    if (GetGasName(exc.label.substr(0, pos)) != gasname) continue;
    exc.prob = r;
    exc.rms = lambda;
    ++nLevelsFound; 
  }
  if (nLevelsFound > 0) {
    std::cout << m_className << "::EnablePenningTransfer:\n"
              << "    Updated transfer probabilities for " << nLevelsFound 
              << " excitation rates.\n";
    AdjustTownsendCoefficient();
  } else {
    std::cerr << m_className << "::EnablePenningTransfer:\n    Warning: present"
              << " gas table has no eligible excitation rates.\n    Ignore this"
              << " message if you are using microscopic tracking only.\n";
  }
  return true;
}

void MediumGas::DisablePenningTransfer() {

  m_rPenningGlobal = 0.;
  m_lambdaPenningGlobal = 0.;

  m_rPenningGas.fill(0.);
  m_lambdaPenningGas.fill(0.);

  if (m_excLevels.empty()) return;
  for (auto& exc : m_excLevels) {
    exc.prob = 0.; 
  }
  AdjustTownsendCoefficient();
}

bool MediumGas::DisablePenningTransfer(std::string gasname) {

  // Get the "standard" name of this gas.
  gasname = GetGasName(gasname);
  if (gasname.empty()) {
    std::cerr << m_className << "::DisablePenningTransfer: Unknown gas name.\n";
    return false;
  }

  // Look for this gas in the present gas mixture.
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      m_rPenningGas[i] = 0.;
      m_lambdaPenningGas[i] = 0.;
      iGas = i;
      break;
    }
  }

  if (iGas < 0) {
    std::cerr << m_className << "::DisablePenningTransfer:\n"
              << "    Requested gas (" << gasname
              << ") is not part of the present gas mixture.\n";
    return false;
  }

  if (m_excLevels.empty()) return true;
  for (auto& exc : m_excLevels) {
    // Try to extract the gas name from the label. 
    const auto pos = exc.label.find('-');
    if (pos == std::string::npos) continue;
    if (GetGasName(exc.label.substr(0, pos)) != gasname) continue;
    exc.prob = 0.; 
  }
  AdjustTownsendCoefficient();
  return true;
}


bool MediumGas::AdjustTownsendCoefficient() {

  // GASSPT

  // Make sure there are Townsend coefficients.
  if (m_eAlp.empty() || m_eAlpNoPenning.empty()) {
    std::cerr << m_className << "::AdjustTownsendCoefficient:\n    "
              << "Present gas table does not include Townsend coefficients.\n";
    return false;
  }
  // Make sure there are excitation and ionisation rates.
  if (m_excLevels.empty() || m_excRates.empty()) {
    std::cerr << m_className << "::AdjustTownsendCoefficient:\n    "
              << "Present gas table does not include excitation rates.\n";
    return false;
  }
  if (m_ionLevels.empty() || m_ionRates.empty()) {
    std::cerr << m_className << "::AdjustTownsendCoefficient:\n    "
              << "Present gas table does not include ionisation rates.\n";
    return false;
  }
  const unsigned int nE = m_eFields.size();
  const unsigned int nB = m_bFields.size();
  const unsigned int nA = m_bAngles.size();
  if (m_debug) {
    std::cout << m_className << "::AdjustTownsendCoefficient:\n"
              << "   Entry         Exc.      Ion.\n";
  } 
  for (unsigned int i = 0; i < nE; ++i) {
    for (unsigned int j = 0; j < nA; ++j) {
      for (unsigned int k = 0; k < nB; ++k) {
        // Compute total ionisation rate.
        double rion = 0.;
        for (const auto& ion : m_ionRates) {
          rion += ion[j][k][i];
        }
        // Compute rate of Penning ionisations.
        double rexc = 0.;
        const unsigned int nexc = m_excLevels.size();
        for (unsigned int ie = 0; ie < nexc; ++ie) {
          rexc += m_excLevels[ie].prob * m_excRates[ie][j][k][i]; 
        }
        if (m_debug) {
          std::cout << FmtInt(i, 4) << FmtInt(j, 4) << FmtInt(k, 4) 
                    << FmtFloat(rexc, 12, 5) << FmtFloat(rion, 12, 5) << "\n"; 
        }
        // Adjust the Townsend coefficient.
        double alpha0 = m_eAlpNoPenning[j][k][i];
        if (alpha0 < -20.) {
          alpha0 = 0.;
        } else {
          alpha0 = m_pressure * exp(alpha0);
        }
        double alpha1 = alpha0;
        if (rion > 0.) alpha1 *= (rexc + rion) / rion;
        m_eAlp[j][k][i] = alpha1 > 0. ? log(alpha1 / m_pressure) : -30.;
      }
    }
  }
  return true;
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
  } else if (gasname == "H2" || gasname == "paraH2") {
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
  } else if (gasname == "C3H7OH" || gasname == "nC3H7OH") {
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

std::string MediumGas::GetGasName(const int gasnumber, const int version) const {

  switch (gasnumber) {
    case 1:
      return "CF4";
      break;
    case 2:
      return "Ar";
      break;
    case 3:
      return "He";
      break;
    case 4:
      return "He-3";
      break;
    case 5:
      return "Ne";
      break;
    case 6:
      return "Kr";
      break;
    case 7:
      return "Xe";
      break;
    case 8:
      return "CH4";
      break;
    case 9:
      return "C2H6";
      break;
    case 10:
      return "C3H8";
      break;
    case 11:
      return "iC4H10";
      break;
    case 12:
      return "CO2";
      break;
    case 13:
      return "neoC5H12";
      break;
    case 14:
      return "H2O";
      break;
    case 15:
      return "O2";
      break;
    case 16:
      return "N2";
      break;
    case 17:
      return "NO";
      break;
    case 18:
      return "N2O";
      break;
    case 19:
      return "C2H4";
      break;
    case 20:
      return "C2H2";
      break;
    case 21:
      return "H2";
      break;
    case 22:
      return "D2";
      break;
    case 23:
      return "CO";
      break;
    case 24:
      return "Methylal";
      break;
    case 25:
      return "DME";
      break;
    case 26:
      return "Reid-Step";
      break;
    case 27:
      return "Maxwell-Model";
      break;
    case 28:
      return "Reid-Ramp";
      break;
    case 29:
      return "C2F6";
      break;
    case 30:
      return "SF6";
      break;
    case 31:
      return "NH3";
      break;
    case 32:
      return "C3H6";
      break;
    case 33:
      return "cC3H6";
      break;
    case 34:
      return "CH3OH";
      break;
    case 35:
      return "C2H5OH";
      break;
    case 36:
      return "C3H7OH";
      break;
    case 37:
      return "Cs";
      break;
    case 38:
      return "F2";
      break;
    case 39:
      return "CS2";
      break;
    case 40:
      return "COS";
      break;
    case 41:
      return "CD4";
      break;
    case 42:
      return "BF3";
      break;
    case 43:
      return "C2H2F4";
      break;
    case 44:
      return version <= 11 ? "He-3" : "TMA";
      break;
    case 45:
      return version <= 11 ? "He" : "paraH2";
      break;
    case 46:
      return version <= 11 ? "Ne" : "nC3H7OH";
      break;
    case 47:
      return "Ar";
      break;
    case 48:
      return "Kr";
      break;
    case 49:
      return "Xe";
      break;
    case 50:
      return "CHF3";
      break;
    case 51:
      return "CF3Br";
      break;
    case 52:
      return "C3F8";
      break;
    case 53:
      return "O3";
      break;
    case 54:
      return "Hg";
      break;
    case 55:
      return "H2S";
      break;
    case 56:
      return "nC4H10";
      break;
    case 57:
      return "nC5H12";
      break;
    case 58:
      return "N2";
      break;
    case 59:
      return "GeH4";
      break;
    case 60:
      return "SiH4";
      break;
    default:
      return "";
      break;
  }
  return "";
}

std::string MediumGas::GetGasName(std::string input) const {
  // Convert to upper-case
  for (unsigned int i = 0; i < input.length(); ++i) {
    input[i] = toupper(input[i]);
  }

  if (input.empty()) return "";

  if (input == "CF4" || input == "FREON" || input == "FREON-14" ||
      input == "TETRAFLUOROMETHANE") {
    return "CF4";
  } else if (input == "AR" || input == "ARGON") {
    return "Ar";
  } else if (input == "HE" || input == "HELIUM" || input == "HE-4" ||
             input == "HE 4" || input == "HE4" || input == "4-HE" || 
             input == "4 HE" || input == "4HE" || input == "HELIUM-4" || 
             input == "HELIUM 4" || input == "HELIUM4") {
    return "He";
  } else if (input == "HE-3" || input == "HE3" || input == "HELIUM-3" ||
      input == "HELIUM 3" || input == "HELIUM3") {
    return "He-3";
  } else if (input == "NE" || input == "NEON") {
    return "Ne";
  } else if (input == "KR" || input == "KRYPTON") {
    return "Kr";
  } else if (input == "XE" || input == "XENON") {
    return "Xe";
  } else if (input == "CH4" || input == "METHANE") {
    return "CH4";
  } else if (input == "C2H6" || input == "ETHANE") {
    return "C2H6";
  } else if (input == "C3H8" || input == "PROPANE") {
    return "C3H8";
  } else  if (input == "C4H10" || input == "ISOBUTANE" || input == "ISO" ||
              input == "IC4H10" || input == "ISO-C4H10" || input == "ISOC4H10") {
    return "iC4H10";
  } else if (input == "CO2" || input == "CARBON-DIOXIDE" ||
      input == "CARBON DIOXIDE" || input == "CARBONDIOXIDE") {
    return "CO2";
  } else if (input == "NEOPENTANE" || input == "NEO-PENTANE" || 
             input == "NEO-C5H12" || input == "NEOC5H12" || 
             input == "DIMETHYLPROPANE" || input == "C5H12") {
    return "neoC5H12";
  } else if (input == "H2O" || input == "WATER" || input == "WATER-VAPOUR" ||
      input == "WATER VAPOUR") {
    return "H2O";
  } else if (input == "O2" || input == "OXYGEN") {
    return "O2";
  } else if (input == "NI" || input == "NITRO" || input == "N2" ||
      input == "NITROGEN") {
    return "N2";
  } else if (input == "NO" || input == "NITRIC-OXIDE" || input == "NITRIC OXIDE" ||
      input == "NITROGEN-MONOXIDE" || input == "NITROGEN MONOXIDE") {
    return "NO";
  } else if (input == "N2O" || input == "NITROUS-OXIDE" || input == "NITROUS OXIDE" ||
      input == "DINITROGEN-MONOXIDE" || input == "LAUGHING-GAS") {
    return "N2O";
  } else if (input == "C2H4" || input == "ETHENE" || input == "ETHYLENE") {
    return "C2H4";
  } else if (input == "C2H2" || input == "ACETYL" || input == "ACETYLENE" ||
      input == "ETHYNE") {
    return "C2H2";
  } else if (input == "H2" || input == "HYDROGEN") {
    return "H2";
  } else if (input == "PARA H2" || input == "PARA-H2" ||
             input == "PARAH2" || input == "PARA HYDROGEN" ||
             input == "PARA-HYDROGEN" || input == "PARAHYDROGEN") {
    return "paraH2";
  } else if (input == "D2" || input == "DEUTERIUM") {
    return "D2";
  } else if (input == "CO" || input == "CARBON-MONOXIDE" ||
      input == "CARBON MONOXIDE") {
    return "CO";
  } else if (input == "METHYLAL" || input == "METHYLAL-HOT" || input == "DMM" ||
      input == "DIMETHOXYMETHANE" || input == "FORMAL" || input == "C3H8O2") {
    // Methylal (dimethoxymethane, CH3-O-CH2-O-CH3, "hot" version)
    return "Methylal";
  } else if (input == "DME" || input == "DIMETHYL-ETHER" || input == "DIMETHYLETHER" ||
      input == "DIMETHYL ETHER" || input == "METHYL ETHER" ||
      input == "METHYL-ETHER" || input == "METHYLETHER" ||
      input == "WOOD-ETHER" || input == "WOODETHER" || input == "WOOD ETHER" ||
      input == "DIMETHYL OXIDE" || input == "DIMETHYL-OXIDE" ||
      input == "DEMEON" || input == "METHOXYMETHANE" || input == "C4H10O2") {
    return "DME";
  } else if (input == "REID-STEP") {
    return "Reid-Step";
  } else if (input == "MAXWELL-MODEL") {
    return "Maxwell-Model";
  } else if (input == "REID-RAMP") {
    return "Reid-Ramp";
  } else if (input == "C2F6" || input == "FREON-116" || input == "ZYRON-116" ||
      input == "ZYRON-116-N5" || input == "HEXAFLUOROETHANE") {
    return "C2F6";
  } else if (input == "SF6" || input == "SULPHUR-HEXAFLUORIDE" ||
      input == "SULFUR-HEXAFLUORIDE" || input == "SULPHUR HEXAFLUORIDE" ||
      input == "SULFUR HEXAFLUORIDE") {
    return "SF6";
  } else if (input == "NH3" || input == "AMMONIA") {
    return "NH3";
  } else if (input == "C3H6" || input == "PROPENE" || input == "PROPYLENE") {
    return "C3H6";
  } else if (input == "C-PROPANE" || input == "CYCLO-PROPANE" ||
      input == "CYCLO PROPANE" || input == "CYCLOPROPANE" ||
      input == "C-C3H6" || input == "CC3H6" || input == "CYCLO-C3H6") {
    return "cC3H6";
  } else if (input == "METHANOL" || input == "METHYL-ALCOHOL" ||
      input == "METHYL ALCOHOL" || input == "WOOD ALCOHOL" ||
      input == "WOOD-ALCOHOL" || input == "CH3OH") {
    return "CH3OH";
  } else if (input == "ETHANOL" || input == "ETHYL-ALCOHOL" ||
      input == "ETHYL ALCOHOL" || input == "GRAIN ALCOHOL" ||
      input == "GRAIN-ALCOHOL" || input == "C2H5OH") {
    return "C2H5OH";
  } else if (input == "PROPANOL" || input == "2-PROPANOL" || input == "ISOPROPYL" ||
      input == "ISO-PROPANOL" || input == "ISOPROPANOL" ||
      input == "ISOPROPYL ALCOHOL" || input == "ISOPROPYL-ALCOHOL" ||
      input == "C3H7OH") {
    return "C3H7OH";
  } else if (input == "NPROPANOL" || input == "N-PROPANOL" || 
             input == "1-PROPANOL" || input == "PROPYL ALCOHOL" ||
             input == "PROPYL-ALCOHOL" || input == "N-PROPYL ALCOHOL" ||
             input == "NC3H7OH" || input == "N-C3H7OH") {
    return "nC3H7OH"; 
  } else if (input == "CS" || input == "CESIUM" || input == "CAESIUM") {
    return "Cs";
  } else if (input == "F2" || input == "FLUOR" || input == "FLUORINE") {
    return "F2";
  } else if (input == "CS2" || input == "CARBON-DISULPHIDE" ||
      input == "CARBON-DISULFIDE" || input == "CARBON DISULPHIDE" ||
      input == "CARBON DISULFIDE") {
    return "CS2";
  } else if (input == "COS" || input == "CARBONYL-SULPHIDE" ||
      input == "CARBONYL-SULFIDE" || input == "CARBONYL SULFIDE") {
    return "COS";
  } else if (input == "DEUT-METHANE" || input == "DEUTERIUM-METHANE" ||
      input == "DEUTERATED-METHANE" || input == "DEUTERATED METHANE" ||
      input == "DEUTERIUM METHANE" || input == "CD4") {
    return "CD4";
  } else if (input == "BF3" || input == "BORON-TRIFLUORIDE" ||
      input == "BORON TRIFLUORIDE") {
    return "BF3";
  } else if (input == "C2HF5" || input == "C2H2F4" || input == "C2F5H" ||
      input == "C2F4H2" || input == "FREON 134" || input == "FREON 134A" ||
      input == "FREON-134" || input == "FREON-134-A" || input == "FREON 125" ||
      input == "ZYRON 125" || input == "FREON-125" || input == "ZYRON-125" ||
      input == "TETRAFLUOROETHANE" || input == "PENTAFLUOROETHANE") {
    // C2H2F4 (and C2HF5).
    return "C2H2F4";
  } else if (input == "TMA" || input == "TRIMETHYLAMINE" || input == "N(CH3)3" ||
      input == "N-(CH3)3") {
    return "TMA";
  } else if (input == "CHF3" || input == "FREON-23" || input == "TRIFLUOROMETHANE" ||
      input == "FLUOROFORM") {
    return "CHF3";
  } else if (input == "CF3BR" || input == "TRIFLUOROBROMOMETHANE" ||
      input == "BROMOTRIFLUOROMETHANE" || input == "HALON-1301" ||
      input == "HALON 1301" || input == "FREON-13B1" || input == "FREON 13BI") {
    return "CF3Br";
  } else if (input == "C3F8" || input == "OCTAFLUOROPROPANE" || input == "R218" ||
      input == "R-218" || input == "FREON 218" || input == "FREON-218" ||
      input == "PERFLUOROPROPANE" || input == "RC 218" || input == "PFC 218" ||
      input == "RC-218" || input == "PFC-218" || input == "FLUTEC PP30" ||
      input == "GENETRON 218") {
    return "C3F8";
  } else if (input == "OZONE" || input == "O3") {
    return "O3";
  } else if (input == "MERCURY" || input == "HG" || input == "HG2") {
    return "Hg";
  } else if (input == "H2S" || input == "HYDROGEN SULPHIDE" || input == "SEWER GAS" ||
      input == "HYDROGEN-SULPHIDE" || input == "SEWER-GAS" ||
      input == "HYDROGEN SULFIDE" || input == "HEPATIC ACID" ||
      input == "HYDROGEN-SULFIDE" || input == "HEPATIC-ACID" ||
      input == "SULFUR HYDRIDE" || input == "DIHYDROGEN MONOSULFIDE" ||
      input == "SULFUR-HYDRIDE" || input == "DIHYDROGEN-MONOSULFIDE" ||
      input == "DIHYDROGEN MONOSULPHIDE" || input == "SULPHUR HYDRIDE" ||
      input == "DIHYDROGEN-MONOSULPHIDE" || input == "SULPHUR-HYDRIDE" ||
      input == "STINK DAMP" || input == "SULFURATED HYDROGEN" ||
      input == "STINK-DAMP" || input == "SULFURATED-HYDROGEN") {
    return "H2S";
  } else if (input == "N-BUTANE" || input == "N-C4H10" || input == "NBUTANE" ||
      input == "NC4H10") {
    return "nC4H10";
  } else if (input == "N-PENTANE" || input == "N-C5H12" || input == "NPENTANE" ||
      input == "NC5H12") {
    return "nC5H12";
  } else if (input == "NI-PHELPS" || input == "NI PHELPS" ||
      input == "NITROGEN-PHELPS" || input == "NITROGEN PHELPHS" ||
      input == "N2-PHELPS" || input == "N2 PHELPS" || input == "N2 (PHELPS)") {
    // Nitrogen
    return "N2 (Phelps)";
  } else if (input == "GERMANE" || input == "GERM" || input == "GERMANIUM-HYDRIDE" ||
      input == "GERMANIUM HYDRIDE" || input == "GERMANIUM TETRAHYDRIDE" ||
      input == "GERMANIUM-TETRAHYDRIDE" || input == "GERMANOMETHANE" ||
      input == "MONOGERMANE" || input == "GEH4") {
    return "GeH4";
  } else if (input == "SILANE" || input == "SIL" || input == "SILICON-HYDRIDE" ||
             input == "SILICON HYDRIDE" || input == "SILICON-TETRAHYDRIDE" ||
             input == "SILICANE" || input == "MONOSILANE" || input == "SIH4") {
    return "SiH4";
  }

  std::cerr << m_className << "::GetGasName:\n"
            << "    Gas " << input << " is not recognized.\n";
  return "";
}

int MediumGas::GetGasNumberGasFile(const std::string& input) const {

  if (input.empty()) return 0;

  if (input == "CF4") {
    return 1;
  } else if (input == "Ar") {
    return 2;
  } else if (input == "He" || input == "He-4") {
    return 3;
  } else if (input == "He-3") {
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
    // Nitric oxide
    return 17;
  } else if (input == "N2O") {
    // Nitrous oxide
    return 18;
  } else if (input == "C2H4") {
    // Ethene
    return 19;
  } else if (input == "C2H2") {
    // Acetylene
    return 20;
  } else if (input == "H2") {
    return 21;
  } else if (input == "D2") {
    // Deuterium
    return 22;
  } else if (input == "CO") {
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
    return 45;
  } else if (input == "nC3H7OH") {
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
    return 54;
  } else if (input == "H2S") {
    return 55;
  } else if (input == "nC4H10") {
    // n-butane
    return 56;
  } else if (input == "nC5H12") {
    // n-pentane
    return 57;
  } else if (input == "N2 (Phelps)") {
    return 58;
  } else if (input == "GeH4") {
    // Germane
    return 59;
  } else if (input == "SiH4") {
    // Silane
    return 60;
  }

  std::cerr << m_className << "::GetGasNumberGasFile:\n"
            << "    Gas " << input << " not found.\n";
  return 0;
}

bool MediumGas::GetPhotoAbsorptionCrossSection(const double e, double& sigma,
                                               const unsigned int i) {
  if (i >= m_nMaxGases) {
    std::cerr << m_className 
              << "::GetPhotoAbsorptionCrossSection: Index out of range.\n";
    return false;
  }

  OpticalData optData;
  if (!optData.IsAvailable(m_gas[i])) return false;
  double eta = 0.;
  return optData.GetPhotoabsorptionCrossSection(m_gas[i], e, sigma, eta);
}
}
