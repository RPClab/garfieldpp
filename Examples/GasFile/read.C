#include <iostream>
#include <cstdlib>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "MediumMagboltz.hh"
#include "MediumSilicon.hh"
#include "FundamentalConstants.hh"
#include "ViewMedium.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
 
  // Setup the gas.
  MediumMagboltz gas;
  gas.LoadGasFile("ar_80_co2_20_2T.gas");
  const std::string path = std::getenv("GARFIELD_HOME");
  gas.LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
  gas.PrintGas();

  ViewMedium view;
  view.SetMedium(&gas);
  view.SetMagneticField(2.);
  
  TCanvas cV("cV", "", 600, 600);
  view.SetCanvas(&cV);
  view.PlotElectronVelocity('e');

  TCanvas cD("cD", "", 600, 600);
  view.SetCanvas(&cD);
  view.PlotElectronDiffusion('e');

  TCanvas cT("cT", "", 600, 600);
  view.SetCanvas(&cT);
  view.PlotElectronTownsend('e');

  TCanvas cA("cA", "", 600, 600);
  view.SetCanvas(&cA);
  view.PlotElectronAttachment('e');

  app.Run(kTRUE);

}
