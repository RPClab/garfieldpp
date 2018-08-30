#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "ViewField.hh"
#include "ViewCell.hh"
#include "ComponentAnalyticField.hh"
#include "ComponentConstant.hh"
#include "MediumMagboltz.hh"
#include "SolidBox.hh" 
#include "SolidTube.hh"
#include "GeometrySimple.hh"
#include "Sensor.hh"
#include "ViewDrift.hh"
#include "FundamentalConstants.hh"
#include "DriftLineRKF.hh"
#include "ViewMedium.hh"
#include "ViewSignal.hh"
#include "Random.hh"
#include "TrackSrim.hh" 
#include "AvalancheMicroscopic.hh"

using namespace Garfield;

  // Prepare structure and list for the panels (move to ComponentNebem3D)
  struct Panel {
    // Perpendicular vector
    double a, b, c;
    // Vertices
    std::vector<double> xv;
    std::vector<double> yv;
    std::vector<double> zv;
    // Colour index
    double colour;
    // Reference to solid to which the panel belongs
    int volume;
  };

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Set up gas
  MediumMagboltz* gas = new MediumMagboltz();
  const double pressure = 760.; 
  const double temperature = 293.15;
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("He", 87.5, "CF4", 12.5);

  // Create and clean up the panels structure
  std::vector<Panel> m_panels;
  m_panels.clear();
  // Geometry and a few boxes and cylinders
  GeometrySimple* geo = new GeometrySimple();

  // A box has 23 parameters:
  // 3 half-lengths
  // 3 centres
  // 3 direction vector
  // 4 cos/sin (derived from the direction vectors)
  // 1 voltage bc
  // 1 eps bc
  // 1 boundary type (see below)
  // 1 charge bc
  // 6 discretisation size

  // Boundary types (see neBEM.c)
  // 1: Conducting surfaces
  // 2: Conducting surfaces with known charge
  // 3: Floating conducting surfaces
  // 4: Dielectric interfaces
  // 5: Dielectric interfaces with known charge
  // 6: E parallel symmetry boundary
  // 7: E perpendicular symmetry boundary

  const double x0 = 1, y0 = 2, z0 = 3, lx = 1, ly = 1, lz = 1;
  SolidBox* box1 = new SolidBox(x0, y0, z0, lx, ly, lz);
  // Add the solid to the geometry, together with the medium inside.
  geo->AddSolid(box1, gas);

  const double x1 = 2, y1 = 7, z1 = 4;
  SolidBox* box2 = new SolidBox(x1, y1, z1, lx, ly, lz);
  geo->AddSolid(box2, gas);

  const double x2 = 3, y2 = 4, z2 = 5, rmin = 1, rmax = 2, l = 5;
  SolidTube* tube1 = new SolidTube(x2, y2, z2, rmin, rmax, l);
  geo->AddSolid(tube1, gas);

  // See whether we can recover the solids
  unsigned int nSolids = geo->GetNumberOfSolids();
  printf("Found %d solids\n", nSolids);
  for (int i=0; i<nSolids; i++) {
    int j;
    if(geo->GetSolid(i)->IsBox()) {
      printf("Solid %d is a box\n", i);
      // Check the correct transmission
      double xc, yc, zc, l1, l2, l3, ctheta, stheta, cphi, sphi;
      geo->GetSolid(i)->GetCentre(xc, yc, zc);
      printf("     Centre:       %g %g %g\n", xc, yc, zc);
      l1 = geo->GetSolid(i)->GetHalfLengthX();
      l2 = geo->GetSolid(i)->GetHalfLengthY();
      l3 = geo->GetSolid(i)->GetHalfLengthZ();
      printf("     Half-lengths: %g %g %g\n", l1, l2, l3);
      geo->GetSolid(i)->GetOrientation(ctheta, stheta, cphi, sphi);
      printf("     Orientation:  %g %g %g %g\n", ctheta, stheta, cphi, sphi);
      // Create panels
      geo->GetSolid(i)->SolidPanels();
    } else if (geo->GetSolid(i)->IsTube()) {
      printf("Solid %d is a tube\n", i);
    }
  }
  // See how many penals we got

  printf("Number of panels: %d\n", int(m_panels.size()));

  // Start waiting loop
  app.Run(kTRUE);
}








