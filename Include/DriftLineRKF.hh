#ifndef G_DRIFTLINE_RKF_H
#define G_DRIFTLINE_RKF_H

#include <vector>
#include <string>

#include "Sensor.hh"
#include "ViewDrift.hh"
#include "Medium.hh"
#include "GeometryBase.hh"


namespace Garfield {
  
  class DriftLineRKF {
    
  public:
 
    DriftLineRKF();
    ~DriftLineRKF() {}

    void SetSensor(Sensor* s);
    void SetIntegrationAccuracy(double intAct);
    void SetMaximumStepSize(double maxStep);

    void EnablePlotting(ViewDrift* view);
    void DisablePlotting();

    void DriftLine(double x0, double y0, double z0, double t0,  
		   double& meanTime, double& rmsTime,                   
		   std::string particleType = "e-");
    
    void SetMaxSteps(int max){maxSteps = max;};

    void EnableDebugging()  {debug = true;}
    void DisableDebugging() {debug = false;}

    void EnableVerbose()  {verbose = true;}
    void DisableVerbose() {verbose = false;}

  private:

    Sensor* sensor;
    Medium* medium;

    double maxStepSize;
    double intAccuracy;

    int maxSteps;

    bool usePlotting;
    int iLine;
    ViewDrift* viewer;

    bool debug;
    
    bool verbose;

    // Used to drift a particle to the edge of a boundary.
    void EndDriftLine();
    // Used to drift a particle to a wire
    void DriftToWire(double x0, double y0, double z0, int iWire);
    // Used by DriftToWire to find the distance to the wires edge
    double DistanceToWire(double x, double y, double z);
    // Used to determine the diffussion over the drift length
    double IntegrateDiffusion(const double x,  const double y,  const double z,
			      const double xe, const double ye, const double ze);
    
    // These variables store the position and radius ofa trapping wire
    double xWire, yWire, rWire;

    struct step{
      //position (initial and final)
      double xi,xf;
      double yi,yf;
      double zi,zf;
      //time (initial and final)
      double ti,tf;
      std::string status;
    };
    std::vector<step> path;

  };

}


#endif
