#ifndef G_COMPONENT_NEBEM_3D_H
#define G_COMPONENT_NEBEM_3D_H
#include "ComponentBase.hh"

namespace Garfield {

class ComponentNeBem3d : public ComponentBase {
 public:
  /// Constructor
  ComponentNeBem3d();
  /// Destructor
  ~ComponentNeBem3d() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  bool GetVoltageRange(double& vmin, double& vmax) override;

  unsigned int GetNumberOfPanels() const { return m_panels.size(); }
  bool GetPanel(unsigned int i, double& a, double& b, double& c,
                std::vector<double>& xv, std::vector<double>& yv,
                std::vector<double>& zv);

  unsigned int GetNumberOfElements() const { return m_elements.size(); }

  void Reset() override;
  void UpdatePeriodicity() override;

  /// Retrieve surface panels, remove contacts and cut polygons to rectangles
  /// and right-angle triangles.
  bool Initialise();

 private:
  /// List of surface panels.
  std::vector<Panel> m_panels;

  // Array of boundary elements (version for nebem 3d)
  struct Element {
    // Geometric type
    int geoType;
    // Center point
    double cX, cY;
    // Half length
    double len;
    // Rotation angle
    double phi;
    // Charge density
    double solution;
    // Boundary condition type and value - have to come from the solids
    //  unsigned int bcType;
    //  double bcValue;
    // Ratio of relative dielectric permittivities
    double lambda;
  };
  std::vector<Element> m_elements;

  /// Isolate the parts of polygon 1 that are not hidden by 2 and vice versa.
  bool EliminateOverlaps(const Panel& panel1, const Panel& panel2,
                         std::vector<Panel>& panelsOut,
                         std::vector<int>& itypo);

  bool TraceEnclosed(const std::vector<double>& xl1,
                     const std::vector<double>& yl1,
                     const std::vector<double>& xl2,
                     const std::vector<double>& yl2, const Panel& originalPanel,
                     std::vector<Panel>& newPanels) const;

  void TraceNonOverlap(
      const std::vector<double>& xp1, const std::vector<double>& yp1,
      const std::vector<double>& xl1, const std::vector<double>& yl1,
      const std::vector<double>& xl2, const std::vector<double>& yl2,
      const std::vector<int>& flags1, const std::vector<int>& flags2,
      const std::vector<int>& links1, const std::vector<int>& links2,
      std::vector<bool>& mark1, int ip1, const Panel& originalPanel,
      std::vector<Panel>& newPanels) const;

  void TraceOverlap(
      const std::vector<double>& xp1, const std::vector<double>& yp1,
      const std::vector<double>& xp2, const std::vector<double>& yp2,
      const std::vector<double>& xl1, const std::vector<double>& yl1,
      const std::vector<double>& xl2, const std::vector<double>& yl2,
      const std::vector<int>& flags1, const std::vector<int>& links1,
      const std::vector<int>& links2, std::vector<bool>& mark1, int ip1,
      int ip2, const Panel& originalPanel, std::vector<Panel>& newPanels) const;

  /// Split a polygon into rectangles and right-angled triangles.
  bool MakePrimitives(const Panel& panelIn,
                      std::vector<Panel>& panelsOut) const;

  /// Check whether a polygon contains parallel lines.
  /// If it does, split it in rectangular and non-rectangular parts.
  bool SplitTrapezium(const Panel panelIn, std::vector<Panel>& stack,
                      std::vector<Panel>& panelsOut, const double epsang) const;
};
}

#endif
