#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include "ComponentNeBem3d.hh"
#include "FundamentalConstants.hh"
#include "GarfieldConstants.hh"
#include "Random.hh"

namespace {

unsigned int NextPoint(const unsigned int i, const unsigned int n) {
  const unsigned int j = i + 1;
  return j < n ? j : 0;
}

unsigned int PrevPoint(const unsigned int i, const unsigned int n) {
  return i > 0 ? i - 1 : n - 1;
}

/// Compute lambda for a point on a line (0 = start, 1 = end) .
double Lambda(const double x1, const double x0, const double x2,
              const double y1, const double y0, const double y2) {
  // Segment of zero length.
  if ((x1 - x2) == 0. && (y1 - y2) == 0.) {
    std::cerr << "ComponentNeBem3d::Lambda: Zero length segment.\n";
    return 2.;
  }

  double xl = 0.;
  const double dx1 = x0 - x1;
  const double dy1 = y0 - y1;
  const double dx2 = x0 - x2;
  const double dy2 = y0 - y2;
  if (dx1 * dx1 + dy1 * dy1 < dx2 * dx2 + dy2 * dy2) {
    // Point nearer to (x1, y1).
    if (fabs(y1 - y2) > fabs(x1 - x2)) {
      xl = dy1 / (y2 - y1);
    } else {
      xl = dx1 / (x2 - x1);
    }
  } else {
    // Point nearer to (x2, y2).
    if (fabs(y1 - y2) > fabs(x1 - x2)) {
      xl = 1. - dy2 / (y1 - y2);
    } else {
      xl = 1. - dx2 / (x1 - x2);
    }
  }
  return xl;
}

/// Determine whether a point (u, v) lies on the straight lines
/// (x1, y1) to (x2, y2).
bool OnLine(const double x1, const double y1, const double x2, const double y2,
            const double u, const double v) {
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v)});
  epsx = std::max(1.e-10, epsx);
  epsy = std::max(1.e-10, epsy);

  if ((fabs(x1 - u) <= epsx && fabs(y1 - v) <= epsy) ||
      (fabs(x2 - u) <= epsx && fabs(y2 - v) <= epsy)) {
    // Point to be examined coincides with start or end.
    return true;
  } else if (fabs(x1 - x2) <= epsx && fabs(y1 - y2) <= epsy) {
    // The line (x1, y1) to (x2, y2) is in fact a point.
    return false;
  }
  double xc = 0., yc = 0.;
  if (fabs(u - x1) + fabs(v - y1) < fabs(u - x2) + fabs(v - y2)) {
    // (u, v) is nearer to (x1, y1).
    const double dx = (x2 - x1);
    const double dy = (y2 - y1);
    const double xl = ((u - x1) * dx + (v - y1) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x1;
      yc = y1;
    } else if (xl > 1.) {
      xc = x2;
      yc = y2;
    } else {
      xc = x1 + xl * dx;
      yc = y1 + xl * dy;
    }
  } else {
    // (u, v) is nearer to (x2, y2).
    const double dx = (x1 - x2);
    const double dy = (y1 - y2);
    const double xl = ((u - x2) * dx + (v - y2) * dy) / (dx * dx + dy * dy);
    if (xl < 0.) {
      xc = x2;
      yc = y2;
    } else if (xl > 1.) {
      xc = x1;
      yc = y1;
    } else {
      xc = x2 + xl * dx;
      yc = y2 + xl * dy;
    }
  }
  // See whether the point is on the line.
  if (fabs(u - xc) < epsx && fabs(v - yc) < epsy) {
    return true;
  }
  return false;
}

/// Determine whether the 2 straight lines (x1, y1) to (x2, y2)
/// and (u1, u2) to (v1, v2) cross at an intermediate point for both lines.
bool Crossing(const double x1, const double y1, const double x2,
              const double y2, const double u1, const double v1,
              const double u2, const double v2, double& xc, double& yc) {
  /// Matrix to compute the crossing point.
  std::array<std::array<double, 2>, 2> a;
  a[0][0] = y2 - y1;
  a[0][1] = v2 - v1;
  a[1][0] = x1 - x2;
  a[1][1] = u1 - u2;
  const double det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  // Initial values.
  xc = 0.;
  yc = 0.;
  // Set tolerances.
  double epsx = 1.e-10 * std::max({fabs(x1), fabs(x2), fabs(u1), fabs(u2)});
  double epsy = 1.e-10 * std::max({fabs(y1), fabs(y2), fabs(v1), fabs(v2)});
  epsx = std::max(epsx, 1.e-10);
  epsy = std::max(epsy, 1.e-10);
  // Check for a point of one line located on the other line.
  if (OnLine(x1, y1, x2, y2, u1, v1)) {
    xc = u1;
    yc = v1;
    return true;
  } else if (OnLine(x1, y1, x2, y2, u2, v2)) {
    xc = u2;
    yc = v2;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x1, y1)) {
    xc = x1;
    yc = y1;
    return true;
  } else if (OnLine(u1, v1, u2, v2, x2, y2)) {
    xc = x2;
    yc = y2;
    return true;
  } else if (fabs(det) < epsx * epsy) {
    // Parallel, non-touching.
    return false;
  }
  // Crossing, non-trivial lines: solve crossing equations.
  const double aux = a[1][1];
  a[1][1] = a[0][0] / det;
  a[0][0] = aux / det;
  a[1][0] = -a[1][0] / det;
  a[0][1] = -a[0][1] / det;
  // Compute crossing point.
  xc = a[0][0] * (x1 * y2 - x2 * y1) + a[1][0] * (u1 * v2 - u2 * v1);
  yc = a[0][1] * (x1 * y2 - x2 * y1) + a[1][1] * (u1 * v2 - u2 * v1);
  // See whether the crossing point is on both lines.
  if (OnLine(x1, y1, x2, y2, xc, yc) && OnLine(u1, v1, u2, v2, xc, yc)) {
    // Intersecting lines.
    return true;
  }
  // Crossing point not on both lines.
  return false;
}

void AddPoints(const std::vector<double>& xp1, const std::vector<double>& yp1,
               const std::vector<double>& xp2, const std::vector<double>& yp2,
               std::vector<double>& xl, std::vector<double>& yl,
               std::vector<int>& flags, std::vector<double>& qs,
               const double epsx, const double epsy) {
  struct Point {
    double x;
    double y;
    int flag;
    double q;
  };

  std::vector<Point> points;

  const unsigned int np1 = xp1.size();
  const unsigned int np2 = xp2.size();
  for (unsigned int i = 0; i < np1; ++i) {
    const double xi0 = xp1[i];
    const double yi0 = yp1[i];
    const double xi1 = xp1[NextPoint(i, np1)];
    const double yi1 = yp1[NextPoint(i, np1)];
    // Add the vertex.
    Point p1;
    p1.x = xi0;
    p1.y = yi0;
    p1.flag = 1;
    p1.q = 0.;
    // If also on 2 or vertex of 2, flag it as crossing or foreign.
    for (unsigned int j = 0; j < np2; ++j) {
      const double xj0 = xp2[j];
      const double yj0 = yp2[j];
      if (fabs(xj0 - xi0) < epsx && fabs(yj0 - yi0) < epsy) {
        p1.flag = 2;
      }
      const double xj1 = xp2[NextPoint(j, np2)];
      const double yj1 = yp2[NextPoint(j, np2)];
      if (OnLine(xj0, yj0, xj1, yj1, xi0, yi0) &&
          (fabs(xj0 - xi0) > epsx || fabs(yj0 - yi0) > epsy) &&
          (fabs(xj1 - xi0) > epsx || fabs(yj1 - yi0) > epsy)) {
        p1.flag = 3;
      }
    }
    points.push_back(std::move(p1));
    // Go over the line segments of the other polygon.
    std::vector<Point> pointsOther;
    for (unsigned int j = 0; j < np2; ++j) {
      const double xj0 = xp2[j];
      const double yj0 = yp2[j];
      // Add vertices of 2 that are on this line.
      if (OnLine(xi0, yi0, xi1, yi1, xj0, yj0) &&
          (fabs(xi0 - xj0) > epsx || fabs(yi0 - yj0) > epsy) &&
          (fabs(xi1 - xj0) > epsx || fabs(yi1 - yj0) > epsy)) {
        Point p2;
        p2.x = xj0;
        p2.y = yj0;
        p2.flag = 2;
        pointsOther.push_back(std::move(p2));
      }
      const double xj1 = xp2[NextPoint(j, np2)];
      const double yj1 = yp2[NextPoint(j, np2)];
      // Add crossing points.
      double xc = 0., yc = 0.;
      bool add = Crossing(xi0, yi0, xi1, yi1, xj0, yj0, xj1, yj1, xc, yc);
      if (add) {
        if ((fabs(xi0 - xc) < epsx && fabs(yi0 - yc) < epsy) ||
            (fabs(xi1 - xc) < epsx && fabs(yi1 - yc) < epsy) ||
            (fabs(xj0 - xc) < epsx && fabs(yj0 - yc) < epsy) ||
            (fabs(xj1 - xc) < epsx && fabs(yj1 - yc) < epsy)) {
          add = false;
        }
        if ((fabs(xi0 - xj0) < epsx && fabs(yi0 - yj0) < epsy) ||
            (fabs(xi0 - xj1) < epsx && fabs(yi0 - yj1) < epsy) ||
            (fabs(xi1 - xj0) < epsx && fabs(yi1 - yj0) < epsy) ||
            (fabs(xi1 - xj1) < epsx && fabs(yi1 - yj1) < epsy)) {
          add = false;
        }
      }
      if (add) {
        Point p2;
        p2.x = xc;
        p2.y = yc;
        p2.flag = 3;
        pointsOther.push_back(std::move(p2));
      }
    }
    // Compute the lambdas for these points.
    for (auto& p : pointsOther) {
      p.q = Lambda(xi0, p.x, xi1, yi0, p.y, yi1);
    }
    // Sort the list by using the lambdas.
    std::sort(
        pointsOther.begin(), pointsOther.end(),
        [](const Point& lhs, const Point& rhs) { return (lhs.q < rhs.q); });
    points.insert(points.end(), pointsOther.begin(), pointsOther.end());
  }

  for (const Point& p : points) {
    xl.push_back(p.x);
    yl.push_back(p.y);
    flags.push_back(p.flag);
    qs.push_back(p.q);
  }
}

/// Determine whether the point (x, y) is located inside of the
/// polygon (xpl, ypl).
void Inside(const std::vector<double>& xpl, const std::vector<double>& ypl,
            const double x, const double y, bool& inside, bool& edge) {
  // Initial settings.
  inside = false;
  edge = false;
  const unsigned int npl = xpl.size();
  if (ypl.size() != npl) return;
  // Special treatment for few points.
  if (npl < 2) {
    return;
  } else if (npl == 2) {
    edge = OnLine(xpl[0], ypl[0], xpl[1], ypl[1], x, y);
    return;
  }
  // Determine the range of the data.
  const double xmin = *std::min_element(std::begin(xpl), std::end(xpl));
  const double xmax = *std::max_element(std::begin(xpl), std::end(xpl));
  const double ymin = *std::min_element(std::begin(ypl), std::end(ypl));
  const double ymax = *std::max_element(std::begin(ypl), std::end(ypl));

  // Set tolerances.
  double epsx = 1.e-8 * std::max(fabs(xmin), fabs(xmax));
  double epsy = 1.e-8 * std::max(fabs(ymin), fabs(ymax));
  epsx = std::max(epsx, 1.e-8);
  epsy = std::max(epsy, 1.e-8);

  // Ensure that we have a range.
  if (fabs(xmax - xmin) <= epsx) {
    if (y >= ymin - epsy && y <= ymax + epsy &&
        fabs(xmax + xmin - 2 * x) <= epsx) {
      edge = true;
    } else {
      edge = false;
    }
  } else if (fabs(ymax - ymin) <= epsy) {
    if (x >= xmin - epsx && x <= xmax + epsx &&
        fabs(ymax + ymin - 2 * y) <= epsy) {
      edge = true;
    } else {
      edge = false;
    }
  }
  // Choose a point at "infinity".
  double xinf = xmin - fabs(xmax - xmin);
  double yinf = ymin - fabs(ymax - ymin);

  unsigned int nIter = 0;
  bool ok = false;
  while (!ok && nIter < 100) {
    ok = true;
    // Loop over the edges counting intersections.
    unsigned int nCross = 0;
    for (unsigned int j = 0; j < npl; ++j) {
      const unsigned int jj = NextPoint(j, npl);
      // Flag points located on one of the edges.
      if (OnLine(xpl[j], ypl[j], xpl[jj], ypl[jj], x, y)) {
        edge = true;
        return;
      }
      // Count mid-line intersects.
      double xc = 0., yc = 0.;
      if (Crossing(x, y, xinf, yinf, xpl[j], ypl[j], xpl[jj], ypl[jj], xc,
                   yc)) {
        ++nCross;
      }
      // Ensure that the testing line doesn't cross a corner.
      if (OnLine(x, y, xinf, yinf, xpl[j], ypl[j])) {
        xinf = xmin - Garfield::RndmUniform() * fabs(xmax - xinf);
        yinf = ymin - Garfield::RndmUniform() * fabs(ymax - yinf);
        ok = false;
        break;
      }
    }
    if (ok) {
      // Set the INSIDE flag.
      if (nCross != 2 * (nCross / 2)) inside = true;
      return;
    }
    ++nIter;
  }

  std::cerr << "ComponentNeBem3d::Inside:\n    Warning. Unable to verify "
            << "whether a point is internal; setting to edge.\n";
  inside = false;
  edge = true;
}

/// Determine whether 2 panels are equal.
bool Equal(const Garfield::Panel& panel1, const Garfield::Panel& panel2,
           const double epsx, const double epsy) {
  const auto& xp1 = panel1.xv;
  const auto& yp1 = panel1.yv;
  const auto& xp2 = panel2.xv;
  const auto& yp2 = panel2.yv;
  if (xp1.empty() || xp2.empty()) return false;
  const unsigned int np1 = xp1.size();
  const unsigned int np2 = xp2.size();

  // Compare all points of 1 with all points of 2.
  for (unsigned int i = 0; i < np1; ++i) {
    // Loop over 2 until a match is found.
    bool match = false;
    for (unsigned int j = 0; j < np2; ++j) {
      if (fabs(xp2[j] - xp1[i]) < epsx && fabs(yp2[j] - yp1[i]) < epsy) {
        match = true;
        break;
      }
      const unsigned int jj = NextPoint(j, np2);
      if (OnLine(xp2[j], yp2[j], xp2[jj], yp2[jj], xp1[i], yp1[i])) {
        match = true;
        break;
      }
    }
    if (!match) return false;
  }

  // Compare all points of 2 with all points of 1.
  for (unsigned int i = 0; i < np2; ++i) {
    // Loop over 1 until a match is found.
    bool match = false;
    for (unsigned int j = 0; j < np1; ++j) {
      if (fabs(xp2[i] - xp1[j]) < epsx && fabs(yp2[i] - yp1[j]) < epsy) {
        match = true;
        break;
      }
      const unsigned int jj = NextPoint(j, np1);
      if (OnLine(xp1[j], yp1[j], xp1[jj], yp1[jj], xp2[i], yp2[i])) {
        match = true;
        break;
      }
    }
    if (!match) return false;
  }

  // If we get this far, the curves are the same.
  return true;
}
}

namespace Garfield {

ComponentNeBem3d::ComponentNeBem3d() : ComponentBase() {
  m_className = "ComponentNeBem3d";
}

void ComponentNeBem3d::ElectricField(const double x, const double y,
                                     const double z, double& ex, double& ey,
                                     double& ez, double& v, Medium*& m,
                                     int& status) {
  ex = ey = ez = v = 0.;
  status = 0;
  // Check if the requested point is inside a medium
  m = GetMedium(x, y, z);
  if (!m) {
    status = -6;
    return;
  }

  if (!m_ready) {
    if (!Initialise()) {
      std::cerr << m_className << "::ElectricField: Initialisation failed.\n";
      status = -11;
      return;
    }
    m_ready = true;
  }
}

void ComponentNeBem3d::ElectricField(const double x, const double y,
                                     const double z, double& ex, double& ey,
                                     double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}

bool ComponentNeBem3d::GetVoltageRange(double& vmin, double& vmax) {
  // Voltage and other bc have to come from the solids
  vmin = vmax = 0;
  return true;
}

bool ComponentNeBem3d::Initialise() {
  // Reset the panel list.
  m_panels.clear();

  if (!m_geometry) {
    std::cerr << m_className << "::Initialise: Geometry not set.\n";
    return false;
  }
  // Be sure we won't have intersections with the bounding box.
  // TODO! Loop over the solids and call PLACYE, PLACHE, PLABXE, PLASPE, ...
  // Loop over the solids.
  std::map<int, Solid::BoundaryCondition> bc;
  std::map<int, double> volt;
  std::map<int, double> eps;
  const unsigned int nSolids = m_geometry->GetNumberOfSolids();
  for (unsigned int i = 0; i < nSolids; ++i) {
    Medium* medium = nullptr;
    const auto solid = m_geometry->GetSolid(i, medium);
    if (!solid) continue;
    // Get the panels.
    solid->SolidPanels(m_panels);
    // Get the boundary condition.
    const auto id = solid->GetId();
    bc[id] = solid->GetBoundaryConditionType();
    volt[id] = solid->GetBoundaryPotential();
    if (!medium) {
      eps[id] = 1.;
    } else {
      eps[id] = medium->GetDielectricConstant();
    }
  }
  // Apply cuts.
  // CALL CELSCT('APPLY')
  // Reduce to basic periodic copy.
  // CALL BEMBAS

  // Find contact panels and split into primitives.

  // *---------------------------------------------------------------------
  // * PLABEM - Prepares panels for BEM applications: removes the contacts
  // *          and cuts polygons to rectangles and right-angle triangles.
  // *---------------------------------------------------------------------

  // Establish tolerances.
  const double epsang = 1.e-6;  // BEMEPA
  const double epsxyz = 1.e-6;  // BEMEPD
  // CALL EPSSET('SET',EPSXYZ,EPSXYZ,EPSXYZ)
  const unsigned int nRef = m_panels.size();

  std::vector<Panel> oldPanels;
  oldPanels.swap(m_panels);
  const unsigned int nOld = oldPanels.size();
  // Keep track of which panels have been processed.
  std::vector<bool> mark(nOld, false);
  // Pick up panels which coincide potentially.
  for (unsigned int i = 0; i < nOld; ++i) {
    // Skip panels already done.
    if (mark[i]) continue;
    // Fetch panel parameters.
    const double a1 = oldPanels[i].a;
    const double b1 = oldPanels[i].b;
    const double c1 = oldPanels[i].c;
    const auto& xp1 = oldPanels[i].xv;
    const auto& yp1 = oldPanels[i].yv;
    const auto& zp1 = oldPanels[i].zv;
    const unsigned int np1 = xp1.size();
    // Establish its norm and offset.
    const double d1 = a1 * xp1[0] + b1 * yp1[0] + c1 * zp1[0];
    if (m_debug) {
      std::cout << "Panel " << i << "\n"
                << "  Norm vector: " << a1 << ", " << b1 << ", " << c1 << ", "
                << d1 << "\n";
    }
    // Rotation matrix.
    std::array<std::array<double, 3>, 3> rot;
    if (fabs(c1) <= fabs(a1) && fabs(c1) <= fabs(b1)) {
      // Rotation: removing C
      rot[0][0] = b1 / sqrt(a1 * a1 + b1 * b1);
      rot[0][1] = -a1 / sqrt(a1 * a1 + b1 * b1);
      rot[0][2] = 0.;
    } else if (fabs(b1) <= fabs(a1) && fabs(b1) <= fabs(c1)) {
      // Rotation: removing B
      rot[0][0] = c1 / sqrt(a1 * a1 + c1 * c1);
      rot[0][1] = 0.;
      rot[0][2] = -a1 / sqrt(a1 * a1 + c1 * c1);
    } else {
      // Rotation: removing A
      rot[0][0] = 0.;
      rot[0][1] = c1 / sqrt(b1 * b1 + c1 * c1);
      rot[0][2] = -b1 / sqrt(b1 * b1 + c1 * c1);
    }
    rot[2][0] = a1;
    rot[2][1] = b1;
    rot[2][2] = c1;
    rot[1][0] = rot[2][1] * rot[0][2] - rot[2][2] * rot[0][1];
    rot[1][1] = rot[2][2] * rot[0][0] - rot[2][0] * rot[0][2];
    rot[1][2] = rot[2][0] * rot[0][1] - rot[2][1] * rot[0][0];
    // Rotate to the x, y plane.
    std::vector<double> xp(np1, 0.);
    std::vector<double> yp(np1, 0.);
    std::vector<double> zp(np1, 0.);
    double zm = 0.;
    for (unsigned int k = 0; k < np1; ++k) {
      xp[k] = rot[0][0] * xp1[k] + rot[0][1] * yp1[k] + rot[0][2] * zp1[k];
      yp[k] = rot[1][0] * xp1[k] + rot[1][1] * yp1[k] + rot[1][2] * zp1[k];
      zp[k] = rot[2][0] * xp1[k] + rot[2][1] * yp1[k] + rot[2][2] * zp1[k];
      zm += zp[k];
    }
    zm /= np1;
    // Store it.
    std::vector<Panel> newPanels;
    std::vector<int> vol1;
    std::vector<int> vol2;
    Panel panel1 = oldPanels[i];
    panel1.xv = xp;
    panel1.yv = yp;
    panel1.zv = zp;
    vol1.push_back(panel1.volume);
    vol2.push_back(-1);
    newPanels.push_back(std::move(panel1));
    // Pick up all matching planes.
    for (unsigned int j = i + 1; j < nRef; ++j) {
      if (mark[j]) continue;
      const double a2 = oldPanels[j].a;
      const double b2 = oldPanels[j].b;
      const double c2 = oldPanels[j].c;
      const auto& xp2 = oldPanels[j].xv;
      const auto& yp2 = oldPanels[j].yv;
      const auto& zp2 = oldPanels[j].zv;
      const unsigned int np2 = xp2.size();
      // See whether this matches the first.
      const double d2 = a2 * xp2[0] + b2 * yp2[0] + c2 * zp2[0];
      // Inner product.
      const double dot = a1 * a2 + b1 * b2 + c1 * c2;
      // Offset between the two planes.
      const double offset = d1 - d2 * dot;
      if (fabs(fabs(dot) - 1.) > epsang || fabs(offset) > epsxyz) continue;
      // Found a match.
      mark[j] = true;
      if (m_debug) std::cout << "Match with panel " << j << "\n";
      // Rotate this plane too.
      xp.assign(np2, 0.);
      yp.assign(np2, 0.);
      zp.assign(np2, 0.);
      zm = 0.;
      for (unsigned int k = 0; k < np2; ++k) {
        xp[k] = rot[0][0] * xp2[k] + rot[0][1] * yp2[k] + rot[0][2] * zp2[k];
        yp[k] = rot[1][0] * xp2[k] + rot[1][1] * yp2[k] + rot[1][2] * zp2[k];
        zp[k] = rot[2][0] * xp2[k] + rot[2][1] * yp2[k] + rot[2][2] * zp2[k];
        zm += zp[k];
      }
      zm /= np2;
      // Store it.
      Panel panel2 = oldPanels[j];
      panel2.xv = xp;
      panel2.yv = yp;
      panel2.zv = zp;
      vol1.push_back(panel2.volume);
      vol2.push_back(-1);
      newPanels.push_back(std::move(panel2));
    }
    std::vector<bool> obsolete(newPanels.size(), false);
    // Cut them as long as needed till no contacts remain.
    unsigned int jmin = 0;
    bool change = true;
    while (change) {
      change = false;
      const unsigned int n = newPanels.size();
      for (unsigned int j = 0; j < n; ++j) {
        if (obsolete[j] || j < jmin) continue;
        if (vol1[j] >= 0 && vol2[j] >= 0) continue;
        const auto& panelj = newPanels[j];
        for (unsigned int k = j + 1; k < n; ++k) {
          if (obsolete[k]) continue;
          if (vol1[k] >= 0 && vol2[k] >= 0) continue;
          const auto& panelk = newPanels[k];
          if (m_debug) std::cout << "Cutting " << j << ", " << k << ".\n";
          // Separate contact and non-contact areas.
          std::vector<Panel> panelsOut;
          std::vector<int> itypo;
          EliminateOverlaps(panelj, panelk, panelsOut, itypo);
          const unsigned int nOut = panelsOut.size();
          if (nOut == 2) {
            // TODO: retrieve epsx, epsy from overlap finding?
            const double epsx = epsxyz;
            const double epsy = epsxyz;
            // If there are just 2 panels, see whether there is a new one.
            const bool equal1 = Equal(panelj, panelsOut[0], epsx, epsy);
            const bool equal2 = Equal(panelj, panelsOut[1], epsx, epsy);
            const bool equal3 = Equal(panelk, panelsOut[0], epsx, epsy);
            const bool equal4 = Equal(panelk, panelsOut[1], epsx, epsy);
            if ((equal1 || equal3) && (equal2 || equal4)) {
              if (m_debug) {
                std::cout << "Original and new panels are identical.\n";
              }
            } else {
              change = true;
            }
          } else {
            change = true;
          }
          if (m_debug) std::cout << "Change flag: " << change << "\n";
          // If there is no change, keep the two panels and proceed.
          if (!change) continue;
          // Flag the existing panels as inactive.
          obsolete[j] = true;
          obsolete[k] = true;

          // Add the new panels.
          for (unsigned int l = 0; l < nOut; ++l) {
            if (itypo[l] == 1) {
              vol1.push_back(std::max(vol1[j], vol2[j]));
              vol2.push_back(-1);
            } else if (itypo[l] == 2) {
              vol1.push_back(std::max(vol1[k], vol2[k]));
              vol2.push_back(-1);
            } else {
              vol1.push_back(std::max(vol1[j], vol2[j]));
              vol2.push_back(std::max(vol1[k], vol2[k]));
            }
            newPanels.push_back(std::move(panelsOut[l]));
            obsolete.push_back(false);
          }
          jmin = j + 1;
          // Restart the loops.
          break;
        }
        if (change) break;
      }
    }
    // And rotate the panels back in place.
    const unsigned int nNew = newPanels.size();
    for (unsigned int j = 0; j < nNew; ++j) {
      if (obsolete[j]) continue;
      std::vector<Panel> panelsOut;
      // Reduce to rectangles and right-angle triangles.
      MakePrimitives(newPanels[j], panelsOut);
      // Loop over the rectangles and triangles.
      for (auto& panel : panelsOut) {
        const auto& up = panel.xv;
        const auto& vp = panel.yv;
        const auto& wp = panel.zv;
        const unsigned int np = up.size();
        // Rotate.
        xp.assign(np, 0.);
        yp.assign(np, 0.);
        zp.assign(np, 0.);
        for (unsigned int k = 0; k < np; ++k) {
          xp[k] = rot[0][0] * up[k] + rot[1][0] * vp[k] + rot[2][0] * wp[k];
          yp[k] = rot[0][1] * up[k] + rot[1][1] * vp[k] + rot[2][1] * wp[k];
          zp[k] = rot[0][2] * up[k] + rot[1][2] * vp[k] + rot[2][2] * wp[k];
        }
        panel.xv = xp;
        panel.yv = yp;
        panel.zv = zp;
        m_panels.push_back(std::move(panel));
      }
    }
  }
  return true;
}

bool ComponentNeBem3d::EliminateOverlaps(const Panel& panel1,
                                         const Panel& panel2,
                                         std::vector<Panel>& panelsOut,
                                         std::vector<int>& itypo) {
  // *-----------------------------------------------------------------------
  // *   PLAOVL - Isolates the parts of plane 1 that are not hidden by 2.
  // *-----------------------------------------------------------------------

  const auto& xp1 = panel1.xv;
  const auto& yp1 = panel1.yv;
  const auto& zp1 = panel1.zv;
  const auto& xp2 = panel2.xv;
  const auto& yp2 = panel2.yv;
  const auto& zp2 = panel2.zv;
  // If the size of either is less than 3, simply return.
  if (xp1.size() <= 2 || xp2.size() <= 2) {
    return true;
  }
  // Compute the various tolerances.
  const double xmin1 = *std::min_element(std::begin(xp1), std::end(xp1));
  const double ymin1 = *std::min_element(std::begin(yp1), std::end(yp1));
  const double xmax1 = *std::max_element(std::begin(xp1), std::end(xp1));
  const double ymax1 = *std::max_element(std::begin(yp1), std::end(yp1));

  const double xmin2 = *std::min_element(std::begin(xp2), std::end(xp2));
  const double ymin2 = *std::min_element(std::begin(yp2), std::end(yp2));
  const double xmax2 = *std::max_element(std::begin(xp2), std::end(xp2));
  const double ymax2 = *std::max_element(std::begin(yp2), std::end(yp2));

  const double xmin = std::min(xmin1, xmin2);
  const double ymin = std::min(ymin1, ymin2);
  const double xmax = std::max(xmax1, xmax2);
  const double ymax = std::max(ymax1, ymax2);

  const double epsx = 1.e-6 * std::max(std::abs(xmax), std::abs(xmin));
  const double epsy = 1.e-6 * std::max(std::abs(ymax), std::abs(ymin));

  const double zsum1 = std::accumulate(std::begin(zp1), std::end(zp1), 0.);
  const double zsum2 = std::accumulate(std::begin(zp2), std::end(zp2), 0.);
  const double zmean = (zsum1 + zsum2) / (zp1.size() + zp2.size());

  std::array<std::vector<double>, 2> xl;
  std::array<std::vector<double>, 2> yl;
  std::array<std::vector<int>, 2> flags;
  std::array<std::vector<double>, 2> qs;
  // Establish the list of special points around polygon 1.
  AddPoints(xp1, yp1, xp2, yp2, xl[0], yl[0], flags[0], qs[0], epsx, epsy);
  // Establish the list of special points around polygon 2.
  AddPoints(xp2, yp2, xp1, yp1, xl[1], yl[1], flags[1], qs[1], epsx, epsy);

  bool ok = true;
  // Look up the cross-links: from plane 1 (2) to plane 2 (1).
  std::array<std::vector<int>, 2> links;
  for (unsigned int ic = 0; ic < 2; ++ic) {
    const unsigned int n1 = xl[ic].size();
    links[ic].assign(n1, -1);
    const unsigned int jc = ic == 0 ? 1 : 0;
    const unsigned int n2 = xl[jc].size();
    for (unsigned int i = 0; i < n1; ++i) {
      unsigned int nFound = 0;
      for (unsigned int j = 0; j < n2; ++j) {
        if (fabs(xl[ic][i] - xl[jc][j]) < epsx &&
            fabs(yl[ic][i] - yl[jc][j]) < epsy) {
          ++nFound;
          links[ic][i] = j;
        }
      }
      if (nFound == 0 && (flags[ic][i] == 2 || flags[ic][i] == 3)) {
        std::cerr << m_className << "::EliminateOverlaps: "
                  << "Warning. Expected match not found (" << ic + 1 << "-"
                  << jc + 1 << ").\n";
        links[ic][i] = -1;
        ok = false;
      } else if (nFound > 1) {
        std::cerr << m_className << "::EliminateOverlaps: "
                  << "Warning. More than 1 match found (" << ic + 1 << "-"
                  << jc + 1 << ").\n";
        links[ic][i] = -1;
        ok = false;
      }
    }
  }

  // List the points for debugging.
  if (m_debug) {
    for (unsigned int j = 0; j < 2; ++j) {
      std::cout << "Polygon " << j << "\n";
      std::cout << " No Type            x            y        Q   links\n";
      const unsigned int n = xl[j].size();
      for (unsigned int i = 0; i < n; ++i) {
        printf("  %3d %5d %13.6f %13.6f %5.3f %3d\n", i, flags[j][i], xl[j][i],
               yl[j][i], qs[j][i], links[j][i]);
      }
    }
  }
  if (!ok) return false;

  for (unsigned int ic = 0; ic < 2; ++ic) {
    // See whether all of 1 (2) is inside 2 (1).
    bool allInside = true;
    const unsigned int np = xl[ic].size();
    for (unsigned int i = 0; i < np; ++i) {
      if (flags[ic][i] != 1) {
        allInside = false;
        break;
      }
      bool inside = false, edge = false;
      if (ic == 0) {
        Inside(xp2, yp2, xl[ic][i], yl[ic][i], inside, edge);
      } else {
        Inside(xp1, yp1, xl[ic][i], yl[ic][i], inside, edge);
      }
      if (!(inside || edge)) {
        allInside = false;
        break;
      }
    }
    if (allInside) {
      // Apparently 1 (2) really is fully inside 2 (1).
      if (ic == 0) {
        if (m_debug) std::cout << "Curve 1 fully inside 2.\n";
        // Write out curve 1.
        panelsOut.push_back(panel1);
      } else {
        if (m_debug) std::cout << "Curve 2 fully inside 1.\n";
        // Write out curve 2.
        panelsOut.push_back(panel2);
      }
      panelsOut.back().zv.assign(panelsOut.back().xv.size(), zmean);
      itypo.push_back(3);
      std::vector<Panel> newPanels;
      if (ic == 0) {
        if (!TraceEnclosed(xl[0], yl[0], xl[1], yl[1], panel2, newPanels)) {
          return false;
        }
      } else {
        if (!TraceEnclosed(xl[1], yl[1], xl[0], yl[0], panel1, newPanels)) {
          return false;
        }
      }
      for (auto& panel : newPanels) {
        panel.zv.assign(panel.xv.size(), zmean);
        if (ic == 0) {
          itypo.push_back(2);
        } else {
          itypo.push_back(1);
        }
      }
      panelsOut.insert(panelsOut.end(), newPanels.begin(), newPanels.end());
      return true;
    }
  }

  for (unsigned int ic = 0; ic < 2; ++ic) {
    std::vector<Panel> newPanels;
    const unsigned int n = xl[ic].size();
    // Identify the parts of 1 (2) that are not overlapped, first mark.
    std::vector<bool> mark(n, false);
    bool done = false;
    while (!done) {
      if (m_debug) {
        std::cout << "Searching for starting point on " << ic + 1 << ".\n";
      }
      done = true;
      // Try and find a new starting point
      for (unsigned int i = 0; i < n; ++i) {
        const unsigned int ii = NextPoint(i, n);
        // Skip parts already processed.
        if (mark[i] || mark[ii]) continue;
        // Skip if mid point is inside other volume.
        bool inside = false, edge = false;
        const double xm = 0.5 * (xl[ic][i] + xl[ic][ii]);
        const double ym = 0.5 * (yl[ic][i] + yl[ic][ii]);
        if (ic == 0) {
          Inside(xp2, yp2, xm, ym, inside, edge);
        } else {
          Inside(xp1, yp1, xm, ym, inside, edge);
        }
        if (inside || edge) continue;
        // Found one.
        done = false;

        if (ic == 0) {
          // Trace this part of 1 outside 2.
          TraceNonOverlap(xp1, yp1, xl[0], yl[0], xl[1], yl[1], flags[0],
                          flags[1], links[0], links[1], mark, i, panel1,
                          newPanels);
        } else {
          // Trace this part of 2 outside 1.
          TraceNonOverlap(xp2, yp2, xl[1], yl[1], xl[0], yl[0], flags[1],
                          flags[0], links[1], links[0], mark, i, panel2,
                          newPanels);
        }
        break;
      }
    }
    for (auto& panel : newPanels) {
      panel.zv.assign(panel.xv.size(), zmean);
      itypo.push_back(ic + 1);
    }
    panelsOut.insert(panelsOut.end(), newPanels.begin(), newPanels.end());
    if (m_debug) {
      std::cout << "No further non-overlapped areas of " << ic + 1 << ".\n";
    }
  }

  // Look for the overlapped parts.
  std::vector<Panel> newPanels;
  const unsigned int n1 = xl[0].size();
  std::vector<bool> mark1(n1, false);
  bool done = false;
  while (!done) {
    done = true;
    if (m_debug) std::cout << "Searching for starting point on overlap.\n";
    for (unsigned int i = 0; i < n1; ++i) {
      // Skip points already processed.
      if (mark1[i]) continue;
      // Skip if not an edge point on both 1 and 2 or internal in 2.
      int ip1 = i;
      int ip2 = links[0][ip1];
      if (ip2 < 0 || flags[0][ip1] == 1) {
        bool inside = false, edge = false;
        Inside(xp2, yp2, xl[0][ip1], yl[0][ip1], inside, edge);
        if (!(inside || edge)) continue;
      } else if (flags[1][ip2] == 1) {
        continue;
      }
      // Found one.
      done = false;
      TraceOverlap(xp1, yp1, xp2, yp2, xl[0], yl[0], xl[1], yl[1], flags[0],
                   links[0], links[1], mark1, ip1, ip2, panel1, newPanels);
      break;
    }
  }
  for (auto& panel : newPanels) {
    panel.zv.assign(panel.xv.size(), zmean);
    itypo.push_back(3);
  }
  panelsOut.insert(panelsOut.end(), newPanels.begin(), newPanels.end());
  // Finished
  if (m_debug) std::cout << "No further overlapped areas.\n";
  return true;
}

bool ComponentNeBem3d::TraceEnclosed(const std::vector<double>& xl1,
                                     const std::vector<double>& yl1,
                                     const std::vector<double>& xl2,
                                     const std::vector<double>& yl2,
                                     const Panel& panel2,
                                     std::vector<Panel>& panelsOut) const {
  const int n1 = xl1.size();
  const int n2 = xl2.size();
  // Find 2 non-crossing connections: JP1-JP2 and KP1-KP2.
  unsigned int nFound = 0;
  int jp1 = 0, jp2 = 0;
  int kp1 = 0, kp2 = 0;
  for (int ip1 = 0; ip1 < n1; ++ip1) {
    const double x1 = xl1[ip1];
    const double y1 = yl1[ip1];
    for (int ip2 = 0; ip2 < n2; ++ip2) {
      if (nFound > 0 && ip2 == jp2) continue;
      const double x2 = xl2[ip2];
      const double y2 = yl2[ip2];
      bool cross = false;
      for (int k = 0; k < n1; ++k) {
        const int kk = NextPoint(k, n1);
        if (k == ip1 || kk == ip1) continue;
        double xc = 0., yc = 0.;
        cross =
            Crossing(x1, y1, x2, y2, xl1[k], yl1[k], xl1[kk], yl1[kk], xc, yc);
        if (cross) break;
      }
      if (cross) continue;
      if (m_debug) std::cout << "No crossing with 1.\n";
      for (int k = 0; k < n2; ++k) {
        const int kk = NextPoint(k, n2);
        if (k == ip2 || kk == ip2) continue;
        double xc = 0., yc = 0.;
        cross =
            Crossing(x1, y1, x2, y2, xl2[k], yl2[k], xl2[kk], yl2[kk], xc, yc);
        if (cross) break;
      }
      if (cross) continue;
      if (nFound == 0) {
        jp1 = ip1;
        jp2 = ip2;
        if (m_debug) {
          std::cout << "First junction: " << jp1 << ", " << jp2 << ".\n";
        }
        ++nFound;
        break;
      } else {
        kp1 = ip1;
        kp2 = ip2;
        double xc = 0., yc = 0.;
        cross = Crossing(x1, y1, x2, y2, xl1[jp1], yl1[jp1], xl2[jp2], yl2[jp2],
                         xc, yc);
        if (!cross) {
          if (m_debug) {
            std::cout << "Second junction: " << kp1 << ", " << kp2 << ".\n";
          }
          ++nFound;
          break;
        }
      }
    }
    if (nFound > 1) break;
  }
  if (nFound < 2) {
    std::cerr << m_className << "::TraceEnclosed: Found no cut-out.\n";
    return false;
  }

  // Create part 1 of area 2.
  std::vector<double> xpl;
  std::vector<double> ypl;
  if (m_debug) std::cout << "Creating part 1 of area 2.\n";
  for (int ip1 = jp1; ip1 <= kp1; ++ip1) {
    if (m_debug) std::cout << "Adding " << ip1 << " on 1.\n";
    xpl.push_back(xl1[ip1]);
    ypl.push_back(yl1[ip1]);
  }
  // Try one way.
  int imax = jp2 < kp2 ? jp2 + n2 : jp2;
  int dir = +1;
  for (int i = kp2; i <= imax; ++i) {
    int ip2 = i % n2;
    if (m_debug) std::cout << "Adding " << ip2 << " on 2.\n";
    xpl.push_back(xl2[ip2]);
    ypl.push_back(yl2[ip2]);
  }
  // Check for undesirable crossings.
  bool ok = true;
  for (int ip1 = 0; ip1 < n1; ++ip1) {
    if (ip1 == jp1 || ip1 == kp1) continue;
    bool inside = false, edge = false;
    Inside(xpl, ypl, xl1[ip1], yl1[ip1], inside, edge);
    if (inside) {
      ok = false;
      break;
    }
  }
  if (!ok) {
    // Use the other way if this failed
    if (m_debug) std::cout << "Trying the other direction.\n";
    xpl.resize(kp1 - jp1 + 1);
    ypl.resize(kp1 - jp1 + 1);
    imax = jp2 < kp2 ? kp2 : kp2 + n2;
    dir = -1;
    for (int i = imax; i >= jp2; --i) {
      const int ip2 = i % n2;
      if (m_debug) std::cout << "Adding " << ip2 << " on 2.\n";
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
    }
  }

  // Save this part.
  Panel newPanel1 = panel2;
  newPanel1.xv = xpl;
  newPanel1.yv = ypl;
  panelsOut.push_back(std::move(newPanel1));
  if (m_debug) std::cout << "Part 1 has " << xpl.size() << " nodes.\n";

  // Create part 2 of area 2.
  xpl.clear();
  ypl.clear();
  if (m_debug) std::cout << "Creating part 2 of area 2.\n";
  imax = jp1 + n1;
  for (int i = kp1; i <= imax; ++i) {
    const int ip1 = i % n1;
    if (m_debug) std::cout << "Adding " << ip1 << " on 1.\n";
    xpl.push_back(xl1[ip1]);
    ypl.push_back(yl1[ip1]);
  }
  // Add the part over area 2.
  if (dir == -1) {
    imax = jp2 > kp2 ? jp2 : jp2 + n2;
    for (int i = imax; i >= kp2; --i) {
      const int ip2 = i % n2;
      if (m_debug) std::cout << "Adding " << ip2 << " on 2.\n";
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
    }
  } else {
    imax = jp2 > kp2 ? kp2 + n2 : kp2;
    for (int i = jp2; i <= imax; ++i) {
      const int ip2 = i % n2;
      if (m_debug) std::cout << "Adding " << ip2 << " on 2.\n";
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
    }
  }
  // Save this part.
  Panel newPanel2 = panel2;
  newPanel2.xv = xpl;
  newPanel2.yv = ypl;
  panelsOut.push_back(std::move(newPanel2));
  if (m_debug) std::cout << "Part 1 has " << xpl.size() << " nodes.\n";
  return true;
}

void ComponentNeBem3d::TraceNonOverlap(
    const std::vector<double>& xp1, const std::vector<double>& yp1,
    const std::vector<double>& xl1, const std::vector<double>& yl1,
    const std::vector<double>& xl2, const std::vector<double>& yl2,
    const std::vector<int>& flags1, const std::vector<int>& flags2,
    const std::vector<int>& links1, const std::vector<int>& links2,
    std::vector<bool>& mark1, int ip1, const Panel& panel1,
    std::vector<Panel>& panelsOut) const {
  const unsigned int n1 = xl1.size();
  const unsigned int n2 = xl2.size();

  // Remember the starting point.
  const int is1 = ip1;
  std::vector<double> xpl;
  std::vector<double> ypl;
  // Add the starting point.
  xpl.push_back(xl1[ip1]);
  ypl.push_back(yl1[ip1]);
  mark1[ip1] = true;
  if (m_debug) std::cout << "Start from point " << ip1 << " on curve 1.\n";
  // Next point.
  ip1 = NextPoint(ip1, n1);
  xpl.push_back(xl1[ip1]);
  ypl.push_back(yl1[ip1]);
  mark1[ip1] = true;
  if (m_debug) std::cout << "Next point is " << ip1 << " on curve 1.\n";

  // Keep track of the curve we are currently following.
  unsigned int il = 1;
  // Direction flag (-1: backward, 0: not set, 1: forward).
  int dir = 0;
  // End-of-curve flag.
  bool eoc = false;
  int ip2 = 0;
  while (!eoc) {
    if (il == 1 && flags1[std::max(ip1, 0)] == 1) {
      // On curve 1 and not on the edge of curve 2?
      ip1 = NextPoint(ip1, n1);
      if (ip1 == is1) {
        eoc = true;
        continue;
      }
      mark1[ip1] = true;
      xpl.push_back(xl1[ip1]);
      ypl.push_back(yl1[ip1]);
      if (m_debug) std::cout << "Went to point " << ip1 << " on curve 1.\n";
    } else if (il == 1) {
      // On curve 1 and on the edge of curve 2?
      ip2 = links1[ip1];
      bool added = false;
      if (dir == +1 || dir == 0) {
        const double xm = 0.5 * (xl2[ip2] + xl2[NextPoint(ip2, n2)]);
        const double ym = 0.5 * (yl2[ip2] + yl2[NextPoint(ip2, n2)]);
        bool inside = false, edge = false;
        Inside(xp1, yp1, xm, ym, inside, edge);
        if (inside) {
          ip2 = NextPoint(ip2, n2);
          il = 2;
          dir = +1;
          ip1 = links2[ip2];
          if (ip1 == is1) {
            eoc = true;
            continue;
          } else if (ip1 >= 0) {
            mark1[ip1] = true;
          }
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          added = true;
          if (m_debug) std::cout << "Added point " << ip2 << " along 2 +.\n";
        }
      }
      if (dir == -1 || dir == 0) {
        const double xm = 0.5 * (xl2[ip2] + xl2[PrevPoint(ip2, n2)]);
        const double ym = 0.5 * (yl2[ip2] + yl2[PrevPoint(ip2, n2)]);
        bool inside = false, edge = false;
        Inside(xp1, yp1, xm, ym, inside, edge);
        if (inside) {
          ip2 = PrevPoint(ip2, n2);
          il = 2;
          dir = -1;
          ip1 = links2[ip2];
          if (ip1 == is1) {
            eoc = true;
            continue;
          } else if (ip1 >= 0) {
            mark1[ip1] = true;
          }
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          added = true;
          if (m_debug) std::cout << "Added point " << ip2 << " along 2 -\n";
        }
      }
      if (!added) {
        ip1 = NextPoint(ip1, n1);
        if (ip1 == is1) {
          eoc = true;
          continue;
        } else if (ip1 >= 0) {
          mark1[ip1] = true;
        }
        xpl.push_back(xl1[ip1]);
        ypl.push_back(yl1[ip1]);
        if (m_debug) std::cout << "Continued over 1.\n";
      }
    } else if (il == 2 && flags2[std::max(ip2, 0)] == 1) {
      // On curve 2 normal vertex (outside 1 hopefully).
      ip2 = dir > 0 ? NextPoint(ip2, n2) : PrevPoint(ip2, n2);
      ip1 = links2[ip2];
      if (ip1 == is1) {
        eoc = true;
        continue;
      } else if (ip1 >= 0) {
        mark1[ip1] = true;
      }
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
      if (m_debug) std::cout << "Went to point " << ip2 << " on 2.\n";
    } else if (il == 2) {
      // On curve 2 and on edge of 1.
      ip1 = links2[ip2];
      ip1 = NextPoint(ip1, n1);
      il = 1;
      if (ip1 == is1) {
        eoc = true;
        continue;
      }
      xpl.push_back(xl1[ip1]);
      ypl.push_back(yl1[ip1]);
      if (m_debug) std::cout << "Resumed 1 at point " << ip1 << ".\n";
    } else {
      // Other cases should not occur.
      std::cerr << m_className << "::TraceNonOverlap: Unexpected case.\n";
      return;
    }
  }

  Panel newPanel = panel1;
  newPanel.xv = xpl;
  newPanel.yv = ypl;
  panelsOut.push_back(std::move(newPanel));
  if (m_debug) {
    std::cout << "End of curve reached, " << xpl.size() << " points.\n";
  }
}

void ComponentNeBem3d::TraceOverlap(
    const std::vector<double>& xp1, const std::vector<double>& yp1,
    const std::vector<double>& xp2, const std::vector<double>& yp2,
    const std::vector<double>& xl1, const std::vector<double>& yl1,
    const std::vector<double>& xl2, const std::vector<double>& yl2,
    const std::vector<int>& flags1, const std::vector<int>& links1,
    const std::vector<int>& links2, std::vector<bool>& mark1, int ip1, int ip2,
    const Panel& panel1, std::vector<Panel>& panelsOut) const {
  int ip1L = -1;
  int ip1LL = -1;

  const unsigned int n1 = xl1.size();
  const unsigned int n2 = xl2.size();

  // Remember the starting points.
  const int is1 = ip1;
  const int is2 = ip2;
  std::vector<double> xpl;
  std::vector<double> ypl;
  xpl.push_back(xl1[ip1]);
  ypl.push_back(yl1[ip1]);
  mark1[ip1] = true;

  // Keep track of the curve we are currently following.
  unsigned int il = 1;
  // Direction flag (-1: backward, 0: not set, 1: forward).
  int dir = 0;
  // End-of-curve flag.
  bool eoc = false;

  if (m_debug) {
    std::cout << "Start from point " << ip1 << " on curve " << il << "\n";
  }
  while (!eoc) {
    ip1LL = ip1L;
    ip1L = ip1;
    if (il == 1) {
      // We are on curve 1. Move to the next point.
      const int ii = NextPoint(ip1, n1);
      // Maybe finished over line 1?
      if (ii == is1) {
        eoc = true;
        continue;
      }
      // See whether the next point of 1 is on the edge or inside of 2.
      bool inside = false, edge = false;
      if (links1[ii] >= 0) {
        edge = true;
      } else if (flags1[ii] == 1) {
        Inside(xp2, yp2, xl1[ii], yl1[ii], inside, edge);
      }
      // If it is, check that it doesn't leave 2 at any stage.
      if (inside || edge) {
        const double xm = 0.5 * (xl1[ip1] + xl1[ii]);
        const double ym = 0.5 * (yl1[ip1] + yl1[ii]);
        Inside(xp2, yp2, xm, ym, inside, edge);
      }
      // If it is, continue over 1.
      if (inside || edge) {
        ip1 = ii;
        if (m_debug) {
          std::cout << "Continued to point " << ip1 << " on " << il << "\n";
        }
        xpl.push_back(xl1[ip1]);
        ypl.push_back(yl1[ip1]);
        mark1[ip1] = true;
        continue;
      }
      // Else we have to continue over 2, ensure we really are on curve 2.
      ip2 = links1[ip1];
      if (ip2 < 0) {
        std::cerr << m_className << "::TraceOverlap: "
                  << "No point 2 reference found; abandoned.\n";
        return;
      }
      // Impose a direction on 2 to avoid returning.
      if (dir == 0) {
        if (links2[NextPoint(ip2, n2)] == ip1LL &&
            links2[PrevPoint(ip2, n2)] == ip1LL) {
          std::cerr << m_className << "::TraceOverlap: "
                    << "Both 2+ and 2- return on 1; not stored.\n";
          return;
        } else if (links2[NextPoint(ip2, n2)] == ip1LL) {
          if (m_debug) std::cout << "2+ is a return to previous point on 1.\n";
          dir = -1;
        } else if (links2[PrevPoint(ip2, n2)] == ip1LL) {
          if (m_debug) std::cout << "2- is a return to previous point on 1.\n";
          dir = +1;
        } else {
          if (m_debug) std::cout << "Both ways are OK.\n";
        }
      }
      // If not, try to continue over 2 in the + direction.
      if (dir == +1 || dir == 0) {
        ip2 = NextPoint(ip2, n2);
        if (ip2 == is2) {
          if (m_debug) std::cout << "Return to start over 2+.\n";
          eoc = true;
          continue;
        }
        Inside(xp1, yp1, xl2[ip2], yl2[ip2], inside, edge);
        if (inside || edge) {
          if (m_debug) std::cout << "Going to 2+ (point " << ip2 << " of 2).\n";
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          dir = +1;
          if (links2[ip2] >= 0) {
            ip1 = links2[ip2];
            mark1[ip1] = true;
            il = 1;
            if (m_debug) {
              std::cout << "This point is also on curve 1: " << ip1 << "\n";
            }
          } else {
            il = 2;
          }
          continue;
        }
        // Continuing in the + direction didn't work so go back a step.
        ip2 = PrevPoint(ip2, n2);
      }
      // Or if this still fails, try 2 in the - direction.
      if (dir == -1 || dir == 0) {
        ip2 = PrevPoint(ip2, n2);
        if (ip2 == is2) {
          if (m_debug) std::cout << "Return to start over 2-\n";
          eoc = true;
          continue;
        }
        Inside(xp1, yp1, xl2[ip2], yl2[ip2], inside, edge);
        if (inside || edge) {
          if (m_debug) std::cout << "Going to 2- (point " << ip2 << " of 2).\n";
          xpl.push_back(xl2[ip2]);
          ypl.push_back(yl2[ip2]);
          dir = -1;
          if (links2[ip2] >= 0) {
            ip1 = links2[ip2];
            mark1[ip1] = true;
            il = 1;
            if (m_debug) {
              std::cout << "This point is also on 1: " << ip1 << ".\n";
            }
          } else {
            il = 2;
          }
          continue;
        }
      }
      // Should not get here.
      if (m_debug) std::cout << "Dead end.\n";
      return;
    } else if (il == 2) {
      // We are on curve 2. Ensure the direction is set.
      if (dir == 0) {
        std::cerr << m_className << "::TraceOverlap: "
                  << "Direction not set; abandoned.\n";
        return;
      }
      // Move to the next point.
      ip2 = dir > 0 ? NextPoint(ip2, n2) : PrevPoint(ip2, n2);
      // Maybe finished over line 2?
      if (ip2 == is2) {
        // Reached the end.
        eoc = true;
        continue;
      }
      // Next step over 2.
      if (m_debug) {
        std::cout << "Stepped over 2 to point " << ip2 << " of 2.\n";
      }
      xpl.push_back(xl2[ip2]);
      ypl.push_back(yl2[ip2]);
      if (links2[ip2] >= 0) {
        ip1 = links2[ip2];
        mark1[ip1] = true;
        il = 1;
        if (m_debug) {
          std::cout << "This point is also on curve 1: " << ip1 << ".\n";
        }
      } else {
        il = 2;
      }
    }
  }

  if (xpl.size() <= 2) {
    if (m_debug) std::cout << "Too few points.\n";
  } else {
    Panel newPanel = panel1;
    newPanel.xv = xpl;
    newPanel.yv = ypl;
    panelsOut.push_back(std::move(newPanel));
  }

  if (m_debug) {
    std::cout << "End of curve reached, " << xpl.size() << " points.\n";
  }
}

bool ComponentNeBem3d::MakePrimitives(const Panel& panelIn,
                                      std::vector<Panel>& panelsOut) const {
  // *-----------------------------------------------------------------------
  // *   PLATRC - Cuts a polygon into right-angled triangles.
  // *-----------------------------------------------------------------------

  // Establish tolerances.
  // TODO! Class member?
  const double epsang = 1.e-6;  // BEMEPA

  if (panelIn.xv.empty() || panelIn.yv.empty() || panelIn.zv.empty()) {
    return false;
  }
  // Determine the mean z value.
  const double zsum =
      std::accumulate(std::begin(panelIn.zv), std::end(panelIn.zv), 0.);
  const double zmean = zsum / panelIn.zv.size();

  std::vector<Panel> stack;
  stack.push_back(panelIn);
  stack.back().zv.clear();
  for (unsigned int k = 0; k < stack.size(); ++k) {
    // Next polygon.
    const auto& xp1 = stack[k].xv;
    const auto& yp1 = stack[k].yv;
    const unsigned int np = xp1.size();
    if (m_debug) std::cout << "Polygon " << k << " with " << np << " nodes.\n";
    if (np <= 2) {
      // Too few nodes.
      if (m_debug) std::cout << "Too few points.\n";
      continue;
    }

    // See whether this is a right-angled triangle.
    if (np == 3) {
      const double x12 = xp1[0] - xp1[1];
      const double y12 = yp1[0] - yp1[1];
      const double x23 = xp1[1] - xp1[2];
      const double y23 = yp1[1] - yp1[2];
      const double x31 = xp1[2] - xp1[0];
      const double y31 = yp1[2] - yp1[0];
      const double x32 = xp1[2] - xp1[1];
      const double y32 = yp1[2] - yp1[1];
      const double x13 = xp1[0] - xp1[2];
      const double y13 = yp1[0] - yp1[2];
      const double x21 = xp1[1] - xp1[0];
      const double y21 = yp1[1] - yp1[0];
      const double s12 = x12 * x12 + y12 * y12;
      const double s32 = x32 * x32 + y32 * y32;
      const double s13 = x13 * x13 + y13 * y13;
      const double s23 = x23 * x23 + y23 * y23;
      const double s31 = x31 * x31 + y31 * y31;
      const double s21 = x21 * x21 + y21 * y21;
      if (fabs(x12 * x32 + y12 * y32) < epsang * sqrt(s12 * s32)) {
        if (m_debug) std::cout << "Right-angled triangle node 2 - done.\n";
        panelsOut.push_back(stack[k]);
        continue;
      } else if (fabs(x13 * x23 + y13 * y23) < epsang * sqrt(s13 * s23)) {
        if (m_debug) std::cout << "Right-angled triangle node 3 - rearrange.\n";
        Panel panel = stack[k];
        panel.xv = {xp1[1], xp1[2], xp1[0]};
        panel.yv = {yp1[1], yp1[2], yp1[0]};
        panelsOut.push_back(std::move(panel));
        continue;
      } else if (fabs(x31 * x21 + y31 * y21) < epsang * sqrt(s31 * s21)) {
        if (m_debug) std::cout << "Right-angled triangle node 1 - rearrange.\n";
        Panel panel = stack[k];
        panel.xv = {xp1[2], xp1[0], xp1[1]};
        panel.yv = {yp1[2], yp1[0], yp1[1]};
        panelsOut.push_back(std::move(panel));
        continue;
      }
    }
    // See whether this is a rectangle.
    if (np == 4) {
      const double x12 = xp1[0] - xp1[1];
      const double y12 = yp1[0] - yp1[1];
      const double x23 = xp1[1] - xp1[2];
      const double y23 = yp1[1] - yp1[2];
      const double x34 = xp1[2] - xp1[3];
      const double y34 = yp1[2] - yp1[3];
      const double x43 = xp1[3] - xp1[2];
      const double y43 = yp1[3] - yp1[2];
      const double x14 = xp1[0] - xp1[3];
      const double y14 = yp1[0] - yp1[3];
      const double x32 = xp1[2] - xp1[1];
      const double y32 = yp1[2] - yp1[1];

      const double s12 = x12 * x12 + y12 * y12;
      const double s23 = x23 * x23 + y23 * y23;
      const double s32 = x32 * x32 + y32 * y32;
      const double s34 = x34 * x34 + y34 * y34;
      const double s43 = x43 * x43 + y43 * y43;
      const double s14 = x14 * x14 + y14 * y14;
      if (fabs(x12 * x32 + y12 * y32) < epsang * sqrt(s12 * s32) &&
          fabs(x23 * x43 + y23 * y43) < epsang * sqrt(s23 * s43) &&
          fabs(x14 * x34 + y14 * y34) < epsang * sqrt(s14 * s34)) {
        if (m_debug) std::cout << "Rectangle.\n";
        panelsOut.push_back(stack[k]);
        continue;
      }
    }
    // See whether there are parallel sides, e.g. a trapezium (UK English).
    if (np >= 4 && SplitTrapezium(stack[k], stack, panelsOut, epsang)) continue;

    // Find a right-angled corner we can cut off.
    if (m_debug) std::cout << "Trying to find a right-angle\n";
    bool corner = false;
    for (unsigned int ip = 0; ip < np; ++ip) {
      // Take only right angles.
      const unsigned int inext = NextPoint(ip, np);
      const unsigned int iprev = PrevPoint(ip, np);

      const double dxprev = xp1[iprev] - xp1[ip];
      const double dyprev = yp1[iprev] - yp1[ip];
      const double dxnext = xp1[inext] - xp1[ip];
      const double dynext = yp1[inext] - yp1[ip];
      if (fabs(dxprev * dxnext + dyprev * dynext) >
          epsang * sqrt((dxprev * dxprev + dyprev * dyprev) *
                        (dxnext * dxnext + dynext * dynext))) {
        continue;
      }
      // Ensure the midpoint is internal.
      if (np > 3) {
        const double xm = 0.5 * (xp1[iprev] + xp1[inext]);
        const double ym = 0.5 * (yp1[iprev] + yp1[inext]);
        bool inside = false, edge = false;
        Inside(xp1, yp1, xm, ym, inside, edge);
        if (!inside) continue;
      }
      // Check all vertex crossings.
      bool cross = false;
      for (unsigned int jp = 0; jp < np; ++jp) {
        // Accept immediate contact.
        const unsigned int jnext = NextPoint(jp, np);
        if (jp == iprev || jp == ip || jp == inext || jnext == iprev ||
            jnext == ip || jnext == inext)
          continue;
        // Check crossing.
        double xc = 0., yc = 0.;
        if (Crossing(xp1[iprev], yp1[iprev], xp1[inext], yp1[inext], xp1[jp],
                     yp1[jp], xp1[jnext], yp1[jnext], xc, yc)) {
          cross = true;
          break;
        }
      }
      if (cross) continue;
      // Found a triangle, introduce shorthand node references.
      if (m_debug) {
        std::cout << "Cutting at right-angled corner " << ip << "\n";
      }
      corner = true;
      Panel panel = stack[k];
      panel.xv = {xp1[iprev], xp1[ip], xp1[inext]};
      panel.yv = {yp1[iprev], yp1[ip], yp1[inext]};
      stack.push_back(std::move(panel));
      // Eliminate this node from the polygon.
      stack[k].xv.erase(stack[k].xv.begin() + ip);
      stack[k].yv.erase(stack[k].yv.begin() + ip);
      stack.push_back(std::move(stack[k]));
      if (m_debug) {
        std::cout << "Going for another pass, NP = " << np << "\n";
      }
      break;
    }
    if (corner) continue;

    // Find any corner we can cut off.
    if (m_debug) std::cout << "Trying to find a corner\n";
    corner = false;
    for (unsigned int ip = 0; ip < np; ++ip) {  // 20
      const unsigned int iprev = PrevPoint(ip, np);
      const unsigned int inext = NextPoint(ip, np);
      // Ensure the midpoint is internal.
      if (np > 3) {
        const double xm = 0.5 * (xp1[iprev] + xp1[inext]);
        const double ym = 0.5 * (yp1[iprev] + yp1[inext]);
        bool inside = false, edge = false;
        Inside(xp1, yp1, xm, ym, inside, edge);
        if (!inside) continue;
      }
      // Check all vertex crossings.
      bool cross = false;
      for (unsigned int jp = 0; jp < np; ++jp) {
        const unsigned int jj = NextPoint(jp, np);
        // Accept immediate contact.
        if (jp == iprev || jp == ip || jp == inext || jj == iprev || jj == ip ||
            jj == inext)
          continue;
        // Check crossing.
        double xc = 0., yc = 0.;
        if (Crossing(xp1[iprev], yp1[iprev], xp1[inext], yp1[inext], xp1[jp],
                     yp1[jp], xp1[jj], yp1[jj], xc, yc)) {
          cross = true;
          break;
        }
      }
      if (cross) continue;
      // Found a triangle, introduce shorthand node references.
      if (m_debug) std::cout << "Cutting at corner " << ip << "\n";
      corner = true;
      const double x1 = xp1[iprev];
      const double x2 = xp1[ip];
      const double x3 = xp1[inext];
      const double y1 = yp1[iprev];
      const double y2 = yp1[ip];
      const double y3 = yp1[inext];
      const double s21 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
      const double s31 = (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1);
      const double s32 = (x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2);
      const double s12 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
      const double s13 = (x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3);
      const double s23 = (x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3);
      // Find the biggest opening angle.
      const double a1 =
          ((x2 - x1) * (x3 - x1) + (y2 - y1) * (y3 - y1)) / sqrt(s21 * s31);
      const double a2 =
          ((x3 - x2) * (x1 - x2) + (y3 - y2) * (y1 - y2)) / sqrt(s32 * s12);
      const double a3 =
          ((x1 - x3) * (x2 - x3) + (y1 - y3) * (y2 - y3)) / sqrt(s13 * s23);
      if (m_debug) {
        const double phi1 = acos(a1);
        const double phi2 = acos(a2);
        const double phi3 = acos(a3);
        std::cout << "Angles: " << RadToDegree * phi1 << ", "
                  << RadToDegree * phi2 << ", " << RadToDegree * phi3 << "\n";
        std::cout << "Sum = " << RadToDegree * (phi1 + phi2 + phi3) << "\n";
      }
      // See whether one angle is more or less right-angled.
      if (fabs(a1) < epsang || fabs(a2) < epsang || fabs(a3) < epsang) {
        if (m_debug) std::cout << "Right-angled corner cut off.\n";
        Panel panel = stack[k];
        if (fabs(a1) < epsang) {
          panel.xv = {x3, x1, x2};
          panel.yv = {y3, y1, y2};
        } else if (fabs(a2) < epsang) {
          panel.xv = {x1, x2, x3};
          panel.yv = {y1, y2, y3};
        } else {
          panel.xv = {x2, x3, x1};
          panel.yv = {y2, y3, y1};
        }
        stack.push_back(std::move(panel));
      } else if (a1 <= a2 && a1 <= a3) {
        if (m_debug) std::cout << "A1 < A2, A3 - adding 2 triangles.\n";
        const double xc = x2 + a2 * (x3 - x2) * sqrt(s12 / s32);
        const double yc = y2 + a2 * (y3 - y2) * sqrt(s12 / s32);
        Panel panel1 = stack[k];
        panel1.xv = {x3, xc, x1};
        panel1.yv = {y3, yc, y1};
        stack.push_back(std::move(panel1));
        Panel panel2 = stack[k];
        panel2.xv = {x2, xc, x1};
        panel2.yv = {y2, yc, y1};
        stack.push_back(std::move(panel2));
      } else if (a2 <= a1 && a2 <= a3) {
        if (m_debug) std::cout << "A2 < A1, A3 - adding 2 triangles.\n";
        const double xc = x3 + a3 * (x1 - x3) * sqrt(s23 / s13);
        const double yc = y3 + a3 * (y1 - y3) * sqrt(s23 / s13);
        Panel panel1 = stack[k];
        panel1.xv = {x1, xc, x2};
        panel1.yv = {y1, yc, y2};
        stack.push_back(std::move(panel1));
        Panel panel2 = stack[k];
        panel2.xv = {x3, xc, x2};
        panel2.yv = {y3, yc, y2};
        stack.push_back(std::move(panel2));
      } else {
        if (m_debug) std::cout << "A3 < A1, A2 - adding 2 triangles.\n";
        const double xc = x1 + a1 * (x2 - x1) * sqrt(s31 / s21);
        const double yc = y1 + a1 * (y2 - y1) * sqrt(s31 / s21);
        Panel panel1 = stack[k];
        panel1.xv = {x1, xc, x3};
        panel1.yv = {y1, yc, y3};
        stack.push_back(std::move(panel1));
        Panel panel2 = stack[k];
        panel2.xv = {x2, xc, x3};
        panel2.yv = {y2, yc, y3};
        stack.push_back(std::move(panel2));
      }
      // Eliminate this node from the polygon.
      stack[k].xv.erase(stack[k].xv.begin() + ip);
      stack[k].yv.erase(stack[k].yv.begin() + ip);
      stack.push_back(std::move(stack[k]));
      if (m_debug) {
        std::cout << "Going for another pass, NP = " << np << "\n";
      }
      break;
    }
    if (corner) continue;
    std::cerr << m_className << "::MakePrimitives:\n    "
              << "Unable to identify a corner to cut, probably a degenerate "
                 "polygon.\n";
    // Next stack element.
  }
  for (auto& panel : panelsOut) {
    panel.zv.assign(panel.xv.size(), zmean);
  }
  return true;
}

bool ComponentNeBem3d::SplitTrapezium(const Panel panelIn,
                                      std::vector<Panel>& stack,
                                      std::vector<Panel>& panelsOut,
                                      const double epsang) const {
  const auto xp1 = panelIn.xv;
  const auto yp1 = panelIn.yv;
  const unsigned int np = xp1.size();
  for (unsigned int ip = 0; ip < np; ++ip) {
    const unsigned int inext = NextPoint(ip, np);
    const double xi0 = xp1[ip];
    const double yi0 = yp1[ip];
    const double xi1 = xp1[inext];
    const double yi1 = yp1[inext];
    const double dxi = xi0 - xi1;
    const double dyi = yi0 - yi1;
    const double si2 = dxi * dxi + dyi * dyi;
    for (unsigned int jp = ip + 2; jp < np; ++jp) {
      const unsigned int jnext = NextPoint(jp, np);
      // Skip adjacent segments.
      if (ip == jp || ip == jnext || inext == jp || inext == jnext) {
        continue;
      }
      const double xj0 = xp1[jp];
      const double yj0 = yp1[jp];
      const double xj1 = xp1[jnext];
      const double yj1 = yp1[jnext];
      const double dxj = xj0 - xj1;
      const double dyj = yj0 - yj1;
      const double sj2 = dxj * dxj + dyj * dyj;
      // Require parallelism.
      const double sij = sqrt(si2 * sj2);
      if (fabs(dxi * dxj + dyi * dyj + sij) > epsang * sij) continue;
      if (m_debug) {
        std::cout << "Found parallel sections: " << ip << ", " << jp << "\n";
      }
      // Avoid division by zero.
      if (sj2 <= 0 || si2 <= 0) {
        std::cerr << m_className << "::SplitTrapezium:\n    "
                  << "Zero norm segment found; skipped.\n";
        continue;
      }
      // Establish the cutting lines.
      const double xl1 =
          ((xi0 - xj0) * (xj1 - xj0) + (yi0 - yj0) * (yj1 - yj0)) / sj2;
      const double xl2 =
          ((xi1 - xj0) * (xj1 - xj0) + (yi1 - yj0) * (yj1 - yj0)) / sj2;
      const double xl3 =
          ((xj0 - xi0) * (xi1 - xi0) + (yj0 - yi0) * (yi1 - yi0)) / si2;
      const double xl4 =
          ((xj1 - xi0) * (xi1 - xi0) + (yj1 - yi0) * (yi1 - yi0)) / si2;
      if (m_debug) {
        std::cout << "xl1 = " << xl1 << ", xl2 = " << xl2 << ", "
                  << "xl3 = " << xl3 << ", xl4 = " << xl4 << "\n";
      }
      // Check that there is at all a rectangle.
      const double r1 = (xl1 + epsang) * (1. + epsang - xl1);
      const double r2 = (xl2 + epsang) * (1. + epsang - xl2);
      const double r3 = (xl3 + epsang) * (1. + epsang - xl3);
      const double r4 = (xl4 + epsang) * (1. + epsang - xl4);
      if ((r1 < 0 && r4 < 0) || (r2 < 0 && r3 < 0)) {
        if (m_debug) std::cout << "No rectangle.\n";
        continue;
      }
      // Determine the rectangular part.
      std::vector<double> xpl(4, 0.);
      std::vector<double> ypl(4, 0.);
      if (r1 >= 0) {
        xpl[0] = xi0;
        ypl[0] = yi0;
        xpl[1] = xj0 + xl1 * (xj1 - xj0);
        ypl[1] = yj0 + xl1 * (yj1 - yj0);
      } else if (r4 >= 0) {
        xpl[0] = xi0 + xl4 * (xi1 - xi0);
        ypl[0] = yi0 + xl4 * (yi1 - yi0);
        xpl[1] = xj1;
        ypl[1] = yj1;
      }
      if (r2 >= 0) {
        xpl[2] = xj0 + xl2 * (xj1 - xj0);
        ypl[2] = yj0 + xl2 * (yj1 - yj0);
        xpl[3] = xi1;
        ypl[3] = yi1;
      } else if (r3 >= 0) {
        xpl[2] = xj0;
        ypl[2] = yj0;
        xpl[3] = xi0 + xl3 * (xi1 - xi0);
        ypl[3] = yi0 + xl3 * (yi1 - yi0);
      }
      // Verify that the midpoints of these lines are internal.
      double xm = 0.5 * (xpl[0] + xpl[1]);
      double ym = 0.5 * (ypl[0] + ypl[1]);
      bool inside = false, edge = false;
      Inside(xp1, yp1, xm, ym, inside, edge);
      if (!(inside || edge)) {
        if (m_debug) std::cout << "Midpoint 1 not internal.\n";
        continue;
      }
      xm = 0.5 * (xpl[2] + xpl[3]);
      ym = 0.5 * (ypl[2] + ypl[3]);
      Inside(xp1, yp1, xm, ym, inside, edge);
      if (!(inside || edge)) {
        if (m_debug) std::cout << "Midpoint 2 not internal.\n";
        continue;
      }

      const unsigned int iprev = PrevPoint(ip, np);
      const unsigned int jprev = PrevPoint(jp, np);
      // Ensure there are no crossings, accepting contact.
      bool cross = false;
      for (unsigned int i = 0; i < np; ++i) {
        if ((i == iprev && r1 >= 0) || i == ip || (i == inext && r2 >= 0) ||
            (i == jprev && r3 >= 0) || i == jp || (i == jnext && r4 >= 0))
          continue;
        const unsigned int ii = NextPoint(i, np);
        double xc = 0., yc = 0.;
        if (Crossing(xp1[i], yp1[i], xp1[ii], yp1[ii], xpl[0], ypl[0], xpl[1],
                     ypl[1], xc, yc) ||
            Crossing(xp1[i], yp1[i], xp1[ii], yp1[ii], xpl[2], ypl[2], xpl[3],
                     ypl[3], xc, yc)) {
          if (m_debug) {
            std::cout << "Crossing (edge " << i << ", " << ii
                      << "), ip = " << ip << ", jp = " << jp << "\n";
          }
          cross = true;
          break;
        }
      }
      if (cross) continue;
      // Add the rectangular part.
      if ((fabs(xl1) < epsang && fabs(xl3) < epsang) ||
          (fabs(1. - xl2) < epsang && fabs(1. - xl4) < epsang)) {
        if (m_debug) std::cout << "Not stored, degenerate.\n";
      } else {
        if (m_debug) std::cout << "Adding rectangle.\n";
        Panel panel = panelIn;
        panel.xv = xpl;
        panel.yv = ypl;
        panelsOut.push_back(std::move(panel));
      }
      // First non-rectangular section.
      xpl.clear();
      ypl.clear();
      for (unsigned int i = jp + 1; i <= ip + np; ++i) {
        const unsigned int ii = i % np;
        xpl.push_back(xp1[ii]);
        ypl.push_back(yp1[ii]);
      }
      if (r1 >= 0 && r4 >= 0) {
        if (m_debug) std::cout << "1-4 degenerate\n";
      } else if (r1 >= 0) {
        if (m_debug) std::cout << "Using 1\n";
        xpl.push_back(xj0 + xl1 * (xj1 - xj0));
        ypl.push_back(yj0 + xl1 * (yj1 - yj0));
      } else if (r4 >= 0) {
        if (m_debug) std::cout << "Using 4\n";
        xpl.push_back(xi0 + xl4 * (xi1 - xi0));
        ypl.push_back(yi0 + xl4 * (yi1 - yi0));
      } else {
        if (m_debug) std::cout << "Neither 1 nor 4, should not happen\n";
      }
      if (xpl.size() < 3) {
        if (m_debug) {
          std::cout << "Not stored, only " << xpl.size() << " vertices.\n";
        }
      } else {
        if (m_debug) {
          std::cout << "Adding non-rectangular part with " << xpl.size()
                    << " vertices.\n";
        }
        Panel panel = panelIn;
        panel.xv = xpl;
        panel.yv = ypl;
        stack.push_back(std::move(panel));
      }
      // Second non-rectangular section.
      xpl.clear();
      ypl.clear();
      for (unsigned int i = ip + 1; i <= jp; ++i) {
        const unsigned int ii = i % np;
        xpl.push_back(xp1[ii]);
        ypl.push_back(yp1[ii]);
      }
      if (r2 >= 0 && r3 >= 0) {
        if (m_debug) std::cout << "2-3 degenerate\n";
      } else if (r2 >= 0) {
        if (m_debug) std::cout << "Using 2\n";
        xpl.push_back(xj0 + xl2 * (xj1 - xj0));
        ypl.push_back(yj0 + xl2 * (yj1 - yj0));
      } else if (r3 >= 0) {
        if (m_debug) std::cout << "Using 3\n";
        xpl.push_back(xi0 + xl3 * (xi1 - xi0));
        ypl.push_back(yi0 + xl3 * (yi1 - yi0));
      } else {
        if (m_debug) std::cout << "Neither 2 nor 3, should not happen\n";
      }
      if (xpl.size() < 3) {
        if (m_debug) {
          std::cout << "Not stored, only " << xpl.size() << " vertices.\n";
        }
      } else {
        if (m_debug) {
          std::cout << "Adding non-rectangular part with " << xpl.size()
                    << " vertices.\n";
        }
        Panel panel = panelIn;
        panel.xv = xpl;
        panel.yv = ypl;
        stack.push_back(std::move(panel));
      }
      return true;
    }
  }
  return false;
}

bool ComponentNeBem3d::GetPanel(const unsigned int i, double& a, double& b,
                                double& c, std::vector<double>& xv,
                                std::vector<double>& yv,
                                std::vector<double>& zv) {
  if (i >= m_panels.size()) {
    std::cerr << m_className << ":: Incorrect panel number.\n";
    return false;
  }
  const auto& panel = m_panels[i];
  a = panel.a;
  b = panel.b;
  c = panel.c;
  xv = panel.xv;
  yv = panel.yv;
  zv = panel.zv;
  return true;
}

void ComponentNeBem3d::Reset() {
  m_panels.clear();
  m_elements.clear();
  m_ready = false;
}

void ComponentNeBem3d::UpdatePeriodicity() {
  std::cerr << m_className << "::UpdatePeriodicity:\n"
            << "    Periodicities are not supported.\n";
}
}
