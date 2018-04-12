#ifndef TRACER_H
#define TRACER_H

#include <vector>
#include "surface.h"

// Point is close if function value less than PTCLOSE
#define PTCLOSE 0.00001

class Point
{
public:
  Point() { x=0.0; y=0.0; z=0.0; }
  Point( double _x, double _y, double _z ) { x=_x; y=_y; z=_z; }
  void set( double _x, double _y, double _z ) { x=_x; y=_y; z=_z; }
  double distsq( Point &p )
    { double dx=p.x-x, dy=p.y-y, dz=p.z-z; return dx*dx+dy*dy+dz*dz; }
  double X() { return x; }
  double Y() { return y; }
  double Z() { return z; }
private:
  double x, y, z;
};

//############################################################################
//  ContourData contains a single contour
class ContourData
{
public:
  ContourData() { computebb(); }
  ~ContourData() { clear(); }
  int  size() { return points.size(); }
  void addpoint( Point &p );
  bool isClose( Point &p, double eps );
  void clear() { points.clear(); }
  void computebb();
  Point getpoint( int i ) { return points[i]; }
private:
  void movetocg( Point &p );
  // contour elements
  vector<Point> points;
  // bounding box
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
};

//############################################################################
class Tracer
{
public:
  Tracer() { epsilon = 0.01; }
  ~Tracer() { clearContours(); }
  void setEpsilon( double eps ) { epsilon = eps; }
  void clearContours() { contours.clear(); }
  void addContour( Surface *s, Point &p );
  int  getNumContours() { return contours.size(); }
  int  getContourSize( int num ) { return contours[num].size(); }
  void getContourPoint( int contour, int point, Point &p ) 
       { p = contours[contour].getpoint(point); }
private:
  void movetocg( Surface *s, Point &p );
  bool pointOnContour( Point &p );
  vector<ContourData> contours;
  double epsilon;
};

#endif  // TRACER_H

