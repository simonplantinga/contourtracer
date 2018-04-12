//############################################################################
#ifndef METABALLS_H
#define METABALLS_H

#include <vector>

#include "surface.h"

//############################################################################
// Single metaball
//############################################################################
class MBall
{
private:
  MBall(double _x, double _y, double _z, double _r, double _w);
  double x;         // centre
  double y;
  double z;
  double r;         // radius
  double w;         // weight

  double rsq;       // squared radius
  double distsq;    // squared distance for last computed point
  double dist;      // distance for last computed point

  friend class MetaBalls;
};

//############################################################################
// Metaball object
//############################################################################
class MetaBalls : public Surface
{
public:
  MetaBalls();
  ~MetaBalls(); 

  double F(double x, double y, double z);
  double Fx(double x, double y, double z);
  double Fy(double x, double y, double z);
  double Fz(double x, double y, double z);
  double Fxx(double x, double y, double z); 
  double Fyy(double x, double y, double z); 
  double Fzz(double x, double y, double z); 
  double Fxy(double x, double y, double z); 
  double Fxz(double x, double y, double z); 
  double Fyz(double x, double y, double z); 

  Interval F(Interval x, Interval y, Interval z);
  Interval Fx(Interval x, Interval y, Interval z);
  Interval Fy(Interval x, Interval y, Interval z);
  Interval Fz(Interval x, Interval y, Interval z);
  Interval Fxx(Interval x, Interval y, Interval z); 
  Interval Fyy(Interval x, Interval y, Interval z); 
  Interval Fzz(Interval x, Interval y, Interval z); 
  Interval Fxy(Interval x, Interval y, Interval z); 
  Interval Fxz(Interval x, Interval y, Interval z); 
  Interval Fyz(Interval x, Interval y, Interval z); 

  double MinX();
  double MaxX();
  double MinY();
  double MaxY();
  double MinZ();
  double MaxZ();

  int NumBalls();
  void AddBall( double x, double y, double z, double r, double w );
  void SetBallCentre( int i, double x, double y, double z );
  void SetBallRadius( int i, double r );
  void SetBallWeight( int i, double w );
  void SetThreshold( double t );

  void GetBallCentre( int i, double &x, double &y, double &z );
  double GetBallRadius( int i );
  double GetBallWeight( int i );
  double GetThreshold();

private:
  void ComputeDistances( double _x, double _y, double _z );
  void ComputeBounds();

  vector<MBall> mballs;
  double threshold;

  double x, y, z;   // last computed point

  double minx, maxx, miny, maxy, minz, maxz;
};

#endif  // METABALLS_H
//############################################################################

