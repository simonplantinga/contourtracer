//############################################################################
#ifndef SURFACE_H
#define SURFACE_H

#include "Interval.h"

//############################################################################
// Abstract class for an implicit surface
//############################################################################
class Surface
{
public:
  Surface() {}    
  virtual ~Surface() {} 

  virtual double F(double x, double y, double z) = 0;
  virtual double Fx(double x, double y, double z) = 0;
  virtual double Fy(double x, double y, double z) = 0;
  virtual double Fz(double x, double y, double z) = 0;
  virtual double Fxx(double x, double y, double z) = 0; 
  virtual double Fyy(double x, double y, double z) = 0; 
  virtual double Fzz(double x, double y, double z) = 0; 
  virtual double Fxy(double x, double y, double z) = 0; 
  virtual double Fxz(double x, double y, double z) = 0; 
  virtual double Fyz(double x, double y, double z) = 0; 

  virtual Interval F(Interval x, Interval y, Interval z) = 0;
  virtual Interval Fx(Interval x, Interval y, Interval z) = 0;
  virtual Interval Fy(Interval x, Interval y, Interval z) = 0;
  virtual Interval Fz(Interval x, Interval y, Interval z) = 0;
  virtual Interval Fxx(Interval x, Interval y, Interval z) = 0; 
  virtual Interval Fyy(Interval x, Interval y, Interval z) = 0; 
  virtual Interval Fzz(Interval x, Interval y, Interval z) = 0; 
  virtual Interval Fxy(Interval x, Interval y, Interval z) = 0; 
  virtual Interval Fxz(Interval x, Interval y, Interval z) = 0; 
  virtual Interval Fyz(Interval x, Interval y, Interval z) = 0; 

  virtual double MinX() = 0;
  virtual double MaxX() = 0;
  virtual double MinY() = 0;
  virtual double MaxY() = 0;
  virtual double MinZ() = 0;
  virtual double MaxZ() = 0;
};

#endif  // SURFACE_H
//############################################################################

