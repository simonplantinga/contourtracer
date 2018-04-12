#include "tracer.h"
#include <math.h>
#define K 0.5   // constant for contour corrector


//#define REPORT  // print results
#ifdef REPORT
#include <stdio.h>
#endif

//############################################################################
//  ContourData members
void
ContourData::addpoint( Point &p )
{
  points.push_back( p );
  if ( p.X() < minx )   minx = p.X();
  if ( p.X() > maxx )   maxx = p.X();
  if ( p.Y() < miny )   miny = p.Y();
  if ( p.Y() > maxy )   maxy = p.Y();
  if ( p.Z() < minz )   minz = p.Z();
  if ( p.Z() > maxz )   maxz = p.Z();
}

bool
ContourData::isClose( Point &p, double eps )
{
  if ( p.X() < minx - eps )     return false;
  if ( p.X() > maxx + eps )     return false;
  if ( p.Y() < miny - eps )     return false;
  if ( p.Y() > maxy + eps )     return false;
  if ( p.Z() < minz - eps )     return false;
  if ( p.Z() > maxz + eps )     return false;

  for ( int i=0; i<points.size(); i++ )
    {
/*
      double dx = p.X() - points[i].X();
      double dy = p.Y() - points[i].Y();
      double dz = p.Z() - points[i].Z();
      if ( dx*dx + dy*dy + dz*dz < eps*eps )
        return true;
*/
      if ( points[i].distsq(p) < eps*eps )
        return true;
    }

  return false; 
}

void
ContourData::computebb()
{
  minx = maxx = 0.0;
  miny = maxy = 0.0;
  minz = maxz = 0.0;

  for ( int i=0; i<points.size(); i++)
    {
      if ( points[i].X() < minx )   minx = points[i].X();
      if ( points[i].X() > maxx )   maxx = points[i].X();
      if ( points[i].Y() < miny )   miny = points[i].Y();
      if ( points[i].Y() > maxy )   maxy = points[i].Y();
      if ( points[i].Z() < minz )   minz = points[i].Z();
      if ( points[i].Z() > maxz )   maxz = points[i].Z();
    }
}

//############################################################################
//  Tracer members
bool
Tracer::pointOnContour( Point &p )
{
  for ( int i=0; i<contours.size(); i++ )
    if ( contours[i].isClose( p, epsilon ) )
      return true;
  return false;
}

void
Tracer::addContour( Surface *s, Point &p )
{
  double delta = 0.1;

  movetocg( s, p );

  if ( pointOnContour( p ) )
    return;

  ContourData contour;
  Point p0 = p;

  contour.addpoint( p0 );


  while ( 1 )
    {
      //##########################
      if ( contour.size() > 15000 )
        break;
      //##########################

      // advance p
      delta = delta * 1.2;
      double x,y,z;

      while (1)
        {
          // compute new position (x,y,z)
          double fx = s->Fx( p.X(), p.Y(), p.Z() );
          double fy = s->Fy( p.X(), p.Y(), p.Z() );
          double fz = s->Fz( p.X(), p.Y(), p.Z() );
          double fxz = s->Fxz( p.X(), p.Y(), p.Z() );
          double fyz = s->Fyz( p.X(), p.Y(), p.Z() );
          double fzz = s->Fzz( p.X(), p.Y(), p.Z() );
          x = fy*fzz - fz*fyz;
          y = fz*fxz - fx*fzz;
          z = fx*fyz - fy*fxz;
          double gradsq = x*x+y*y+z*z;
          x = p.X() + delta * ( x / gradsq );
          y = p.Y() + delta * ( y / gradsq );
          z = p.Z() + delta * ( z / gradsq );
          
          // compute bounding box (bx,by,bz)
          double dist = sqrt( (x-p.X())*(x-p.X())
                            + (y-p.Y())*(y-p.Y())
                            + (z-p.Z())*(z-p.Z()) ) * 1.1;
          Interval bx(p.X()-dist, p.X()+dist);
          Interval by(p.Y()-dist, p.Y()+dist);
          Interval bz(p.Z()-dist, p.Z()+dist);

          // compute w at p
          double wx = fy*fzz - fz*fyz;
          double wy = fz*fxz - fx*fzz;
          double wz = fx*fyz - fy*fxz;
          double nw = sqrt(wx*wx+wy*wy+wz*wz);
          wx = wx / nw;
          wy = wy / nw;
          wz = wz / nw;

          // compute w over box
          Interval ifx = s->Fx(bx,by,bz);
          Interval ify = s->Fy(bx,by,bz);
          Interval ifz = s->Fz(bx,by,bz);
          Interval ifxz = s->Fxz(bx,by,bz);
          Interval ifyz = s->Fyz(bx,by,bz);
          Interval ifzz = s->Fzz(bx,by,bz);

          Interval iwx = ify*ifzz - ifz*ifyz;
          Interval iwy = ifz*ifxz - ifx*ifzz;
          Interval iwz = ifx*ifyz - ify*ifxz;
          Interval inw = sqrt(sqr(iwx)+sqr(iwy)+sqr(iwz));

          iwx = iwx / inw;
          iwy = iwy / inw;
          iwz = iwz / inw;

          Interval I = wx*iwx + wy*iwy + wz*iwz;

          // perform the interval test
          if ( I.inf() > 0.71 )
            break;
          else
            delta = delta * 0.8;

          if (delta < 0.00001)    // JUST TO BE SAFE!!!
            break;
        }
      p.set( x, y, z );

      movetocg( s, p );

      if ( contour.size() > 10 )
        if ( p.distsq(p0) < delta*delta )
          break;

#ifdef REPORT
  printf("added (%f, %f, %f); %i points on contour\n",x,y,z,contour.size());
#endif
      // add to contour
      contour.addpoint( p );
    } 

  cout << "Contour size " << contour.size() << endl;

  // add contour to contour components
  contours.push_back( contour );
}

void 
Tracer::movetocg( Surface *s, Point &p )
{
  double f;
  double fx,fy,fz,fxz,fyz,fzz;
  double x,y,z;
  double gs;

  x = p.X();
  y = p.Y();
  z = p.Z();
  while ( fabs(s->F(x,y,z)) > PTCLOSE || fabs(s->Fz(x,y,z)) > PTCLOSE )
    {
      f = s->F(x,y,z);
      if ( fabs(f) > PTCLOSE )
        {
          fx = s->Fx(x,y,z);
          fy = s->Fy(x,y,z);
          fz = s->Fz(x,y,z);
          gs = fx*fx + fy*fy + fz*fz;
          x -= f * fx / gs;
          y -= f * fy / gs;
          z -= f * fz / gs;
        }

      fz = s->Fz(x,y,z);
      if ( fabs(fz) > PTCLOSE )
        {
          fxz = s->Fxz(x,y,z);
          fyz = s->Fyz(x,y,z);
          fzz = s->Fzz(x,y,z);
          gs = fxz*fxz + fyz*fyz + fzz*fzz;
          x -= fz * fxz / gs;
          y -= fz * fyz / gs;
          z -= fz * fzz / gs;
        }
    }
  p.set( x, y, z);
}

