//############################################################################
#include "metaballs.h"

//############################################################################
// Metaball implementation using the C2 kernel 
// k(r) = -6*r^5 + 15*r^4 - 10*r^3 + 1
//############################################################################

//############################################################################
// MBall
//############################################################################
MBall::MBall( double _x, double _y, double _z, double _r, double _w )
{
  x=_x; y=_y; z=_z;
  r=_r; w=_w;
  rsq = r * r;    // precompute squared radius
}

//############################################################################
// Metaball object
//############################################################################
//----------------------------------------------------------------------------
MetaBalls::MetaBalls()
  : Surface()
{
  x = 1.0;    // force recomputation
  ComputeDistances( 0.0, 0.0, 0.0 );

  SetThreshold( 1.0 );
  ComputeBounds();
}

MetaBalls::~MetaBalls()
{
  mballs.clear();
}

//----------------------------------------------------------------------------
double 
MetaBalls::F(double x, double y, double z)
{
  double sum = 0.0;
  double t;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );
          double r4 = r2 * r2;

          sum += mballs[i].w * ( -6 * r4 * r + 15 * r4 - 10 * r2 * r + 1 );
        }
    }
  return sum - threshold;
}

double 
MetaBalls::Fx(double x, double y, double z)
{
  double sum = 0.0;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w 
                  * ( -30*r2*r + 60*r2 - 30*r)
                  * ( x - mballs[i].x ) / ( mballs[i].rsq );
        }
    }
  return sum;
}

double 
MetaBalls::Fy(double x, double y, double z)
{
  double sum = 0.0;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w 
                  * ( -30*r2*r + 60*r2 - 30*r)
                  * ( y - mballs[i].y ) / ( mballs[i].rsq );
        }
    }
  return sum;
}

double 
MetaBalls::Fz(double x, double y, double z)
{
  double sum = 0.0;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w 
                  * ( -30*r2*r + 60*r2 - 30*r)
                  * ( z - mballs[i].z ) / ( mballs[i].rsq );
        }
    }
  return sum;
}

double 
MetaBalls::Fxx(double x, double y, double z) 
{
  double sum = 0.0;
  double t;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w
                  * ( (-90*r2 + 120*r - 30 )
                        * ( x - mballs[i].x ) * ( x - mballs[i].x ) 
                        / ( mballs[i].rsq * r ) 
                      + ( -30*r2*r + 60*r2 - 30*r ) )
                  / mballs[i].rsq;
        }
    }
  return sum;
}

double 
MetaBalls::Fyy(double x, double y, double z) 
{
  double sum = 0.0;
  double t;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w
                  * ( (-90*r2 + 120*r - 30 )
                        * ( y - mballs[i].y ) * ( y - mballs[i].y ) 
                        / ( mballs[i].rsq * r ) 
                      + ( -30*r2*r + 60*r2 - 30*r ) )
                  / mballs[i].rsq;
        }
    }
  return sum;
}

double 
MetaBalls::Fzz(double x, double y, double z) 
{
  double sum = 0.0;
  double t;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w
                  * ( (-90*r2 + 120*r - 30 )
                        * ( z - mballs[i].z ) * ( z - mballs[i].z ) 
                        / ( mballs[i].rsq * r ) 
                      + ( -30*r2*r + 60*r2 - 30*r ) )
                  / mballs[i].rsq;
        }
    }
  return sum;
}

double 
MetaBalls::Fxy(double x, double y, double z) 
{
  double sum = 0.0;
  double t;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w 
                  * ( -90*r2 + 120*r - 30 )
                  * ( y - mballs[i].y ) * ( x - mballs[i].x )
                  / ( mballs[i].rsq * mballs[i].rsq * r );
        }
    }
  return sum;
}

double 
MetaBalls::Fxz(double x, double y, double z) 
{
  double sum = 0.0;
  double t;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w 
                  * ( -90*r2 + 120*r - 30 )
                  * ( z - mballs[i].z ) * ( x - mballs[i].x )
                  / ( mballs[i].rsq * mballs[i].rsq * r );
        }
    }
  return sum;
}

double 
MetaBalls::Fyz(double x, double y, double z) 
{
  double sum = 0.0;
  double t;

  ComputeDistances( x, y, z );
  for ( int i=0; i<mballs.size(); i++ )
    {
      if ( mballs[i].distsq < mballs[i].rsq )
        {
          double r2 = mballs[i].distsq / mballs[i].rsq;
          double r = sqrt( r2 );

          sum += mballs[i].w 
                  * ( -90*r2 + 120*r - 30 )
                  * ( z - mballs[i].z ) * ( y - mballs[i].y )
                  / ( mballs[i].rsq * mballs[i].rsq * r );
        }
    }
  return sum;
}

//----------------------------------------------------------------------------
Interval 
MetaBalls::F(Interval x, Interval y, Interval z)
{
  Interval result(0);
  Interval dx, dy, dz, distsq;
  double rsq;
  
  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx)+sqr(dy)+sqr(dz);
      rsq = mballs[i].rsq;
      // only process intervals intersecting the influence region
      if ( distsq.inf() < rsq )
        {
          // ignore part of interval outside influence region
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          Interval r4 = sqr( r2 );
          result += mballs[i].w * ( -6*r4*r + 15*r4 - 10*r2*r + 1 );
        }
    }
  return result - threshold;
}

Interval 
MetaBalls::Fx(Interval x, Interval y, Interval z)
{
  Interval result(0);
  Interval dx, dy, dz, distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx)+sqr(dy)+sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dx.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * ( -30*r2*r + 60*r2 - 30*r ) * dx / rsq;
        }
    }
  return result;
}

Interval 
MetaBalls::Fy(Interval x, Interval y, Interval z)
{
  Interval result(0);
  Interval dx, dy, dz, distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx)+sqr(dy)+sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dy.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * ( -30*r2*r + 60*r2 - 30*r ) * dy / rsq;
        }
    }
  return result;
}

Interval 
MetaBalls::Fz(Interval x, Interval y, Interval z)
{
  Interval result(0);
  Interval dx, dy, dz, distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx)+sqr(dy)+sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dz.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * ( -30*r2*r + 60*r2 - 30*r ) * dz / rsq;
        }
    }
  return result;
}

Interval 
MetaBalls::Fxx(Interval x, Interval y, Interval z) 
{
  Interval result(0);
  Interval dx, dy, dz, distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx)+sqr(dy)+sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dx.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * ( (-90*r2+120*r-30)*sqr(dx) / (rsq*r)
                                    + (-30*r2*r+60*r2-30*r) ) / rsq;
        }
    }
  return result;
}

Interval 
MetaBalls::Fyy(Interval x, Interval y, Interval z) 
{
  Interval result(0);
  Interval dx, dy, dz, distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx)+sqr(dy)+sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dy.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * ( (-90*r2+120*r-30)*sqr(dy) / (rsq*r)
                                    + (-30*r2*r+60*r2-30*r) ) / rsq;
        }
    }
  return result;
}

Interval 
MetaBalls::Fzz(Interval x, Interval y, Interval z) 
{
  Interval result(0);
  Interval dx, dy, dz, distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx)+sqr(dy)+sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dz.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * ( (-90*r2+120*r-30)*sqr(dz) / (rsq*r)
                                    + (-30*r2*r+60*r2-30*r) ) / rsq;
        }
    }
  return result;
}

Interval 
MetaBalls::Fxy(Interval x, Interval y, Interval z) 
{
  Interval result(0);
  Interval dx, dy, dz;
  Interval distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++)
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx) + sqr(dy) + sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dx.intersect( Interval(-mballs[i].r, mballs[i].r) );
          dy.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * (-90*r2+120*r-30) * dy * dx / (rsq*rsq*r);
        }
    }
  return result;
}

Interval 
MetaBalls::Fxz(Interval x, Interval y, Interval z) 
{
  Interval result(0);
  Interval dx, dy, dz;
  Interval distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++)
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx) + sqr(dy) + sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dx.intersect( Interval(-mballs[i].r, mballs[i].r) );
          dz.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * (-90*r2+120*r-30) * dz * dx / (rsq*rsq*r);
        }
    }
  return result;
}

Interval 
MetaBalls::Fyz(Interval x, Interval y, Interval z) 
{
  Interval result(0);
  Interval dx, dy, dz;
  Interval distsq;
  double rsq;

  for ( int i=0; i<mballs.size(); i++)
    {
      dx = x - mballs[i].x;
      dy = y - mballs[i].y;
      dz = z - mballs[i].z;
      distsq = sqr(dx) + sqr(dy) + sqr(dz);
      rsq = mballs[i].rsq;
      if ( distsq.inf() < rsq )
        {
          distsq.intersect( Interval(0, rsq) );
          Interval r2 = distsq / rsq;
          Interval r = sqrt( r2 );
          dy.intersect( Interval(-mballs[i].r, mballs[i].r) );
          dz.intersect( Interval(-mballs[i].r, mballs[i].r) );
          result += mballs[i].w * (-90*r2+120*r-30) * dz * dy / (rsq*rsq*r);
        }
    }
  return result;
}

//----------------------------------------------------------------------------
double 
MetaBalls::MinX()
{
  return minx;
}

double 
MetaBalls::MaxX()
{
  return maxx;
}

double 
MetaBalls::MinY()
{
  return miny;
}

double 
MetaBalls::MaxY()
{
  return maxy;
}

double 
MetaBalls::MinZ()
{
  return minz;
}

double 
MetaBalls::MaxZ()
{
  return maxz;
}

//----------------------------------------------------------------------------
// Functions to modify the metaballs
//
int
MetaBalls::NumBalls()
{
  return mballs.size();
}

void
MetaBalls::AddBall( double x, double y, double z, double r, double w )
{
  mballs.push_back( MBall(x,y,z,r,w) );
  if ( x-r < minx ) minx = x-r;
  if ( x+r > maxx ) maxx = x+r;
  if ( y-r < miny ) miny = y-r;
  if ( y+r > maxy ) maxy = y+r;
  if ( z-r < minz ) minz = z-r;
  if ( z+r > maxz ) maxz = z+r;
}

void 
MetaBalls::SetBallCentre( int i, double x, double y, double z )
{
  mballs[i].x = x;
  mballs[i].y = y;
  mballs[i].z = z;
  ComputeBounds();
}

void 
MetaBalls::SetBallRadius( int i, double r )
{
  mballs[i].r = r;
  mballs[i].rsq = r * r;
}

void 
MetaBalls::SetBallWeight( int i, double w )
{
  mballs[i].w = w;
}

void 
MetaBalls::SetThreshold( double t )
{
  threshold = t;
}

void 
MetaBalls::GetBallCentre( int i, double &x, double &y, double &z )
{
  x = mballs[i].x;
  y = mballs[i].y;
  z = mballs[i].z;
}

double
MetaBalls::GetBallRadius( int i )
{
  return mballs[i].r;
}

double 
MetaBalls::GetBallWeight( int i )
{
  return mballs[i].w;
}

double
MetaBalls::GetThreshold()
{
  return threshold;
}

//----------------------------------------------------------------------------
// Internal functions
void 
MetaBalls::ComputeDistances( double _x, double _y, double _z )
{
  double dx, dy, dz, square;

  if ( x==_x && y==_y && z==_z )
    return;
  
  for ( int i=0; i<mballs.size(); i++ )
    {
      dx = _x - mballs[i].x;
      dy = _y - mballs[i].y;
      dz = _z - mballs[i].z;
      square = dx*dx + dy*dy + dz*dz;
      mballs[i].distsq = square;
      mballs[i].dist = sqrt( square );
    }
  x = _x;
  y = _y;
  z = _z;
}

void
MetaBalls::ComputeBounds()
{
  if ( mballs.size() == 0 )
    {
      minx = maxx = 0.0;
      miny = maxy = 0.0;
      minz = maxz = 0.0;
    }
  else
    {
      minx = maxx = mballs[0].x;
      miny = maxy = mballs[0].y;
      minz = maxz = mballs[0].z;
      
      for ( int i=0; i<mballs.size(); i++ )
        {
          if ( mballs[i].x - mballs[i].r < minx ) 
            minx = mballs[i].x - mballs[i].r;
          if ( mballs[i].x + mballs[i].r > maxx ) 
            maxx = mballs[i].x + mballs[i].r;
          if ( mballs[i].y - mballs[i].r < miny ) 
            miny = mballs[i].y - mballs[i].r;
          if ( mballs[i].y + mballs[i].r > maxy ) 
            maxy = mballs[i].y + mballs[i].r;
          if ( mballs[i].z - mballs[i].r < minz ) 
            minz = mballs[i].z - mballs[i].r;
          if ( mballs[i].z + mballs[i].r > maxz ) 
            maxz = mballs[i].z + mballs[i].r;
        }
    }
}

