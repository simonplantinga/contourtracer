//############################################################################
#include "solver.h"

//############################################################################
// Solve F=0, Fz=0, (Fy=0 OR Fzz=0)
void
Solver::Solve3()
{
  Interval ix( surface->MinX(), surface->MaxX() );
  Interval iy( surface->MinY(), surface->MaxY() );
  Interval iz( surface->MinZ(), surface->MaxZ() );

  numSolutions = 0;
  DoSolve3( ix, iy, iz, 1 );    // F=Fz=Fy=0
  DoSolve3( ix, iy, iz, 2 );    // F=Fz=Fzz=0
}

void 
Solver::Solve4()
{
  Interval ix( surface->MinX(), surface->MaxX() );
  Interval iy( surface->MinY(), surface->MaxY() );
  Interval iz( surface->MinZ(), surface->MaxZ() );
  Interval it( -1, 1 );

  numSolutions = 0;
  DoSolve4( ix, iy, iz, it );
}

void
Solver::DoSolve3( Interval x, Interval y, Interval z, int pass )
{
  // is there any room left for more solutions?
  if ( numSolutions >= MAX_SOLUTIONS )
    return;

  // can this box contain a solution?
  if ( !surface->F(x, y, z).contains(0.0) )
    return;
  if ( !surface->Fz(x, y, z).contains(0.0) )
    return;
  if ( pass == 1 )
    {
      if ( !surface->Fy(x, y, z).contains(0.0) )
        return;
    }
  else
    {
      if ( !surface->Fzz(x, y, z).contains(0.0) )
        return;
    }
  
  // box large enough to split?
  if ( (x.diam() > eps) || (y.diam() > eps) || (z.diam() > eps) )
    {
      if ( ( x.diam() > y.diam() ) && ( x.diam() > z.diam() ) )
        {
          DoSolve3( Interval( x.inf(), x.mid() ), y, z, pass );
          DoSolve3( Interval( x.mid(), x.sup() ), y, z, pass );
        }
      else if ( y.diam() > z.diam() )
        {
          DoSolve3( x, Interval( y.inf(), y.mid() ), z, pass );
          DoSolve3( x, Interval( y.mid(), y.sup() ), z, pass );
        }
      else
        {
          DoSolve3( x, y, Interval( z.inf(), z.mid() ), pass );
          DoSolve3( x, y, Interval( z.mid(), z.sup() ), pass );
        }
    }
  else
    {
      // test if this is a solution
      int a;
      
      a = FSign( surface->F( x.inf(), y.inf(), z.inf() ) );
      if (   a == FSign( surface->F( x.inf(), y.inf(), z.sup() ) )
          && a == FSign( surface->F( x.inf(), y.sup(), z.inf() ) )
          && a == FSign( surface->F( x.inf(), y.sup(), z.sup() ) )
          && a == FSign( surface->F( x.sup(), y.inf(), z.inf() ) )
          && a == FSign( surface->F( x.sup(), y.inf(), z.sup() ) )
          && a == FSign( surface->F( x.sup(), y.sup(), z.inf() ) )
          && a == FSign( surface->F( x.sup(), y.sup(), z.sup() ) ) )
        return;

      a = FSign( surface->Fz( x.inf(), y.inf(), z.inf() ) );
      if (   a == FSign( surface->Fz( x.inf(), y.inf(), z.sup() ) )
          && a == FSign( surface->Fz( x.inf(), y.sup(), z.inf() ) )
          && a == FSign( surface->Fz( x.inf(), y.sup(), z.sup() ) )
          && a == FSign( surface->Fz( x.sup(), y.inf(), z.inf() ) )
          && a == FSign( surface->Fz( x.sup(), y.inf(), z.sup() ) )
          && a == FSign( surface->Fz( x.sup(), y.sup(), z.inf() ) )
          && a == FSign( surface->Fz( x.sup(), y.sup(), z.sup() ) ) )
        return;

      if ( pass == 1 )
        {
          a = FSign( surface->Fy( x.inf(), y.inf(), z.inf() ) );
          if (   a == FSign( surface->Fy( x.inf(), y.inf(), z.sup() ) )
              && a == FSign( surface->Fy( x.inf(), y.sup(), z.inf() ) )
              && a == FSign( surface->Fy( x.inf(), y.sup(), z.sup() ) )
              && a == FSign( surface->Fy( x.sup(), y.inf(), z.inf() ) )
              && a == FSign( surface->Fy( x.sup(), y.inf(), z.sup() ) )
              && a == FSign( surface->Fy( x.sup(), y.sup(), z.inf() ) )
              && a == FSign( surface->Fy( x.sup(), y.sup(), z.sup() ) ) )
            return;
        }
      else
        {
          a = FSign( surface->Fzz( x.inf(), y.inf(), z.inf() ) );
          if (   a == FSign( surface->Fzz( x.inf(), y.inf(), z.sup() ) )
              && a == FSign( surface->Fzz( x.inf(), y.sup(), z.inf() ) )
              && a == FSign( surface->Fzz( x.inf(), y.sup(), z.sup() ) )
              && a == FSign( surface->Fzz( x.sup(), y.inf(), z.inf() ) )
              && a == FSign( surface->Fzz( x.sup(), y.inf(), z.sup() ) )
              && a == FSign( surface->Fzz( x.sup(), y.sup(), z.inf() ) )
              && a == FSign( surface->Fzz( x.sup(), y.sup(), z.sup() ) ) )
            return;
        }

      // add solution to list

      // 1. move point to cg
      double px = x.mid();
      double py = y.mid();
      double pz = z.mid();
      while ( fabs(surface->F(px,py,pz)) > CLOSETOSURFACE
            ||fabs(surface->Fz(px,py,pz)) > CLOSETOSURFACE )
        {
          if ( fabs(surface->F(px,py,pz)) > CLOSETOSURFACE )
            {
              double f = surface->F(px,py,pz);
              double fx = surface->Fx(px,py,pz);
              double fy = surface->Fy(px,py,pz);
              double fz = surface->Fz(px,py,pz);
              double gs = fx*fx + fy*fy + fz*fz;
              px -= f * fx / gs;
              py -= f * fy / gs;
              pz -= f * fz / gs;
            }
          if ( fabs(surface->Fz(px,py,pz)) > CLOSETOSURFACE )
            {
              double fz = surface->Fz(px,py,pz);
              double fxz = surface->Fxz(px,py,pz);
              double fyz = surface->Fyz(px,py,pz);
              double fzz = surface->Fzz(px,py,pz);
              double gs = fxz*fxz + fyz*fyz + fzz*fzz;
              px -= fz * fxz / gs;
              py -= fz * fyz / gs;
              pz -= fz * fzz / gs;
            }

        }
      // 2. test if point already found
      for ( int i=0; i<numSolutions; i++ )
        {
          if ( (sx[i]-px)*(sx[i]-px)
              +(sy[i]-py)*(sy[i]-py)
              +(sz[i]-pz)*(sz[i]-pz) < eps*eps )
            {
              return;
            }
        }
      sx[numSolutions] = px;
      sy[numSolutions] = py;
      sz[numSolutions] = pz;
      numSolutions++;
    }
}

void
Solver::DoSolve4( Interval x, Interval y, Interval z, Interval t )
{
  Interval f, fx, fy, fz, fxz, fyz, fzz, grad;

  // can we store more solutions?
  if ( numSolutions >= MAX_SOLUTIONS )
    return;

  f = surface->F(x, y, z);

  // can box contain a solution?

  // 1. F(x,y,z) = 0
  if ( !f.contains(0.0) )
    return;

  fx = surface->Fx(x, y, z);
  fy = surface->Fy(x, y, z);
  fz = surface->Fz(x, y, z);

  grad = sqrt( sqr(fx) + sqr(fy) + sqr(fz) );
  
  // 2. Fz(x,y,z) = theta * || nabla F ||
  if ( !(fz - t * grad).contains( 0.0 ) )
    return;

  fxz = surface->Fxz(x, y, z);
  fyz = surface->Fyz(x, y, z);

  // 3. Fx * Fyz = Fy * Fxz
  if ( !( fx * fyz - fy * fxz ).contains( 0.0 ) )
    return;

  fzz = surface->Fzz(x, y, z);

  // 4. Fy * Fzz = theta * Fyz * || nabla F ||
  if ( !( fy * fzz - t * fyz * grad ).contains( 0.0 ) )
    return;

  // 5. Fx * Fzz = theta * Fxz * || nabla F ||
  if ( !( fx * fzz - t * fxz * grad ).contains( 0.0 ) )
    return;

  // box large enough to split?
  if ( (x.diam()>eps) || (y.diam()>eps) || (z.diam()>eps) || (t.diam()>eps) )
    {
      if ( (x.diam()>y.diam()) && (x.diam()>z.diam()) && (x.diam()>t.diam()) )
        {
          DoSolve4( Interval( x.inf(), x.mid() ), y, z, t );
          DoSolve4( Interval( x.mid(), x.sup() ), y, z, t );
        }
      else if ( (y.diam()>z.diam()) && (y.diam()>t.diam()) )
        {
          DoSolve4( x, Interval( y.inf(), y.mid() ), z, t );
          DoSolve4( x, Interval( y.mid(), y.sup() ), z, t );
        }
      else if ( z.diam() > t.diam() )
        {
          DoSolve4( x, y, Interval( z.inf(), z.mid() ), t );
          DoSolve4( x, y, Interval( z.mid(), z.sup() ), t );
        }
      else
        {
          DoSolve4( x, y, z, Interval( t.inf(), t.mid() ) );
          DoSolve4( x, y, z, Interval( t.mid(), t.sup() ) );
        }
    }
  else
    {
      // test for a solution

        // IS THIS REALLY A SOLUTION?

      // add solution to list
      sx[numSolutions] = x.mid();
      sy[numSolutions] = y.mid();
      sz[numSolutions] = z.mid();
      st[numSolutions] = t.mid();
      numSolutions++;
    }
}

int
Solver::FSign( double d )
{
  if ( d < 0 )
    return -1;
  else if ( d > 0 )
    return 1;
  return 0;
}

