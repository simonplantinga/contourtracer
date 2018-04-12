//############################################################################
#ifndef SOLVER_H
#define SOLVER_H

#include "surface.h"

#define MAX_SOLUTIONS   1000
#define CLOSETOSURFACE 0.0000001

class Solver
{
public:
  Solver( Surface *s ) { surface = s; numSolutions = 0; eps = 0.01; }
  void Solve3();
  void Solve4();
  int GetNumSolutions() { return numSolutions; }
  void GetSolution3( int i, double &x, double &y, double &z )
    { x = sx[i]; y = sy[i]; z = sz[i]; }
  void GetSolution4( int i, double &x, double &y, double &z, double &t )
    { x = sx[i]; y = sy[i]; z = sz[i]; t = st[i]; }
  void SetEpsilon( double e ) { eps = e; }
  double GetEpsilon() { return eps; }
private:
  void DoSolve3( Interval x, Interval y, Interval z, int pass );
  void DoSolve4( Interval x, Interval y, Interval z, Interval t );
  int  FSign( double d );

  Surface *surface;
  int numSolutions;
  double sx[MAX_SOLUTIONS];
  double sy[MAX_SOLUTIONS];
  double sz[MAX_SOLUTIONS];
  double st[MAX_SOLUTIONS];
  double eps;
};

#endif  // SOLVER_H

