//############################################################################
//  main.cpp

#define EPS_OUTPUT
#define POV_OUTPUT

#include <qapplication.h>
#include <qwidget.h>
#include <qpainter.h>

#include "metaballs.h"
#include "solver.h"
#include "tracer.h"

#ifdef EPS_OUTPUT
#include <stdio.h>
#endif

#ifdef POV_OUTPUT
#include <stdio.h>
#endif

#define EPSILON 0.05
//#define ROT_ANGLE 0.2
#define ROT_ANGLE 0.15

//############################################################################
//  view widget
class QuickView : public QWidget
{
public:
  QuickView( QWidget *parent = 0, const char *name = 0 )
      : QWidget( parent, name ) { mousePressed = false; }
  void setSurface( Surface *s ) { surface = s; }
protected:
  void paintEvent( QPaintEvent * );
  void mousePressEvent( QMouseEvent * );
  void mouseMoveEvent( QMouseEvent * );
  void mouseReleaseEvent( QMouseEvent * );
private:
  void realToScreen( double x, double y, double z, int &sx, int &sy );
  void screenToReal( int sx, int sy, double &x, double &y, double &z );
  void plotSolutions( QPainter *p );
  void plotContours( QPainter *p );
public:
  void RotateXPos();
  void RotateXNeg();
  void RotateYPos();
  void RotateYNeg();
private:
  Surface *surface;
  bool mousePressed;
  QPoint mousePos;
};

void
QuickView::realToScreen( double x, double y, double z, int &sx, int &sy )
{
  sx = (int)(width() * (x-surface->MinX())/(surface->MaxX()-surface->MinX()));
  sy = (int)(height()* (y-surface->MinY())/(surface->MaxY()-surface->MinY()));
}

void
QuickView::screenToReal( int sx, int sy, double &x, double &y, double &z )
{
  x=surface->MinX()+((double)sx/width())*(surface->MaxX()-surface->MinX());
  y=surface->MinY()+((double)sy/height())*(surface->MaxY()-surface->MinY());
  z=surface->MinZ();
}

void 
QuickView::paintEvent( QPaintEvent *e )
{
  QPainter p( this );

  // plot solutions to equations
  plotSolutions( &p );
}

void
QuickView::mousePressEvent( QMouseEvent *e )
{
  if ( e->button() != LeftButton )
    return;
  mousePressed = true;
  mousePos = e->pos();
}

void
QuickView::mouseMoveEvent( QMouseEvent *e )
{
  if ( !mousePressed )
    return;
  QPoint pnt = e->pos();

  if ( pnt.x() > mousePos.x() )
    RotateYPos();
  if ( pnt.x() < mousePos.x() )
    RotateYNeg();
  if ( pnt.y() > mousePos.y() )
    RotateXNeg();
  if ( pnt.y() < mousePos.y() )
    RotateXPos();
  mousePos = pnt;
  repaint();
}

void
QuickView::mouseReleaseEvent( QMouseEvent *e )
{
  if ( e->button() == LeftButton )
    mousePressed = false;
}

void
QuickView::plotSolutions( QPainter *p )
{
  Solver solver( surface );
  Tracer tracer;
  double x, y, z;
  int sx, sy;

  // solve equations
  solver.SetEpsilon( EPSILON );
  solver.Solve3();

  // plot solutions to equations
  for ( int i=0; i<solver.GetNumSolutions(); i++ )
    {
      solver.GetSolution3( i, x, y, z );
      realToScreen( x, y, z, sx, sy );
      p->drawRect( sx-2, sy-2, 5, 5 );
    }



  // compute contours
//  for ( int i=20; i<21; i++ )
  for ( int i=0; i<solver.GetNumSolutions(); i++ )
    {
      tracer.setEpsilon( EPSILON );
//      tracer.setEpsilon(1);
      solver.GetSolution3( i, x, y, z );
      
      Point pt(x,y,z);
      tracer.addContour( surface, pt );
    }
  
  cout << "Number of contours: " << tracer.getNumContours() << endl;

  // plot contours
  for ( int i=0; i<tracer.getNumContours(); i++ )
    {
      Point pt;
      double x, y, z;
      double xn, yn, zn;
      tracer.getContourPoint( i, 0, pt );
      xn = pt.X(); yn = pt.Y(); zn = pt.Z();
      for ( int c=1; c<tracer.getContourSize( i ); c++ )
        {
          x = xn; y = yn; z = zn;
          tracer.getContourPoint( i, c, pt );
          xn = pt.X(); yn = pt.Y(); zn = pt.Z();
          int x1, y1, x2, y2;
          realToScreen( x, y, z, x1, y1 );
          realToScreen( xn, yn, zn, x2, y2 );
          p->drawLine( x1, y1, x2, y2 );
        }
    }

//  cout << tracer.getNumContours() << " components\n";
#ifdef EPS_OUTPUT
{
  FILE *fp = fopen("contour.eps","wt");
  fprintf(fp, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(fp, "%%%%CreationDate: 00/00/00\n");
  fprintf(fp, "%%%%BoundingBox: 0 0 1000 1000\n");
  fprintf(fp, "%%%%EndComments\n");
  // plot solutions
  for ( int i=0; i<solver.GetNumSolutions(); i++ )
    {
      solver.GetSolution3( i, x, y, z );
      sx = (int)(1000 * (x-surface->MinX())/(surface->MaxX()-surface->MinX()));
      sy = (int)(1000 * (y-surface->MinY())/(surface->MaxY()-surface->MinY()));
      fprintf(fp, "newpath\n");      
      fprintf(fp, "%i %i moveto\n", sx-5, sy-5);
      fprintf(fp, "%i %i lineto\n", sx+5, sy-5);
      fprintf(fp, "%i %i lineto\n", sx+5, sy+5);
      fprintf(fp, "%i %i lineto\n", sx-5, sy+5);
      fprintf(fp, "closepath\nstroke\n");
    }
  // plot contours
  for ( int i=0; i<tracer.getNumContours(); i++ )
    {
      Point pt;
      double x, y, z;
      double xn, yn, zn;
      int x1, y1, x2, y2;
      tracer.getContourPoint( i, 0, pt );
      xn = pt.X(); yn = pt.Y(); zn = pt.Z();
      x1 = (int)(1000 
              * (xn-surface->MinX())/(surface->MaxX()-surface->MinX()));
      y1 = (int)(1000 
              * (yn-surface->MinY())/(surface->MaxY()-surface->MinY()));
      fprintf( fp, "newpath\n%i %i moveto\n", x1, y1 );

      for ( int c=1; c<tracer.getContourSize( i ); c++ )
        {
          x = xn; y = yn; z = zn;
          tracer.getContourPoint( i, c, pt );
          xn = pt.X(); yn = pt.Y(); zn = pt.Z();
          x1 = (int)(1000 
              * (xn-surface->MinX())/(surface->MaxX()-surface->MinX()));
          y1 = (int)(1000 
              * (yn-surface->MinY())/(surface->MaxY()-surface->MinY()));
          x2 = (int)(1000 
                  * (xn-surface->MinX())/(surface->MaxX()-surface->MinX()));
          y2 = (int)(1000 
                  * (yn-surface->MinY())/(surface->MaxY()-surface->MinY()));
          p->drawLine( x1, y1, x2, y2 );
          fprintf( fp, "%i %i lineto\n", x2, y2 );
        }
      fprintf( fp, "stroke\n" );
    }
  fprintf(fp, "showpage\n");
  fclose(fp);  
}
#endif

#ifdef POV_OUTPUT
{
  FILE *fp = fopen("contour.inc","wt");
  Point Pt;
  for (int i=0; i<tracer.getNumContours(); i++)
    {
      for (int c=1; c<tracer.getContourSize(i); c++)
        {
          tracer.getContourPoint(i, c, Pt);
          fprintf(fp, "sphere {<%f,%f,%f> RADIUS}\n",Pt.X(),Pt.Y(),Pt.Z());
        }
    }
  fclose(fp);
  exit(0);
}
#endif
#ifdef EPS_OUTPUT
  exit(0);
#endif
}

void
QuickView::RotateXPos()
{
  MetaBalls *m = (MetaBalls *)surface;
  double x, y, z;
  for ( int i=0; i<m->NumBalls(); i++ )
    {
      m->GetBallCentre( i, x, y, z );
      y = cos(ROT_ANGLE) * y - sin(ROT_ANGLE) * z;
      z = sin(ROT_ANGLE) * y + cos(ROT_ANGLE) * z;
      m->SetBallCentre( i, x, y, z );
    }
}

void
QuickView::RotateXNeg()
{
  MetaBalls *m = (MetaBalls *)surface;
  double x, y, z;
  for ( int i=0; i<m->NumBalls(); i++ )
    {
      m->GetBallCentre( i, x, y, z );
      y = cos(-ROT_ANGLE) * y - sin(-ROT_ANGLE) * z;
      z = sin(-ROT_ANGLE) * y + cos(-ROT_ANGLE) * z;
      m->SetBallCentre( i, x, y, z );
    }
}

void
QuickView::RotateYPos()
{
  MetaBalls *m = (MetaBalls *)surface;
  double x, y, z;
  for ( int i=0; i<m->NumBalls(); i++ )
    {
      m->GetBallCentre( i, x, y, z );
      x = cos(ROT_ANGLE) * x - sin(ROT_ANGLE) * z;
      z = sin(ROT_ANGLE) * x + cos(ROT_ANGLE) * z;
      m->SetBallCentre( i, x, y, z );
    }
}

void
QuickView::RotateYNeg()
{
  MetaBalls *m = (MetaBalls *)surface;
  double x, y, z;
  for ( int i=0; i<m->NumBalls(); i++ )
    {
      m->GetBallCentre( i, x, y, z );
      x = cos(-ROT_ANGLE) * x - sin(-ROT_ANGLE) * z;
      z = sin(-ROT_ANGLE) * x + cos(-ROT_ANGLE) * z;
      m->SetBallCentre( i, x, y, z );
    }
}

//############################################################################

int
main( int argc, char *argv[] )
{
  MetaBalls metaballs;

  QApplication app( argc, argv );

  QuickView qv;
  
  metaballs.AddBall( 1, 1, 1, 1.8, 1 );
  metaballs.AddBall( 1, 1,-1, 1.8, 1 );
  metaballs.AddBall( 1,-1, 1, 1.8, 1 );
  metaballs.AddBall( 1,-1,-1, 1.8, 1 );
  metaballs.AddBall(-1, 1, 1, 1.8, 1 );
  metaballs.AddBall(-1, 1,-1, 1.8, 1 );
  metaballs.AddBall(-1,-1, 1, 1.8, 1 );
  metaballs.AddBall(-1,-1,-1, 1.8, 1 );

  metaballs.SetThreshold( 0.5 );

  qv.setSurface( &metaballs );
  qv.RotateXPos();
  qv.RotateXPos();
//  qv.RotateXPos();
//  qv.RotateYPos();
  qv.RotateYPos();
  qv.setGeometry( 100, 100, 600, 600 );
  app.setMainWidget( &qv );
  qv.show();
  return app.exec();
}

