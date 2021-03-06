/* ParticleExample
 no recenter on particle files
 */

#include <slsimlib.h>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
#include "gridmap.h"
#include <math.h>
#include "particle_halo.h"
#include "utilities.h"
#include "concave_hull.h"
#include "uniform_lens.h"
//#include "gadget.hh"

using namespace std;

int main(int arg,char **argv){
  
    time_t t;
    time(&t);
  
  COSMOLOGY cosmo(Planck18);
  const long Nstars = 5000;//00;
  const double zl = 0.5;  // redshift of lens
  const double zs = 1.0;  // redshift of source
  const int Ninit = 1000;  // initioal size of grid of rays
  const double kappa_star = 0.5;  // optical depth
  const double mass = 1;         /// star mass in solar masses
  
  
  double Re = sqrt( 1.0 / cosmo.SigmaCrit(zl,zs)/ PI );  // Einstein radius of stars
  double Rphys = Re * sqrt( mass * Nstars / kappa_star);      // radius of region to but stars in
  
  //long seed = -11920;
  long seed = time(&t);
  Utilities::RandomNumbers_NR ran(seed);
  
  std::vector<StarType> stars(Nstars);
  Point_2d center = {0,0};
  Point_2d rotation = {0,0};   // no rotation of stars
  for(StarType &s : stars){
    double r = Rphys*sqrt(ran()) ;
    double theta = 2*PI*ran();
    s.x[0] = r*cos(theta) + center[0];
    s.x[1] = r*sin(theta) + center[1];
    s.x[2] = 0;  // in a plane perpendicular to the line of sight
    s.Mass = mass; 
  }
  
  LensHaloParticles<StarType> lens_halo(stars,zl,cosmo,rotation,false,0,false);
  
  // put int a compinsating negative density
  LensHaloKappaDisk disk(zl,-kappa_star*cosmo.SigmaCrit(zl, zs),Rphys,Point_2d(0,0), cosmo);
  
  
  Lens lens(&seed,zs,cosmo);
  lens.moveinMainHalo(lens_halo,true);
  lens.insertMainHalo(disk,false);
  
  //GridMap grid(&lens,Ninit,center.x, Ninit * Re / cosmo.angDist(zl) / 10 );
  GridMap grid(&lens,Ninit,center.x, 3*Rphys/cosmo.angDist(zl) );

  // This makes a nice image of the inverse magnification
  // if you cut off the pixel values < -500 and use a histogram
  // scaling
  PixelMap map = grid.writePixelMapUniform(INVMAG);
  map.printFITS("!inverse_mag.fits");
  grid.writeFitsUniform(KAPPA,"!kappa.fits");
  grid.writeFitsUniform(GAMMA,"!gamma.fits");
  grid.writeFitsUniform(GAMMA1,"!gamma1.fits");
  grid.writeFitsUniform(ALPHA1,"!alpha1.fits");
  grid.writeFitsUniform(ALPHA,"!alpha.fits");

  double Dl = cosmo.angDist(zl);
  
  // make a source with a Gaussian profile
  SourceGaussian source(center,  Re / Dl,zl);

  // print an image of the source
  double brightness = grid.RefreshSurfaceBrightnesses(&source);
  map = grid.getPixelMap(1);
  
  map.printFITS("!source_image.fits");
  
  // Let's add a little background shear and convergence
  Point_2d shear(0.01,0);
  double kappa = 0.001;
  LensHaloUniform background(zl,kappa,shear,cosmo);
  lens.insertMainHalo(background,false);
  
  // make a light curve
  double velocity = 1000*kmpersecTOmpcperday/Dl; // relative speed of source
  Point_2d direction(1,1);  // direction the source moves in
  direction.unitize();
  
  for(int day = 0 ; day < 200 ; ++day){
    source.setTheta( center + direction*velocity*day);
    brightness = grid.RefreshSurfaceBrightnesses(&source);
    std::cout << day << " " << brightness << std::endl;
  }
}



