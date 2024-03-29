/* -*- c++ -*- ----------------------------------------------------------
   Modified by R
   Reflectie particle method: 10.1103/PhysRevE.57.7259
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/rpm,FixWallRPM)

#else

#ifndef LMP_FIX_WALL_RPM_H
#define LMP_FIX_WALL_RPM_H

#include "fix.h"
#include "random_park.h"

namespace LAMMPS_NS {

class FixWallRPM : public Fix {
 public:
  FixWallRPM(class LAMMPS *, int, char **);
  virtual ~FixWallRPM();
  int setmask();
  void init();
  void post_integrate();

 protected:
  int nwall;
  int wallwhich[6],wallstyle[6];
  double coord0[6];
  char *varstr[6];
  int varindex[6];
  int varflag;
  double xscale,yscale,zscale;
  double p;
  RanPark rng;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Wall defined twice in fix wall/reflect command

Self-explanatory.

E: Cannot use fix wall/reflect in periodic dimension

Self-explanatory.

E: Cannot use fix wall/reflect zlo/zhi for a 2d simulation

Self-explanatory.

E: Variable name for fix wall/reflect does not exist

Self-explanatory.

E: Variable for fix wall/reflect is invalid style

Only equal-style variables can be used.

W: Should not allow rigid bodies to bounce off relecting walls

LAMMPS allows this, but their dynamics are not computed correctly.

*/
