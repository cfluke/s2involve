#ifndef __ROUTINES_H__
#define __ROUTINES_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "s2plot.h"


typedef struct {
   fitsfile *fptr;                      /* Pointer to FITS file */
   long naxes[3];
   float crv[3], crp[3], cde[3];
   float dmin, dmax;
   float low, high;
   float rdmin, rdmax;
   float obsfreq;
   char **label;

   float ***array;

   XYZ min, max;
   XYZ rflag;
   XYZ range, mp;
   long axmax;
   XYZ vp;
} FITSCube;


typedef struct {
   int tx, ty, tz;
   int sx, sy, sz;
   int fx, fy, fz;
} ScaleAxis;



float ***initVolume(int nx, int ny, int nz, float val);

int errorFITS(int status);
FITSCube readFITScubeHeader(char *fname, int debug);
void readFITScube(FITSCube *cube, int debug);
FITSCube copyCubeHeader(FITSCube c);
void averageCube(FITSCube *avg, FITSCube c, ScaleAxis sa);
void extractCube(FITSCube *avg, FITSCube c, ScaleAxis sa, int mid[3], int delta[3]);
void averExtractCube(FITSCube *avg, FITSCube c, ScaleAxis sa, XYZ roi0, XYZ roi1);
ScaleAxis setTarget(int target, long naxes[3]);





#endif
