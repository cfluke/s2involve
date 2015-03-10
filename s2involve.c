#include <stdio.h>
#include <stdlib.h>
#include "routines.h" 
#include "s2plot.h"

#define EXITERROR 1

#define ALPHA 0.7
#define MAXHANDLE 4
#define LABELLEN 32
#define _CMAP (char *)"iron"
#define _STRETCH 1.0
#define _CFIRST 100
#define _CMAX 255
#define _DMIN -1E9		/* A very negative number */
#define _DMAX +1E4		/* A very positive number */

#define EMPTY  0
#define TEX3D  1
#define TEX2D  2
#define TEXPDF 3

typedef struct {
   float amin, amax;           /* Alpha channel minimum and maximum */
   float dmin, dmax;           /* Current plotted data minimum and maximum */
   int c1, cr;                 /* Base colour index and colour index range */
   char cmap[64];              /* Current colour map name */
} VisParameters;

typedef struct {
   int N;
   float *bin;			/* Needs to be allocated */
   float dmin, dmax;		/* Data min and data max */
   float bmin, bmax;		/* Min/Max data values within the bins */
   float delta;			/* Width of an individual bin */
   float scale;			/* Unit scaling factor */
   float sum;			/* Sum over all bin values */
   int visible;
   int log;
} Histogram;

typedef struct {
   float min, max;
   float mid, range;
   double mean;
   double sd;
} DataValues;

typedef struct {
   long Naxis;					/* NAXIS */
   long naxis[4];				/* NAXISn: n=1..4 */
#ifdef NOTNOW
   char object[128];				/* OBJECT */
   float bscale;				/* BSCALE */
   float bzero;					/* BZERO */
#endif
} Header;


typedef struct {
   unsigned int tid[3];
   long N[3];
   float ***vr;
   VisParameters vp;
} CubeTexture; 

typedef struct {
   long N[3];
   float ***vol;
   DataValues dv;
} RawCube;

typedef struct {
   int input[3];
   int target[3];
   int scale[3];
   int f[3];
} Axes;

typdef struct {
   char **axis;
} Labels;

typedef struct {
   int pid;
   int ctype;
   int cubeID;
   XYZ min, max; 
   CubeTexture ct;
   Histogram h;
   int handles[MAXHANDLE];
   DataValues dv;
   Axes axes;
   Labels labels;
   float stretch;
   long target[3];
} Panel;

 
typedef struct {
   int Nx, Ny;			/* Number of panels */
   Panel **panel;
   int Ncube;
   char cmap[64];
   Header  *header;
   RawCube *cube;
} Configuration;

/* Global Variables */

Configuration config;
int debug = 1;

/* Function prototypes */
Histogram initHistogram(int Nbin);
DataValues characteriseData(float ***vol, long N[3]);


void ds2dvrXXX(int, int, int);
int kcb(unsigned char *key);
void cb(double *t, int *kc, int *value);


Panel initPanel(int cubeID, int ctype, int i, int j, int Nx, int Ny, int Nbin)
{
   Panel p;
   p.cubeID = cubeID;
   p.ctype = ctype;
   p.min.x = (float)i/(float)Nx;
   p.max.x = (float)(i+1)/(float)Nx;
   p.min.y = 1.0-(float)(j+1)/(float)Ny;
   p.max.y = 1.0-(float)j/(float)Ny;
   p.min.z = p.max.z = 0.0;
   if ((i==0) && (j==0)) {
      p.pid = 0; 
      xs2mp(0, p.min.x,p.min.y,p.max.x,p.max.y);
   } else {
      p.pid = xs2ap(p.min.x,p.min.y,p.max.x,p.max.y);
   }
   p.h  = initHistogram(Nbin);
   for (i=0;i<MAXHANDLE;i++) p.handles[i] = 0;
 
   for (i=0;i<3;i++) {
      p.axes.input[i]  = 1;
      p.axes.target[i] = 1;
      p.axes.scale[i]  = 1;
      p.axes.f[i]      = 1;
   }

   p.labels.axis = (char **)calloc(3, sizeof(char *));
   for (i=0;i<3;i++) {
      p.labels.axis[i] = (char *)calloc(LABELLEN, sizeof(char));
   }

   xs2cp(p.pid);
   xs2lpc(0, p.pid);
   s2swin(0,1,0,1,0,1);
   COLOUR active = { 0,0,0 };
   COLOUR inactive = { 0.0,0.0,0.0 };
   xs2spp(active, inactive, 2);

   xs2cp(0);
   return p;
}


Configuration readConfiguration(FILE *cfp, int Nbin)
{
   Configuration c;
   char string[128];
   int ctype  = TEX3D;
   int cubeID = 0;
   int i,j;
   float stretch = 1.0;
   char cmap[64];  
   sprintf(cmap,"%s",_CMAP);
   long NN[3]= { 128, 128, 128 };

/* Read the number of columns from the config file */
   fgets(string,128,cfp);

   while (strncmp(string,"END",3) != 0) {
      switch (string[0]) {
         case 'a' : /* Axes information */
                    if (sscanf(string,"axi %ld %ld %ld", &NN[0], &NN[1], &NN[2]) != 3) {
                       fprintf(stderr,"Missing axis values in config.txt\n");
                       fclose(cfp);
                       exit(EXITERROR);
   		    }
                    break;
         case 'c' : /* Column information */
                    if (sscanf(string,"col %d", &c.Nx) != 1) {
                       fprintf(stderr,"Missing col value in config.txt\n");
                       fclose(cfp);
                       exit(EXITERROR);
   		    }
                    break;
         case 'm' : /* Colour map information */
                    if (sscanf(string,"map %s", cmap) != 1) {
                       fprintf(stderr,"Missing colour map name in config.txt\n");
                       fclose(cfp);
                       exit(EXITERROR);
   		    }
                    break;
         case 'r' : /* Row information */
                    if (sscanf(string,"row %d", &c.Ny) != 1) {
                       fprintf(stderr,"Missing row value in config.txt\n");
                       fclose(cfp);
                       exit(EXITERROR);
   		    }
                    break;
         case 's' : /* Stretch factor information */
                    if (sscanf(string,"str %f", &stretch) != 1) {
                       fprintf(stderr,"Missing stretch factor value in config.txt\n");
                       fclose(cfp);
                       exit(EXITERROR);
   		    }
                    break;
 
      }
      fgets(string,128,cfp);
   }

/* Allocate memory for panel array */
   c.panel = (Panel **)calloc(c.Nx, sizeof(Panel *));
   for (i=0;i<c.Nx;i++) {
      c.panel[i] = (Panel *)calloc(c.Ny, sizeof(Panel));
      for (j=0;j<c.Ny;j++) {
         c.panel[i][j].pid = 0;
      }
   }

/* CJF: Better to allow multiple colour maps */
/* Copy colour map */
   sprintf(c.cmap,"%s",cmap);

/* Read the configuration information for the each panel */
   char stype[16];
   int *cidx = (int *)calloc(c.Nx*c.Ny,sizeof(int));
   int pidx = 0;

   fgets(string,128,cfp);
   while (!feof(cfp)) {
      if (sscanf(string,"%d %d %s %d",&i, &j, stype, &cubeID) !=4) {
         fprintf(stderr,"Bad format: %s\n",string);
         fclose(cfp);
         exit(EXITERROR);
      }
      if (strncmp(stype,"3DT",3) == 0) { 
         ctype = TEX3D;
      } else if (strncmp(stype,"2DT",3) == 0) { 
         ctype = TEX2D;
      } else if (strncmp(stype,"PDF",3) == 0) { 
         ctype = TEXPDF;
      } else {
         fprintf(stderr,"Unknown texture type: %s\n",stype);
         fclose(cfp);
         exit(EXITERROR);
      }

/* Configuration file uses [1..Nx], [1..Ny] - convert to C-array values */
      i--;
      j--;
      cubeID--;
      c.panel[i][j] = initPanel(cubeID, ctype, i,j, c.Nx,c.Ny, Nbin);
      c.panel[i][j].target[0] = NN[0];
      c.panel[i][j].target[1] = NN[1];
      c.panel[i][j].target[2] = NN[2];
      c.panel[i][j].stretch = stretch;

/* Attach the callbacks to the new panel */
      xs2cp(c.panel[i][j].pid);
      cidx[pidx] = pidx;
      cs2scbx(cb, &cidx[pidx]);		/* Attach dynamic callback to this panel */
      cs2skcb(kcb);			/* Attach keypress callback to this panel */
      xs2cp(0);

      pidx++;
      fgets(string,128,cfp);
   }
   xs2cp(0);
   return c;
}

Configuration initConfiguration(int Nx, int Ny, int Nbin, char *cmap)
{
   Configuration c;
   c.Nx = Nx;
   c.Ny = Ny;
   int ctype = TEX3D;
   int i,j;
   int pidx = 0;
   int *cidx = (int *)calloc(c.Nx*c.Ny,sizeof(int));
   long NN[3] = { 128, 128, 128 }; 

   c.panel = (Panel **)calloc(c.Nx, sizeof(Panel *));
   for (i=0;i<c.Nx;i++) {
      c.panel[i] = (Panel *)calloc(c.Ny, sizeof(Panel));
      for (j=0;j<c.Ny;j++) {
         ctype = i%3;
         c.panel[i][j] = initPanel(i, ctype, i,j, c.Nx,c.Ny, Nbin);
	 c.panel[i][j].target[0] = NN[0];
	 c.panel[i][j].target[1] = NN[1];
	 c.panel[i][j].target[2] = NN[2];
         c.panel[i][j].stretch = _STRETCH;
         xs2cp(c.panel[i][j].pid);
         cidx[pidx] = pidx;
         cs2scbx(cb, &cidx[pidx]);	/* Attach callback to this panel */
         cs2skcb(kcb);			/* Attach keypress callback to this panel */
         pidx++;
      }
   }

   sprintf(c.cmap,"%s",cmap);

   xs2cp(0);
   return c;
}



Histogram initHistogram(int Nbin)
{
   Histogram h;
   h.N    = Nbin;
   h.bin  = (float *)calloc(h.N, sizeof(float));
   h.dmin = 0.0;
   h.dmax = 1.0;
   h.bmin = +1E9;
   h.bmax = -1e9;
   h.scale = 1.0;
   h.delta = (h.dmax-h.dmin)/(float)(h.N-1);
   h.sum   = 0.0;
   h.visible = 0;
   h.log = 0;
   return h;
}

VisParameters initVisParameters(float dmin, float dmax, int c1, int cr, char *cmap)
{
   VisParameters p;
   p.dmin = dmin;
   p.dmax = dmax;
   p.amin = 0.0;
   p.amax = ALPHA;
   sprintf(p.cmap,"%s",cmap);
   p.c1 = c1;
   p.cr = cr;
   return p;
}

CubeTexture initCubeTexture(Axes a, int Nbin)
{
   int i;
   CubeTexture ct;
   for (i=0;i<3;i++) {
     ct.tid[i] = 0;
     ct.N[i] = a.target[i];
   }

   ct.vr = initVolume(ct.N[0],ct.N[1],ct.N[2],-1E99);

/* VisParameters should be initialised outside of this function */

   return ct;
}


void textureVolumeRender(unsigned int texid, long Naxis[3])
{
  // These would be args to the vol renderer once it's built into S2PLOT 
  int adim = Naxis[0];
  int bdim = Naxis[1];
  int cdim = Naxis[2];
  int a1 = 0;
  int a2 = adim-1;
  int b1 = 0;
  int b2 = bdim-1;
  int c1 = 0;
  int c2 = cdim-1;
  float tr[] = {0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 0.};
  tr[1] = 1./(float)(a2);
  tr[6] = 1./(float)(b2);
  tr[11] = 1./(float)(c2);

  char itrans = 's';
  float ialpha = 1.00;

  // these are the vertices of a unit cube
  XYZ unitverts[] = {{0, 0, 0},
		     {1, 0, 0},
		     {0, 1, 0},
		     {1, 1, 0},
		     {0, 0, 1},
		     {1, 0, 1},
		     {0, 1, 1},
		     {1, 1, 1}};
  // the vertices of the data being displayed
  XYZ dataverts[8];
  // the world vertices of the data being displayed (= dataverts * tr)
  XYZ worldverts[8];

  // these are the edges of the cube, made up by joining these unitverts:
  int edges[12][2] = {{0,1}, {0,2}, {1,3}, {2,3},
		      {4,5}, {4,6}, {5,7}, {6,7},
		      {0,4}, {1,5}, {2,6}, {3,7}};

  // 1. get camera position and view direction in world coords
  XYZ campos, upvec, viewdir, right;
  ss2qc(&campos, &upvec, &viewdir, 1);
  Normalise(&viewdir);
  right = CrossProduct(viewdir, upvec);
  
  // 2. find indices of first and last vertices: first is that vertex
  //    which a plane normal to the viewdir crosses, travelling towards
  //    the centre of the cube, from the camera position.
  int i;
  int near_vtx, far_vtx;
  float near_dist, far_dist;
  XYZ tmp;
  float thisdist;  
  near_vtx = far_vtx = -1;
  near_dist = 9e30;
  far_dist = -9e30;
  for (i = 0; i < 8; i++) {
    // calculate this data vertex
    dataverts[i].x = a1 + unitverts[i].x * (a2-a1);
    dataverts[i].y = b1 + unitverts[i].y * (b2-b1);
    dataverts[i].z = c1 + unitverts[i].z * (c2-c1);
    
    // and this world vertex position
    worldverts[i].x = tr[0] + tr[1] * dataverts[i].x 
      + tr[2] * dataverts[i].y + tr[3] * dataverts[i].z;
    worldverts[i].y = tr[4] + tr[5] * dataverts[i].x 
      + tr[6] * dataverts[i].y + tr[7] * dataverts[i].z;
    worldverts[i].z = tr[8] + tr[9] * dataverts[i].x 
      + tr[10] * dataverts[i].y + tr[11] * dataverts[i].z;

    // and now its distance from the camera position measured along the
    // view direction
    tmp = VectorSub(campos, worldverts[i]);
    thisdist = DotProduct(tmp, viewdir);
    if (thisdist < near_dist) {
      near_vtx = i;
      near_dist = thisdist;
    }
    if (thisdist > far_dist) {
      far_vtx = i;
      far_dist = thisdist;
    }

  }


  // 3. step from near distance to far distance, and calculate the 
  //    bounds of each polygon slice (intersection of cube and plane).
  XYZ p1, p2;
  int plidx; // plane index
  float fracdist; // 0 to 1 (near to far)
  XYZ pip; // point-in-plane
  XYZ pipvd; // point-in-plane, but along viewdir: should be centred!
  PLANE theplane;
  double mu;
  XYZ pt, pt2;

  // and here we place up to 6 vertices for a sliced polygon
  int npolyverts;
  int polyverts[6]; // which edge?
  float polyfracs[6]; // how far along edge?
  // and this is the position angle of the vertex in the viewplane
  float polyangs[6];

  int j,k;
  float ang;

  float xpts[7], ypts[7], zpts[7];
  XYZ iP[6], iTC[6];

  int NPL; // the number of planes we will draw

  // now scale NPL by dot product of (nearvtx - farvtx) . viewdir
  // because this says what is the "depth" of planes...
  pt2 = VectorSub(worldverts[near_vtx], worldverts[far_vtx]);
  ang = DotProduct(pt2, viewdir);
  pt2 = VectorSub(worldverts[0], worldverts[7]); // diagonal
  ang /= Modulus(pt2);
  NPL = (int)(2.0 * ang * sqrt(adim*adim+bdim*bdim+cdim*cdim) );

/*** NOTE: WHAT WAS THIS FOR? */
/*
  int subsample = (*kc % 5) + 1;
  NPL /= (float)subsample;
*/

  // loop in reverse order so farthest planes added to list (and then
  // drawn) first.
  for (plidx = NPL; plidx > 0; plidx--) {

    fracdist = (float)plidx / (float)(NPL+1);
    thisdist = near_dist + fracdist * (far_dist-near_dist);
    
    // point-in-plane along near_vtx to far_vtx line
    pip = VectorSub(worldverts[near_vtx], worldverts[far_vtx]);
    pip = VectorMul(pip, fracdist);
    pip = VectorAdd(worldverts[near_vtx], pip);

    // plane equation: for n={a,b,c}, the plane is n.p=-d, giving
    // ax+by+cz+d = 0
    // So all we do is calculate what d is:
    theplane.a = viewdir.x;
    theplane.b = viewdir.y;
    theplane.c = viewdir.z;
    theplane.d = -1. * DotProduct(viewdir, pip);

    // point-in-plane along viewdir
    p2 = VectorAdd(campos, viewdir);
    if (!LinePlane(campos, p2, theplane, &mu, &pipvd)) {
      fprintf(stderr, "Viewdir doesn't intersect plane: impossible!!!\n");
      exit(EXITERROR);
    }

    npolyverts = 0;

    for (i = 0; i < 12; i++) {
      p1 = worldverts[edges[i][0]];
      p2 = worldverts[edges[i][1]];
      if (LinePlane(p1, p2, theplane, &mu, &pt)) {
	if ((mu >= 0) && (mu <= 1.)) {
	  // get position angle of vertex
	  pt2 = VectorSub(pipvd, pt);
	  Normalise(&pt2);
	  ang = atan2(DotProduct(pt2, upvec), DotProduct(pt2, right));
	  // and insert in list
	  j = 0;
	  while ((j < npolyverts) && (polyangs[j] < ang)) {
	    j++;
	  }
	  k = npolyverts - 1;
	  while (k >= j) {
	    polyverts[k+1] = polyverts[k];
	    polyfracs[k+1] = polyfracs[k];
	    polyangs[k+1] = polyangs[k];
	    k--;
	  }
	  k++;
	  polyverts[k] = i;
	  polyfracs[k] = mu;
	  polyangs[k] = ang;
	  npolyverts++;
	}
      }
    }

    // ok, we have the edges, fraction along those edges, and the 
    // position angle in the view plane of each vertex of this poylgon.
    // Now we need to draw the polygon in eg. clockwise order...
    for (i = 0; i < npolyverts; i++) {
      p1 = worldverts[edges[polyverts[i]][0]];
      p2 = worldverts[edges[polyverts[i]][1]];
      mu = polyfracs[i];
      pt = VectorAdd(p1, VectorMul(VectorSub(p1, p2), mu));
      xpts[i] = pt.x;
      ypts[i] = pt.y;
      zpts[i] = pt.z;

      if (i == 0) {
	xpts[npolyverts] = pt.x;
	ypts[npolyverts] = pt.y;
	zpts[npolyverts] = pt.z;
      }

      // here are the XYZ arrays for 3d texturing via ns2texpoly3d...
      iP[i] = pt;
      p1 = unitverts[edges[polyverts[i]][0]];
      p2 = unitverts[edges[polyverts[i]][1]];
      iTC[i] = VectorAdd(p1, VectorMul(VectorSub(p1, p2), mu));

    }

    if (ss2qrm() == WIREFRAME) {
      s2sci(2 + (plidx % 12));
      s2line(npolyverts+1, xpts, ypts, zpts);
    } else {
      ns2texpoly3d(iP, iTC, npolyverts, texid, itrans, ialpha);
    }

  }
}

void plotHistogram(Histogram h, DataValues dv, VisParameters vp)
{
   ss2tsc("clr");
   float z = 0.001, ybase = 0.02, yheight = 0.06; 
   float x1 = 0.10, x2 = 0.90, dx = x2 - x1;

/* Draw a line at the bottom of the histogram */
   ns2line(x1,ybase,z, x2,ybase,z, 0,0,0);

/* Install the current colour map */
   s2icm(vp.cmap,vp.c1,vp.c1+vp.cr);

   int i;
   float x,y,yy,v;
   long idx;
   float r,g,b;
   XYZ P[4];

   P[0].x = P[1].x = P[2].x = P[3].x = 0.0;
   P[0].y = P[1].y = P[2].y = P[3].y = ybase;
   P[0].z = P[1].z = P[2].z = P[3].z = 0.002;
   COLOUR col;
   float lbmin = log10(h.bmin);
   float lbfac = yheight/(log10(h.bmax)-log10(h.bmin));
   
/* Draw and colour each histogram bin.  If the bin is outside the plotted */
/* data range, draw outline otherwise colour by colour map value for bin. */

   for (i=0;i<h.N;i++) {
      if (h.bin[i] > 0) {
         x = ((float)i/(float)(h.N-1))*dx + x1; 
         if (h.log) {
            y = (log10(h.bin[i])-lbmin)*lbfac + ybase;
         } else {
            y = yheight*h.bin[i]/h.bmax + ybase;
         }
         P[0].x = P[3].x = x;
         P[1].x = P[2].x = x + dx/(float)(h.N-1);
         P[0].y = P[1].y = y;

	 v = dv.min + i/(h.delta);
         if ((v > vp.dmin) && (v < vp.dmax)) {
            idx = ((v-vp.dmin)/(vp.dmax-vp.dmin))*vp.cr + vp.c1;
            s2qcr(idx, &r,&g,&b);
            col.r = r;
            col.g = g;
            col.b = b;
            ns2vf4(P, col);
         } else {
            col.r = 0.0;
            col.g = 0.0;
            col.b = 0.0;
	    ns2vline(P[0],P[1],col);
	    if (i>0) {
               if (h.bin[i-1] > 0) {
                  if (h.log) {
                     yy = (log10(h.bin[i-1])-lbmin)*lbfac + ybase;
                  } else {
                     yy = yheight*h.bin[i-1]/h.bmax + ybase;
                  }
	          ns2line(P[0].x,P[0].y,P[0].z, P[0].x, yy,P[0].z, 0,0,0); 
               } else {
	          ns2line(P[0].x,P[0].y,P[0].z, P[0].x, ybase,P[0].z, 0,0,0); 
               }
            }

            if (i<(h.N-1)) {
               if (h.bin[i+1] > 0) {
                  if (h.log) {
                     yy = (log10(h.bin[i+1])-lbmin)*lbfac + ybase;
                  } else {
                     yy = yheight*h.bin[i+1]/h.bmax + ybase;
                  }
	          ns2line(P[1].x,P[1].y,P[1].z, P[1].x,yy,P[1].z, 0,0,0); 
               } else {
	          ns2line(P[1].x,P[1].y,P[1].z, P[1].x,ybase,P[1].z, 0,0,0); 
               }
            }
         }
      }
   }

/* Put a coloured background behind the histogram */
   P[0].x = P[3].x = 0.0;
   P[1].x = P[2].x = 1.0;
   P[0].y = P[1].y = ybase+1.20*yheight;
   P[2].y = P[3].y = 0.0;
   P[0].z = P[1].z = P[2].z = P[3].z = 0.005;
   col.r = col.g = col.b = 1.0;
   ns2vf4(P, col);

/*
   ns2thline(P[0].x,P[0].y,P[0].z, P[1].x,P[0].y,P[0].z, 1,1,0, 2);
*/

/* Draw a short vertical line at the location of the histogram mean */
   x = (dv.mean-h.dmin)/(h.dmax-h.dmin)*dx + x1;
   ns2thline(x,ybase+1.05*yheight,z, x,ybase+1.15*yheight,z, 0,0,0,1);

/* Draw a horizontal line spanning from mean-std.dev to mean+std.dev */
   x = (dv.mean-dv.sd-h.dmin)/(h.dmax-h.dmin)*dx + x1;
   y = (dv.mean+dv.sd-h.dmin)/(h.dmax-h.dmin)*dx + x1;
   ns2thline(x,ybase+1.15*yheight,z, y,ybase+1.15*yheight,z, 0,0,0,1);

   char string[32];
   float pad = 1.0;
   s2sch(1.0);
   float xx1, xx2, yy1, yy2, dyy, dxx;
   s2sci(S2_PG_BLACK);

/* Write label for overall data minimum */
   sprintf(string,"%.2f",h.dmin*h.scale);	/* CJF: Need scale factor */
   s2qtxtxy(&xx1,&xx2,&yy1,&yy2, 0,0,0, string, pad);
   dyy = yy2-yy1;
   dxx = xx2-xx1;
   s2textxy(x1-1.15*dxx,ybase+0.5*dyy,z,string);

/* Write label for overall data maximum */
   sprintf(string,"%.2f",h.dmax*h.scale);	/* CJF: Need scale factor */
   s2qtxtxy(&xx1,&xx2,&yy1,&yy2, 0,0,0, string, pad);
   dyy = yy2-yy1;
   dxx = xx2-xx1;
   s2textxy(x2+0.15*dxx,ybase+0.5*dyy,z,string);

/* Write label for plotted data minimum */
   sprintf(string,"%.2f",vp.dmin*h.scale);	/* CJF: Need scale factor */
   s2qtxtxy(&xx1,&xx2,&yy1,&yy2, 0,0,0, string, pad);
   dyy = yy2-yy1;
   dxx = xx2-xx1;
   x = (vp.dmin-h.dmin)/(h.dmax-h.dmin)*dx + x1;
   s2textxy(x-1.10*dxx,ybase-0.5*dyy,z,string);
/* Draw vertical line at plotted data min */
   ns2thline(x,ybase,z, x,ybase-0.01,z, 0,0,0, 1);

/* Write label for plotted data maximum */
   sprintf(string,"%.2f",vp.dmax*h.scale);	/* CJF: Need scale factor */
   s2qtxtxy(&xx1,&xx2,&yy1,&yy2, 0,0,0, string, pad);
   dyy = yy2-yy1;
   dxx = xx2-xx1;
   x = (vp.dmax-h.dmin)/(h.dmax-h.dmin)*dx + x1;
   s2textxy(x+0.15*dxx,ybase-0.5*dyy,z,string);
/* Draw vertical line at plotted data max */
   ns2thline(x,ybase,z, x,ybase-0.01,z, 0,0,0, 1);

   ss2tsc("");
}


void printMetaData(FILE *stream, Histogram h, DataValues dv, VisParameters vp)
{
   fprintf(stream,"Scale &  $\\times$%.3f & Data scale factor\\\\\n",h.scale);
   fprintf(stream,"dmin & %.2f & Minimum data value \\\\\n",h.dmin*h.scale);
   fprintf(stream,"dmax & %.2f & Maximum data value \\\\\n",h.dmax*h.scale);
   fprintf(stream,"avg &  %.2f & Sample mean\\\\\n",dv.mean*h.scale);
   fprintf(stream,"$\\sigma$ &  %.2f & Sample standard deviation\\\\\n",dv.sd*h.scale);

   fprintf(stream,"Pdmin & %.2f & Minimum plotted data value \\\\\n",vp.dmin*h.scale);
   fprintf(stream,"Pdmax & %.2f & Maximum plotted data value \\\\\n",vp.dmax*h.scale);
   fprintf(stream,"$\\Delta$ &  %.3f & Histogram bin width\\\\\n",1.0*h.delta*h.scale);

   fprintf(stream,"Cval &  %d & Number of unique colour map values\\\\\n",vp.cr);
   fprintf(stream,"Nbin &  %d & Number of histogram bins\\\\\n",h.N);
   fprintf(stream,"Log & %c & Axis type for histogram \\\\\n",((h.log==1)?'Y':'N'));
/*
   fprintf(stream,"mode &  %.2f & Histogram mode\\\\\n",s.mode);
   fprintf(stream,"Nmode & %.2d & Number of histogram bin containing mode\\\\ \n",s.Nmode);
*/


}


int plotit = 0;
int kcb(unsigned char *key)
{
   if (*key == 'P') {
      plotit = config.Nx*config.Ny;
   }
   return 0;
}


void cb(double *t, int *kc, int *value)
{
/* Reset coordinates of current window and set default drawing colour */
   s2swin(0,1,0,1,0,1);
   s2sci(S2_PG_WHITE);

   pushVRMLname("ANON");
   s2box("BCDE",0,0,"BCDE",0,0,"BCDE",0,0);

/* Row and column for the current panel */
   int row = *value%config.Ny;
   int col = *value/config.Ny;

/* Determine the type of texture object to display */
   if (config.panel[col][row].ctype == TEX3D) { 
      textureVolumeRender(config.panel[col][row].ct.tid[0], config.panel[col][row].ct.N); 
   } else if (config.panel[col][row].ctype == TEX2D) { 
       ds2dvr(config.panel[col][row].ct.tid[0], 1);
   } else if (config.panel[col][row].ctype == TEXPDF) { 
       pushVRMLname((char *)"VRSET1");
       ds2dvrXXX(config.panel[col][row].ct.tid[0],1,1);
       pushVRMLname((char *)"VRSET2");
       ds2dvrXXX(config.panel[col][row].ct.tid[1],1,2);
       pushVRMLname((char *)"VRSET3");
       ds2dvrXXX(config.panel[col][row].ct.tid[2],1,3);
   }
   
/* CJF: Break this from a key-press? */
   config.panel[col][row].h.log = ((*kc+1)%2);


   labelaxes(config.panel[col][row]);

   if (plotit > 0) {
      char pname[128];
      sprintf(pname,"plot.%d.%d",col,row);
      
      FILE *fp = fopen(pname,"w");
      printMetaData(fp,config.panel[col][row].h, config.panel[col][row].dv, 
			config.panel[col][row].ct.vp);
      if (plotit == 1) {
         ss2wtga(pname);
	 FILE *pfp = fopen("packagehead.csh","w");
/* CJF: HACK...*/
         float ybase = 0.02, yheight = ybase+1.2*0.06;
         fprintf(pfp,"%s %f",pname,yheight);
         fclose(pfp);
      }
      plotit--;
   }

   if (config.panel[col][row].h.visible) 
      plotHistogram(config.panel[col][row].h, config.panel[col][row].dv, config.panel[col][row].ct.vp);

   if (plotit < 0) plotit = 0;
}


#ifdef NOTYET
float transfer(float *dval)             /* Alpha transfer function */
{
   if (*dval < gdmin) return 0.0;
   if (*dval > gdmax) return 0.0;
   else return (*dval-gdmin)*gdrange;
}
#endif

void cube2pdf(unsigned int tid[3], float ***array, long N[3], VisParameters p)
{

   int i;
   float tr[12];

/* CJF: Box centred or on edges? */
   for (i=0;i<12;i++) tr[i] = 0.0;
   tr[ 0]  = 0; tr[ 1]  = 1.0/(float)(N[0]-1);
   tr[ 4]  = 0; tr[ 6]  = 1.0/(float)(N[1]-1);
   tr[ 8]  = 0; tr[11]  = 1.0/(float)(N[2]-1);

   char itrans = 's';
   int c1 = p.c1;
   int cr = p.cr;
   s2scir(c1,c1+cr);
   s2icm(p.cmap,c1,c1+cr);

/*CJF: Look at this step */
   ns2sevas(0,0,1);
   for (i=0;i<3;i++) {
/* ns2cvra */
      tid[i] = ns2cvr(array, N[0],N[1],N[2], 0,N[0]-1, 0,N[1]-1, 0,N[2]-1,
                        tr, itrans, p.dmin, p.dmax, p.amin, p.amax);
   }
}

unsigned int cube2texture(float ***array, long N[3], VisParameters p, int tid)
{
   float dmin = p.dmin;
   float dmax = p.dmax;
   float alpha = p.amax;
   int c1 = p.c1;
   int cr = p.cr;
   int w,h,d;
   s2scir(c1,c1+cr);
   s2icm(p.cmap,c1,c1+cr);

   unsigned int id;
   if (tid < 0) {
      id = ss2c3dt(N[0], N[1], N[2]);
   } else {
      id = tid;
   }
   unsigned char *bits = (unsigned char *)ss2g3dt(id, &w, &h, &d);
   memset(bits, (unsigned char)0, w*h*d*4);

   float r,g,b;
   int d_idx, t_idx;
   int i,j,k;
   unsigned char x;
   float cmax   = (float)(cr-1.0);
   float scale  = (cmax)/(dmax-dmin);
   int   c2 = cr-1;

/* Arbitrary scaling factor */
   float factor = (cmax*alpha)/powf(cmax, 2.0);

   for (i=0;i<N[0];i++) {
      for (j=0;j<N[1];j++) {
         for (k=0;k<N[2];k++) {
            d_idx = i * N[1] * N[2] + j * N[2] + k;
            t_idx = k * N[1] * N[0] + j * N[0] + i;
            if ((array[i][j][k] > dmin) && (array[i][j][k] < dmax)) {
               x = (unsigned char)((array[i][j][k] - dmin)*scale);
               s2qcr(c1 + x, &r, &g, &b);
               bits[t_idx*4 + 3] = (int)(powf(x,2.0)*factor);
            } else {
               r = 0;
               g = 0;
               b = 0;
               bits[t_idx*4 + 3] = 0;
            }
            bits[t_idx*4 + 0] = (int)(r * c2);
            bits[t_idx*4 + 1] = (int)(g * c2);
            bits[t_idx*4 + 2] = (int)(b * c2);
         }
      }
   }
   ss2pt(id);
   return id;
}




Header initHeader(void)
{
   Header h;
   h.Naxis = 0;
   int i;
   for (i=0;i<4;i++) { h.naxis[i] = 0; }
#ifdef NOTNOW
   h.bscale = 1.0;
   h.bzero  = 0.0;
   sprintf(h.object,"--");
#endif
   return h;
}

float keyword2float(char *key, char **keywords, int Nkey)
{
   float keyval = -9e99;
   int i;
   int len = strlen(key);
   char dummy[32];
   for (i=0;i<Nkey;i++) {
      if (strncmp(key, keywords[i], len) == 0) {
         sscanf(keywords[i],"%s = %f",dummy, &keyval);
         return keyval;
      } 
   }
   return keyval;
}

long keyword2long(char *key, char **keywords, int Nkey)
{
   long keyval = -1;
   int i;
   int len = strlen(key);
   char dummy[32];
   for (i=0;i<Nkey;i++) {
      if (strncmp(key, keywords[i], len) == 0) {
         sscanf(keywords[i],"%s = %ld",dummy, &keyval);
         return keyval;
      } 
   }
   return keyval;
}


Header parseHeader(char *headerfile)
{
   Header h = initHeader();
   FILE *fp = fopen(headerfile,"r");
   if (fp == NULL) {
      fprintf(stderr,"No header file: %s\n", headerfile);
      exit(EXITERROR);
   }

   char string[256];
   fgets(string,256,fp); 
   int Nkey = 0;
   while (!feof(fp)) {  
      fgets(string,256,fp); 
      Nkey++; 
   }
   rewind(fp);
   char **keywords = (char **)calloc(Nkey,sizeof(char *));
   int i;
   for (i=0;i<Nkey;i++) {
      keywords[i] = (char *)calloc(256,sizeof(char));
      fgets(keywords[i],256,fp); 
   }
   
   h.Naxis = keyword2long("NAXIS",keywords,Nkey);
   if (h.Naxis < 3) {
      fprintf(stderr,"Not a data cube\n");
      exit(EXITERROR);
   }

   char nstring[10];
   for (i=0;i<h.Naxis;i++) {
      sprintf(nstring,"NAXIS%d",i+1); 
      h.naxis[i] = keyword2long(nstring,keywords,Nkey);
      if (h.naxis[i] < 0) {
         fprintf(stderr,"Negative dimension: %ld\n",h.naxis[i]);
         exit(EXITERROR);
      }
   }
   
#ifdef NOTNOW
   h.bscale = keyword2float("BSCALE",keywords,Nkey);
   h.bzero  = keyword2float("BZERO",keywords,Nkey);
#endif
   

   fclose(fp);

   return h;
}


RawCube readRawCube(long N[3], char *fname)
{
   RawCube c;
   FILE *fp = fopen(fname,"rb");
   long i,j,k;
   for (i=0;i<3;i++) c.N[i] = N[i];

   c.vol = initVolume(N[0], N[1], N[2],-1E99); 

   c.dv.min = +1e9; 
   c.dv.max = -1e9;
   c.dv.mean = 0;
   c.dv.sd = 0;
   for (i=0;i<N[0];i++) {
      for (j=0;j<N[1];j++) {
         for (k=0;k<N[2];k++) {
             fread(&c.vol[i][j][k],1,sizeof(float),fp);
         }
      }
   }
   fclose(fp);

   c.dv = characteriseData(c.vol, N);

   return c;
}

void averageData(float ***vr, float ***vol, Axes a, float dmin, float dmax)
{
   int i, j, k;
   int ii=0, jj=0, kk=0;

   int mx = a.input[0]/2;
   int my = a.input[1]/2;
   int mz = a.input[2]/2;

   int ix1 = mx - a.f[0]/2;
   int ix2 = mx + a.f[0]/2;

   int iy1 = my - a.f[1]/2;
   int iy2 = my + a.f[1]/2;

   int iz1 = mz - a.f[2]/2;
   int iz2 = mz + a.f[2]/2;

   int ix, iy, iz;
   double sum = 0;
   float tmp = 0;

   ix = 0;
   for (i=ix1;i<ix2;i+=a.scale[0]) {
      iy = 0;
      for (j=iy1;j<iy2;j+=a.scale[1]) {
          iz = 0;
          for (k=iz1;k<iz2;k+=a.scale[2]) {
	     sum = 0;
	     for (ii=0;ii<a.scale[0];ii++) {
	        for (jj=0;jj<a.scale[1];jj++) {
	            for (kk=0;kk<a.scale[2];kk++) {
	               if (!isinf(tmp = vol[i+ii][j+jj][k+kk]) && (tmp > dmin) && (tmp < dmax)) {
		          sum += tmp;
                       }
                    }
                }
             }
             vr[ix][iy][iz] = sum;
	     iz++;
          }
          iy++;
      }
      ix++;
   }
}

void strideData(float ***vr, float ***vol, Axes a, float dmin, float dmax)
{
   long i, j, k;
   long ii=0, jj=0, kk=0;

   int mx = a.input[0]/2;
   int my = a.input[1]/2;
   int mz = a.input[2]/2;
   int ix1 = mx - a.f[0]/2;
   int ix2 = mx + a.f[0]/2;
   int iy1 = my - a.f[1]/2;
   int iy2 = my + a.f[1]/2;
   int iz1 = mz - a.f[2]/2;
   int iz2 = mz + a.f[2]/2;

/* Needs to be fixed for extraction or averaging */
   for (i=ix1;i<ix2;i+=a.scale[0]) {
      for (j=iy1;j<iy2;j+=a.scale[1]) {
          for (k=iz1;k<iz2;k+=a.scale[2]) {
	     if (vr[ii][jj][kk] > dmin) {
	        if (vr[ii][jj][kk] < dmax) {
                   vr[ii][jj][kk] = vol[i][j][k];
                } else {
                   vr[i][j][k] = (dmin*10);			
                }
             } else {
/* A CJF: bit of trickery - how to deal with out of range values */
                vr[i][j][k] = (dmin*10);			
             } 
	     kk++;
          }
          kk = 0;
          jj++;
      }
      jj = 0;
      ii++;
   }
}

DataValues characteriseData(float ***vol, long N[3])
{
   DataValues dv;
   dv.min   = +1E99;
   dv.max   = -1E99;
   dv.mean  = 0.0;
   dv.sd    = 0.0;
   dv.range = 0.0;
   dv.mid   = 0.0;

   long i, j, k;
   for (i=0;i<N[0];i++) {
      for (j=0;j<N[1];j++) {
         for (k=0;k<N[2];k++) {
	    if (!isinf(vol[i][j][k])) {
	       if (vol[i][j][k] < dv.min) dv.min = vol[i][j][k];
	       if (vol[i][j][k] > dv.max) dv.max = vol[i][j][k];
	       dv.mean += vol[i][j][k];
	       dv.sd   += (vol[i][j][k]*vol[i][j][k]);;
            }
         }
      }
   }

   long Ntot = N[0]*N[1]*N[2];			/* Total number of voxels */
   dv.mean   = dv.mean/(float)Ntot;		/* Sample mean */
   dv.sd     = sqrt(dv.sd/(float)Ntot-(dv.mean*dv.mean))*sqrt((float)(Ntot)/(float)(Ntot-1));
						/* Sample standar deviation */
   dv.range  = dv.max-dv.min;			/* Range */
   dv.mid    = 0.5*(dv.max+dv.min);		/* Mid-point */

   return dv;
}

void generateHistogram(Histogram *h, DataValues dv, float ***vr, long *N)
{
   int i, j, k;					/* Loop variables */

/* Create N bins and set to zero */
   for (i=0;i<h->N;i++) { 
      h->bin[i] = 0;
   }
   
   h->dmin = dv.min;				/* Raw cube data minimum and data maximum */
   h->dmax = dv.max;
   h->delta = (float)(h->N-1)/(h->dmax-h->dmin);/* Step size */

   long idx;					/* Index of histogram bin */
   for (i=0;i<N[0];i++) {			/* Loop over axes */
      for (j=0;j<N[1];j++) {
         for (k=0;k<N[2];k++) {
	    if (!isinf(vr[i][j][k])) {		/* Is this a usable data values? */
	       idx = (floor)((vr[i][j][k] - h->dmin)*h->delta);
						/* Calculate bin for current data value */
	       h->bin[idx] += 1.0;		/* Incremement count in this bin */ 
            }
         }
      }
   }

   h->sum = 0;					/* Total number of items in bin */
   for (i=0;i<h->N;i++) { 
      h->sum = h->sum + h->bin[i];
      if ((h->bin[i] > 0) && (h->bin[i] < h->bmin)) h->bmin = h->bin[i];
      if (h->bin[i] > h->bmax)  h->bmax = h->bin[i]; 
   }
   h->visible = 1;
   h->log     = 0;

   while ((h->scale) * fabs(h->dmin) < 1) {
      h->scale *= 10.0; 
   }
}


Axes scaleAxes(long target[3], long input[3])
{
   Axes a;
   int i;
   for (i=0;i<3;i++) {
      a.input[i]  = input[i];
      a.target[i] = target[i];
      if (input[i] < a.target[i]) a.target[i] = input[i];
      a.scale[i]  = input[i]/a.target[i];
      a.f[i]      = a.scale[i]*a.target[i];

   }
   return a;
}


void generateTextures(Configuration c)
{
/* CJF: Improvement = unique colour map per cube */

   int c1 = _CFIRST, cr = _CMAX, c2 = cr-1;	/* Indices for colour map */
   RawCube *cube = c.cube;			/* Convenience pointer to RawCube */
   int i,j;					/* Loop variables */

   for (i=0;i<config.Nx;i++) {
      for (j=0;j<config.Ny;j++) {
         Panel *cp   = &config.panel[i][j];	/* Convenience pointer to current Panel */
	 RawCube *cb = &cube[cp->cubeID];	/* Convenience pointer to current Cube */

/* Set the viewport for the curent panel */
         xs2cp(cp->pid);
         s2svp(-cb->N[0],cb->N[0],-cb->N[1],cb->N[1],-cb->N[2]*cp->stretch,cb->N[2]*cp->stretch);

/* Modify the data to make it suitable for visualising */
         cp->axes = scaleAxes(cp->target, cb->N);
         cp->ct   = initCubeTexture(cp->axes, _CMAX);
         averageData(cp->ct.vr, cb->vol, cp->axes, _DMIN, _DMAX); 

/* CJF: At this point, the original data cube can be freed */
         cp->dv = characteriseData(cp->ct.vr, cp->ct.N);
         generateHistogram(&cp->h, cp->dv, cp->ct.vr, cp->ct.N); 

/* Make a reasonable guess at the data range to display */
         cp->ct.vp = initVisParameters(cp->dv.mean, cp->dv.mean+8.0*cp->dv.sd, c1, cr, c.cmap); 

/* Create the textures depending on the input type */
         switch (cp->ctype) {
            case TEX3D  : cp->ct.tid[0] = cube2texture(cp->ct.vr, cp->ct.N, cp->ct.vp, -1); 
		          break;
            case TEX2D  : cube2pdf(cp->ct.tid, cp->ct.vr, cp->ct.N, cp->ct.vp);
		          break;
            case TEXPDF : cube2pdf(cp->ct.tid, cp->ct.vr, cp->ct.N, cp->ct.vp);
		          break;
         }
         c1 += c2;				/* Next segement of colour space */
      }
   }

   xs2cp(0);
}

void vrlights(void)
{
   ss2spt(PERSPECTIVE);
   ss2srm(SHADE_FLAT);
   COLOUR amb = {0.97, 0.97, 0.97};
   ss2sl(amb, 0, NULL, NULL, 0);
}

Configuration configure(int argc, char *argv[])
{
   Configuration c;
   int i, off = 1, Ncube = 1;
   char cfile[256];

/* Look for a configuration file */
   if (argv[1][0] == '-') {
      int len = strlen(argv[1]);
      for (i=1;i<len;i++) {
         cfile[i-1] = argv[1][i];
      }
      cfile[i-1] = '\0';
      off = 2;
      FILE *cfp = fopen(cfile,"r");
      if (cfp != NULL) {
         c = readConfiguration(cfp, _CMAX);
         fclose(cfp);
      } else {
/* CJF: Needs some testing - what happens if no configuration file? */
         Ncube  = (argc-off)/2;
         c = initConfiguration(Ncube, Ncube, _CMAX, _CMAP);
      }
   } else {
      off = 1;
      Ncube  = (argc-off)/2;
      c = initConfiguration(Ncube, Ncube, _CMAX, _CMAP);
   }

   c.Ncube  = (argc-off)/2;
   c.header = (Header  *)calloc(c.Ncube, sizeof(Header));
   c.cube   = (RawCube *)calloc(c.Ncube, sizeof(RawCube));

   for (i=0;i<c.Ncube;i++) {
      c.header[i]    = parseHeader(argv[off]);
      c.cube[i].N[0] = c.header[i].naxis[0];
      c.cube[i].N[1] = c.header[i].naxis[1];
      c.cube[i].N[2] = c.header[i].naxis[2];
      c.cube[i]      = readRawCube(c.cube[i].N, argv[off+1]);
      off+=2;
   }

   return c;
}

int main(int argc, char *argv[])
/* Usage: s2raw -config.txt <header1.fits> <raw1.bin> [headier2.fits> <raw2.bin> ...] */
/* Assumes: Pairs of data files and header files. Named configuration file */
{

/* Start the S2PLOT display environment */
   s2opend("/?",argc,argv);

/* Read the configuration files: what will be displayed? */
   config = configure(argc, argv);

/* Generate the S2PLOT texture items depending on display mode */
   generateTextures(config);   

/* Draw and label a box: required to get camera locations right */
   s2swin(-1,1,-1,1,-1,1);
   s2box("BCDET",0,0,"BCDET",0,0,"BCDET",0,0);

/* Set up the S2PLOT lights for volume rendering */
   vrlights();
   
/* Hand control to S2PLOT display loop */
   s2show(1);

   return 0;
}
