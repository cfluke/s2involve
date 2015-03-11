#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "routines.h"

#define RAWFILE    "raw.bin"
#define HEADERFILE "header.txt"


int listhead(char *fname, char *oname)
/* Based on the listhead example: */
/* https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html#listhead */

{
    FILE *ofp = fopen(oname,"w");	
    if (ofp == NULL) {
       fprintf(stderr,"EXIT: Could not open output file: %s\n",oname);
       exit(1);
    }

    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int single = 0, hdupos, nkeys, ii;

    if (!fits_open_file(&fptr, fname, READONLY, &status))
    {
      fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

      /* List only a single header if a specific extension was given */ 
      if (hdupos != 1 || strchr(fname, '[')) single = 1;

      for (; !status; hdupos++)  /* Main loop through each extension */
      {
        fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */

        fprintf(ofp, "# Header listing for HDU #%d:\n", hdupos);

        for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */

           if (fits_read_record(fptr, ii, card, &status))break;
           fprintf(ofp,"%s\n", card);
        }
        fprintf(ofp,"END\n");  /* terminate listing with END */

        if (single) break;  /* quit if only listing a single header */

        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
      }

      if (status == END_OF_FILE)  status = 0; /* Reset after normal error */

      fits_close_file(fptr, &status);
    }

    fclose(ofp);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}


FITSCube readFITS(char *fname)
{
   FILE *fp = fopen(fname, "r");
   if (fp == NULL) {            /* Could not open the specified file */
      fprintf(stderr,"EXIT: Could not open input FITS file: %s\n",fname);
      exit(1);
   }
   fclose(fp);

   FITSCube c = readFITScubeHeader(fname, 1);

   readFITScube(&c, 1);

   return c;
}

void fits2raw(char *fname, long naxes[3], float ***array)
{
   long i, j, k;
   FILE *ofp = fopen(fname,"wb");
   if (ofp == NULL) {
      fprintf(stderr,"EXIT: Could not open output RAW file: %s\n",fname);
      exit(1);
   }
   for (i=0;i<naxes[0];i++) {
      for (j=0;j<naxes[1];j++) {
          for (k=0;k<naxes[2];k++) {
             fwrite(&array[i][j][k],1,sizeof(float),ofp);
          }
      }
   }
   fclose(ofp);
}

int main(int argc, char *argv[])
{
/* Extract the header from the FITS file */
   listhead(argv[1], HEADERFILE);		

/* Read the FITS file primary data cube */
   FITSCube cube = readFITS(argv[1]);  		

/* Convert the primary data cube to raw format */
   fits2raw(RAWFILE, cube.naxes, cube.array);	
  
   return 0;
}
