/*
 Copyright (C)	2014-2016 Vigeesh Gangadharan <vigeesh@kis.uni-freiburg.de>,
	Kiepenheuer-Institut fuer Sonnenphysik, Freiburg, Germany.

  Permission is hereby granted, free of charge, to any person obtaining a
  copy of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.
*/

/*
 *  postel_remap.c
 *
 *  Rempas a Postel projected data to Plate-Carree for the FLD pipeline
 *  The module takes in a 3-dimensional real dataset from the mtrack
 *
 *  Responsible:
 *      Vigeesh Gangadharan                             vigeesh@kis.uni-freiburg.de
 *
 *  Usage:
 *    postel_remap [-ts] in= out= ...
 *
 *  Arguments: (type		default         description)
 *	in	DataSet none            Input dataset
 *			A series of records containing 3-d data cubes as one
 *			of their segments. It is assumed that all records in
 *			the dataset are from the same dataseries, or at least
 *			share a common segment structure
 *          segment string	-		Name of the input data segment (only
 *			required if the input series has multiple 3-dimensional
 *			segments)
 *	out	DataSeries	-		Output data series;  must
 *
 *  Flags
 *
 *  Bugs:
 *
 *
 *  Future Updates
 *
 *
 *  Revision history
 *      
 */

/**
   @defgroup
   @ingroup

   @brief A brief module description goes here

   @par Synopsis:
   @code
   ModuleName [-bdfxGEN_FLAGS] in=<record-set query> out=<out series>
   @endcode

   The first paragraph of an in-depth description goes here.

   This is the second paragraph of the description

   @par Flags:
   @c -b: Some flag <br>
   @c -d: Another flag <br>
   @c -f: Yet another flag <br>
   @c -x: A fourth flag <br>

   @par GEN_FLAGS:
   Ubiquitous flags present in every module.
   @ref jsoc_main

   @param in A record-set query that identifies input records.
   @param out The output series.

   @par Exit_Status:
   Brief description of abnormal, non-zero, exit values.

   @par Example:
   Brief description of the first example goes here
   @code
   ModuleName -bx in=<hmi.lev0> out=<hmi.lev03>
   @endcode

   @par Example:
   Brief description of the second example goes here
   @code
   ModuleName -f in=<hmi.lev0> out=<hmi.lev05>ModuleName -f in=<hmi.lev0> out=<hmi.lev05>
   @endcode

   @bug
   A description of any bugs goes here.

   @par Code:
   The doxygen code that makes this page is here:
   @verbinclude  ./doxygen_moduletemplate.txt
*/


// All the nessasary include files
#include "drms.h"
#include "jsoc.h"
#include "jsoc_main.h"
#include "fitsio.h"

#include "keystuff.c"	/* Required to propagate keys*/


//#include "fitsio.h"
//#include "fstats.h"

// The name of the module
char *module_name = "jsoc_postel_stretch";
char *module_desc = "Stretch Postel maps for FLD";
char *version_id = "0.1";

/*  global declaration of missing to be initialized as NaN  */
float missing_val=0.0 / 0.0;;
#define INTERP_NEAREST_NEIGHBOR	(1)
#define INTERP_BILINEAR	(2)


#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

/* All the input arguments
 * input series as 'in'
 * output series as 'out'
 * latitude of the center of the map as 'center_lat'
 * longitude of the center of the map as 'center_lon'
 * width in longitudes "width_lon'
 * height in latitudes "height_lat'
 * flag -t for tests, no writing to series done
 */

ModuleArgs_t module_args[] =
{
		{ARG_STRING, "in", "NOT SPECIFIED",  "Input data series."},
		{ARG_STRING, "out", "NOT SPECIFIED",  "Output data series."},
		{ARG_FLAG, "t", "0", "testing only, no files written"},
		{ARG_INT, "inseg", "0", "Input segment number"},
		{ARG_INT, "outseg", "0", "Output segment number"},
		{ARG_FLOAT, "rhow", "0", "Height/Radial of the map section in degrees"},
		{ARG_FLOAT, "phiw", "360", "Width/Azimuthal of the map section in degrees"},
		{ARG_INT, "rows", "1024", "Number of pixels in radial"},
		{ARG_INT, "cols", "1024", "Number of pixels in azhimuthal"},
		{ARG_STRING, "requestid", "NA", "RequestID if called as an export processing step."},
		{ARG_END}
};

#define     Deg2Rad    (M_PI/180.0)
#define     Rad2arcsec (3600.0/Deg2Rad)
#define     arcsec2Rad (Deg2Rad/3600.0)
#define     Rad2Deg    (180.0/M_PI)


char *propagate[] = {"CarrRot", "CMLon", "CDELT1","CDELT2","CDELT3","CRPIX1","CRPIX2",
		"MidTime", "Duration", "LonSpan", "T_START", "T_STOP", "Coverage", "Quality",
		"MapScale", "Width", "Height"};



struct ObsInfo_struct
{
	// from observation info
	TIME  t_obs;
	double rsun_obs, obs_vr, obs_vw, obs_vn;
	double crpix1, crpix2, cdelt1, cdelt2, crota2;
	double crval1, crval2;
	double cosa, sina;
	double obs_b0;
	// observed point
	// parameters for observed point
	double x,y,r;
	double rho;
	double lon;
	double lat;
	double sinlat, coslat;
	double sig;
	double mu;
	double chi;
	double obs_v;
};

typedef struct ObsInfo_struct ObsInfo_t;


int DoIt(void)
{

	/*
	 * DRMS definitions.
	 * Setting up the input/output records, segments and arrays
	 */
	int status = DRMS_SUCCESS;
	DRMS_RecordSet_t *inRS, *outRS;
	DRMS_Record_t *inRec,*outRec;
	DRMS_Segment_t *inSeg, *outSeg, *logSeg;
	DRMS_Array_t *inArray, *outArray;
	CmdParams_t *params = &cmdparams;

	/* Variable definitions */
	/* */

	fitsfile *fptr;			/* for the test fits file	*/

	const char *filename;
	char *buf;
	const char *series = NULL;
	char history[4096];
	char *mapproj;
	char yesno[10];
	char *inStr;
	char *dir;				/* path of the final map 	*/
	char module_ident[64];	/* for drms initialization 	*/

	char *inSer;			// DRMS series names
	char logfilename[DRMS_MAXPATHLEN];	//drms logs
	char *outsegname;					//drms name of the output segment
	char source[DRMS_MAXQUERYLEN];

	int verbose_logs;				//drms logs
	int nframes;
	int usestdin = 0;
	int irec, nrecs;
	int img_cnt;
	int quality=0;
	int azhpixsiz=200;//roi_wpix; //roi_width/maplonstp;
	int rhopixsiz=200;//roi_hpix; //roi_height/maplatstp;
	int map_cnt;
	int ijk;		// index
	int ni,iframe;	// indices
	int inx, iny;	// indicies
	int in_nx, in_ny;	// indicies
	int segs,isegnum,osegnum,found	;		//for drms initialization

	int keyct = sizeof (propagate) / sizeof (char *); 		/* Number of keywrods to propagate */

	long im_xaxis, im_yaxis, im_taxis, nelements;			// for the test fits file

	float *inData,*outData; 	// to store the input and output data array
	float *img_xhp, *img_yhp;
	float *map_xhp, *map_yhp, *map_dhp;
	float mapazhstp= 0.125; //roi_width/roi_wpix ; //0.175955;
	float maprhostp = 0.125; //roi_height/roi_hpix ; //0.175955;
	float phistep;
	float rhostep;

	float y00,y01,y02,y03,t,u;		/* Variables for map interpolation	*/
	int i,j,ij,nii1,nij1,niij1;
	float * mapData;
	int il,iu,jl,ju,im,jm;

	float width,height,extent;

	double dsun,rsun;
	double b_0, l_0;
	double map_azi,map_rho;
	double xdegstp,ydegstp;
	double x0,y0;

	struct stat file_stat;

    float *readData;
	int fstatus = 0;
	float nullval=0.0;
	int anynul;

	const char *inQuery = params_get_str(&cmdparams, "in");
	const char *outSer = params_get_str(&cmdparams, "out");
	const char *requestid = params_get_str(&cmdparams, "requestid");
	int crop = params_get_int(&cmdparams, "c");
	int as_is = params_get_int(&cmdparams, "u");
	float roi_rhow = params_get_float(&cmdparams,"rhow");		/*radial from 0 to rho_max degrees*/
	float roi_phiw = params_get_float(&cmdparams,"phiw");		/*phi=0,360 degrees*/
	int roi_rows = params_get_int(&cmdparams,"rows");
	int roi_cols = params_get_int(&cmdparams,"cols");
	int test = params_get_int(&cmdparams, "t");
	int inseg = params_get_int(&cmdparams, "inseg");
	int outseg = params_get_int(&cmdparams, "outseg");
	int full_header = params_get_int(&cmdparams, "h") || strcmp(requestid, "NA");
	int verbose = params_isflagset (params, "v");
	char *insegname; //= strdup (params_get_str (params, "segment"));

	/**
	 * Initializing the Input record
	 *
	 */

	/* Get the record from the input query string */
	inRS = drms_open_records (drms_env, inQuery, &status);
	if (status)	{
		fprintf (stderr, "Error (%s): drms_open_records() returned %d for "
				"dataset:\n  %s\n", module_ident, status, inQuery);
		DIE("Check if the series exists\nQuitting");
	}

	/* Check how many records are there in the input series */
	nrecs = inRS->n;
	if (nrecs < 1)	fprintf (stderr, "No records found in input data set %s\n", inQuery);

	/* Select the record */
	/* and read the number of data segments in the selected record */
	inRec = inRS->records[0];
	inSer = strdup (inQuery);
	if ((segs = drms_record_numsegments (inRec)) < 1) {
		fprintf (stderr, "Error: no data segments in input data series:\n");
		fprintf (stderr, "  %s\n", inSer);
	}

	found = 0;
	/* Select the segment with the data, i.e that contains 3D array */
	/* the series has two data segments 'Log' and 'V' */
	for (int n = 0; n < segs; n++) {
		inSeg = drms_segment_lookupnum (inRec, n);
		if (inSeg->info->naxis != 3) continue;
		if (!found) isegnum = n;
		found++;
	}
	/* Error if there were no 3d array (rank=3 data segment) */
	if (!found) {
		fprintf (stderr, "Error: no segment of rank 3 in input data series:\n");
		fprintf (stderr, "  %s\n", inSer);
	}

	/* If an array of rank=3 is found, and there are more than one rank 3 array,
	 * please specify the name of the segment to be used */
	if (found > 1) {
		if (strcmp (insegname, "Not Specified")) {
			inSeg = drms_segment_lookup (inRec, insegname);
			if (!inSeg) {
				fprintf (stderr,
						"Warning: requested segment %s not found in input data series:\n",
						insegname);
				fprintf (stderr, "  %s\n", inSer);
				inSeg = drms_segment_lookupnum (inRec, isegnum);
				fprintf (stderr, "  using segement %s\n", inSeg->info->name);
			} else if (inSeg->info->naxis != 3) {
				fprintf (stderr,
						"Warning: requested segment %s in input data series:\n", insegname);
				fprintf (stderr, "  %s is not 3-dimensional", inSer);
				inSeg = drms_segment_lookupnum (inRec, isegnum);
				fprintf (stderr, " using segment %s\n", inSeg->info->name);
			} else isegnum = inSeg->info->segnum;
		} else {
			fprintf (stderr,
					"Warning: multiple segments of rank 3 in input data series:\n");
			fprintf (stderr, "  %s\n", inSer);
			fprintf (stderr, "  using %s\n", inSeg->info->name);
		}
	}


	/**
	 * Initializing the Output record
	 *
	 */
	/* Create the record in the DRMS*/
	if (!(outRS = drms_create_records (drms_env, nrecs, (char *) outSer, DRMS_PERMANENT,
			&status))) {
		fprintf (stderr, "Error: unable to create %d records in series %s\n",
				nrecs, outSer);
		fprintf (stderr, "       drms_create_records() returned status %d\n", status);
		return 1;
	}
	if (verbose) printf ("creating %d record(s) in series %s\n", nrecs, outSer);

	/*  Check output data series structure  */
	outRec = drms_recordset_getrec (outRS, 0);
	if (!outRec) {
		fprintf (stderr, "Error accessing record %d in series %s\n", 0, outSer);
		drms_close_records (outRS, DRMS_FREE_RECORD);
		return 1;
	}

	/*Check how many data segments are there with rank=3 i.e 3d array*/
	segs = drms_record_numsegments (outRec);
	found = 0;
	for (int n = 0; n < segs; n++) {
		outSeg = drms_segment_lookupnum (outRec, n);
		if (outSeg->info->naxis != 3) continue;
		if (!found) outsegname = strdup (outSeg->info->name);
		found++;
	}
	if (found < 1) {
		fprintf (stderr,
				"Error: no data segment of dimension 3 and appropriate size in output series %s\n", outSer);
		drms_close_records (outRS, DRMS_FREE_RECORD);
		return 1;
	}
	/* If an array of rank=3 is found, and there are more than one rank 3 array, please specify the name of the segment to be used */
	if (found > 1) {
		if (strcmp (outsegname, "Not Specified")) {
			outSeg = drms_segment_lookup (outRec, outsegname);
			if (!inSeg) {
				fprintf (stderr,
						"Warning: requested segment %s not found in input data series:\n",
						outsegname);
				fprintf (stderr, "  %s\n", outSer);
				inSeg = drms_segment_lookupnum (outRec, osegnum);
				fprintf (stderr, "  using segement %s\n", outSeg->info->name);
			} else if (outSeg->info->naxis != 3) {
				fprintf (stderr,
						"Warning: requested segment %s in input data series:\n", outsegname);
				fprintf (stderr, "  %s is not 3-dimensional", inSer);
				outSeg = drms_segment_lookupnum (outRec, isegnum);
				fprintf (stderr, " using segment %s\n", outSeg->info->name);
			} else isegnum = outSeg->info->segnum;
		} else {

			fprintf (stderr,
					"Warning: multiple data segments of dimension 3 and appropriate size in output series %s\n", outSer);
			fprintf (stderr, "       using \"%s\"\n", outsegname);
		}
	}

	logSeg = drms_segment_lookup (outRec, "Log");
	if (logSeg) drms_segment_filename (logSeg, logfilename);
	else if (verbose) {
		fprintf (stderr,
				"Warning: segment \"Log\" not present in output series %s\n", outSer);
		fprintf (stderr, "         verbose logging turned off\n");
		verbose_logs = 0;
	}

	if (verbose_logs) {
		logSeg = drms_segment_lookup (outRec, "Log");
		if (logSeg) drms_segment_filename (logSeg, logfilename);
		//fopen (logfilename, "w");
	}



	/* Get the relevant keyword values from the DRMS*/
	quality = drms_getkey_int(inRec, "QUALITY", &status);
	mapproj=drms_getkey_string(inRec,"MapProj",&status);
	width = drms_getkey_float(inRec, "Width", &status);
	height = drms_getkey_float(inRec, "Height", &status);

	/* Assuming that the crpix1/2 arein the center of the image*/
	extent = (((width) < (height)) ? (width/2.) : (height/2.));


	if (strcmp(mapproj,"Postel"))
				DIE("-----------\n\n\n\nUsed only with \"Postel\" projected data, read the manual\n\n\n\n\n-----------\n");


	/*
	 * Setting uo the arrays
	 */

	azhpixsiz=roi_rows;//roi_wpix; //roi_width/maplonstp;
	rhopixsiz=roi_cols;//roi_hpix; //roi_height/maplatstp;
	phistep=roi_phiw/azhpixsiz;

	if (roi_rhow == 0) roi_rhow=extent;
	rhostep=roi_rhow/rhopixsiz;

	irec=0;	//assuming that there is only one record
	inRec = inRS->records[irec];


	/*
	 * This part of the code reds in test data to check if the projections are OKay
	 */

	im_xaxis=160;
	im_yaxis=160;
	im_taxis=2;

	nelements=im_xaxis*im_yaxis*im_taxis;

	readData = (float *) calloc(nelements,sizeof(float));			//for read individual slice
	char fpname[DRMS_MAXPATHLEN];
	//strcpy(fpname,"/dat/seismo/vigeesh/mapping/Circle.fits");
	strcpy(fpname,"/dat/seismo/vigeesh/mapping/spot_grid.fits");

	fits_open_file(&fptr,fpname, READONLY, &fstatus);
	/* Write the array of integers to the image */
	fits_read_img(fptr, TFLOAT, 1, nelements, &nullval,readData, &anynul , &fstatus);
	fits_close_file(fptr, &fstatus);            /* close the file */

	/*
	 * End of test data read
	 */


	if (status || (!status && quality >= 0))
	{
		inSeg = drms_segment_lookupnum(inRec, inseg);
		inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);

		in_nx 	= inArray->axis[0];
		in_ny 	= inArray->axis[1];
		nframes = inArray->axis[2];

		int outDims[3] = {azhpixsiz,rhopixsiz,nframes};
		outArray = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
		outData = (float *) calloc(azhpixsiz*rhopixsiz*nframes,sizeof(float));

		inData = (float *)inArray->data;

		xdegstp=drms_getkey_double(inRec, "CDELT1",&status);
		ydegstp=drms_getkey_double(inRec, "CDELT2",&status);
		x0=drms_getkey_double(inRec, "CRPIX1",&status);
		y0=drms_getkey_double(inRec, "CRPIX2",&status);
		rsun=drms_getkey_double(inRec, "RSUN_REF",&status);
		dsun=drms_getkey_double(inRec, "DSUN_REF",&status);
		b_0 =drms_getkey_double(inRec, "CRLT_OBS",&status)*Deg2Rad;


		/*
		 * The img_?hp contains the Postel images that is grabbed from the mtrack
		 */

		img_cnt	= in_nx*in_ny*nframes;
		img_xhp = (float *)malloc (img_cnt * sizeof (float));
		img_yhp = (float *)malloc (img_cnt * sizeof (float));

		map_cnt=azhpixsiz*rhopixsiz*nframes;

		//long naxis2=3;
		long naxes2[3] = {azhpixsiz, rhopixsiz,nframes};
		long naxes[3] = {in_nx , in_ny, nframes};

		mapData=(float *)malloc (map_cnt * sizeof (float));
		map_xhp=(float *)malloc (map_cnt * sizeof (float)); 	// these arrays contain the x,y values in terms of Postel coordinates
		map_yhp=(float *)malloc (map_cnt * sizeof (float));		//

		/*
		 * This here calculates the axis values
		 * The step size of the grid is taken from the crdelt keywords
		 * the axis value of the image is then calculated from this values
		 * img_  - 10 rad to +10 rad
		 *
		 *
		 */
		for (iframe=0; iframe<nframes; iframe++)
		{
			ni=0;
			for (iny = 0; iny < in_ny; iny += 1)
				for (inx = 0; inx < in_nx;  inx += 1) {
					ni = iframe*in_nx*in_ny+in_nx*iny + inx;
					//Compute the HeliProj axis value at each pixel loc
					img_xhp[ni]=xdegstp*(inx-x0)*Deg2Rad;
					img_yhp[ni]=ydegstp*(iny-y0)*Deg2Rad;
				}

			/*
			 * The map_?hp contains the co-ordinates of the stretched images
			 * This is grid on which the final map will be interpolated
			 * the map goes from map_rho = 0-rho_max and map_azi = 0 to 2pi
			 */


			ni=0;
			// This here calculates the axis values
			for (iny = 0; iny < rhopixsiz; iny += 1)
				for (inx = 0; inx < azhpixsiz;  inx += 1) {
					ni = iframe*azhpixsiz*rhopixsiz+azhpixsiz*iny + inx;
					//the axis value in proper units
					//map_rho=0.125*iny*Deg2Rad;				/* Compute the axis for rho */
					//map_azi=2.0*M_PI/azhpixsiz*inx;			/* Compute the axis for phi */
					map_rho=rhostep*iny*Deg2Rad;				/* Compute the axis for rho */
					map_azi=phistep*inx*Deg2Rad;			/* Compute the axis for phi */

					//this contains the computed value of the img pixels
					//theta = atan2(map_hglt,map_hgln);
					map_xhp[ni]= map_rho * cos(map_azi) ;
					map_yhp[ni]= map_rho * sin(map_azi) ;
				}


			/*
			 * The interpolation begins here
			 */

			ni=0;
			for (j=0; j < naxes2[1]; j++)
			{
				for (i=0; i < naxes2[0]; i++)
				{

					ijk=iframe*naxes2[0]*naxes2[1]+j*naxes2[0]+i;	// ijk coz, if there are more than one final image (map_
					//printf("%d,%d,%d,%d\n",i,j,irec,ijk);
					ij=ijk;							// ij is the index for thesingle image (map_
					il=0;										//
					iu=naxes[0]-1;								//	These are the extremas of the index in the original (img_
					jl=0;										//
					ju=naxes[1]-1;								//
					while (iu-il > 1 && ju-jl > 1)
					{
						im = (iu+il) >> 1; //locate the mid point
						jm = (ju+jl) >> 1;

						if (map_xhp[ij] >= img_xhp[iframe*naxes[0]*naxes[1]+jm*naxes[0]+im])
						{
							if (map_yhp[ij] >= img_yhp[iframe*naxes[0]*naxes[1]+jm*naxes[1]+im])
							{
								il=im;
								jl=jm;
							}
							else
							{
								il=im;
								ju=jm;
							}
						}
						else
						{
							if (map_yhp[ij] >= img_yhp[iframe*naxes[0]*naxes[1]+jm*naxes[1]+im])
							{
								iu=im;
								jl=jm;
							}
							else
							{
								iu=im;
								ju=jm;
							}
						}

					}

					ni = iframe*naxes[0]*naxes[1]+jl*naxes[1] + il;
					nii1 = iframe*naxes[0]*naxes[1]+jl*naxes[1] + iu;
					nij1 = iframe*naxes[0]*naxes[1]+ju*naxes[1] +il;
					niij1 = iframe*naxes[0]*naxes[1]+ju*naxes[1] + iu;

					y00 = inData[ni];
					y01 = inData[nii1];
					y02 = inData[niij1];
					y03 = inData[nij1];

					t = (map_xhp[ij] - img_xhp[ni])/(img_xhp[nii1] - img_xhp[ni]);
					u = (map_yhp[ij] - img_yhp[ni])/(img_yhp[nij1] - img_yhp[ni]);

					if ((y00 == missing_val) | (y01 == missing_val) | (y02 ==missing_val) | (y03 ==missing_val))
					{
						outData[ijk] = missing_val;
					}
					else
					{
						if (test) {
							outData[ijk] = readData[ni];
						}
						else {
							outData[ijk] = inData[ni]; //(1-t)*(1-u)*y00 + t*(1-u)*y01 + t*u*y02 + (1-t)*u*y03;
						}
					}
				} // i-loop
			} // j-loop
            
            //progress-bar
            int barwidth=70;
            //printf("%d, %d, %f",iframe,nframes, (float) iframe/(float)nframes);
            printf("%3d%% [", (int)((float)iframe/(float)nframes * 100) );
            for (int ip=0; ip < (int)(float)iframe/(float)barwidth; ip++)
                printf("=");
            for (int ip=(int)(float)iframe/(float)barwidth; ip<barwidth; ip++)
                printf(" ");
            //printf("]\n\033[F\033[J");
            printf("]\r");
            fflush(stdout);
            //DIE("GG");
		}	// iframe-loop
	}

	/* Copy the data into the output array structure */
	memcpy (outArray->data, outData, map_cnt * sizeof (float));

	/* Some hard coded stuff */
	float bzero=0.0;
	float bscale=0.5;
	outArray->bzero = bzero;
	outArray->bscale = bscale;

	outRec = outRS->records[0];
	outSeg = drms_segment_lookupnum(outRec, outseg);

    //the directory path in which the file is availabe is now in variable dir
    drms_record_directory(outRec, &dir, 1);
    printf("The final map is available at:\n %s/%s.fits\n",&dir,outSeg->info->name);

    /* Propagate the keys */
    propagate_keys (outRec, inRec, propagate, keyct);

    drms_sprint_rec_query (source, inRec);
    check_and_set_key_str   (outRec, "MapProj", "PlateCarre");
    check_and_set_key_str   (outRec, "Source", source);
    check_and_set_key_str   (outRec, "Module", module_ident);
    check_and_set_key_str   (outRec, "BLD_VERS", jsoc_version);
    check_and_set_key_str   (outRec, "Input", inSer);
    check_and_set_key_time  (outRec, "Created", CURRENT_SYSTEM_TIME);
    //for sunspots
    check_and_set_key_double(outRec, "LatHG", 0.0);
    check_and_set_key_double(outRec, "LonHG", 0.0);

	/* Writing the record to the drms */
	if(drms_segment_write(outSeg,outArray,1))
		fprintf(stderr, "Can not write file for the segment.\n");

	/* Close the input record */
	if(drms_close_records(inRS, DRMS_FREE_RECORD))
		fprintf(stderr, "Error closing Record for input series\n");

	/* Close the output record */
	if(drms_close_records(outRS, DRMS_INSERT_RECORD))
		fprintf(stderr, "Error closing Record for output series\n");


	drms_free_array (inArray);
	drms_free_array (outArray);
	//free_all (clat, clon, mai, delta_rot, maps, last, maplat, maplon, map_sinlat,
	//		offsun, orecs, log, reject_list);
	printf("\n---END---\n");
	return (DRMS_SUCCESS);

} // end of DoIt
