/*
  Copyright (C) 2014-2016 Vigeesh Gangadharan <vigeesh@kis.uni-freiburg.de>,
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
 *  fld.c
 *
 *  Calculate the FLD Coefficients of an input time series Legendre coefficients
 *  The module takes in a 3-dimensional real dataset from the mtrack
 *
 *  Responsible:
 *      Vigeesh Gangadharan                             vigeesh@kis.uni-freiburg.de
 *
 *  Usage:
 *    fld [-lvx] in= out= ...
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
 *			contain at least one 3-dimensional segment
 *    mask_in	Float		0.9375		Inner apodization radius
 *    mask_ex	Float		1.0		Outer apodization radius
 *    apodize	Float		0.96875		Temporal apodization edge
 *    fbin	int		0		output frequency binning
 *
 *  Flags
 *	-l	output direct tpower spectrum rather than scaled log
 *	-n	do not save results (diagnostic only)
 *	-v	run verbose
 *	-x	use double-precision calculation
 *
 *  Bugs:
 *
 *
 *  Future Updates
 *    Use drms_keyword_lookup()
 *    Create new data series as needed
 *    Parallelize over multiple records
 *    Add recreation option for target output records
 *    Add parameters for overriding essential input keywords
 *    Add parameter for overriding keyword list for propagation
 *    Modify to also work with a series of two-dimensional images - maybe
 *	not such a good idea
 *
 *  Revision history is at end of file.
 */



/**
   @defgroup jsoc_fld prepares map for the fld code
   @ingroup su_util

   @brief .

   @par Synopsis:
   @code
   jsoc_fld in=input_data out=output_data
   @endcode

   This is a general purpose module that takes a series of input 
   data.

   The image is not registered to solar center by this module.  The image will be rotated by
   a flip-flip procedure is the CROTA2 parameter is near +-180.0 unless the -u flag (unchanged) is
   present.  If the -c flag is present the image will be cropped at the solar limb before scaling.
   If the -h flag or requestid parameter is present the output segments will have full headers.


   @par Flags:
   @c
   -c  Crop before scaling.  Use rsun_obs/cdelt1 for limb radius.
   -h  Write full FITS headers.
   -u  Leave unchanged, do NOT rotate by 180 degrees if CROTA2 is near 180.  default is to do flip-flip method so
       image is norths up and no pixel values are changed.

   @par GEN_FLAGS:
   Ubiquitous flags present in every module.
   @ref jsoc_main

   @param in  The input data series.
   @param out The output series.

   @par Exit_Status:

   @par Example:

   @code
   jsoc_map4fld -c in=hmi.v_45s[2013.12.18_00] out='kis_vigeesh.v_bin' width=120 height=120
   @endcode

   @bug
   None known so far.

 */

// All the nessasary include files
#include "drms.h"
#include "jsoc.h"
#include "jsoc_main.h"
#include "keystuff.c"


#include <fftw3.h>
#include <cfitsio.h> //Angular Brackets for library files
#include <hdf5.h>
#include <mpi.h>


#include "fldjob.h"
#include "decomp.h"
#include "backend.h"

#define HFILE "freq_lmg.h5"


/*http://stackoverflow.com/questions/24883461/hdf5-updating-a-cell-in-a-table-of-integers */
#define VERIFY(a) do { if((a)<0) { fprintf(stderr,"Failure line %d.\n",__LINE__); exit(-1);}}while(0)


// The name of the module
char *module_name = "fld";
char *module_desc = "HMI Dopplergrams to FLD Coefficients";
char *version_id = "0.1";

/*  global declaration of missing to be initialized as NaN  */
float missing_val=0.0 / 0.0;;
#define INTERP_NEAREST_NEIGHBOR	(1)
#define INTERP_BILINEAR	(2)
#define INTERP_CUBIC_CONVOLUTION (3)
#define INTERP_RESAMPLE (4)

//
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
		{ARG_STRING,	"segment", "Not Specified",
				"input data series segment; ignored if series only has one 3-d segment"},
				{ARG_FLAG, "c", "0", "Crop at rsun_obs."},
				{ARG_FLAG, "u", "0", "do not rotate by 180 if needed."},
				{ARG_FLAG, "t", "0", "testing only, no files written"},
				{ARG_INT, "inseg", "0", "Input segment number"},
				{ARG_INT, "outseg", "0", "Output segment number"},
				{ARG_FLOAT, "width", "0", "Width of the map section in degrees"},
				{ARG_FLOAT, "height", "0", "Height of the map section in degrees"},
				{ARG_INT, "wpix", "0", "Width of the map section in degrees"},
				{ARG_INT, "hpix", "0", "Height of the map section in degrees"},
				{ARG_FLOAT, "clat", "0", "Latitude of the center"},
				{ARG_FLOAT, "clon", "0", "Longitude of the center"},
				{ARG_STRING, "decomp", "FLD", "Decomposition method: SHD, FHD, FLD"},
				{ARG_STRING, "apodf", "Tukey", "Apodization function: Tukey"},
				{ARG_INT, "apodw", "8", "Apodization width"},
				{ARG_INT, "apodm", "0", "Apodization margin"},
				{ARG_INT, "lpfwhm", "0", "Low pass FWHM"},
				{ARG_INT, "m", "0", "Absolute m range needed"},
				{ARG_INT, "l_max", "0", "Maximum l range needed"},
				{ARG_INT, "interp", "2", "Map Interpolation: 1-Nearest, 2-Bilinear, 3-Cubic-Conv, 4-AST"},
				{ARG_STRING, "requestid", "NA", "RequestID if called as an export processing step."},
				{ARG_END}
};

#define     Deg2Rad    (M_PI/180.0)
#define     Rad2arcsec (3600.0/Deg2Rad)
#define     arcsec2Rad (Deg2Rad/3600.0)
#define     Rad2Deg    (180.0/M_PI)


char *propagate[] = {"CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
		"MidTime", "Duration", "LonSpan", "T_START", "T_STOP", "Coverage", "Quality",
		"MapScale", "MapProj","Width", "Height", "lMax","mAbs", "Cadence"};

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


//trying to get fldcomplex_t into the hd5
typedef struct {
     float re;   /*real part */
     float im;   /*imaginary part */
} complex_t;
complex_t tmp;  /*used only to compute offsets */


// fft part comlex 1d dft
//need to change thhis toa more elegent gcc complex.h format



void rebinArraySF(DRMS_Array_t *out, DRMS_Array_t *in);
int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);
int crop_image(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);
const char *get_input_recset(DRMS_Env_t *drms_env, const char *in);
int map4fld();

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus);
ObsInfo_t *GetMinObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus);

static int cleanup (int error, DRMS_RecordSet_t *irecs, DRMS_Record_t *orec,
		DRMS_Array_t *orig, DRMS_Array_t *pspec, int dispose) {
	if (orig) drms_free_array (orig);
	if (pspec) drms_free_array (pspec);
	if (irecs) drms_close_records (irecs, DRMS_FREE_RECORD);
	if (orec) {
		if (error) drms_close_record (orec, DRMS_FREE_RECORD);
		else drms_close_record (orec, dispose);
	}
	return error;
}

int DoIt(void)
{
	//DRMS Related stuff
	int status = DRMS_SUCCESS;
	DRMS_RecordSet_t *inRS, *outRS;
	DRMS_Record_t *irec, *orec = NULL;
	DRMS_Segment_t *iseg, *record_segment, *oseg, *logseg;
	DRMS_Array_t *orig = NULL, *pspec = NULL;
	CmdParams_t *params = &cmdparams;

	char *inStr;
	//Input arguments related stuff
	const char *inQuery = params_get_str(&cmdparams, "in");
	const char *outSeries = params_get_str(&cmdparams, "out");
	const char *requestid = params_get_str(&cmdparams, "requestid");
	int crop = params_get_int(&cmdparams, "c");
	int as_is = params_get_int(&cmdparams, "u");
	int inseg = params_get_int(&cmdparams, "inseg");
	int outseg = params_get_int(&cmdparams, "outseg");
	int full_header = params_get_int(&cmdparams, "h") || strcmp(requestid, "NA");
	int dp_calc = params_isflagset (params, "x");
	int verbose = params_isflagset (params, "v");
	int no_save = params_isflagset (params, "n");
	char *seg_name = strdup (params_get_str (params, "segment"));
	double apode_edge = params_get_double (params, "apodize");

	double x0,y0;

	float roi_lon,roi_lat,roi_width,roi_height;

	const float apod_width = params_get_float(&cmdparams,"apodw");
	const float apod_margin = params_get_float(&cmdparams,"apodm");
	const char *apod_func = params_get_str(&cmdparams,"apodf");
	const float lowpass_fwhm = params_get_float(&cmdparams,"lpfwhm");

	int l_max;
	int m_abs;
	int key_n, kstat, propct;
	char **copykeylist;
	int keyct = sizeof (propagate) / sizeof (char *);

	char *osegname;
	int verbose_logs;
	char logfilename[DRMS_MAXPATHLEN];
	double dval;
	char module_ident[64];

	int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;
	int  rgnct, segs, segct, isegnum,  found;

	char *inser;
	int nrecs;

	// TODO: Change this to a more elegant intput stream
	int nmaps;
	float maplatstp;
	float maplonstp;


	TIME t_start;
	TIME t_end;
	TIME t_stop;


	//input fits cube file related stuff
	char fname[DRMS_MAXPATHLEN];

	char source[DRMS_MAXQUERYLEN];

	char *dir;
	char *fpname;

	char outhdf[DRMS_MAXPATHLEN];

	/* Initializing the Input record */
	inRS = drms_open_records (drms_env, inQuery, &status);
	if (status) {
		fprintf (stderr, "Error (%s): drms_open_records() returned %d for dataset:\n",
				module_ident, status);
		fprintf (stderr, "  %s\n", inQuery);
		return cleanup (1, inRS, orec, orig, pspec, dispose);
	}

	rgnct = inRS->n;
	if (rgnct < 1) {
		fprintf (stderr, "No records found in input data set %s\n", inQuery);
		return cleanup (1, inRS, orec, orig, pspec, dispose);
	}

	irec = inRS->records[0];
	inser = strdup (inQuery);
	if ((segs = drms_record_numsegments (irec)) < 1) {
		fprintf (stderr, "Error: no data segments in input data series:\n");
		fprintf (stderr, "  %s\n", inser);
		return cleanup (1, inRS, orec, orig, pspec, dispose);
	}
	found = 0;
	//the series has two data segments 'Log' and 'V'

	/* Getting the full file path of the input series */
	iseg = drms_segment_lookupnum (irec, 0);
	drms_record_directory(irec, &fpname, 1);
	//strcpy(&fpname,&fname);
	strcat(&fpname,"/");
	strcat(&fpname,iseg->info->name);
	strcat(&fpname,".h5");
	printf("Input file\t%s\n",&fpname);


	roi_lon=drms_getkey_float(irec, "LonHG",&status);
	roi_lat=drms_getkey_float(irec, "LatHG",&status);
	roi_width=drms_getkey_float(irec, "Width",&status);
	roi_height=drms_getkey_float(irec, "Height",&status);
	maplonstp=drms_getkey_double(irec, "CDELT1",&status);
	maplatstp=drms_getkey_double(irec, "CDELT2",&status);
	x0=drms_getkey_double(irec, "CRPIX1",&status);
	y0=drms_getkey_double(irec, "CRPIX2",&status);

	t_start=drms_getkey_time(irec, "T_START",&status);
	t_end =drms_getkey_time(irec, "T_STOP",&status);

	l_max =drms_getkey_int(irec, "lMax",&status);
	m_abs =drms_getkey_int(irec, "mAbs",&status);

	char *mapproj;
	mapproj=drms_getkey_string(irec,"MapProj",&status);

	//TODO: Needs to be propagated
	if (strcmp(mapproj,"PlateCarree"))
			DIE("-----------\n\n\n\nRight now, only works for \"PlateCarree\" Projection\n\n\n\n\n-----------\n");


	/* Initializing the Output record */
	if (!(outRS = drms_create_records (drms_env, rgnct, outSeries, DRMS_PERMANENT,
			&status))) {
		fprintf (stderr, "Error: unable to create %d records in series %s\n",
				rgnct, outSeries);
		fprintf (stderr, "       drms_create_records() returned status %d\n", status);
		return 1;
	}
	if (verbose) printf ("creating %d record(s) in series %s\n", rgnct, outSeries);
	if (verbose && dispose == DRMS_FREE_RECORD)
		printf ("experimental run, output records will not be saved\n");
	/*  check output data series struct  */
	orec = drms_recordset_getrec (outRS, 0);
	if (!orec) {
		fprintf (stderr, "Error accessing record %d in series %s\n", 0, outSeries);
		drms_close_records (outRS, DRMS_FREE_RECORD);
		return 1;
	}

	segct = drms_record_numsegments (orec);
	found = 0;
	for (int n = 0; n < segct; n++) {
		record_segment = drms_segment_lookupnum (orec, n);
		if (!found) osegname = strdup (record_segment->info->name);
		found++;
	}

	if (found < 1) {
		fprintf (stderr,
				"Error: no data segment of dimension 3 and appropriate size in output series %s\n", outSeries);
		drms_close_records (outRS, DRMS_FREE_RECORD);
		return 1;
	}
	record_segment = drms_segment_lookup (orec, osegname);
	if (found > 1) {
		fprintf (stderr,
				"Warning: multiple data segments of dimension 3 and appropriate size in output series %s\n", outSeries);
		fprintf (stderr, "       using \"%s\"\n", osegname);
	}

	logseg = drms_segment_lookup (orec, "Log");
	if (logseg) drms_segment_filename (logseg, logfilename);
	else if (verbose) {
		fprintf (stderr,
				"Warning: segment \"Log\" not present in output series %s\n", outSeries);
		fprintf (stderr, "         verbose logging turned off\n");
		verbose_logs = 0;
	}

	if (verbose_logs) {
		logseg = drms_segment_lookup (orec, "Log");
		if (logseg) drms_segment_filename (logseg, logfilename);
		//fopen (logfilename, "w");
	}


	/* This is where the reading and fft begins */
	float duration;
	float cadence;
	duration = drms_getkey_float(irec, "Duration",&status);
	cadence = drms_getkey_double(irec, "Cadence",&status);
	int nframes = (int) duration/cadence;

	long nl = l_max+1; 					//TODO: needs to be propagated
	long nm = 1 + (2 * m_abs);          //TODO: needs to be propagated
	//int maps_per_chunk =2;
	//int nchunks=3;                      //TODO: needs to be propagated
	long fullframe = nframes; //TODO: needs to be propagated
	complex_t *timedata;
	complex_t *freqdata;

	long nElemTot = nm * nl * fullframe;


    long naxes = 3;
    long naxis[3] = {nm, nl, fullframe};

    long mynaxes = 1;
    long mynaxis[1] = {fullframe};

    long timeElemTot = nm * nl * fullframe;

	/* Allocate memory for the time and frequency hyperslabs */
	timedata = (struct complex_t*) malloc (fullframe * sizeof (complex_t));
	freqdata = (struct complex_t*) malloc (fullframe * sizeof (complex_t));

	/* FFT Initialize (single precision only) */
	fftwf_plan fplan;
	fftwf_complex *inff, *outff;
	inff = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fullframe);
	outff = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fullframe);
	fplan = fftwf_plan_dft_1d(fullframe, inff, outff, FFTW_FORWARD, FFTW_ESTIMATE); //single precision

	/* Output filename */
	//TODO: Should come from the jsd file definiton
	strcpy(&outhdf,"freq_lmg.h5");

	if (access( &fpname, F_OK ) == 0 ) {

		/* HDF5 variable definitions */
		hid_t       hfile_id, hdataset_id, hdataspace_id, hdatatype, grp;
		hid_t       houtfile_id, houtset_id, houtspace_id, outgrp;
		hid_t		hmemspace_id;
		herr_t      hstatus;
		hid_t 		complex_id;

		/* Define the hyperslab in the input and the output dataset */
		hsize_t     count[3];              /* size of subset in the file */
		hsize_t     offset[3];             /* subset offset in the file */
		/* Define the hyperslab in the input and the output dataset */
		hsize_t     count_slab[1];              /* size of subset in the file */
		hsize_t     offset_slab[1];             /* subset offset in the file */

		/* Open the file */
		hfile_id = H5Fopen(&fpname, H5F_ACC_RDONLY, H5P_DEFAULT);
	    /* Open the group in the file */
	    grp  = H5Gopen2(hfile_id, "/lmg", H5P_DEFAULT);
	    /* Open the dataset in the file */
	    hdataset_id = H5Dopen2(hfile_id,"/lmg/Alm", H5P_DEFAULT);
	    /* Get Dataspace and allocate memory for read buffer */
		hdataspace_id = H5Dget_space(hdataset_id);
		/* Create the Datatype */
		complex_id = H5Tcreate (H5T_COMPOUND, sizeof(tmp));
		if (H5Tinsert (complex_id, "real", HOFFSET(complex_t, re), H5T_NATIVE_FLOAT))
			fprintf(stderr, "Unable to add a real part to the complex datatype.\n");
		if (H5Tinsert (complex_id, "imaginary", HOFFSET(complex_t, im), H5T_NATIVE_FLOAT))
			fprintf(stderr, "Unable to add a imaginary part to the complex datatype.\n");
		/* Create a memory space for hmemsapce_id to store the time slice*/
		hmemspace_id = H5Screate_simple (mynaxes, mynaxis, NULL);
		/* Step through the l's and m's to extract the time slab */

		/* Create a new hdf5 file */
		houtfile_id = H5Fcreate(HFILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		printf("File Created");
		//create a group in the file
		outgrp = H5Gcreate2(houtfile_id,"/lmg", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		/* Create a memory space for houtsapce_id */
		houtspace_id = H5Screate_simple (naxes, naxis, NULL);
		//Dataset for Alm's
		houtset_id = H5Dcreate2(houtfile_id, "/lmg/Alm", complex_id, houtspace_id,
				H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);

	    for (int step_m=0;step_m<nm;step_m++)
	    {
	    	for (int step_l=0;step_l<nl;step_l++)
	    	{
	    		/* Hyperslab selection */
	    		offset[0] = step_m;  // m
	    		offset[1] = step_l;  // l
	    		offset[2] = 0;

	    		count[0]  = 1;
	    		count[1]  = 1;						/* this has to be 1, coz. we need 1 element in that direction */
	    		count[2]  = fullframe;
	    	    /* Select the specific time hyperslab from the input data */
	    		if (H5Sselect_hyperslab (hdataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL))
	    			fprintf(stderr, "Unable to select the memory hyperslab.\n");
	    		/* Use the same offset and count as the input data */
	    		if (H5Sselect_hyperslab (houtspace_id, H5S_SELECT_SET, offset, NULL, count, NULL))
	    			fprintf(stderr, "Unable to select the memory hyperslab.\n");
	    		/* No offset for the time hyperslab */
	    		offset_slab[0]=0;
	    		count_slab[0]=fullframe;
	    		if (H5Sselect_hyperslab (hmemspace_id, H5S_SELECT_SET, offset_slab, NULL, count_slab, NULL))
	    			fprintf(stderr, "Unable to select the memory hyperslab.\n");
	    		/* Read the hyperslab (needle) */
	    		if (H5Dread(hdataset_id, complex_id, hmemspace_id, hdataspace_id, H5P_DEFAULT, timedata))
	    			fprintf(stderr, "Unable to read the time slab from the input dataspace.\n");
	        	//DIE("-----------\n\n\n\nDone Alm\n\n\n\n\n-----------\n");
	    	    /* compute FFT */
	    		/* copy the time axis to fft inff */
	    	    //TODO: This should be modified, such that the hdf5 can write native ansi c complex
	    	    for (int i=0; i<fullframe; i++)
	    	    {
	    		  	 	inff[i] = timedata[i].re + timedata[i].im * I;
	    	    }
	    	    fftwf_execute(fplan); /* repeat as needed */

	    		for (int i=0; i<fullframe; i++)
	    		{
	    			freqdata[i].re = crealf(outff[i]);
	    			freqdata[i].im = cimagf(outff[i]);
	    		}

	    	    /* Write the freq hyperslab to the output dataspace */
	    		//TODO: Check the frequency axis
	    		if (H5Dwrite(houtset_id, complex_id, hmemspace_id, houtspace_id, H5P_DEFAULT, freqdata))
	    			fprintf(stderr, "Unable to write the frequency slab to the output dataspace.\n");
	    	}
	    }

		/* Close the datatype */
		if (H5Tclose(complex_id))
			fprintf(stderr, "Unable to close complex datatype.\n");
    	/* Terminate access to the mem space. */
    	if (H5Sclose(hmemspace_id))
    		fprintf(stderr, "Unable to terminate access to memory dataspace.\n");
		/* Terminate access to the data space. */
    	if (H5Sclose(hdataspace_id))
    		fprintf(stderr, "Unable to terminate access to input dataspace.\n");
		if (H5Sclose(houtspace_id))
			fprintf(stderr, "Unable to terminate access to output dataspace.\n");
    	/* End access to the dataset and release resources used by it. */
    	if (H5Dclose(hdataset_id))
    		fprintf(stderr, "Unable to terminate access to input dataset.\n");
    	if (H5Dclose(houtset_id))
    		fprintf(stderr, "Unable to terminate access to output dataset.\n");




    	/* Do for the Blm */
    	hdataset_id = H5Dopen2(hfile_id,"/lmg/Blm", H5P_DEFAULT);
    	/* Get Dataspace and allocate memory for read buffer */
    	hdataspace_id = H5Dget_space(hdataset_id);
    	/* Create the Datatype */
    	complex_id = H5Tcreate (H5T_COMPOUND, sizeof(tmp));
    	if (H5Tinsert (complex_id, "real", HOFFSET(complex_t, re), H5T_NATIVE_FLOAT))
    		fprintf(stderr, "Unable to add a real part to the complex datatype.\n");
    	if (H5Tinsert (complex_id, "imaginary", HOFFSET(complex_t, im), H5T_NATIVE_FLOAT))
    		fprintf(stderr, "Unable to add a imaginary part to the complex datatype.\n");
    	/* Create a memory space for hmemsapce_id to store the time slice*/
    	hmemspace_id = H5Screate_simple (mynaxes, mynaxis, NULL);
    	/* Step through the l's and m's to extract the time slab */
    	houtspace_id = H5Screate_simple (naxes, naxis, NULL);
    	//Dataset for Alm's
    	houtset_id = H5Dcreate2(houtfile_id, "/lmg/Blm", complex_id, houtspace_id,
    			H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);

    	for (int step_m=0;step_m<nm;step_m++)
    	{
    		for (int step_l=0;step_l<nl;step_l++)
    		{
    			/* Hyperslab selection */
    			offset[0] = step_m;  // m
    			offset[1] = step_l;  // l
    			offset[2] = 0;

    			count[0]  = 1;
    			count[1]  = 1;						/* this has to be 1, coz. we need 1 element in that direction */
    			count[2]  = fullframe;
    			/* Select the specific time hyperslab from the input data */
    			if (H5Sselect_hyperslab (hdataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL))
    				fprintf(stderr, "Unable to select the memory hyperslab.\n");
    			/* Use the same offset and count as the input data */
    			if (H5Sselect_hyperslab (houtspace_id, H5S_SELECT_SET, offset, NULL, count, NULL))
    				fprintf(stderr, "Unable to select the memory hyperslab.\n");
    			/* No offset for the time hyperslab */
    			offset_slab[0]=0;
    			count_slab[0]=fullframe;
    			if (H5Sselect_hyperslab (hmemspace_id, H5S_SELECT_SET, offset_slab, NULL, count_slab, NULL))
    				fprintf(stderr, "Unable to select the memory hyperslab.\n");
    			/* Read the hyperslab (needle) */
    			if (H5Dread(hdataset_id, complex_id, hmemspace_id, hdataspace_id, H5P_DEFAULT, timedata))
    				fprintf(stderr, "Unable to read the time slab from the input dataspace.\n");

    			/* compute FFT */
    			/* copy the time axis to fft inff */
    			//TODO: This should be modified, such that the hdf5 can write native ansi c complex
    			for (int i=0; i<fullframe; i++)
    			{
    				inff[i] = timedata[i].re + timedata[i].im * I;
    			}
    			fftwf_execute(fplan); /* repeat as needed */

    			for (int i=0; i<fullframe; i++)
    			{
    				freqdata[i].re = crealf(outff[i]);
    				freqdata[i].im = cimagf(outff[i]);
    			}

    			/* Write the freq hyperslab to the output dataspace */
    			//TODO: Check the frequency axis
    			if (H5Dwrite(houtset_id, complex_id, hmemspace_id, houtspace_id, H5P_DEFAULT, freqdata))
    				fprintf(stderr, "Unable to write the frequency slab to the output dataspace.\n");
    		}
    	}

    	/* Close the datatype */
    	if (H5Tclose(complex_id))
    		fprintf(stderr, "Unable to close complex datatype.\n");
    	/* Terminate access to the mem space. */
    	if (H5Sclose(hmemspace_id))
    		fprintf(stderr, "Unable to terminate access to memory dataspace.\n");
    	/* Terminate access to the data space. */
    	if (H5Sclose(hdataspace_id))
    		fprintf(stderr, "Unable to terminate access to input dataspace.\n");
    	if (H5Sclose(houtspace_id))
    		fprintf(stderr, "Unable to terminate access to output dataspace.\n");
    	/* End access to the dataset and release resources used by it. */
    	if (H5Dclose(hdataset_id))
    		fprintf(stderr, "Unable to terminate access to input dataset.\n");
    	if (H5Dclose(houtset_id))
    		fprintf(stderr, "Unable to terminate access to output dataset.\n");
    	if (H5Gclose(outgrp))
    		fprintf(stderr, "Unable to close output group.\n");

    	/* Close the group */
    	if (H5Gclose(grp))
    		fprintf(stderr, "Unable to close input group.\n");

    	/* Close the file. */
    	if (H5Fclose(hfile_id))
    		fprintf(stderr, "Unable to close input hdf5 file.\n");
    	if (H5Fclose(houtfile_id))
    		fprintf(stderr, "Unable to close output hdf5 file.\n");
	}

	/* FFTW3 Finilize */
	fftwf_destroy_plan(fplan);
	fftwf_free(inff);
	fftwf_free(outff);



    hid_t dataspace_id,file_id,grp, attribute_id;
    hid_t atype;
	/* Writing the Attribute for the HDF5 file */
    file_id = H5Fopen(HFILE, H5F_ACC_RDWR, H5P_DEFAULT);

    //open the group in the file
    grp  = H5Gopen2(file_id, "/lmg", H5P_DEFAULT);

    long nrank = 1;
    long ndims[1] = {1};
    dataspace_id = H5Screate_simple(nrank, ndims, NULL);
    float	attr_flt[1];
    double	attr_dbl[1];
    int		attr_int[1];
    char	attr_char[1];
    double result;

    for (int n = 0; n < keyct; n++)
    {
    	if (drms_keyword_type(drms_keyword_lookup(irec,propagate[n],1))==DRMS_TYPE_FLOAT)
    	{
    		dval = drms_getkey_float (irec, propagate[n], &status);
    		attr_flt[0]=dval;
    		attribute_id = H5Acreate(grp, propagate[n], H5T_IEEE_F32BE, dataspace_id,
    				H5P_DEFAULT, H5P_DEFAULT);
    		status=H5Awrite(attribute_id,H5T_NATIVE_FLOAT,attr_flt);
    		status=H5Aclose(attribute_id);
    	}
    	if (drms_keyword_type(drms_keyword_lookup(irec,propagate[n],1))==DRMS_TYPE_DOUBLE)
    	{
    		dval = drms_getkey_double (irec, propagate[n], &status);
    		attr_dbl[0]=dval;
    		attribute_id = H5Acreate(grp, propagate[n], H5T_IEEE_F64LE, dataspace_id,
    				H5P_DEFAULT, H5P_DEFAULT);
    		status=H5Awrite(attribute_id,H5T_NATIVE_DOUBLE,attr_dbl);
    		status=H5Aclose(attribute_id);
    	}
    	if (drms_keyword_type(drms_keyword_lookup(irec,propagate[n],1))==DRMS_TYPE_INT)
    	{
    		dval = drms_getkey_int (irec, propagate[n], &status);
    		attr_int[0]=dval;
    		attribute_id = H5Acreate(grp, propagate[n], H5T_STD_I32BE, dataspace_id,
    		   				H5P_DEFAULT, H5P_DEFAULT);
    		status=H5Awrite(attribute_id,H5T_NATIVE_INT,attr_int);
    		status=H5Aclose(attribute_id);
    	}
    	if (drms_keyword_type(drms_keyword_lookup(irec,propagate[n],1))==DRMS_TYPE_TIME)
    	{
    		dval = drms_getkey_time (irec, propagate[n], &status);
    		result = drms2time(DRMS_TYPE_TIME, &dval, status);
    		attr_dbl[0]=result;
    		attribute_id = H5Acreate(grp, propagate[n], H5T_IEEE_F64LE, dataspace_id,
    				H5P_DEFAULT, H5P_DEFAULT);
    		status=H5Awrite(attribute_id,H5T_NATIVE_DOUBLE,attr_dbl);
    		status=H5Aclose(attribute_id);
    	}
    }


    status=H5Sclose(dataspace_id);
    status=H5Gclose(grp);
    status=H5Fclose(file_id);


	/* Writing the record to the drms */
	if(drms_segment_write_from_file(record_segment,&outhdf))
		fprintf(stderr, "Can not write hdf5 file for the segment.\n");

    propagate_keys (orec, irec, propagate, keyct);

    drms_sprint_rec_query (source, irec);
	inser = strdup (inQuery);
    check_and_set_key_str   (orec, "Source", source);
    check_and_set_key_str   (orec, "Module", module_ident);
    check_and_set_key_str   (orec, "BLD_VERS", jsoc_version);
    check_and_set_key_str   (orec, "Input", inser);
    check_and_set_key_time  (orec, "Created", CURRENT_SYSTEM_TIME);



	/* Close the input record */
	if(drms_close_records(inRS, DRMS_FREE_RECORD))
		fprintf(stderr, "Error closing Record for input series\n");

	/* Close the output record */
	if(drms_close_records(outRS, DRMS_INSERT_RECORD))
		fprintf(stderr, "Error closing Record for output series\n");
	printf("\n---END---\n");

	return (DRMS_SUCCESS);
} // end of DoIt

// ----------------------------------------------------------------------

#define CHECK(keyname) {if (status) {fprintf(stderr,"Keyword failure to find: %s, status=%d\n",keyname,status); *rstatus=status; return(NULL);}}

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus)
{
	TIME t_prev;
	DRMS_Record_t *rec;
	TIME t_obs;
	//double dv;
	ObsInfo_t *ObsLoc;
	int status;

	if (!seg || !(rec = seg->record))
	{ *rstatus = 1; return(NULL); }

	ObsLoc = (pObsLoc ? pObsLoc : (ObsInfo_t *)malloc(sizeof(ObsInfo_t)));
	if (!pObsLoc)
		memset(ObsLoc, 0, sizeof(ObsInfo_t));

	t_prev = ObsLoc->t_obs;
	t_obs = drms_getkey_time(rec, "T_OBS", &status); CHECK("T_OBS");

	if (t_obs <= 0.0)
	{ *rstatus = 2; return(NULL); }

	if (t_obs != t_prev)
	{
		ObsLoc->crpix1 = drms_getkey_double(rec, "CRPIX1", &status); CHECK("CRPIX1");
		ObsLoc->crpix2 = drms_getkey_double(rec, "CRPIX2", &status); CHECK("CRPIX2");
		ObsLoc->crval1 = drms_getkey_double(rec, "CRVAL1", &status); CHECK("CRVAL1");
		ObsLoc->crval2 = drms_getkey_double(rec, "CRVAL2", &status); CHECK("CRVAL2");
		ObsLoc->cdelt1 = drms_getkey_double(rec, "CDELT1", &status); CHECK("CDELT1");
		ObsLoc->cdelt2 = drms_getkey_double(rec, "CDELT2", &status); CHECK("CDELT1");
		ObsLoc->crota2 = drms_getkey_double(rec, "CROTA2", &status); if (status) ObsLoc->crota2 = 0.0; // WCS default
		ObsLoc->sina = sin(ObsLoc->crota2*Deg2Rad);
		ObsLoc->cosa = sqrt (1.0 - ObsLoc->sina*ObsLoc->sina);
		ObsLoc->rsun_obs = drms_getkey_double(rec, "RSUN_OBS", &status);
		if (status)
		{
			double dsun_obs = drms_getkey_double(rec, "DSUN_OBS", &status); CHECK("DSUN_OBS");
			ObsLoc->rsun_obs = asin(696000000.0/dsun_obs)/arcsec2Rad;
		}
		ObsLoc->obs_vr = drms_getkey_double(rec, "OBS_VR", &status); CHECK("OBS_VR");
		ObsLoc->obs_vw = drms_getkey_double(rec, "OBS_VW", &status); CHECK("OBS_VW");
		ObsLoc->obs_vn = drms_getkey_double(rec, "OBS_VN", &status); CHECK("OBS_VN");
		ObsLoc->obs_b0 = drms_getkey_double(rec, "CRLT_OBS", &status); CHECK("CRLT_OBS");
		ObsLoc->t_obs = t_obs;
	}
	*rstatus = 0;
	return(ObsLoc);
}

/* GetMinObsInfo - gets minimum standard WCS keywords for e.g. heliographic mapped data */
ObsInfo_t *GetMinObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus)
{
	TIME t_prev;
	DRMS_Record_t *rec;
	TIME t_obs;
	//double dv;
	ObsInfo_t *ObsLoc;
	int status;

	if (!seg || !(rec = seg->record))
	{ *rstatus = 1; return(NULL); }

	ObsLoc = (pObsLoc ? pObsLoc : (ObsInfo_t *)malloc(sizeof(ObsInfo_t)));
	if (!pObsLoc)
		memset(ObsLoc, 0, sizeof(ObsInfo_t));

	t_prev = ObsLoc->t_obs;
	t_obs = drms_getkey_time(rec, "T_OBS", &status); CHECK("T_OBS");

	if (t_obs <= 0.0)
	{ *rstatus = 2; return(NULL); }

	if (t_obs != t_prev)
	{
		ObsLoc->crpix1 = drms_getkey_double(rec, "CRPIX1", &status); CHECK("CRPIX1");
		ObsLoc->crpix2 = drms_getkey_double(rec, "CRPIX2", &status); CHECK("CRPIX2");
		ObsLoc->crval1 = drms_getkey_double(rec, "CRVAL1", &status); CHECK("CRVAL1");
		ObsLoc->crval2 = drms_getkey_double(rec, "CRVAL2", &status); CHECK("CRVAL2");
		ObsLoc->cdelt1 = drms_getkey_double(rec, "CDELT1", &status); CHECK("CDELT1");
		ObsLoc->cdelt2 = drms_getkey_double(rec, "CDELT2", &status); CHECK("CDELT1");
		ObsLoc->crota2 = drms_getkey_double(rec, "CROTA2", &status); if (status) ObsLoc->crota2 = 0.0; // WCS default
		ObsLoc->sina = sin(ObsLoc->crota2*Deg2Rad);
		ObsLoc->cosa = sqrt (1.0 - ObsLoc->sina*ObsLoc->sina);
		ObsLoc->rsun_obs = drms_getkey_double(rec, "RSUN_OBS", &status);
		ObsLoc->obs_vr = drms_getkey_double(rec, "OBS_VR", &status);
		ObsLoc->obs_vw = drms_getkey_double(rec, "OBS_VW", &status);
		ObsLoc->obs_vn = drms_getkey_double(rec, "OBS_VN", &status);
		ObsLoc->obs_b0 = drms_getkey_double(rec, "CRLT_OBS", &status);
		ObsLoc->t_obs = t_obs;
	}
	*rstatus = 0;
	return(ObsLoc);
}

// ----------------------------------------------------------------------

// In cases known to not have compact slotted series and cadence is specified
// generate explicit recordset list of closest good record to desired grid
// First get vector of times and quality
// Then if vector is not OK, quit.
// then: make temp file to hold recordset list
//       start with first time to define desired grid,
//       make array of desired times.
//       make empty array of recnums
//       search vector for good images nearest desired times
//       for each found time, write record query


const char *get_input_recset(DRMS_Env_t *drms_env, const char *inQuery)
{
	static char newInQuery[102];
	TIME epoch = (cmdparams_exists(&cmdparams, "epoch")) ? params_get_time(&cmdparams, "epoch") : 0;
	DRMS_Array_t *data;
	TIME t_start, t_stop, t_now, t_want, t_diff, this_t_diff;
	int status = 1;
	int nrecs, irec;
	int nslots, islot;
	long long *recnums;
	TIME *t_this, half;
	TIME cadence;
	double *drecnum, *dquality;
	int quality;
	long long recnum;
	char keylist[DRMS_MAXQUERYLEN];
	static char filename[100];
	char *tmpdir;
	FILE *tmpfile;
	char newIn[DRMS_MAXQUERYLEN];
	char seriesname[DRMS_MAXQUERYLEN];
	char *lbracket;
	char *at = index(inQuery, '@');
	if (at && *at && (strncmp(inQuery,"aia.lev1[", 9)==0 ||
			strncmp(inQuery,"hmi.lev1[", 9)==0 ||
			strncmp(inQuery,"aia.lev1_nrt2[",14)==0 ||
			strncmp(inQuery,"hmi.lev1_nrt[", 13)==0 ))
	{
		char *ip=(char *)inQuery, *op=newIn, *p;
		long n, mul;
		while ( *ip && ip<at )
			*op++ = *ip++;
		ip++; // skip the '@'
		n = strtol(ip, &p, 10); // get digits only
		if (*p == 's') mul = 1;
		else if (*p == 'm') mul = 60;
		else if (*p == 'h') mul = 3600;
		else if (*p == 'd') mul = 86400;
		else
		{
			fprintf(stderr,"cant make sense of @xx cadence spec for aia or hmi lev1 data");
			return(NULL);
		}
		cadence = n * mul;
		ip = ++p;  // skip cadence multiplier
		while ( *ip )
			*op++ = *ip++;
		*op = '\0';
		half = cadence/2.0;
		sprintf(keylist, "T_OBS,QUALITY,recnum");
		data = drms_record_getvector(drms_env, newIn, keylist, DRMS_TYPE_DOUBLE, 0, &status);
		if (!data || status)
		{
			fprintf(stderr,"getkey_vector failed status=%d\n", status);
			return(NULL);
		}
		nrecs = data->axis[1];
		irec = 0;
		t_this = (TIME *)data->data;
		dquality = (double *)data->data + 1*nrecs;
		drecnum = (double *)data->data + 2*nrecs;
		if (epoch > 0.0)
		{
			int s0 = (t_this[0] - epoch)/cadence;
			TIME t0 = s0*cadence + epoch;
			t_start = (t0 < t_this[0] ? t0 + cadence : t0);
		}
		else
			t_start = t_this[0];
		t_stop = t_this[nrecs-1];
		nslots = (t_stop - t_start + cadence/2)/cadence;
		recnums = (long long *)malloc(nslots*sizeof(long long));
		for (islot=0; islot<nslots; islot++)
			recnums[islot] = 0;
		islot = 0;
		t_want = t_start;
		t_diff = 1.0e9;
		for (irec = 0; irec<nrecs; irec++)
		{
			t_now = t_this[irec];
			quality = (int)dquality[irec] & 0xFFFFFFFF;
			recnum = (long long)drecnum[irec];
			this_t_diff = fabs(t_now - t_want);
			if (quality < 0)
				continue;
			if (t_now <= (t_want-half))
				continue;
			while (t_now > (t_want+half))
			{
				islot++;
				if (islot >= nslots)
					break;
				t_want = t_start + cadence * islot;
				this_t_diff = fabs(t_now - t_want);
				t_diff = 1.0e8;
			}
			if (this_t_diff <= t_diff)
				recnums[islot] = recnum;
			t_diff = fabs(t_now - t_want);
		}
		if (islot+1 < nslots)
			nslots = islot+1;  // take what we got.
		strcpy(seriesname, inQuery);
		lbracket = index(seriesname,'[');
		if (lbracket) *lbracket = '\0';
		tmpdir = getenv("TMPDIR");
		if (!tmpdir) tmpdir = "/tmp";
		sprintf(filename, "%s/hg_patchXXXXXX", tmpdir);
		mkstemp(filename);
		tmpfile = fopen(filename,"w");
		for (islot=0; islot<nslots; islot++)
			if (recnums[islot])
				fprintf(tmpfile, "%s[:#%lld]\n", seriesname, recnums[islot]);
		fclose(tmpfile);
		free(recnums);
		drms_free_array(data);
		sprintf(newInQuery,"@%s", filename);
		return(newInQuery);
	}
	else
		return(inQuery);
}



int drms_segment_write_hdf5(DRMS_Segment_t *seg, char *infile) {
  char *filename;            /* filename without path */
  char outfile[DRMS_MAXPATHLEN];
  FILE *in, *out;            /* input and output file stream */
  size_t read_size;          /* number of bytes on last read */
  const unsigned int bufsize = 16*1024;
  char *buf = malloc(bufsize*sizeof(char)); /* buffer for data */

  if (seg->info->scope == DRMS_CONSTANT &&
      seg->info->cseg_recnum) {
    fprintf(stderr, "ERROR in drms_segment_write: constant segment has already"
	    " been initialized. Series = %s.\n",  seg->record->seriesinfo->seriesname);
    return DRMS_ERROR_INVALIDACTION;
  }

  if (seg->record->readonly) {
    fprintf(stderr, "ERROR in drms_segment_write_from_file: Can't use "
	    "on readonly segment\n");
    return DRMS_ERROR_RECORDREADONLY;
  }

  // check protocol
  if (seg->info->protocol != DRMS_GENERIC)  {
    fprintf(stderr, "ERROR in drms_segment_write_from_file: Can't use "
	    "on non-DRMS_GENERIC segment.  Series = %s.\n", seg->record->seriesinfo->seriesname);
    return DRMS_ERROR_INVALIDACTION;
  }


  // strip path from infile
  filename = rindex(infile, '/');
  if (filename)
    filename++;
  else
    filename = infile;

  CHECKSNPRINTF(snprintf(outfile, DRMS_MAXPATHLEN, "%s/" DRMS_SLOTDIR_FORMAT "/%s",
			 seg->record->su->sudir, seg->record->slotnum, filename), DRMS_MAXPATHLEN);

  printf("the outfile is = %s",&outfile);
  hid_t       hfile_id, hdataset_id, hdataspace_id;
  hsize_t     hdims[2];
  herr_t      hstatus;

  //

  hfile_id = H5Fcreate(&outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Create the data space for the dataset. */
  hdims[0] = 4;
  hdims[1] = 6;
  hdataspace_id = H5Screate_simple(2, hdims, NULL);

  /* Create the dataset. */
  hdataset_id = H5Dcreate(hfile_id, "/dset", H5T_STD_I32BE, hdataspace_id,
		  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* End access to the dataset and release resources used by it. */
  hstatus = H5Dclose(hdataset_id);

  /* Terminate access to the data space. */
  hstatus = H5Sclose(hdataspace_id);

  /* Close the file. */
  hstatus = H5Fclose(hfile_id);

  CHECKSNPRINTF(snprintf(seg->filename, DRMS_MAXSEGFILENAME, "%s", filename), DRMS_MAXSEGFILENAME);
  free(buf);
  buf = NULL;

  if (seg->info->scope == DRMS_CONSTANT &&
      !seg->info->cseg_recnum) {

    if (seg->record->lifetime == DRMS_TRANSIENT) {
      fprintf(stderr, "Error: cannot set constant segment in a transient record\n");
      goto bailout;
   }
    return drms_segment_set_const(seg);
  }

  return DRMS_SUCCESS;

 bailout1:
  unlink(outfile);
 bailout:
  if (buf)
    free(buf);
  return 1;
}


int drms_segment_set_const(DRMS_Segment_t *seg) {
  XASSERT(seg->info->scope == DRMS_CONSTANT &&
	  !seg->info->cseg_recnum);
  DRMS_Record_t *rec = seg->record;
  seg->info->cseg_recnum = rec->recnum;

  // write back to drms_segment table
  char stmt[DRMS_MAXQUERYLEN];
  char *namespace = ns(rec->seriesinfo->seriesname);
  sprintf(stmt, "UPDATE %s." DRMS_MASTER_SEGMENT_TABLE
	  " SET cseg_recnum = %lld WHERE seriesname = '%s' and segmentname = '%s'",
	  namespace, rec->recnum, rec->seriesinfo->seriesname, seg->info->name);
  free(namespace);
  if(drms_dms(rec->env->session, NULL, stmt)) {
    fprintf(stderr, "Failed to update drms_segment table for constant segment\n");
    return DRMS_ERROR_QUERYFAILED;
  }
  return DRMS_SUCCESS;
}

