/*
 Copyright (C)	2014-2016 Vigeesh Gangadharan <vigeesh@kis.uni-freiburg.de>,
	Kiepenheuer-Institut fuer Sonnenphysik, Freiburg, Germany.
 Copyright (C) 2009-2010 Hans-Peter Doerr <doerr@kis.uni-freiburg.de>,
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
 *  tld.c
 *
 *  Calculate the Time series of Legendre Coefficients of an input dataset
 *  The module takes in a 3-dimensional real dataset from the mtrack
 *
 *  Responsible:
 *      Vigeesh Gangadharan                             vigeesh@kis.uni-freiburg.de
 *
 *  Usage:
 *    tld [-ts] in= out= ...
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
 *	-t	test runs
 *	-s  save fits
 *
 *  Bugs:
 *
 *
 *  Future Updates
 *
 *
 *  Revision history
 *     2014: VG - Adapted HPs code for drms 
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


// All the necessary include files
#include "drms.h"
#include "jsoc.h"
#include "jsoc_main.h"

#include "keystuff.c"

#include "hdf5.h"
#include "fitsio.h"
#include <mpi.h>
#include "fldjob.h"
#include "decomp.h"
#include "backend.h"

#define HFILE "lmg.h5"					/*TODO: VG- Remove Hardcoded*/


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

//enum DecompMode {NA = 0, FLD, FHD, SHD};

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
		{ARG_INT, "m", "25", "Absolute m range needed"},
		{ARG_INT, "l_max", "180", "Maximum l range needed"},
		{ARG_INT, "interp", "2", "Map Interpolation: 1-Nearest, 2-Bilinear, 3-Cubic-Conv, 4-AST"},
		{ARG_STRING, "copy",  "+", "comma separated list of keys to propagate forward"},
		{ARG_END}
};

#define     Deg2Rad    (M_PI/180.0)
#define     Rad2arcsec (3600.0/Deg2Rad)
#define     arcsec2Rad (Deg2Rad/3600.0)
#define     Rad2Deg    (180.0/M_PI)

char *propagate[] = {"CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
		"MidTime", "Duration", "LonSpan", "T_START", "T_STOP", "Coverage", "Quality",
		"MapScale", "MapProj","Width", "Height"};

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

struct mpi_struct {
	char  hdfname[FLD_MAXPATH+1];
	float   roi_lat;
	float	roi_lon;
	float 	roi_width;
	float 	roi_height;
	float	maplatstp;
	float	maplonstp;
	int	l_max;
	int 	m_abs;
	float	cadence;
	int	nframes;
	float	apod_width;
	float	apod_margin;
	char	apod_func[80];
	float	lowpass_fwhm;
	char	fpname[FLD_MAXPATH+1];
} mpi_send;


typedef struct ObsInfo_struct ObsInfo_t;

//trying to get fldcomplex_t into the hd5
typedef struct {
     float re;   /*real part */
     float im;   /*imaginary part */
} complex_t;

complex_t tmp;  /*used only to compute offsets */

void rebinArraySF(DRMS_Array_t *out, DRMS_Array_t *in);


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
	//DRMS_Array_t *inArray, *outArray;
	CmdParams_t *params = &cmdparams;


	//Input arguments related stuff
	const char *inQuery = params_get_str(&cmdparams, "in");
	const char *outSeries = params_get_str(&cmdparams, "out");
	int verbose = params_isflagset (params, "v");
	int no_save = params_isflagset (params, "n");
	char *seg_name = strdup (params_get_str (params, "segment"));
	const float apod_width = params_get_float(&cmdparams,"apodw");
	const float apod_margin = params_get_float(&cmdparams,"apodm");
	const char *apod_func = params_get_str(&cmdparams,"apodf");
	const char *decomp = params_get_str(&cmdparams,"decomp");
	const float lowpass_fwhm = params_get_float(&cmdparams,"lpfwhm");
	const int l_max = params_get_int(&cmdparams,"l_max");
	const int m_abs = params_get_int(&cmdparams,"m_abs");
	const int interp = params_get_int(&cmdparams,"interp");


	//enum DecompMode decomp;

	char *osegname;
	char logfilename[DRMS_MAXPATHLEN];
	char module_ident[64];
	char *mapproj;
	char *inser;
	char source[DRMS_MAXQUERYLEN];
	char fpname[FLD_MAXPATH+1];
	char hdfname[FLD_MAXPATH+1];

	int keyct = sizeof (propagate) / sizeof (char *);
	int verbose_logs;
	int nframes;
	int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;
	int segs, segct, isegnum,  found;
	int nrecs;
	int numprocs, myrank;

	float duration;
	float cadence;
	float maplatstp;
	float maplonstp;
	float roi_lon,roi_lat,roi_width,roi_height;

	double x0,y0;
	double dval;

	//struct mpi_struct mpi_send;


	TIME t_start;
	TIME t_end;

	int argc;

	//MPI Stuff
	MPI_Status *mpi_stat;
	MPI_Init (&argc,NULL);
	
	int provided;

	int buff;
	
	//MPI_Init_thread(&argc, NULL, MPI_THREAD_SINGLE, &provided);
	//if (provided < MPI_THREAD_MULTIPLE)
	//{
  	// printf("Error: the MPI library doesn't provide the required thread level\n");
   	//MPI_Abort(MPI_COMM_WORLD, 0);
	//}

	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm comm  = MPI_COMM_WORLD;
	//MPI_Info info  = MPI_INFO_NULL;
	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Info_set(info, "romio_cb_write", "enable");
	MPI_Info_set(info, "romio_ds_write", "disable");
	MPI_Info_set(info, "romio_cb_read", "enable");
	MPI_Info_set(info, "romio_ds_read", "disable");
	//MPI_Info_set(info, "cb_buffer_size", "1000000000");


	//creating struct serialization in C and transfer over MPI
	/* create a type for struct car */
    	

	const int nitems=16;
    	int          blocklengths[16] = {FLD_MAXPATH+1,1,1,1,1,1,1,1,1,1,1,1,1,80,1,FLD_MAXPATH+1};
    	MPI_Datatype types[16] = {MPI_CHAR,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_INT,MPI_INT,MPI_FLOAT,MPI_INT,MPI_FLOAT,MPI_FLOAT,MPI_CHAR,MPI_FLOAT,MPI_CHAR};
    	MPI_Datatype mpi_send_type;
    	MPI_Aint     offsets[16];

	//const int nitems=4;
	//int blocklengths[4]={FLD_MAXPATH+1,1,1,FLD_MAXPATH+1};
	//MPI_Datatype types[4]={MPI_CHAR,MPI_FLOAT,MPI_INT,MPI_CHAR};
	//MPI_Datatype mpi_send_type;
	//MPI_Aint     offsets[4];

	offsets[0] = offsetof(struct mpi_struct, hdfname);
        offsets[1] = offsetof(struct mpi_struct, roi_lat);
        offsets[2] = offsetof(struct mpi_struct, roi_lon);
        offsets[3] = offsetof(struct mpi_struct, roi_width);
        offsets[4] = offsetof(struct mpi_struct, roi_height);
        offsets[5] = offsetof(struct mpi_struct, maplatstp);
        offsets[6] = offsetof(struct mpi_struct, maplonstp);
        offsets[7] = offsetof(struct mpi_struct, l_max);
        offsets[8] = offsetof(struct mpi_struct, m_abs);
        offsets[9] = offsetof(struct mpi_struct, cadence);
        offsets[10] = offsetof(struct mpi_struct, nframes);
        offsets[11] = offsetof(struct mpi_struct, apod_width);
        offsets[12] = offsetof(struct mpi_struct, apod_margin);
        offsets[13] = offsetof(struct mpi_struct, apod_func);
        offsets[14] = offsetof(struct mpi_struct, lowpass_fwhm);
        offsets[15] = offsetof(struct mpi_struct, fpname);


	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_send_type);
    	MPI_Type_commit(&mpi_send_type);


	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    	int name_len;
    	MPI_Get_processor_name(processor_name, &name_len);

	//MPI_Barrier(MPI_COMM_WORLD);

	if (myrank ==0)
	{
		/**
		 * Initializing the Input record
		 *
		 */

		/* Get the record from the input query string */
		inRS = drms_open_records (drms_env, inQuery, &status);
		if (status) {
			fprintf (stderr, "Error (%s): drms_open_records() returned %d for dataset:\n",
					module_ident, status);
			fprintf (stderr, "  %s\n", inQuery);
		}

		/* Check how many records are there in the input series */
		nrecs = inRS->n;
		if (nrecs < 1) {
			fprintf (stderr, "No records found in input data set %s\n", inQuery);
		}

		/* Select the record */
		/* and read the number of data segments in the selected record */
		inRec = inRS->records[0];
		inser = strdup (inQuery);
		if ((segs = drms_record_numsegments (inRec)) < 1) {
			fprintf (stderr, "Error: no data segments in input data series:\n");
			fprintf (stderr, "  %s\n", inser);
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
			fprintf (stderr, "  %s\n", inser);
		}
		/* If an array of rank=3 is found, and there are more than one rank 3 array,
		 * please specify the name of the segment to be used */
		if (found > 1) {
			if (strcmp (seg_name, "Not Specified")) {
				inSeg = drms_segment_lookup (inRec, seg_name);
				if (!inSeg) {
					fprintf (stderr,
							"Warning: requested segment %s not found in input data series:\n",
							seg_name);
					fprintf (stderr, "  %s\n", inser);
					inSeg = drms_segment_lookupnum (inRec, isegnum);
					fprintf (stderr, "  using segement %s\n", inSeg->info->name);
				} else if (inSeg->info->naxis != 3) {
					fprintf (stderr,
							"Warning: requested segment %s in input data series:\n", seg_name);
					fprintf (stderr, "  %s is not 3-dimensional", inser);
					inSeg = drms_segment_lookupnum (inRec, isegnum);
					fprintf (stderr, " using segment %s\n", inSeg->info->name);
				} else isegnum = inSeg->info->segnum;
			} else {
				fprintf (stderr,
						"Warning: multiple segments of rank 3 in input data series:\n");
				fprintf (stderr, "  %s\n", inser);
				fprintf (stderr, "  using %s\n", inSeg->info->name);
			}
		}


		/* Get the relevant keyword values from the input record*/
		roi_lon=drms_getkey_float(inRec, "LonHG",&status);
		roi_lat=drms_getkey_float(inRec, "LatHG",&status);
		roi_width=drms_getkey_float(inRec, "Width",&status);
		roi_height=drms_getkey_float(inRec, "Height",&status);
		maplatstp=drms_getkey_double(inRec, "CDELT1",&status);
		maplonstp=drms_getkey_double(inRec, "CDELT2",&status);
		x0=drms_getkey_double(inRec, "CRPIX1",&status);
		y0=drms_getkey_double(inRec, "CRPIX2",&status);
		t_start=drms_getkey_time(inRec, "T_START",&status);
		t_end =drms_getkey_time(inRec, "T_STOP",&status);
		mapproj=drms_getkey_string(inRec,"MapProj",&status);
		/* To get the input number of records */
		duration = drms_getkey_float(inRec, "Duration",&status);
		cadence = drms_getkey_double(inRec, "CDELT3",&status);
		nframes = (int) duration/cadence;

		
		/* Check the number of frames and number of processors assigned */
		if (numprocs >= nframes)
			DIE("-----------\n\nToo many processors for too little frames\n\n\n-----------\n");

		/* Check the Map projection */
		if (!strcmp(mapproj,"Postel"))
            DIE("-----------\n\nUse postel_remap to remap to \"PlateCarree\" Projection\n\n\n-----------\n");

		/**
		 * Initializing the Output record
		 *
		 */
		/* Create the record in the DRMS*/
		if (!(outRS = drms_create_records (drms_env, nrecs, outSeries, DRMS_PERMANENT,
				&status))) {
			fprintf (stderr, "Error: unable to create %d records in series %s\n",
					nrecs, outSeries);
			fprintf (stderr, "       drms_create_records() returned status %d\n", status);
			return 1;
		}
		if (verbose) printf ("creating %d record(s) in series %s\n", nrecs, outSeries);


		if (verbose && dispose == DRMS_FREE_RECORD)
			printf ("experimental run, output records will not be saved\n");

		/*  Check output data series structure  */
		outRec = drms_recordset_getrec (outRS, 0);
		if (!outRec) {
			fprintf (stderr, "Error accessing record %d in series %s\n", 0, outSeries);
			drms_close_records (outRS, DRMS_FREE_RECORD);
			return 1;
		}

		/*Check how many data segments are there */
		segct = drms_record_numsegments (outRec);
		found = 0;
		for (int n = 0; n < segct; n++) {
			outSeg = drms_segment_lookupnum (outRec, n);
			if (!found) osegname = strdup (outSeg->info->name);
			found++;
		}

		if (found < 1) {
			fprintf (stderr,
					"Error: no data segment of dimension 3 and appropriate size in output series %s\n", outSeries);
			drms_close_records (outRS, DRMS_FREE_RECORD);
			return 1;
		}
		outSeg = drms_segment_lookup (outRec, osegname);
		if (found > 1) {
			fprintf (stderr,
					"Warning: multiple data segments of dimension 3 and appropriate size in output series %s\n", outSeries);
			fprintf (stderr, "       using \"%s\"\n", osegname);
		}

		logSeg = drms_segment_lookup (outRec, "Log");
		if (logSeg) drms_segment_filename (logSeg, logfilename);
		else if (verbose) {
			fprintf (stderr,
					"Warning: segment \"Log\" not present in output series %s\n", outSeries);
			fprintf (stderr, "         verbose logging turned off\n");
			verbose_logs = 0;
		}

		if (verbose_logs) {
			logSeg = drms_segment_lookup (outRec, "Log");
			if (logSeg) drms_segment_filename (logSeg, logfilename);
		}

		/* To get the segment name with rank 3 i.e "V" in our case */

		
		for (int n = 0; n < segs; n++) {
			inSeg = drms_segment_lookupnum (inRec, n);
			if (inSeg->info->naxis != 3) continue;
			drms_record_directory(inRec, &fpname, 1);
			//strcpy(&fpname,&fname);
			strcat(&fpname,"/");
			strcat(&fpname,inSeg->info->name);
			strcat(&fpname,".fits");
			//printf("rank=%d, filename=%s\n",myrank,&fpname);
		}

		//printf("rank=%d, filename=%s\n",myrank,&fpname);
		//The HDF file is
		//strcpy(&hdfname,osegname);
		//strcat(&hdfname,".h5");
		//strcpy(mpi_send.hdfname,&hdfname);
		strcpy(&hdfname,osegname);
              	strcat(&hdfname,".h5");
                       
		strcpy(&mpi_send.hdfname,hdfname);
		mpi_send.roi_lat=roi_lat;
                mpi_send.roi_lon=roi_lon;
                mpi_send.roi_width=roi_width;
                mpi_send.roi_height=roi_height;
                mpi_send.maplatstp=maplatstp;
                mpi_send.maplonstp=maplonstp;
                //mpi_send.l_max=l_max;
                //mpi_send.m_abs=m_abs;
                mpi_send.cadence=cadence;
                mpi_send.nframes=nframes;
                //mpi_send.apod_width=apod_width;
                //mpi_send.apod_margin=apod_margin;
                //mpi_send.lowpass_fwhm=lowpass_fwhm;
                //strcpy(mpi_send.apod_func,&apod_func);
               	strcpy(&mpi_send.fpname,fpname);
                  //      printf("Before: rank=%d, filename=%s\n",myrank,&fpname);
                        //MPI_Bcast(&mpi_send, 1, mpi_send_type, 0, MPI_COMM_WORLD);
                        //printf("After: rank=%d, filename=%s\n",myrank,mpi_send.hdfname);


	}
	MPI_Barrier(MPI_COMM_WORLD);
	fflush(stdout);
	//printf("Before: rank=%d, filename=%s,%s\n",myrank,&mpi_send.hdfname,&mpi_send.fpname);	
	//fflush(stdout);
	MPI_Bcast(&mpi_send, 1, mpi_send_type, 0, MPI_COMM_WORLD);
	//fflush(stdout);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("After: rank=%d, filename=%s,%f,%d,%s\n",myrank,&mpi_send.hdfname,mpi_send.roi_lat,mpi_send.l_max,&mpi_send.fpname);

	strcpy(hdfname,mpi_send.hdfname);
        roi_lat=mpi_send.roi_lat;
        roi_lon=mpi_send.roi_lon;
       	roi_width=mpi_send.roi_width;
        roi_height=mpi_send.roi_height;
        maplatstp=mpi_send.maplatstp;
        maplonstp=mpi_send.maplonstp;
        //l_max=&mpi_send.l_max;
        //m_abs=&mpi_send.m_abs;
        cadence=mpi_send.cadence;
        nframes=mpi_send.nframes;
        //apod_width=&mpi_send.apod_width;
        //apod_margin=&mpi_send.apod_margin;
        //lowpass_fwhm=&mpi_send.lowpass_fwhm;
        //strcpy(apod_func,mpi_send.apod_func);
        strcpy(fpname,mpi_send.fpname);
	MPI_Barrier(MPI_COMM_WORLD);
        fflush(stdout);
	//printf("After: rank=%d, filename=%s,%f,%d,%s\n",myrank,&mpi_send.hdfname,roi_lat,l_max,&mpi_send.fpname);
	//printf("After: rank=%d, filename=%s,%f,%d,%s\n",myrank,&hdfname,roi_lat,l_max,&fpname);

	//else {
		//MPI_Bcast(&hdfname, FLD_MAXPATH+1, MPI_CHAR, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&roi_lat, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&roi_lon, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&roi_width, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&roi_height, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&maplatstp, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&maplonstp, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&l_max, 1, MPI_INT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&m_abs, 1, MPI_INT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&cadence, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&nframes, 1, MPI_INT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&apod_width, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&apod_margin, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&apod_func, 80, MPI_CHAR, 0, MPI_COMM_WORLD);
        	////MPI_Bcast(&DecompMode, 80, MPI_INT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&lowpass_fwhm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&fpname, FLD_MAXPATH+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		//printf("rank=%d, filename=%d\n",myrank,mpi_send.l_max);
		//printf("Message from processor %s, rank %d"
	        //   " out of %d processors\n",
        	//   processor_name, myrank, numprocs);

	//}

	/*TODO: Clean this up: The MPI Broadcasting to all the processors */
	fflush(stdout);

	/* FLD part starts here */
	fldjob job = fldjob_alloc ();
	
        /* Some MPI checks */
	fflush(stdout);
 	MPI_Barrier(MPI_COMM_WORLD);	

	bool errorflag = false;

	//DIE("-----------\n\nJust testing....\n\n\n");
	if (job)
	{
		//job parameters

		job->roi_lat = DEG_TO_RAD (roi_lat);
		job->roi_lon = DEG_TO_RAD (roi_lon);
		job->roi_width = DEG_TO_RAD (roi_width);
		job->roi_height = DEG_TO_RAD (roi_height);
		job->mapscale_theta = DEG_TO_RAD (maplatstp);
		job->mapscale_phi = DEG_TO_RAD (maplonstp);
		job->l_min = 0;
		job->l_max = l_max;
		job->nl= l_max+1;
		job->m_abs=m_abs;
		job->nm=1 + (2 * job->m_abs);
		job->nchunks = numprocs-1;							//TODO: numprocs or numprocs-1 Remove this hardcoded stuff
		job->maps_per_chunk = nframes/job->nchunks;			//TODO: Remove this hardcoded stuff
		if (nframes%(job->nchunks) > 0)
			job->maps_per_chunk=job->maps_per_chunk+1;
		if (job->maps_per_chunk/2 > 1)
			job->min_maps_per_chunk=job->maps_per_chunk/2;
		else
			job->min_maps_per_chunk=1;
		job->apod_width=apod_width;
		job->apod_margin=apod_margin;
		//strncpy (job->apod_func, apod_func,80);
		strncpy (job->apod_func, "Tukey",80);
		//strncpy (job->apod_func, &apod_func,80);				//TODO: VG: Change this later
		//job->decomp_mode=DecompMode;					//TODO: VG: Change this later
		job->lowpass_fwhm=lowpass_fwhm;
		job->cadence=cadence;								//TODO: Remove this hardcoded stuff
		job->do_save_maps=false;							//TODO: Is used to save fits file if requested
		job->chunk_dmap_mask=NULL;
		job->data_scale_factor=1;							//TODO: Remove this hardcoded stuff
		job->roi_npix_phi = (int) (job->roi_width / job->mapscale_phi);
		job->roi_npix_theta = (int) (job->roi_height / job->mapscale_theta);
		//job->roi_width = job->roi_npix_phi * job->mapscale_phi;
		//job->roi_height = job->roi_npix_theta * job->mapscale_theta;
		job->do_demean = FALSE;								//TODO: Remove this hardcoded stuff
		MPI_Barrier(MPI_COMM_WORLD);
		strncpy (job->input_file, &fpname, DRMS_MAXPATHLEN);

		//printf("filename=====%s\n\n\n\n",job->input_file);

	
		/* VG: This is used only when do_save_maps=true, this will output lmg and maps in the hardcoded folder */
	    if (job->do_save_maps)
	    {
		getcwd (job->jobdir, FLD_MAXPATH);
		strcat(job->jobdir,"/../fld_testing");

		if (myrank == 0)
		{
			mode_t process_mask1 = umask(0);
			mkdir(job->jobdir, S_IRWXU | S_IRWXG | S_IRWXO);
			umask(process_mask1);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		strcat(job->jobdir,"/");
	    }


		/* HDF5 Processing */
		/* Create the file if it does not exist */
		if( access( hdfname, F_OK ) == -1 )
		{
			/* HDF5 Variable definitions */
			hid_t       hfile_id, hdataset_id, hdataspace_id;
			hid_t		grp;
			hid_t 		complex_id;
			/* Setting profile for Parallel I/O */
			hid_t		plist_id;
			herr_t      hstatus;

			plist_id = H5Pcreate(H5P_FILE_ACCESS);

			if (H5Pset_fapl_mpio(plist_id, comm, info))
				printf('Error reported in setting file access property list with Parellel I/O access');

			long naxes = 3;
			long fullframe= job->maps_per_chunk * job->nchunks; //number of records
			long mynaxis[3] = {job->nm, job->nl, fullframe};

			/* Create the file */
			hfile_id = H5Fcreate(HFILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
			/* close property list */
			hstatus = H5Pclose(plist_id);

			/* Create the group */
			grp = H5Gcreate2(hfile_id,"/lmg", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			/* Create a Dataspace */
			hdataspace_id = H5Screate_simple(naxes,mynaxis,NULL);
			/* Create the Complex Datatype */
			complex_id = H5Tcreate (H5T_COMPOUND, sizeof(tmp));
			if (H5Tinsert (complex_id, "real", HOFFSET(complex_t, re), H5T_NATIVE_FLOAT))
				fprintf(stderr, "Unable to add a real part to the complex datatype.\n");
			if (H5Tinsert (complex_id, "imaginary", HOFFSET(complex_t, im), H5T_NATIVE_FLOAT))
				fprintf(stderr, "Unable to add an imaginary part to the complex datatype.\n");
			/* Create the Dataset for Alm's */
			hdataset_id = H5Dcreate2(hfile_id, "/lmg/Alm", complex_id, hdataspace_id,
					H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);
			/* End access to the dataset and release resources used by it. */
			if (H5Dclose(hdataset_id))
				fprintf(stderr, "Unable to terminate access to dataset.\n");
			/* Create the Dataset for Blm's */
			hdataset_id = H5Dcreate2(hfile_id, "/lmg/Blm", complex_id, hdataspace_id,
					H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);
			/* End access to the dataset and release resources used by it. */
			if (H5Dclose(hdataset_id))
				fprintf(stderr, "Unable to terminate access to dataset.\n");
			/* Close the datatype */
			if (H5Tclose(complex_id))
				fprintf(stderr, "Unable to close complex datatype.\n");
		    /* Terminate access to the data space. */
			if (H5Sclose(hdataspace_id))
				fprintf(stderr, "Unable to terminate access to input dataspace.\n");
			/* Close the group */
			if (H5Gclose(grp))
				fprintf(stderr, "Unable to close input group.\n");
			/* Close the file. */
			if (H5Fclose(hfile_id))
				fprintf(stderr, "Unable to close input hdf5 file.\n");
		}
		

		MPI_Barrier(MPI_COMM_WORLD);

		if (FLD_SUCCESS == fld_decomp_init (job))
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if (FLD_SUCCESS != fld_decomp_schedule_job (job)) {
				logmesg (LERROR, "couldn't process job, aborting.\n");
				errorflag = true;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			fld_decomp_finalize (job);
		}
		else {
			logmesg (LERROR, "couldn't initialize decomposer, aborting.\n");
			errorflag = true;
		}

		
		//printf("Here: rank=%d, filename=%s,%f,%d,%s\n",myrank,&hdfname,roi_lat,l_max,&fpname);
                //printf("Here.................:%s\n",myrank);
                //DIE("-----------\n\nJust testing....\n\n\n");



		MPI_Barrier(MPI_COMM_WORLD);

		fldjob_free (job);
	}
	else {
		logmesg (LERROR, "couldn't initialize job handle, aborting.\n");
		errorflag = true;
	}

	if (errorflag) {

		MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
		return EXIT_FAILURE;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* Decomposition and saving into HDF5 file ends here */

	/* All the attributes to the HDF5 file are assigned here */
	if (myrank ==0)
	{

	    hid_t dataspace_id,file_id,grp, attribute_id;

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
	    //char	attr_char[1];
	    double result;

	    for (int n = 0; n < keyct; n++)
	    {
	    	if (drms_keyword_type(drms_keyword_lookup(inRec,propagate[n],1))==DRMS_TYPE_FLOAT)
	    	{
	    		dval = drms_getkey_float (inRec, propagate[n], &status);
	    		attr_flt[0]=dval;
	    		attribute_id = H5Acreate(grp, propagate[n], H5T_IEEE_F32BE, dataspace_id,
	    				H5P_DEFAULT, H5P_DEFAULT);
	    		status=H5Awrite(attribute_id,H5T_NATIVE_FLOAT,attr_flt);
	    		status=H5Aclose(attribute_id);
	    	}
	    	if (drms_keyword_type(drms_keyword_lookup(inRec,propagate[n],1))==DRMS_TYPE_DOUBLE)
	    	{
	    		dval = drms_getkey_double (inRec, propagate[n], &status);
	    		attr_dbl[0]=dval;
	    		attribute_id = H5Acreate(grp, propagate[n], H5T_IEEE_F64LE, dataspace_id,
	    				H5P_DEFAULT, H5P_DEFAULT);
	    		status=H5Awrite(attribute_id,H5T_NATIVE_DOUBLE,attr_dbl);
	    		status=H5Aclose(attribute_id);
	    	}
	    	if (drms_keyword_type(drms_keyword_lookup(inRec,propagate[n],1))==DRMS_TYPE_INT)
	    	{
	    		dval = drms_getkey_int (inRec, propagate[n], &status);
	    		attr_int[0]=dval;
	    		attribute_id = H5Acreate(grp, propagate[n], H5T_STD_I32BE, dataspace_id,
	    		   				H5P_DEFAULT, H5P_DEFAULT);
	    		status=H5Awrite(attribute_id,H5T_NATIVE_INT,attr_int);
	    		status=H5Aclose(attribute_id);
	    	}
	    	if (drms_keyword_type(drms_keyword_lookup(inRec,propagate[n],1))==DRMS_TYPE_TIME)
	    	{
	    		dval = drms_getkey_time (inRec, propagate[n], &status);
	    		result = drms2time(DRMS_TYPE_TIME, &dval, &status);
	    		attr_dbl[0]=result;
	    		attribute_id = H5Acreate(grp, propagate[n], H5T_IEEE_F64LE, dataspace_id,
	    				H5P_DEFAULT, H5P_DEFAULT);
	    		status=H5Awrite(attribute_id,H5T_NATIVE_DOUBLE,attr_dbl);
	    		status=H5Aclose(attribute_id);
	    	}
	    	/*
	    	if (drms_keyword_type(drms_keyword_lookup(outRec,propagate[n],1))==DRMS_TYPE_STRING)
	    	{
	    		dval = drms_getkey_string (inRec, propagate[n], &status);
	    		attr_char[0]=dval;
	    		atype=H5Tcopy(H5T_C_S1);
	            H5Tset_size(atype, 5);
	            H5Tset_strpad(atype,H5T_STR_NULLTERM);
	    		attribute_id = H5Acreate(grp, propagate[n], atype, dataspace_id,
	    				H5P_DEFAULT, H5P_DEFAULT);
	    		status=H5Awrite(attribute_id,H5T_NATIVE_CHAR,attr_char);
	    		status=H5Aclose(attribute_id);
	    		printf("Testing, no string attribute");
	    	}
	    	*/

	    }


	    status=H5Sclose(dataspace_id);
	    status=H5Gclose(grp);
	    status=H5Fclose(file_id);

	    propagate_keys (outRec, inRec, propagate, keyct);
	    drms_sprint_rec_query (source, inRec);
		inser = strdup (inQuery);
	    check_and_set_key_double (outRec, "lMax", l_max);
	    check_and_set_key_double (outRec, "mAbs", m_abs);
	    check_and_set_key_str   (outRec, "Source", source);
	    check_and_set_key_str   (outRec, "Module", module_ident);
	    check_and_set_key_str   (outRec, "BLD_VERS", jsoc_version);
	    check_and_set_key_str   (outRec, "Input", inser);
	    check_and_set_key_time  (outRec, "Created", CURRENT_SYSTEM_TIME);


		/* Writing the record to the DRMS */
		if(drms_segment_write_from_file(outSeg,hdfname))
			fprintf(stderr, "Cannot write hdf5 file for the segment.\n");

		/* Close the input record */
		if(drms_close_records(inRS, DRMS_FREE_RECORD))
			fprintf(stderr, "Error closing Record for input series\n");

		/* Close the output record */
		if(drms_close_records(outRS, DRMS_INSERT_RECORD))
			fprintf(stderr, "Error closing Record for output series\n");

	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize ();
	remove(HFILE);

	printf("---END---\n");
	return (DRMS_SUCCESS);
	//return EXIT_SUCCESS;
} // end of DoIt



#define CHECK(keyname) {if (status) {fprintf(stderr,"Keyword failure to find: %s, status=%d\n",keyname,status); *rstatus=status; return(NULL);}}
