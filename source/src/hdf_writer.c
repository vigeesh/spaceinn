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

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/utsname.h>

#include <mpi.h>
#include "hdf5.h"
#include "config.h"
#include "fldjob.h"
#include "fldmath.h"
#include "datacube.h"
#include "fldfits.h"
#include "backend.h"


//trying to get fldcomplex_t into the hd5
typedef struct {
     float re;   /*real part */
     float im;   /*imaginary part */
} complex_t;
complex_t tmp;  /*used only to compute offsets */



int fldhdf5_write_cCube (
    fldjob j,
    cCube c,
    char *file,
    fitscard *header,
    int ncards,
    char *comment,
    int chunkid,
    char *aorb
    )
{

	int numprocs, myrank;

	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	//MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info;
	MPI_Info_create(&info);

    //int status = 0;

	#define FILE "lmg.h5"
    size_t nElemTot = c->nFrames * c->nPixX * c->nPixY;

	//size_t nelem_lmgrid = c->nPixX * c->nPixY;
	//cCube needle_lm = cCube_alloc (nelem_lmgrid, c->nFrames, 1);
	//cCube needle_lm = cCube_alloc (c->nFrames, c->nPixX, c->nPixY);

	//cCube needle_Blm = cCube_alloc (nelem_lmgrid, c->nFrames, 1);;

    hsize_t naxes = 3;
    hsize_t naxis[3] = {c->nPixX, c->nPixY, c->nFrames};

    double starttime, endtime;

    fflush(stdout);

	complex_t *mydata;

	mydata = (complex_t*) malloc (nElemTot * sizeof (complex_t*));

//	printf("Size of mydata is %zu\n",sizeof(mydata));

	int thisiscount=0;
	for (int midx = 0; midx < c->nPixX; midx++)
	{
		for (int lidx = 0; lidx < c->nPixY; lidx++)
		{
			for (int mapidx = 0; mapidx < c->nFrames; mapidx++)
				{
					mydata[thisiscount].re = crealf(cCube_getPixel (c, mapidx, midx, lidx));
					mydata[thisiscount].im = cimagf(cCube_getPixel (c, mapidx, midx, lidx));
					thisiscount+=1;
			}
		}
	}

	hid_t       hfile_id, hdataset_id, hdataspace_id;
	hid_t		grp, hmemspace_id;
	herr_t      hstatus;
	//hid_t 		lmg_id,plist;  /* property list */

	hsize_t     count[3];              /* size of subset in the file */
	hsize_t     offset[3];             /* subset offset in the file */

	hid_t 		complex_id;

    /*
     *
     printf("\n\nbefore_list_creation\n\n");
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	printf("\n\nafter_list_creation\n\n");

	printf("\n\nbefore dxpl in\n\n");
    if (H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT))
    	printf('Error reported in MPIO collective transfer property assignement');
	printf("\n\nafter_dxpl\n\n");
	*/

	fflush(stdout);

	//printf("what will happen %d ", myrank);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("here %d :\n", myrank);


	if( access( FILE, F_OK ) != -1 ) {
		//printf("file exists\n");
		/* Specify size and shape of subset to write. */

		offset[0] = 0;
		offset[1] = 0;
		offset[2] = naxis[2]*chunkid;

		count[0]  = naxis[0];
		count[1]  = naxis[1];
		count[2]  = naxis[2];

//		plist_id = H5Pcreate(H5P_FILE_ACCESS);
//		printf("rank=%d\n",myrank);
//		printf("Alm/Blm....in %d\n",myrank);
//		printf("1....\n");
		//H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
//		if (H5Pset_fapl_mpio(plist_id, comm, info))
//			printf('Error reported in setting file access property list with Parellel I/O access');

//		printf("2....\n");
		//printf('Error reported in setting the profiles for Parellel I/O');
	    //open the file

	    hfile_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);

//		printf("3....\n");
//		hstatus = H5Pclose(plist_id);
//		printf("4.....\n");

	    //open the group in the file
	    grp  = H5Gopen2(hfile_id, "/lmg", H5P_DEFAULT);

		if (strcmp(aorb,"Alm"))
		{

			//printf("Alm: %d is writing file %s\n", chunkid, file);
			//printf("Alm: %d Writing on Proc: %i\n",endtime-starttime, myrank);
    	    hdataset_id = H5Dopen2(hfile_id,"/lmg/Alm", H5P_DEFAULT);

    	    /* Create memory space with size of subset. Get file dataspace
    	           and select subset from file dataspace. */
    	    hmemspace_id = H5Screate_simple (naxes, naxis, NULL);

    	    hdataspace_id = H5Dget_space (hdataset_id);
    	    hstatus = H5Sselect_hyperslab (hdataspace_id, H5S_SELECT_SET, offset,
    	                                      NULL, count, NULL);

    	    /* Write a subset of data to the dataset */
			complex_id = H5Tcreate (H5T_COMPOUND, sizeof(tmp));
			hstatus = H5Tinsert (complex_id, "real", HOFFSET(complex_t, re), H5T_NATIVE_FLOAT);
			hstatus = H5Tinsert (complex_id, "imaginary", HOFFSET(complex_t, im), H5T_NATIVE_FLOAT);
    	    //hstatus = H5Dwrite(hdataset_id, complex_id, hmemspace_id, hdataspace_id, H5P_DEFAULT, mydata);

    	}

		if (strcmp(aorb,"Blm"))
		{
			//printf("Blm: %d Writing on Proc: %i\n",endtime-starttime, myrank);
			hdataset_id = H5Dopen2(hfile_id,"/lmg/Blm", H5P_DEFAULT);

			/* Create memory space with size of subset. Get file dataspace
		    	           and select subset from file dataspace. */

			hmemspace_id = H5Screate_simple (naxes, naxis, NULL);

			hdataspace_id = H5Dget_space (hdataset_id);
			hstatus = H5Sselect_hyperslab (hdataspace_id, H5S_SELECT_SET, offset,
					NULL, count, NULL);

			/* Write a subset of data to the dataset, then read the
		    	       entire dataset back from the file.  */
			complex_id = H5Tcreate (H5T_COMPOUND, sizeof(tmp));
			hstatus = H5Tinsert (complex_id, "real", HOFFSET(complex_t, re), H5T_NATIVE_FLOAT);
			hstatus = H5Tinsert (complex_id, "imaginary", HOFFSET(complex_t, im), H5T_NATIVE_FLOAT);
		}

		starttime = MPI_Wtime();
		hid_t	plist_aid;
	    plist_aid = H5Pcreate(H5P_DATASET_XFER);

	    if (H5Pset_dxpl_mpio(plist_aid, H5FD_MPIO_INDEPENDENT))
	    	fprintf(stderr, "Unable to set property list.\n");

	    //MPI_Barrier(MPI_COMM_WORLD);
		if (H5Dwrite(hdataset_id, complex_id, hmemspace_id, hdataspace_id, plist_aid, mydata))
	    {
	    	fprintf(stderr, "Unable to write data.\n");
	    }
		endtime = MPI_Wtime();
		//MPI_Barrier(MPI_COMM_WORLD);
		free(mydata);

		/* close property list */
    	hstatus = H5Pclose(plist_aid);

		printf("It took %f s to write on Proc. %d\n",endtime-starttime, myrank);
		/* Close the datatype */
		hstatus = H5Tclose(complex_id);

    	/* End access to the dataset and release resources used by it. */
    	hstatus = H5Dclose(hdataset_id);

    	/* Terminate access to the data space. */
    	hstatus = H5Sclose(hmemspace_id);

    	/* Close the group */
    	hstatus = H5Gclose(grp);

    	/* Close the file. */
    	hstatus = H5Fclose(hfile_id);

    }

    return FLD_SUCCESS;
}

