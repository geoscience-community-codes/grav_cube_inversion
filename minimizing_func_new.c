/* 
	 File Name:   minimizing function.c

	 Program Name:  grav_parallel        
	 Subroutine Name(s): test_bounds(), init_globals(), get_points(),
                       setup_prisms(), get_prisms(),
                       minimizing_func(),
                       assign_new_params(), init_optimal_params(), 
                       printout_points(), printout_parameters(),
                       printout_model(), _free(), rmse()
                       
	 Release Date:         April 1, 2020
	 Release Version:      1.0

	 VERSION/REVISION HISTORY
	 
	 Date: April 1, 2020, Author: Laura Connor        
	 Initial Release, Version 1.0
	 
	 
	 DISCLAIMER/NOTICE
	 
	 This computer code/material was prepared as an account of work
	 performed by the Center for Nuclear Waste Regulatory Analyses (CNWRA)
	 for the Division of Waste Management of the Nuclear Regulatory
	 Commission (NRC), an independent agency of the United States
	 Government. The developer(s) of the code nor any of their sponsors
	 make any warranty, expressed or implied, or assume any legal
	 liability or responsibility for the accuracy, completeness, or
	 usefulness of any information, apparatus, product or process
	 disclosed, or represent that its use would not infringe on
	 privately-owned rights.
	 
	 IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW WILL THE SPONSORS
	 OR THOSE WHO HAVE WRITTEN OR MODIFIED THIS CODE, BE LIABLE FOR
	 DAMAGES, INCLUDING ANY LOST PROFITS, LOST MONIES, OR OTHER SPECIAL,
	 INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR
	 INABILITY TO USE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
	 BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY THIRD PARTIES OR A
	 FAILURE OF THE PROGRAM TO OPERATE WITH OTHER PROGRAMS) THE PROGRAM,
	 EVEN IF YOU HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES,
	 OR FOR ANY CLAIM BY ANY OTHER PARTY.
	 
	 
	 PURPOSE: 
	 The purpose of each subroutine in this file is described above each
	 individual routine.
	 
	 PROGRAMMING LANGUAGE:  ANSI C 
	 
	 GLOBAL VARIABLES:

	   DEBUG : if zero, no debugging output is produced
                   if > zero, various debugging outputs are written to stderr and the log files
	   NUM_OF_PARAMS : an integer,  the number of prism parameters that will be simultaneously modeled
	   NUM_OF_VERTICES : an integer, the number of vertices of the simplex model
                             (i.e. the number of sets of parameters being simultaneously modeled)
		             this value is always one greater static int num_pts = 0;that the NUM_OF_PARAMS
	   TOLERANCE :  the program runs until the goodness-of-fit values, resulting from a comparison of the calculated
	                with the observed magnetic values, all fall within the range of this value
	   _LO[LAST_PARAM] :  an array of the minimum parameter values
	   _HI[LAST_PARAM] :  an array of the maximum parameter values

	 REFERENCES: Numerical Recipes
	 
	 PROGRAM FLOW:
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <gc.h>
#include "prototypes.h"

/* The maximum line length */
#define MAX_LINE 200
#define README "parameters.README"

static unsigned int SEED = 0;
static POINT *p_all=NULL;
static int total_pts = 0;

/* local node varialbles */
static int procs=-1;
static int my_rank=-1;
static int my_count=-1; /* byte count of node's POINT array */
static int *displ=NULL;
static int *recv_ct=NULL; /* pointers to arrays of integers based on total number of nodes */

static int num_pts = 0;
static POINT *pt=NULL;
static PRISM *pr=NULL;
static PARAMETER P;
static FILE *log_file=NULL;
static double **GRID=NULL;

/****************************************************************
FUNCTION: test_bounds
DESCRIPTION: This function bounds a value if it is 
greater than or less than its preset boundaries.
INPUTS: (IN) int param  (the index of the parameter being tested)
        (IN/OUT) double *try  (the actual value being tested)
        (IN) double bound  (a variable boundary value)
        in this case the bound is the depth_to_top
        The depth_to_bottom should always be greater than the
        depth_to_top (depths are positive down)
OUTPUT: none 
*****************************************************************/
void test_bounds(int param, double *try, double bound) {

  if (param == DEPTH_TO_BOT) {
    //if (*try < bound) *try = bound;
    if ((*try - bound) < 50) *try = bound;
   }
  
  if (*try < LO_PARAM(param)) *try = LO_PARAM(param); 
    
  else if (*try > HI_PARAM(param)) *try = HI_PARAM(param); 
  
}

/****************************************************************
FUNCTION: init_globals
DESCRIPTION: This function reads a configuration file
and sets some global variables.
INPUTS:  (IN) char *  (complete path to the configuration file)
OUTPUTS: int 1=error, 0=no error
 ****************************************************************/
int init_globals(char *config_file, INPUTS *in) {
	
  FILE *conf_file;
  char buf[1][30], **ptr1;
  char line[MAX_LINE];
  char space[4] = "\n\t ";
  char *token;
  int i;
  
  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  
  /* Get my process rank, an integer */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  conf_file = fopen(config_file, "r");
  if (conf_file == NULL) {
    fprintf(stderr, 
	    "[%d-of-%d]\tCannot open configuration file=[%s]:[%s]. Exiting.\n", 
	    my_rank, procs, config_file, strerror(errno)); 
    return 1;
  }
  fprintf(stderr, "[Node %d]Reading config file: %s\n",my_rank, config_file );
  
  ptr1 = (char **)&buf[0];
  while (fgets(line, MAX_LINE, conf_file) != NULL) { 

    if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') continue;
     
    token = strtok_r(line, space, ptr1);
    // fprintf(stderr, "%s\t", token);
    
    if (!strncmp(token, "TOLERANCE", strlen("TOLERANCE"))) {
      token = strtok_r(NULL, space, ptr1);
      TOLERANCE = strtod(token, NULL);
      fprintf(log_file, "TOLERANCE = %f\n", TOLERANCE);
      if (!TOLERANCE) TOLERANCE = 0.00001;
    }
    else if (!strncmp(token, "MIN_ROC_DENSITY", strlen("MIN_ROC_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[DENSITY] = strtod(token, NULL);
      P.density = _LO[DENSITY];
      fprintf(log_file, "MIN ROC DENSITY=%.2f\n", _LO[DENSITY]); 
    }
    else if (!strncmp(token, "MAX_ROC_DENSITY", strlen("MAX_ROC_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      _HI[DENSITY] = strtod(token, NULL);
      fprintf(log_file, "MAX ROC DENSITY=%.2f\n", _HI[DENSITY]);
    }
    else if (!strncmp(token, "MIN_NORTHING", strlen("MIN_NORTHING"))) {
      token = strtok_r(NULL,space,ptr1);
      P.min_northing = strtod(token, NULL);
      fprintf(log_file, "MIN Northing = %f\n",  P.min_northing);
    }
    else if (!strncmp(token, "MAX_NORTHING", strlen("MAX_NORTHING"))) {
      token = strtok_r(NULL,space,ptr1);
      P.max_northing = strtod(token, NULL);
      fprintf(log_file, "MAX Northing = %f\n", P.max_northing);
    }
    else if (!strncmp(token, "MIN_EASTING", strlen("MIN_EASTING"))) {
      token = strtok_r(NULL,space,ptr1);
      P.min_easting = strtod(token, NULL);
      fprintf(log_file, "MIN Easting = %f\n",  P.min_easting);
   }
   else if (!strncmp(token, "MAX_EASTING", strlen("MAX_EASTING"))) {
      token = strtok_r(NULL,space,ptr1);
      P.max_easting = strtod(token, NULL);
      fprintf(log_file, "MAX Easting = %f\n", P.max_easting);
    }
    else if (!strncmp(token, "MIN_DEPTH_TO_TOP", strlen("MIN_DEPTH_TO_TOP"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[DEPTH_TO_TOP] = strtod(token, NULL);
      P.depth_to_top = _LO[DEPTH_TO_TOP]; 
      fprintf(log_file, "MIN DEPTH TO TOP = %f\n", _LO[DEPTH_TO_TOP]);
    }  
    else if (!strncmp(token, "MAX_DEPTH_TO_TOP", strlen("MAX_DEPTH_TO_TOP"))) {
      token = strtok_r(NULL,space,ptr1);
      _HI[DEPTH_TO_TOP] = strtod(token, NULL);
      //P.depth_to_top = _HI[DEPTH_TO_TOP]; 
      fprintf(log_file, "MAX DEPTH TO TOP = %f\n", _HI[DEPTH_TO_TOP]);
    }
    else if (!strncmp(token, "MIN_DEPTH_TO_BOTTOM", strlen("MIN_DEPTH_TO_BOTTOM"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[DEPTH_TO_BOT] = strtod(token, NULL);
      fprintf(log_file, "MIN DEPTH TO BOTTOM = %f\n", _LO[DEPTH_TO_BOT]);
   }
   else if (!strncmp(token, "MAX_DEPTH_TO_BOTTOM", strlen("MAX_DEPTH_TO_BOTTOM"))) {
      token = strtok_r(NULL,space,ptr1);
      _HI[DEPTH_TO_BOT] = strtod(token, NULL); 
      fprintf(log_file, "MAX DEPTH TO BOTTOM = %f\n", _HI[DEPTH_TO_BOT]);   
    }    
    else if (!strncmp(token, "SPACING", strlen("SPACING"))) {
      token = strtok_r(NULL,space,ptr1);
      P.sp = strtod(token, NULL);
      fprintf(log_file, "Spacing = %f\n", P.sp);
    }
    else if (!strncmp(token, "SEED", strlen("SEED"))) {
      token = strtok_r(NULL, space, ptr1);
      SEED = (unsigned int)atoi(token);
      fprintf(log_file, "SEED = %u\n", SEED);
    }
    else if (!strncmp(token, "OBS_GRAV_FILE", strlen("OBS_GRAV_FILE"))) {
    	token = strtok_r(NULL, space, ptr1);
    	in->points_file = (char*) GC_MALLOC(sizeof(char) * (strlen(token)+1));
			if (in->points_file == NULL) 
			{
				fprintf(stderr, 
				        "\n[INITIALIZE] Out of Memory assigning filename!\n");
				return 1;
			}
			memset(in->points_file, '\0', strlen(in->points_file));
			strcpy(in->points_file, token);
			fprintf(log_file, "OBS_GRAV_FILE = %s\n", in->points_file);
    }
    else continue;
  }
  if (in->points_file == NULL) {
  	fprintf(stderr, 
				        "\n[INITIALIZE] No gravity observation file specified!\n");
				return 1;
	}
  fprintf(log_file, "Top Surface from %.2f to %.2f\n", _LO[DEPTH_TO_TOP], _HI[DEPTH_TO_TOP]);
  fprintf(log_file, "Bottom Surface from %.2f to %.2f\n", _LO[DEPTH_TO_BOT], _HI[DEPTH_TO_BOT]);
  
 fprintf(stderr, "[%d]Read complete\n", my_rank); 
 
  NUM_OF_PARAMS = setup_prisms() + 2;
  
  fprintf(log_file, "NUM_OF_PARAMS=%d\n", NUM_OF_PARAMS);
  NUM_OF_VERTICES = NUM_OF_PARAMS + 1;

  GRID = (double **)GC_MALLOC((size_t)P.row * sizeof(double));
  if (GRID == NULL) {
    fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for GRID rows:[%s]\n",
	    my_rank, procs, strerror(errno));
	 (void) fclose(conf_file);
    return -1;
  } 
  
  else {
    for (i=0; i < P.row; i++) {
      GRID[i] = (double *)GC_MALLOC((size_t)P.col * sizeof(double));
      if (GRID[i] == NULL) {
	     fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for grid row %d:[%s]\n",
		    my_rank, procs, i, strerror(errno));
		  (void) fclose(conf_file);
	     return -1;
      }
    }
  }
 
  (void) fclose(conf_file);
  return 0;
} 



/*****************************************************************
FUNCTION:  get_points
DESCRIPTION:  This function reads northing,easting coordinates 
and associated magnetic value (nT) from a file of observed data
into a POINTS array.
The total number of points read are divided up between 
nodes so that each node can calculate the magnetic field value at
its portion of the points read.
INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
 ****************************************************************/
int get_points(FILE *in) {
  
  char line[MAX_LINE]; /*maximum line read */ 
  int i, ret;
  
  int extra = 0; /* remaining points to calculate if total does not divide evenly amount nodes */
  int my_start; /* starting line in points file (local) */
  int pts_read = 0; /* number of points read so far (local) */

 /* if (DEBUG == 2) fprintf(log_file, "ENTER[get_points]\n");*/
  
  while (fgets(line, MAX_LINE, in) != NULL)  {
  	if (line[0] == '#' || line[0] == '\n') continue;
   total_pts++;
  }
  rewind(in);
  fprintf(log_file, "  Total Number of points=%d\n", total_pts);
  
  /* Calculate number of points to calculate and starting line in file */
  /* if total points does not divide equally among nodes let highest numbered node do the remainder */
  extra = total_pts % procs;
  if (extra) {
  	num_pts = (int)(total_pts/procs);
  	my_start = my_rank * num_pts;
	if (my_rank == (procs -1)) num_pts += extra;
  }
  else {
  	num_pts = (int)(total_pts/procs);
  	my_start = my_rank * num_pts;
  }
  
  /* Allocate global storage for POINT structures being calculated. 
   * Only needs to be done on root node. 
   */
  if ( !my_rank ) { /* code for master node */
    p_all = (POINT *)GC_MALLOC((size_t)total_pts * sizeof(POINT));
    if (p_all == NULL) {
      fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for all points:[%s]\n",
              my_rank, procs, strerror(errno));
      fclose(in);
      return -1;
    }
    fprintf(log_file, "  TOTAL BYTE COUNT=%ld\n", 
	    (total_pts * sizeof(POINT)) );
    
    /* The size of these arrays of integers are based on the total number of nodes used. */
    displ = (int *)GC_MALLOC((size_t)procs * sizeof(int));
    recv_ct = (int *)GC_MALLOC((size_t)procs * sizeof(int));
    
    /* For each node, calculate the total number of bytes of point storage needed */
    displ[0] = 0;
    // fprintf(stderr, "Total bytes of data = %ld\n", total_pts * sizeof(POINT));
    for (i=0; i < procs; i++) {
      recv_ct[i] = num_pts * sizeof(POINT);
      if (extra) if (i == (procs-1)) recv_ct[i] = (num_pts + extra) * sizeof(POINT);
      
      /* Calculate each node's total byte displacement into the POINT array of data */
      if (i>0) displ[i] = displ[i-1] + recv_ct[i-1];
    /*  if (DEBUG==5) fprintf(stderr,"RECV_CT[%d]=%d bytes  DISPL[%d]= %d\n", i, recv_ct[i], i, displ[i]); */
    }
  } /* end code for master node */
  
  /* Allocate memory for each node's POINT structures (in bytes). */
  my_count = num_pts * sizeof(POINT);
  fprintf(log_file,"  MY BYTE COUNT=%d\n", my_count);
  
  pt = (POINT *) GC_MALLOC((size_t)my_count);
  if (pt == NULL) {
    fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for points:[%s]\n",
            my_rank, procs, strerror(errno));
    fclose(in);
    return -1;
  }
  fprintf(log_file,"  Number of points to calculate=%d\n\tStarting point=%d\n",
	  num_pts, my_start);
  
  /* Each node reads from the points file  and stores its fraction of points to calculate */
  /* for ( i = 0; i < total_pts; i++) { */
  i=0;
  while (i < total_pts) {
    fgets(line, MAX_LINE, in);
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      while (ret = sscanf(line, "%lf %lf %lf", 
			 &(pt+pts_read)->easting, 
			 &(pt+pts_read)->northing,
			 &(pt+pts_read)->observed), 
		  ret != 3) {
		  	
        if (ret == EOF && errno == EINTR) continue; 
        fprintf(stderr, "[%d-of-%d]\t[line=%d,ret=%d] Did not read in 3 points:[%s]\n", 
              my_rank, procs, i+1,ret, strerror(errno));
        fclose(in);
        return -1;
      }
      if ( i >= my_start ) {
        pts_read++;
        if (pts_read == num_pts) break;
      }
    }
    i++;
 }
  fprintf(log_file,"EXIT[get_points]:[%d-of-%d]Read %d points.\n", 
	  my_rank, procs, pts_read);
  fflush(log_file);
  fclose(in);
  return 0;
}


/**************************************************************
FUNCTION:  setup_prisms
DESCRIPTION:  
INPUTS: 
OUTPUTS: int 
***************************************************************/
int setup_prisms(void) {
  int count, x, y, i;
  double xmin, ymax; /*ymin*/

    /* if (DEBUG == 2) fprintf(log_file,"ENTER[setup_prisms]\n"); */
     
    P.row = (int)ceil((P.max_northing - P.min_northing) / P.sp);
    
    P.col = (int)ceil((P.max_easting - P.min_easting) / P.sp);
     
    P.N_units = P.row * P.col;
  
    fprintf(log_file, 
	    "Number of rows = %d\nNumber of cols = %d\nNumber of Prisms = %d\n", 
	    P.row, P.col, P.N_units);
  
    pr = (PRISM *)GC_MALLOC((size_t)P.N_units * sizeof(PRISM));
    if (pr == NULL) {
      fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for prisms:[%s]\n",
            my_rank, procs, strerror(errno));
      return -1;
    } 
    
    count = 0;
    /* for ( ymin = P.min_northing, y = 0; y < P.row; ymin += P.sp, y++) {*/
    for (ymax = P.max_northing, y = 0; y < P.row; ymax -= P.sp, y++) {
      for ( xmin = P.min_easting, x = 0; x < P.col; xmin += P.sp, x++) {
        /*(pr+count)->south = ymin+.0001;*/
        (pr+count)->south = ymax - P.sp +.0001;
        (pr+count)->north = ymax + .0001;
        /*(pr+count)->north = ymin + P.sp + .0001;*/
        (pr+count)->west = xmin+.0001;
        (pr+count)->east = xmin + P.sp + .0001;
        (pr+count)->depth_to_bottom = 1.0;
        count++;
      }
    }
  
    for (i = 0; i < count; i++) {
      fprintf (log_file, "[%d]: %f to %f,  %f to %f\n", i,
      (pr+i)->west,
	     (pr+i)->east,
	     (pr+i)->south,
	     (pr+i)->north);
    }     		
    fflush(log_file);
   return P.N_units;
}
  
/**************************************************************
FUNCTION:  rmse
DESCRIPTION:  
INPUTS: NONE
OUTPUTS: double rmse 
***************************************************************/
double rmse(void) {
  int i;
  double rmse=0.0, error;
  
  for (i=0; i < total_pts; i++) {
  	 error = (p_all+i)->calculated - (p_all+i)->observed; 
    rmse += (error*error);
  }
  rmse /= (total_pts);
  rmse = sqrt(rmse);
  return rmse;
}

/*****************************************************************
FUNCTION: minimizing_func
DESCRIPTION: this is where the nodes assign new parameter values 
to the prisms, calculate new magnetic values and compare the
calculated value with the observed value using the chi-squared test. 
INPUTS: (IN)  double param[]  (an array of new prism parameters) 
RETURN:  double, the result of the rmse or chi-squared test
 *****************************************************************/
double minimizing_func(double param[]) {

  int i, ret;
  double fit;
  
 /* if (DEBUG == 2) fprintf(log_file, "  ENTER[minimizing_func]node=%d\n", my_rank); */
 // fprintf(stderr, "  ENTER[minimizing_func]node=%d\n", my_rank);
  
  if ( !my_rank ) {
      
    /* Send the updated parameters to the slave nodes */  
    
    for (i = 1; i < procs; i++) {
    	//fprintf(stderr, "  \tSending %d parameters to node %dof%d .. ",(int)NUM_OF_PARAMS, i, procs);
      ret = MPI_Send((void *)param, NUM_OF_PARAMS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    /*  if (DEBUG == 2) fprintf(stderr, "  \tParameters sent to node %d, MPIret=%d\n", i, ret); */
   // fprintf(stderr, "  Parameters sent to node %d, MPIret=%d\n", i, ret);
    }
  }


  /* Every node assigns the new parameters to their copy of the array of PRISM's */
  assign_new_params( param );
    
 /* Every node can now calculate A gbox (gravity) value for each of their subset of POINTs */ 
  for (i = 0;  i < num_pts;  i++) {
      (pt+i)->calculated = gbox(pt+i, pr, &P);  
  }
  
  /* Gather all of the calculated gravity values from each node into a single POINT array.
     This MPI function orders the bytes of data from each node in rank order (0, 1, etc).
     
     fprintf(stderr, "[NODE_%d] my_count=%d, my_data=%d bytes p_all=%d bytes\n",
     my_rank, my_count, (int)sizeof recv_ct, (int)sizeof p_all);
  */
  if (ret = MPI_Gatherv(pt, /* address of each proc's BYTE array of data */
			my_count, /* number of BYTES of data from each node */
			MPI_BYTE, /* data_type of each unit of sent of data */
			p_all, /* address of global point array */
			recv_ct, /* address of array of received data */
			displ, /* array of received displacements */
			MPI_BYTE, /* received datatype */
			0, /* root proc */
			MPI_COMM_WORLD), !ret)	{
				
      if ( !my_rank ) {
	  /* Only the master node calculates a new goodness-of-fit value */
	  
	     fit = rmse();
      } else 
	     fit = 0.0;
  }
  else {
      fprintf(stderr, "ERROR: ret=%d\n", ret);
      return 0.0;
  }
/*  if (DEBUG == 2) fprintf(log_file, "  EXIT[minimizing_func]\t[%d-of-%d] ret=%f\n\n", 
			  my_rank, procs, fit); */
  return fit;
}

/****************************************************************** 
FUNCTION:  assign_new_params
The function assigns updated parameter values to the anomaly being modeled.
This happens before each calulation of the magnetic value.
INPUTS: (IN)  double param[]  (an array of new prism parameters to be tested) 
RETURN:  none
*******************************************************************/
/* void assign_new_params( double param[]) { */
void assign_new_params(double param[]) {

  int num, row, col;
  
 P.density = param[DENSITY]; 
 P.depth_to_top = param[DEPTH_TO_TOP]; 
 
  /* assign the new parameters to the grid and calculate the grid border */
  create_grid(param, GRID, P);
  
  /* Now assign the new GRID values back to the prism model */    
  num = 0;
    for (row = 0; row < P.row; row++) {
      for (col = 0; col < P.col; col++) {
	     (pr+num)->depth_to_bottom = GRID[row][col];
	     num++;       
      }
    }
}

void set_LOG(FILE *log  ) {
  log_file = log;
}

void close_logfile(void) {
  fclose(log_file);
}

/****************************************************************************
FUNCTION: init_optimal_params
This function initially sets values for all sets of parameters. The values
are randomly chosen between the minimum and maximum values specified for that
parameter.
INPUTS: (IN/OUT) double op[][NUM_OF_PARAMS]  the 2-D array of parameter sets
RETURN:  none
******************************************************************************/
void init_optimal_params(double op[][NUM_OF_PARAMS]) { /* init_optimal_params */

  int vert=0, parm=0; 
  
  fprintf(stderr, "ENTER[init_optimal_params]: NUM_OF_PARAMS=%d \n", NUM_OF_PARAMS);
  srand(SEED);
  
  for (vert=0; vert < NUM_OF_VERTICES; vert++) { /* for loop */
      
    /* The first parameter is the surface_to_top parameter of the prisms.
	 For each set of possible parameters randomly select an initial single top value for 
	 the prisms. This random value should fall within the LO_PARAM - HI_PARAM range.
    */
    op[vert][DEPTH_TO_TOP] = 
	 (double)LO_PARAM(DEPTH_TO_TOP) + 
	 ((double)(HI_PARAM(DEPTH_TO_TOP) - LO_PARAM(DEPTH_TO_TOP)) * (double)rand()/(RAND_MAX+1.0));
    /*  if (DEBUG == 3) fprintf(stderr, "  param[%d][%d]=%f ", vert, SURF_TO_TOP, op[vert][SURF_TO_TOP]);*/
    
    /* The second parameter is the rock density. 
	 For each set of prisms randomly select an initial density value 
	 within the LO_INTENSITY - HI_INTENSITY range.
    */
    op[vert][DENSITY] = 
	 (double)LO_PARAM(DENSITY) + 
	 ((double)(HI_PARAM(DENSITY)-LO_PARAM(DENSITY)) * (double)rand()/(RAND_MAX+1.0));
    /* if (DEBUG == 3) fprintf(stderr, "  param[%d][%d]=%f ", vert, DENSITY, op[vert][DENSITY]); */
     
    /* The remaining parameters are the surface-to-bot values for each of the
	 prisms. Each prism initally gets a random value within the LO_PARAM - HI_PARAM
	 range.  This value must not be greater than the value selected for the 
	 surface-to-bottom value for the current parameter set.
    */
    for (parm = DEPTH_TO_BOT; parm < NUM_OF_PARAMS; parm++) { /* for loop */
	   op[vert][parm] = 
	   (double)LO_PARAM(DEPTH_TO_BOT) + 
	   ((double)(HI_PARAM(DEPTH_TO_BOT) - LO_PARAM(DEPTH_TO_BOT)) * (double)rand()/(RAND_MAX+1.0));
	   
	   if (op[vert][parm] < op[vert][DEPTH_TO_TOP]) { op[vert][parm] = op[vert][DEPTH_TO_TOP]; }
	   /*if (DEBUG == 3) fprintf(stderr, "  param[%d][%d]=%f ", vert, parm, op[vert][parm]); */
      } /* END for loop */
    } /* END for loop */
  
  /*
    for ( vert=0; vert < NUM_OF_VERTICES; vert++ ) {    
    printf("\tPrism[%d]: ", param);
    for ( param=0; param < NUM_OF_PARAMS; param++)
    printf("%f ", optimal_param[vert][param]); 
    }
  */
  fprintf(stderr, "\nEXIT[init_optimal_params].\n");
}

/*************************************************************************
FUNCTION:  printout_points
DESCRIPTION:  This function prints out the northing and easting
coordinates along with the stored calculated magnetic value to the
file "points.out".
INPUTS:  none
OUTPUTS:  none
 ************************************************************************/
void printout_points(void) {

  int i;
  FILE *out_pt;
  FILE *out;

  out_pt = fopen(CALCULATED_GRAV, "w");
  if (out_pt == NULL) {
    fprintf(stderr, "Cannot open CALCULATED_GRAV file=[%s]:[%s]. Printing to STDOUT.\n", 
	    CALCULATED_GRAV, strerror(errno)); 
    out = stdout;
  } else 
    out = out_pt;
    
  for (i=0; i < total_pts; i++) 
    fprintf(out, "%f %f %f \n", 
	    (p_all+i)->easting, (p_all+i)->northing, (p_all+i)->calculated);
	    
  if (out == out_pt) fclose(out);
}

/*************************************************************************
FUNCTION:   printout_model
DESCRIPTION:  This function prints out the prism locations and each prism's
set of parameters to the file "prisms.out".
INPUTS:  none
OUTPUTS:  none
 ************************************************************************/
void printout_model(void) {

  int i;  
  FILE *model;
  FILE *out;
  FILE *out2;

  model = fopen(PRISM_BOT_DEPTH, "w");
  out2 = fopen(PRISM_GEOMETRY, "w");
  if (model == NULL) {
    fprintf(stderr, 
	 "Cannot output model to file:[%s]. Printing to STDOUT.\n", strerror(errno)); 
    out = stdout;
  } else {
	  out = model;
  }
  
  if (out2 != NULL) {
	 fprintf(out2, "%f %.2f\n",
	 P.depth_to_top,
    P.density);
    for (i=0; i < P.N_units; i++)
	   fprintf(out2, "%f %f %f %f %.2f\n", 
	   (pr+i)->west,
		(pr+i)->east,
		(pr+i)->south,
		(pr+i)->north,
		(pr+i)->depth_to_bottom);
     (void) fclose(out2);
    }
    
    /*fprintf(stderr, "Printing out prism model\n"); */
    for (i=0; i < P.N_units; i++) {
	   /*fprintf(stderr, "[%d] %d ", P.N_units, i);*/
	   fprintf(out, "%f %f %f\n", 
		((pr+i)->west + (pr+i)->east)/2.0,
		((pr+i)->south + (pr+i)->north)/2.0,
		0.0 - (pr+i)->depth_to_bottom);
    }
    
  if (out == model) fclose(out);
}

/*************************************************************************
FUNCTION:   printout_parameters
DESCRIPTION:  This function prints out a time/date stamp and the 
parameters used and modified during the modeling process.
INPUTS:  double chi  The best chi value.
OUTPUTS:  none
 ************************************************************************/
void printout_parameters(double chi) {

  FILE *out;
  time_t mytime;
  int i;
  double mini;
  
  out = fopen(README, "w");
  if (out == NULL) {
    fprintf(stderr, 
	    "Cannot output parameters to file:[%s]. Printing to STDOUT.\n", 
	    strerror(errno)); 
    out = stdout;
  }
  mytime = time(&mytime);
  mini = P.depth_to_top;
  
  for (i=0; i < P.N_units; i++)
	    if ((pr+i)->depth_to_bottom > mini) mini = (pr+i)->depth_to_bottom;
  
  fprintf(out, "%s\nBest RMSE = %.1f\n\t%f to %f (Easting)\n\t%f to %f (Northing)\n\tSpacing = %.1f (meters)\n\nNumber of Parameters = %d\nNumber of Vertices = %d\nTolerance = %f\nPaP.rameter Ranges:\n\tRock Density: %2.f %2.f\n\tDepth-to-top: %.2f  %.2f\n\tDepth-to-bottom: %.2f  %.2f\n\t\nModeled Anomaly:\n\tBottom of Thickest Prism: %.2f (meters)\n\tTop Depth: %.2f(meters)\n\t\tRock Density: %.2f\n",
	  asctime(localtime(&mytime)),
	  chi,
	  P.min_easting, P.max_easting,
     P.min_northing, P.max_northing,
     P.sp,
	  NUM_OF_PARAMS,
	  NUM_OF_VERTICES,
	  TOLERANCE,
	  LO_PARAM(DENSITY), HI_PARAM(DENSITY),
	  LO_PARAM(DEPTH_TO_TOP), HI_PARAM(DEPTH_TO_TOP),
	  LO_PARAM(DEPTH_TO_BOT), HI_PARAM(DEPTH_TO_BOT),
     mini,
	  P.depth_to_top,
	  P.density);
  
  if (out != stdout) fclose(out);

}
