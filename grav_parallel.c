/* 
	 File Name:   grav_parallel.c

	 Program Name:  mag_parallel        
	 Subroutine Names: _exit_now(), main()
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
	 This computer program:
	 
	 a) Calculates the gravity anomaly due to a set of vertical sided prisms using
	 the forward solution proposed by Rao and Babu (1991, Geophysics).
	 
	 b) Solves for the forward solution at many grid points using parallel
	 programming techniques and MPI (Message Passing Interface).
	 
	 c) Compares observed and calculated gravity anomalies using standard
	 goodness-of-fit measures (i.e., chi-squared, rmse)
	 
	 d) Alters parameters (depth to top of prisms, depth to base of prisms,
	 magnetization) using the downhill simplex method in multiple dimensions (Nelder
	 and Meade, 1965, Computer Journal, 7: 308). The simplex method is modified to
	 incorporate rule sets specific to magnetic data modeling.
	 
	 e) outputs a best-fit model together with  calculated solution at each grid
	 point, which may be subsequently visualized using generic mapping tools or
	 similar software.
	 
	 
	 PROGRAMMING LANGUAGE:  ANSI C 
	 
	 USAGE: mpirun -np <number of processors> grav_parallel <configuration file> <data input file>
	 
	 GLOBAL VARIABLES:
	 NUM_OF_PARAMS : an integer,  the number of prism parameters that will be simultaneously modelled
	 NUM_OF_VERTICES : an integer, the number of vertices of the simplex model
                     (i.e. the number of sets of parameters being simultaneously modeled)
										 this value is always one greater that the NUM_OF_PARAMS
	 TOLERANCE :  the program runs until the goodness-of-fit  values, resulting from a comparison of the calculated
                with the observed gravity values, all fall within the range of this value
	 _LO[LAST_PARAM] :  an array of the minimum parameter values
	 _HI[LAST_PARAM] :  an array of the maximum parameter values

	 REFERENCES: 
	 
	 PROGRAM FLOW:
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <mpi.h>
#include "prototypes.h"

#define LOG_FILE "node_"

/* The following Global Variables are assigned some default values */
int NUM_OF_PARAMS = 0;
int NUM_OF_VERTICES = 1;
double TOLERANCE = 1.0e-2;
/*
int ROWS = 1;
int COLS = 1;
*/
double _LO[LAST_PARAM] = {0.1, 0.1, 0.1};
double _HI[LAST_PARAM] = {1.0, 1.0, 1.0};

/* Input file handles */					
FILE *in_pts = NULL;  
FILE *log_file = NULL;

/******************************************************************************************** 
Program: grav_parallel
Authors: Chuck and Laura Connor
Date: April 1, 2020
Language: Ansi C

Purpose:  This program calculates the gravity field due to a set of vertical -sided prisms
Data points where calculations are made are input from a file. 
Outputs: 3 files,
1: the x,y, point location and calculated gravity field.
2: prism locations and depth to bottom of each prism
3: boundaries of each prism
*/

/*
	Main program
*/

int main(int argc, char *argv[]) {
 
  char log_name[25];
  double quit = 0.0;
  int i;
  int my_rank; /* process rank of each node (local) */
  int procs; /* number of nodes used for processing */
  double chi;
  INPUTS In;

  /* Start up MPI */
  MPI_Init(&argc, &argv);

   /* Get my process rank, an integer */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* Check for correct number of comand line arguments */
  if (argc < 2 || argc > 2) {
    if (!my_rank)
      fprintf(stderr, 
	      " Check comand line arguments,\nUSAGE: %s <config file>\n\n", argv[0]);
    MPI_Finalize();
    return(0);
  }

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  
	/* Each node opens a file for logging */
  sprintf(log_name, "%s%d", LOG_FILE, my_rank);
  // fprintf(stderr, "%s\n", log_name);

  log_file  = fopen(log_name, "w+");
  if (log_file == NULL) {
    fprintf(stderr, "Cannot open LOG file=[%s]:[%s]. Exiting.\n", 
	    log_name, strerror(errno));
	    MPI_Finalize();
    return(0);
  }
  set_LOG(log_file);
  
  /* Initialize */
  In.points_file = NULL;
  
	/* Initialize the global variables with inputs from the configuration file. */
  if ( init_globals(argv[1], &In) ) {
    (void) fclose(log_file);
    MPI_Finalize();
    return(0);
  }
	/* Only the Master node (node_0) will print out some values to the console */
	// fprintf(stderr, "My rank: %d\n", my_rank);
	
  if (!my_rank)
    fprintf(stderr, 
"NUM_OF_PARAMS = %d\nNUM_OF_VERTICES = %d\nTOLERANCE = %f\nROC-Density: low=%.2f high=%.2f\nDepth-to-top: low=%.2f high=%.2f\nDepth-to-bottom: low=%.2f high=%.2f\n\n",
	    NUM_OF_PARAMS, 
	    NUM_OF_VERTICES, 
	    TOLERANCE,
	    LO_PARAM(DENSITY), HI_PARAM(DENSITY),
	    LO_PARAM(DEPTH_TO_TOP), HI_PARAM(DEPTH_TO_TOP),
	    LO_PARAM(DEPTH_TO_BOT), HI_PARAM(DEPTH_TO_BOT));
 
	 /* Read in the data values and locations */
    in_pts = fopen(In.points_file, "r"); 
    if (in_pts == NULL) {
      fprintf(stderr, "Cannot open POINTS file=[%s]:[%s]. Exiting.\n", 
	     In.points_file, strerror(errno));
      (void) fclose(log_file);
      MPI_Finalize();
      return(0);
    }
    
    if (get_points(in_pts) ) {
      (void) fclose(in_pts);
      (void) fclose(log_file);
      MPI_Finalize();
      return(0);
   }
  
    /* finished with input - run the optimization */


   if ( my_rank ) { /* slave */
   /* the slave nodes now wait to be called
    * upon by the master to calculate their portion of the magnetic values
    */
    slave(my_rank, log_file); 
  }
  
  else { /* master */ 
    chi = master(); 

    /* Send all slaves a quitin' time signal (i.e. a single zero value */
    for ( i = 1; i < procs; i++ )
      MPI_Send((void *)&quit, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

		/* The Master node prints out a README file listing some input parameters and changed values */
		printout_parameters(chi);

  } /* end master code */

  (void) fclose(log_file);
  MPI_Finalize();
  return(0);
}



