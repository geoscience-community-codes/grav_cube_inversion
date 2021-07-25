/* 
	 File Name:   master.c

	 Program Name:  grav_parallel        
	 Subroutine Name(s): master()
	 Release Date:       April 1, 2020
	 Release Version:    1.0        
	 
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
	 This subroutine is run by the master node. The master node
   initializes each vertex's parameter set and then adjusts the parameters
   to find the optimal set of parameter values, (i.e. one that gives the
	 smallest goodness-of-fit  value when compared with the observed or measured values. 
   Slave nodes are called upon to help
	 calculate the new field values to test the adjusted prism parameter 
	 sets.
	 
	 PROGRAMMING LANGUAGE:  ANSI C 
	 
	 GLOBAL VARIABLES:

	 DEBUG : if zero, no debugging output is produced
           if > zero, various debugging outputs are written to stderr and the log files
	 NUM_OF_PARAMS : an integer,  the number of prism parameters that will be simultaneously modelled
	 NUM_OF_VERTICES : an integer, the number of vertices of the simplex model
                     (i.e. the number of sets of parameters being simultaneouly modelled)
										 this value is always one greater that the NUM_OF_PARAMS
	 TOLERANCE :  the program runs until the goodness-of-fit  values, resulting from a comparison of the calculated
                with the observed magnetic values, all fall within the range of this value
	 _LO[LAST_PARAM] :  an array of the minimum parameter values
	 _HI[LAST_PARAM] :  an array of the maximumj parameter values

	 REFERENCES: 
	 
	 PROGRAM FLOW:
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <mpi.h>
#include "prototypes.h"

/*********************************************************************
INPUTS:  none
RETURN:  double (the best goodness-of-fit  value)
 ********************************************************************/
double master(void) {

  int vert, param;
  
  /* A table containing values for parameters:
   * (surf_to_bot, south_edge, north_edge, east_edge, west_edge, etc.) 
   * for each vertex of the simplex. 
   * There must be one more vertex than the number of parameters. */
  double optimal_param[NUM_OF_VERTICES][NUM_OF_PARAMS];
  
  /* An array of values returned by the minimizing function, 
     for each vertex of the simplex. */
  double minimizing_func_value[NUM_OF_VERTICES];

  int num_evals; /* the number of function evaluations taken */

  /* the set of parameters we are trying to optimize */
  double param_val[NUM_OF_PARAMS]; 

 if (DEBUG) fprintf(stderr, "ENTER[master]\n");
    
    /* initial parameter guesses : optimal_parameter[vertex][parameter]*/
 
    init_optimal_params(optimal_param); 
    
    /*if (DEBUG == 2) {
      for (vert = 0; vert < ( NUM_OF_PARAMS + 1); vert++)
				for ( param=0; param < NUM_OF_PARAMS; param++) 
					fprintf(stderr, "\t[%d][%d]: %f\n", vert, param, optimal_param[vert][param]);
    }

  
    the number of vertices equals one more than the number of parameters being optimized */
    for (vert=0; vert < NUM_OF_VERTICES; vert++) { 
		//fprintf(stderr, "VERT=%d ",vert);
      for (param=0; param < NUM_OF_PARAMS; param++) {
        //fprintf(stderr, "%d ",param);
		  param_val[param] = optimal_param[vert][param];
		 // fprintf(stderr, "p_val=%f ",param_val[param]);
	   }
	  // fprintf(stderr, "\n");
      minimizing_func_value[vert] = minimizing_func( param_val );
      /*if (DEBUG) fprintf(stderr, "[%d]chi=%f\n", vert, minimizing_func_value[vert]);*/
      //fprintf(stderr, "[%d]chi=%f ", vert, minimizing_func_value[vert]);
    }
	
	 fprintf(stderr, "\n");

    /* the dimension of the simplex equals the number of parameters being optimized */
   // fprintf(stderr, "TOLERANCE = %e\n", (double)TOLERANCE);
    optimize_params(optimal_param, 
		    minimizing_func_value,  
		    TOLERANCE, 
		    minimizing_func, 
		    &num_evals);
    
    for ( vert=0; vert < NUM_OF_VERTICES; vert++ ) {
      fprintf(stderr,"[%d]chi=%f\n", vert, minimizing_func_value[vert]);
    }
    
    fprintf(stderr, "BEST FIT = %f\n", minimizing_func_value[0]); 
    /* for ( param=0; param < NUM_OF_PARAMS; param++) 
	  fprintf(stderr, "\tPrism[%d]: %f\n", param, optimal_param[i][param]); */
    
    printout_points();

    for (param=0; param < NUM_OF_PARAMS; param++)
      param_val[param] =  optimal_param[0][param];

    assign_new_params( param_val );
    printout_model();
    if (DEBUG) fprintf(stderr, "EXIT[master]\n");
    return  minimizing_func_value[0];
}










