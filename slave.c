/* 
	 File Name:   slave.c

	 Program Name:  mag_parallel        
	 Subroutine Name(s): slave(int, FILE *)
	 Release Date:       April 1, 2020
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
	  This function is run by each of the slave nodes.
    The slave nodes wait for new prism parameters from the master node.
    After receiving the new parameters, the node calls  minimizing_func() 
    with the new prism parameters and calculates its part of the field
    values.

	 
	 PROGRAMMING LANGUAGE:  ANSI C 
	 
	 GLOBAL VARIABLES:

	 DEBUG : if zero, no debugging output is produced
           if > zero, various debugging outputs are written to stderr and the log files
	 NUM_OF_PARAMS : an integer,  the number of prism parameters that will be simultaneously modelled
	 NUM_OF_VERTICES : an integer, the number of vertices of the simplex model
                     (i.e. the number of sets of parameters being simultaneouly modelled)
										 this value is always one greater that the NUM_OF_PARAMS
	 TOLERANCE :  the program runs until the CHI values, resulting from a comparison of the calculated
                with the observed magnetic values, all fall within the range of this value
	 _LO[LAST_PARAM] :  an array of the minimum parameter values
	 _HI[LAST_PARAM] :  an array of the maximum parameter values

	 REFERENCES: 
	 
	 PROGRAM FLOW:
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gc.h>
#include "prototypes.h"

#define QUIT 0

static double *recv_buffer=NULL;

/******************************************************************
INPUTS:  (IN)  int my_rank  (this slave node's id)
         (IN)  FILE *log_file  (this node's log_file handle)
RETURN:  none
 *****************************************************************/
void slave(int my_rank, FILE *log_file) {

  int ret;
  MPI_Status status;
  
  fprintf(log_file, "Slave[%d] here, ready ....\n",my_rank);
  recv_buffer = (double *)GC_MALLOC((size_t)NUM_OF_PARAMS * sizeof(double));
  if (recv_buffer == NULL) {
    fprintf(log_file, "No room for receive buffer. Exiting/n");
      return;
  }
  for (;;) {
    ret = MPI_Recv( (void *)recv_buffer, NUM_OF_PARAMS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status );
    if ( recv_buffer[0] == QUIT ) { 
      fprintf(log_file, "\treceived QUIT [%d] . . .", ret);
      break;
    }
   /* fprintf(stderr, "Receiving buffer ..."); */
    ret = minimizing_func(recv_buffer);
   /* fprintf(stderr, "Ret=%d\n",ret); */
  }

  fprintf(log_file, "Slave exiting ret=%d.\n", ret);
}
