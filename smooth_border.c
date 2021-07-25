/* 
	 File Name:   smooth_border.c

	 Program Name:  grav_parallel        
	 Subroutine Name(s): smooth(double *param, PARAMETER P)
	                     create_grid(double *param, double **GRID, PARAMETER P)
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

	 REFERENCES: 
	 
	 PROGRAM FLOW:
*/

#include "parameters.h"
#include "common_structures.h"

/******************************************************************
FUNCTION: smooth
DESCRIPTION: This function smooths the border cells of the grid
             by giving the outer cells the same value as their 
             inner neighboring cells.
INPUTS:  (IN/OUT)  double *param : the array of prism parameters being optimized 
         (IN)  PARAMETER P : structure of model parameters
RETURN:  none
 ****************************************************************
void smooth(double *param, PARAMETER P, FILE *logfile) {

  int parm, ct, index1, index2;

 Now smooth the border cells 
  parm = DEPTH_TO_BOT;
  
  for (ct = 1; ct < P.col - 2; ct++) {
	fprintf(log_file, "ct=%d ", ct);
	  index1 = ct; 
	  index2 = ct + P.col;
	  param[parm+index1] = param[parm+index2];
	  
	  fprintf(logfile, "1=%d 2=%d ", index1, index2);
		
	 
	  index1 = (P.row - 1) * P.col + ct;
	  index2 = (P.row - 2) * P.col + ct;
	  param[parm+index1] = param[parm+index2];
	  
	  fprintf(logfile, "1=%d 2=%d ", index1, index2);   
	  
  }

  for (ct = 0; ct < P.row - 1; ct++) {
	  fprintf(log_file, "ct=%d ", ct);
	  index1 = ct * P.col;
	  index2 = ct * P.col + 1;
	  param[parm+index1] = param[parm+index2]; 
	  
	  fprintf(logfile, "1=%d 2=%d ", index1, index2);	   
	 
	  
	  index1 = ct * P.col + P.col - 1;
	  index2 = ct * P.col + P.col - 2;
	  param[parm+index1] = param[parm+index2];
	  
	  fprintf(logfile, "1=%d 2=%d ", index1, index2);	  
	 
  }
}
*/

/******************************************************************
FUNCTION: create_grid
DESCRIPTION: This function fills in the interior of the grid with
             the optimized parameters.
INPUTS:  (IN/OUT)  double *param : the array of prism parameters being optimized 
         (IN)  PARAMETER P : structure of model parameters
         (IN/OUT) double **GRID : 2-D array of grid cells
RETURN:  none
 *****************************************************************/
void create_grid(double *param, double **GRID, PARAMETER P) {
  
  int x, y, i;
  int parm = DEPTH_TO_BOT;
  
  /* Fill in the interior of the grid with the optimized parameters */
  for (y=1; y < P.row-1; y++)
    for (x=1; x < P.col-1; x++)
      GRID[y][x] = param[parm++];

  /* Now calculate the border cells.
     based on the interior neighbor */

  /* NORTH/SOUTH 
  GRID[1][0] = (GRID[1][1] + GRID[2][1])/2.0;
  GRID[P.col-2][0] = (GRID[P.col-3][1] + GRID[P.col-2][1])/2.0;
  GRID[1][P.row-1] = (GRID[1][P.row-2] + GRID[2][P.row-2])/2.0;
  GRID[P.col-2][P.row-1] = (GRID[P.col-3][P.row-2] + GRID[P.col-2][P.row-2])/2.0;
  
  GRID[1][0] = GRID[P.col-2][0] = GRID[1][P.row-1] = GRID[P.col-2][P.row-1] = param[SURF_TO_BOT];
  */
  i=0;
  while (i < 2) {
    if ( !i ) {
      y = 0;
      for (x = 1; x < P.col-1; x++)
	      /*GRID[y][x] = GRID[y+1][x];*/
	      GRID[y][x] = P.depth_to_top;
    }
    else {
      y = P.row-1;
      for (x = 1; x < P.col-1; x++)
	      /* GRID[y][x] = GRID[y-1][x]; */
	      GRID[y][x] = P.depth_to_top;
    }  
      /* if (x > 1 && x < P.col-2) 
	 GRID[x][y-1] = (GRID[x-1][y] + GRID[x][y] + GRID[x-1][y])/3.0; */      
    i++;
  }

  /* EAST/WEST 
  GRID[0][1] = (GRID[1][1] + GRID[2][1])/2.0;
  GRID[0][P.row-2] = (GRID[P.col-3][1] + GRID[P.col-2][1])/2.0;
  GRID[P.col-1][1] = (GRID[1][P.row-2] + GRID[2][P.row-2])/2.0;
  GRID[P.col-1][P.row-2] = (GRID[P.col-3][P.row-2] + GRID[P.col-2][P.row-2])/2.0;
  GRID[0][1] = GRID[0][P.row-2] = GRID[P.col-1][1] = GRID[P.col-1][P.row-2] = param[SURF_TO_BOT];
  */

  i=0;
  while (i < 2) {
    if (!i){
      x = 0;
      for (y = 1; y < P.row-1; y++)
	      /* GRID[y][x] = GRID[y][x+1]; */
	      GRID[y][x] = P.depth_to_top;
    }
    else {
      x = P.col-1;
      for (y = 1; y < P.row-1; y++)
	      /* GRID[y][x] = GRID[y][x-1]; */
	      GRID[y][x] = P.depth_to_top;
    } 
    i++;
  }

  /* 4 CORNERS 
  GRID[0][0] = GRID[1][1];
  GRID[P.row-1][0] = GRID[P.row-2][1];
  GRID[0][P.col-1] = GRID[1][P.col-2];
  GRID[P.row-1][P.col-1] = GRID[P.row-2][P.col-2];
  */
  GRID[0][0] = P.depth_to_top;
  GRID[P.row-1][0] = P.depth_to_top;
  GRID[0][P.col-1] = P.depth_to_top;
  GRID[P.row-1][P.col-1] = P.depth_to_top;
}
