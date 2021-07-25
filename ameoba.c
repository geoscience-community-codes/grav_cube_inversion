/* 
	 File Name:   ameoba.c

	 Program Name:  grav_parallel        
	 Subroutine Name(s): evaluate(), optimize_params(), smooth_model()
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
	 These subroutines optimize the values of the parameter sets in response to
    goodness-of-fit values returned by the minimizing function. Those parameter sets that return
	 the lowest goodness-of-fit  values are saved. Before a parameter set is tested it is smoothed
	 to screen out high frequency anomalies. When all parameter sets fall within the tolerance
   set by the user the function returns the best parameter set (i.e. the one that
   generated the lowest goodness-of-fit  when compared to the observed data value. If the inversion
	 does not converge within a maximum number of operations, the function returns the
	 best set received so far.

	 PROGRAMMING LANGUAGE:  ANSI C 
	 
	 GLOBAL VARIABLES:

	 DEBUG : if zero, no debugging output is produced
           if > zero, various debugging outputs are written to stderr and the log files
	 NUM_OF_PARAMS : an integer,  the number of prism parameters that will be simultaneously modeled
	 NUM_OF_VERTICES : an integer, the number of vertices of the simplex model
                     (i.e. the number of sets of parameters being simultaneously modeled)
										 this value is always one greater that the NUM_OF_PARAMS

	 REFERENCES: Numerical Recipies
	 
	 PROGRAM FLOW:
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <gc.h>
#include "prototypes.h"

#define TINY 1.0e-10

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

//static double **MODEL_GRID = NULL;
/************************************************************************************
 * INPUTS:
 * double op[][NUM_OF_PARAMS]  :  (in/out) a 2-D array of optimal parameters
 * double mfv[]    :  (in) an array of minimizing function return values
 * double psum[]   :  (in/out) an array of values to be optimized
 * double (*funk)  :  (in) pointer to the minimizing function
 * int worst       :  (in) vertex with the highest value
 * double extrapolation_factor  :  (in)

 * RETURN:  double :  try - an optimized value from the minimizing function 
 ***************************************************************************************/ 

double evaluate(double op[][NUM_OF_PARAMS], double mfv[], double psum[], 
		double (*funk)(double []), int worst, double extrapolation_factor) {
			
  int param;
  double fac1,fac2,try,*ptry;
  int prism_param; 
  
  /* if (DEBUG == 2) fprintf(stderr, "ENTER[evaluate]: worst=%d\n", worst); */
  ptry = (double *)GC_MALLOC((size_t)NUM_OF_PARAMS * sizeof(double));

  fac1 = (1.0 - extrapolation_factor) / NUM_OF_PARAMS;
  fac2 = fac1 - extrapolation_factor;
  
  for (param = 0; param < NUM_OF_PARAMS; param++){
    ptry[param] = psum[param] * fac1 - op[worst][param] * fac2;
    prism_param = param;
    if (prism_param > 1) prism_param = DEPTH_TO_BOT;
    /*if (prism_param >= LAST_PARAM) prism_param = LAST_PARAM - 1;*/
    test_bounds(prism_param, &ptry[param], ptry[0]);
  }
 //smooth_model(&ptry[DEPTH_TO_BOT]);
 
  /* Evaluate the function at the trial vertex. */
  try = (*funk)(ptry);

  /* If <try> value is better than the worse, move the worse vertex. */
  if (try < mfv[worst]) {
    /* if (DEBUG) fprintf(stderr, "TRY=%.0f  ",try); */
    mfv[worst] = try; 
    for (param = 0; param < NUM_OF_PARAMS; param++) {
      psum[param] += ptry[param] - op[worst][param];
      op[worst][param] = ptry[param];
    }
  }
  /* if (DEBUG) fprintf(stderr, "EXIT[evaluate]: try=%f mfv[worst]=%f \n",try, mfv[worst]); */
 /* free(ptry); */
  return try;
}


/***********************************************************************************************
 * INPUTS:
 * double op[][NUM_OF_PARAMS]  :  a 2-D array,
 *                                each row represents a vertice of the simplex,
 *                                each column represents a parameter to be optimized (0..n-1)
 *                                (i.e. the coordinates of a vertex in n-dimensional space)
 *
 * double mfv[]  :  array of minimizing function returns for each vertex
 * double tol    :  tolerance
 * double (*funk):  minimizing_function

 * RETURN: none
 ************************************************************************************************/ 
void optimize_params(double op[][NUM_OF_PARAMS], double mfv[], double tol,
		     double (*funk)(double []), int *num_evals) {
  int param, vert;
  int worst; /* vertex with the highest value */
  int better; /* vertex with the next-highest value */
  int best; /* vertex with the lowest value */
  double rtol, sum, swap, save, try, *psum;
  /* int i; */

  
  fprintf(stderr, "ENTER[optimize_params] ...\n");

  psum=(double *)GC_MALLOC((size_t)NUM_OF_PARAMS * sizeof(double));
  if (psum == NULL) {
    fprintf(stderr, "\t[optimize_params]Cannot malloc memory for psum:[%s]\n",
	    strerror(errno));
    return;
  }
  *num_evals = 0;

/*MODEL_GRID = (double **)GC_MALLOC((size_t)ROWS * sizeof(double));
  if (MODEL_GRID == NULL) {
    fprintf(stderr, "Cannot malloc memory for MODEL_GRID rows:[%s]\n",
	    strerror(errno));
    return;
  } 
  else {
    for (i=0; i < ROWS; i++) {
      MODEL_GRID[i] = (double *)GC_MALLOC((size_t)COLS * sizeof(double));
      if (MODEL_GRID[i] == NULL) {
	fprintf(stderr, "Cannot malloc memory for grid row %d:[%s]\n",
		i, strerror(errno));
	return;
      }
    }
  } 
*/

  /* GET PSUM (i.e. sum up each column of parameter values) */
  for (param = 0; param < NUM_OF_PARAMS; param++) {
    for (sum = 0.0, vert = 0; vert < NUM_OF_VERTICES; vert++) 
      sum += op[vert][param];
    psum[param] = sum;

  } 

  for (;;) {
   
    best = 0;
    worst = (mfv[0] > mfv[1]) ? (better = 1,0) : (better = 0,1);
    
    for (vert=0; vert < NUM_OF_VERTICES; vert++) {
      if (mfv[vert] <= mfv[best]) best = vert;
      if (mfv[vert] > mfv[worst]) {
	     better = worst;
	     worst = vert;
      } else if (mfv[vert] > mfv[better] && vert != worst) better = vert;
    }
    
    rtol = 2.0 * fabs(mfv[worst] - mfv[best]) / (fabs(mfv[worst]) + fabs(mfv[best]) + TINY);
    if (rtol < tol) {
      SWAP(mfv[0], mfv[best])
	   for (param = 0; param < NUM_OF_PARAMS; param++) 
	     SWAP(op[0][param], op[best][param]) 
	   break;
    }
    
    if (*num_evals >= NMAX) {
      fprintf(stderr, "\t[optimize_params]NMAX[%d] exceeded\n",NMAX);
      SWAP(mfv[0], mfv[best])
	   for (param = 0; param < NUM_OF_PARAMS; param++) 
	     SWAP(op[0][param], op[best][param]) 
	   break;
    }
    
    *num_evals += 2;
 
    /* Begin a new iteration 
       First extrapolate by a factor of -1. 
    */

    try = evaluate(op, mfv, psum, funk, worst, -1.0);
    // fprintf(stderr, "%d[%.1f][%.1f]  ", *num_evals, try, mfv[best]);
    
    /* If <try> gives a result better than the best,
       then try an extra extapolation by a factor of 2.
    */
    if (try <= mfv[best]) {
    	// fprintf(stderr, "eval->Best  ");
      fprintf(stderr, "%d[%.4f]  ", *num_evals, try);
      try = evaluate(op, mfv, psum, funk, worst, 2.0);
    }
     
    /* If <try> is worse than the 'better' ,
     *  look for an intermediate 'better' . 
     */
    else if (try >= mfv[better]) {
    	// fprintf(stderr, "%d-BEST[%d]=%.1f  ", *num_evals, better, try);
      save = mfv[worst]; 
      fprintf(stderr, "^");    
      try = evaluate(op, mfv, psum, funk, worst, 0.5);    
      
      /* If <try> is still worse than the worst,
	      contract around the best vertex. */
      if (try >= save) {
        fprintf(stderr, "<>");	
	     for (vert = 0; vert < NUM_OF_VERTICES; vert++) {
	       if (vert != best) {
	         for (param = 0; param < NUM_OF_PARAMS; param++)
	           op[vert][param] = psum[param] = 0.5 *(op[vert][param] + op[best][param]);
	         mfv[vert] = (*funk)(psum);
	       }
	       // else fprintf(stderr, "%d-contract_to_VERT[%d] ", *num_evals, vert);
	     }
	     *num_evals += NUM_OF_PARAMS;
	
	     /* GET PSUM (i.e. sum up each column of parameters) */
	     for (param = 0; param < NUM_OF_PARAMS; param++) {
	       for (sum = 0.0, vert = 0; vert < NUM_OF_VERTICES; vert++) sum += op[vert][param];
	       psum[param] = sum;
	     }
      }    
    } 
    
    else --(*num_evals);
    /* if (DEBUG) fprintf(stderr, "\t[optimize_params]NUM_EVAL=%d VERT=%d CHI=%f\n", *num_evals, best, try); */ 
    
    if (!(*num_evals % 1000)) {
      fprintf(stderr, "model->out ");
      printout_model();
      printout_points();
    }
  }
/*free(psum);
  
 for (i=0; i < ROWS; i++) {
    free(MODEL_GRID[i]);
  }
  free(MODEL_GRID); 
*/
  fprintf(stderr,"EXIT[optimize_params]: NUM_EVAL=%d RMSE=%f\n", *num_evals, mfv[best]);
}

