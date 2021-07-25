/* 
	 File Name:   prototypes.h

	 Program Name:  gbox_parallel        
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
	 
	 
	 PURPOSE: This file lists prototypes of all functions used within the program.
	 
	 PROGRAMMING LANGUAGE:  ANSI C 
*/

#include "parameters.h"
#include "common_structures.h"

void optimize_params(double op[][NUM_OF_PARAMS], double mfv[], double tol,
double (*funk)(double []), int *num_evals);
/*void smooth_model(double *m);*/
double minimizing_func(double param[]);
void test_bounds(int param, double *try, double bound);
void init_optimal_params( double op[][NUM_OF_PARAMS]);
void assign_new_params( double []);
int init_globals(char *config_file, INPUTS *in);
int get_points(FILE *in);
void printout_model(void);
void printout_points(void);
void printout_parameters(double chi);
int setup_prisms(void);
void create_grid(double *param, double **GRID, PARAMETER P);
void slave(int my_rank, FILE *log_file);
double master(void);
void set_LOG(FILE *log_file);
double rmse(void);
double gbox(POINT *pt, PRISM *pr, PARAMETER *pa);
