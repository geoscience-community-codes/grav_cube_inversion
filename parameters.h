/* 
	 File Name:   parameters.h

	 Program Name:  grav_parallel        
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
	 
	 
	 PURPOSE: This file describes defined values and the global values.
	 
	 PROGRAMMING LANGUAGE:  ANSI C 
*/

enum {PRISMS}; /* specifies which method to use for calculating the gravity field */

/* identifies the parameters being modelled */
enum {DEPTH_TO_TOP, DENSITY, DEPTH_TO_BOT, LAST_PARAM}; 
/* enum {SURF_TO_BOT, INTENSITY, INC_ROC, DEC_ROC, SURF_TO_TOP, LAST_PARAM}; */ 
extern int ROWS;
extern int COLS;
extern int NUM_OF_PARAMS;
extern int NUM_OF_VERTICES;
extern double TOLERANCE;
extern double _LO[];
extern double _HI[];
 
#ifndef DEBUG
#define DEBUG 0
#endif

#define NMAX 500000
/*
#define POINTS_OUT "points.out"
#define PRISMS_OUT "prisms.out"
#define PRISMS_FWD "prisms.fwd"
*/
#define CALCULATED_GRAV "calculated_grav.out"
#define PRISM_GEOMETRY "prism_geometry.out"
#define PRISM_BOT_DEPTH "prism_bottoms.out"
#define PRISM_TOP_DEPTH "prism_tops.out"
#define LO_PARAM(p) (double)_LO[(p)]
#define HI_PARAM(p) (double)_HI[(p)]
