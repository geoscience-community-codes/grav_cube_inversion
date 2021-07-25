/* 
	 File Name:   common_structures.h

	 Program Name:  grav_parallel        
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
	 
	 
	 PURPOSE: This file describes common data structures used within mag_parallel.
	 
	 PROGRAMMING LANGUAGE:  ANSI C 
*/

typedef struct point
{
  double easting;  /* easting coordinate of an observed/calculated location */
  double northing; /* northing coordinate of an observed/calculated location */
  double elev;  /* elevation of a location */
  double observed; /* the measured value at this location */
  double calculated; /* the calculated value at this location */
} POINT;

/* properties of a single prism */
typedef struct prism{
  /* boundaries of the prism w.r.t. the coordinate system being used, in meters */
  double south; /*south_edge*/ 
  double north; /*north_edge*/
  double west; /*west_edge*/
  double east; /*east_edge*/
	/* surface to bottom distance for a prism (meters)*/
  double depth_to_bottom; /*depth_to_bottom*/
  int b; /* border cell, b=0 or b=1, b=1 means it is a border cell */
} PRISM;

/* the model parameters */
typedef struct parameter{
	/* the individuals of the model will be created within these boundaries */
	double max_northing; /* maximum northing of the field being modeled */
	double min_northing; /* minimum northing of the field being modeled */
	double max_easting; /* maximum easting of the field being modeled */
	double min_easting; /* minimum easting of the field being modeled */
	double sp; /* the spacing used to create an individual of the model */
	double density; /* rock density */
	double depth_to_top; /* the depth to the top prism (same for all individuals) */
	int N_units; /* the total number of individuals in the model */
	int row; /* the total number of individual rows in the model */
	int col; /* the total number of individual columns in the model */
	int Npoints; /* the total number of points used to create the individuals (row * col)*/
} PARAMETER;

typedef struct inputs {
  char *points_file;
} INPUTS;
  






