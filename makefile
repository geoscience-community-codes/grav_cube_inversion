#	 File Name:   makefile
#
#	 Program Name:  grav_parallel        
#	 Release Date:         April 1, 2020
#	 Release Version:      1.0

#	 VERSION/REVISION HISTORY
#	 
#	 Date: December, 2003, Author: Laura Connor        
#	 Initial Release, Version 1.0
#	 
#	 8-20-2003
#	 Changed location of mpicc to not be hard-coded.

#	 DISCLAIMER/NOTICE
#	 
#	 This computer code/material was prepared as an account of work
#	 performed by the Center for Nuclear Waste Regulatory Analyses (CNWRA)
#	 for the Division of Waste Management of the Nuclear Regulatory
#	 Commission (NRC), an independent agency of the United States
#	 Government. The developer(s) of the code nor any of their sponsors
#	 make any warranty, expressed or implied, or assume any legal
#	 liability or responsibility for the accuracy, completeness, or
#	 usefulness of any information, apparatus, product or process
#	 disclosed, or represent that its use would not infringe on
#	 privately-owned rights.
#	 
#	 IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW WILL THE SPONSORS
#	 OR THOSE WHO HAVE WRITTEN OR MODIFIED THIS CODE, BE LIABLE FOR
#	 DAMAGES, INCLUDING ANY LOST PROFITS, LOST MONIES, OR OTHER SPECIAL,
#	 INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR
#	 INABILITY TO USE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
#	 BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY THIRD PARTIES OR A
#	 FAILURE OF THE PROGRAM TO OPERATE WITH OTHER PROGRAMS) THE PROGRAM,
#	 EVEN IF YOU HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES,
#	 OR FOR ANY CLAIM BY ANY OTHER PARTY.
#	 
#	 
#	 PURPOSE: This file is used by the unix program 'make' to direct
#  the compilation of the program.       

CC=mpicc
DEBUG=0
O=O3 
# W=Wfatal-errors
W=Wall

grav_parallel-bot:	master.o slave.o ameoba.o grav_parallel.o minimizing_func_new.o smooth_border.o gbox.o
		$(CC) -$(O) -$(W) -o grav_parallel-bot\
		master.o\
		slave.o\
		ameoba.o\
		grav_parallel.o\
		minimizing_func_new.o -lm\
		smooth_border.o\
		gbox.o -lgc -ldl

master.o:		master.c parameters.h makefile
			$(CC) -$(O) -$(W) -DDEBUG=$(DEBUG) -c master.c

slave.o:		slave.c parameters.h makefile
			$(CC) -$(O) -$(W) -DDEBUG=$(DEBUG) -c slave.c

ameoba.o:		ameoba.c parameters.h makefile
			$(CC) -$(O) -$(W) -DDEBUG=$(DEBUG) -c ameoba.c

minimizing_func_new.o:	minimizing_func_new.c common_structures.h makefile 
			$(CC) -$(O) -$(W) -DDEBUG=$(DEBUG) -c minimizing_func_new.c

gbox.o:			gbox.c common_structures.h makefile
			$(CC) -$(O) -$(W) -DDEBUG=$(DEBUG) -c gbox.c 

grav_parallel.o:	grav_parallel.c common_structures.h parameters.h prototypes.h makefile
			$(CC) -$(O) -$(W) -DDEBUG=$(DEBUG) -c grav_parallel.c

clean:
	rm *.o grav_parallel-bot
