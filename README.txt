The main code to compile is labelled 
IsingModel.c
and the other script files are used by GNUPlot to output graphs. 

This code was written on Windows, and is known and tested to compile properly on a Windows system. The specific differece surrounds the terminals and commands for GNUPlot which were written to plot on a Windows computer with X11 separately installed.

Modifications to the code have been done to hopefully allow it to compile and run on a Mac. It was written by looking at how Gnuplot_From_C was written as that ran on a Mac and reverse-engineering it. 

The difference involves adding an extra definition in the main code which is then called before the plot. The definition of the X_SCRIPT is also changed, where X refers to specific plots. An example is given below:
-----------------------------------------------------------------------------------------
On Mac:
-----------------------------------------------------------------------------------------
#define GNUPLOT_EXE "gnuplot" //gnuplot defined to output on command line and plot
#define EPLOT_SCRIPT "EPlot.script" //contains energy plotting choices

and later on in the code

snprintf(command, sizeof(command), "%s %s",GNUPLOT_EXE, EPLOT_SCRIPT );
system(command);
----------------------------------------------------------------------------------------
On Windows:
----------------------------------------------------------------------------------------
#define EPLOT_SCRIPT "./EPlot.script" //contains energy plotting choices

and later on in the code

snprintf(command, sizeof(command), "%s", EPLOT_SCRIPT );
system(command);
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------


There is a difference in the .script file too, with the text below added for 
Windows at the beginning of the file:
----------------------------------------------------------------------------------------

#!/bin/gnuplot -persist
set terminal x11

----------------------------------------------------------------------------------------

Hopefully the Mac version will compile and run as it does on Windows, and as is 
demonstrated by the output seen in the report for this assignment.