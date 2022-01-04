/*
------------------------------------------------------------------------------------------
PHYM004
Final Project - Option C - The 2D Ising Model
Date: 25/03/21
Author: Joe Salkeld 670008743
Line Length = 90 characters for readability
------------------------------------------------------------------------------------------

-----------------------------------------WELCOME------------------------------------------

This code was written to perform the 2D Ising model of a magnetic material. It has three 
possible outputs which the user selects at runtime. These show the effect of temperature
and the applied field on the energy, magnetisation, heat capacity and susceptibility of
the material. 

The code is compiled using a command 

$gcc IsingModel.c -Wall -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas 
-lm -o ising.exe

After compilation, the user has three options:

-B x		Applying a constant B field x and varying the temperature

-T x		Applying a constant temperature x and varying the B field

-D x y z 	Applies temperature x, B field y, and number of output images z

The code is then run like:

$ ./ising.exe -T 100

N.b When running -D, the dimensions of the system need to be input into the .script file
DomainPlot.script in order for the colour gradient to be properly scaled. A uniform 
magnetisation will cause the colour palette to be ignored and return a uniform white.

The operations in this code are called using a getoptlong() function in the main, adapted 
from Dr C.D.H Williams' code.
The code will output different graphs using GNUPlot depending on the option selected at
runtime. Listed below are the outputs from different options:

-B outputs:
	A graph of energy vs temperature
	A graph of magnetisation vs temperature
	A graph of heat capacity vs temperature
	A graph of susceptibility vs temperature
	
-T outputs:
	A graph of magnetisations vs applied field
	
-D outputs:
	A 2D graph of the lattice with coloured pixels representing magnetisation orientation
	
These script files are all included with the code.

Below is an example of the variables used for the report, along with their average runtime

-B -> DIM = 20, TEMPMAX = 2500, TEMPSTEP = 10, 			Time = 0m 29 s
-B -> DIM = 30, TEMPMAX = 2500, TEMPSTEP = 10, 			Time = 1m 02 s
-B -> DIM = 40, TEMPMAX = 2500, TEMPSTEP = 10, 			Time = 1m 53 s
-B -> DIM = 50, TEMPMAX = 2500, TEMPSTEP = 10, 			Time = 2m 55 s

-T -> DIM = 20, MAXFIELD = 20000, FIELDSTEP = 100, 		Time = 0m 22 s
-T -> DIM = 30, MAXFIELD = 20000, FIELDSTEP = 100, 		Time = 0m 48 s
-T -> DIM = 40, MAXFIELD = 20000, FIELDSTEP = 100, 		Time = 1m 28 s
-T -> DIM = 50, MAXFIELD = 20000, FIELDSTEP = 100, 		Time = 2m 15 s

-D -> DIM = 200, TEMP = 5, FIELD = 0, No. of Image = 5, Time = 0m 11 s
-D -> DIM = 300, TEMP = 5, FIELD = 0, No. of Image = 5, Time = 0m 26 s
-D -> DIM = 400, TEMP = 5, FIELD = 0, No. of Image = 5, Time = 1m 59 s

------------------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <gsl/gsl_rng.h>

//defining well known constants
#define KB 1.38064852e-23 //Boltzmann Constant J/K
#define MU_B 9.27400999e-24 //Bohr Magneton J/T 

#define J 6.44e-21 //Exchange constant in J. Values from [2] in report
//[Iron = 6.44e-21] [Cobalt = 8.24e-21], [Nickel = 3.24e-21].


//Dimensions of the system to be created (DIM x DIM)
#define DIM 30

//Script Name Definitions
#define EPLOT_SCRIPT "./EPlot.script" //contains energy plotting choices
#define MPLOT_SCRIPT "./MPlot.script" //contains magnetisation plotting choices
#define HCPLOT_SCRIPT "./CPlot.script" //contains heat capacity plotting choices
#define SUSPLOT_SCRIPT "./ChiPlot.script" //contains susceptibility plotting choices
#define BRILLOUIN_SCRIPT "./BrillPlot.script" //contains brillouin plotting choices
#define DOMAINPLOT_SCRIPT "./DomainPlot.script" //contains domain plotting choices


//Temperature Scale Definitions:
#define TEMPSTEP 10 //defines the temperature step in Kelvin
#define MAXTEMP 2500 //defines maximum temperature for constant B plots in Kelvin


//B Field Scale Definitions:
//B-Field needs to be on the order of 10^4 to match exchange energy, see report.
#define FIELDSTEP 100 //defines the B-field step in T
#define MAXFIELD 20000 //defines maximum B-field for Brillouin plots in T

//defines the number of times the system is evolved for statistical values
#define EVOLVE 4000 //a value of 4000 at least is recommended for good averages.

//A Macro which allows access of element (x,y) in the matrix by using the modulo operator
#define S(x,y) (Sysmem->matrix[(DIM+(x))%(DIM)][(DIM+(y))%(DIM)]) 

typedef struct system {
	int matrix[DIM][DIM];
	gsl_rng * Rangen;
} System;

typedef struct properties {
	double M;
	double C;
	double Temp;
	double Energy;
	double Energy_av;
	double Energy_avsquare;
	double Energy_stddev;
	double Mag_av;
	double Mag_avsquare;
	double Mag_stddev;
	double Chi;
	double B;
} Properties;


/*
------------------------------------------------------------------------------------------
This is a legacy function which is not used in the final code. However, it does print the 
matrix that has been generated straight to the console which is very helpful when wanting 
to see small matrices for debugging or checking purposes so I have left it in here for 
that purpose.
------------------------------------------------------------------------------------------
*/
/*
static void printing(struct system *Sysmem){
	//printing the matrix to console
	printf("\n");
	for (int x= 0; x <DIM; x++){
		printf("\n");
		for (int y = 0; y <DIM; y++){
			printf("%d\t", Sysmem->matrix[x][y]);
		}
	}
	printf("\n");
}
*/

/*
------------------------------------------------------------------------------------------
Function which generates a matrix of size DIM x DIM, consisting of 1 and -1 values 
randomly allocated through the use of a Mersenne Twister RNG. For more information
regarding the Mersenne Twister RNG, see [3] in the report. The RNG gives
a random number which is then checked to see if it is even or odd. Even random numbers
are assigned a 1 and odd are assigned a -1, and this is then stored in the matrix.
------------------------------------------------------------------------------------------
*/
static struct system * matgen(struct system *Sysmem){ 
	
	int n;
	
	for(int i = 0; i <DIM; i++){
		for (int j = 0; j<DIM; j++){
			//getting a random number from RNG and checking if it is even or odd
			int x = gsl_rng_get (Sysmem->Rangen);  
			if (x%2 == 0){
				n = 1;
			}
			if (x%2 !=0){
				n = -1;
			}
			Sysmem->matrix[i][j] = n;
		}
	}
	
	return Sysmem;	
}


/*
------------------------------------------------------------------------------------------
Function calculates the energy of the system by using equation 1 from the
report. 
------------------------------------------------------------------------------------------
*/
static struct properties * energy(struct system *Sysmem, struct properties *Prop){ 
	
	Prop->Energy = 0;
	
	for (int x=0; x<DIM; x++){
		for (int y=0; y<DIM; y++){
			Prop->Energy -= S(x,y)*(J *(S(x+1,y) + S(x,y+1)) + (MU_B * Prop->B));
		}
	}
	
	return Prop;	
}


/*
------------------------------------------------------------------------------------------
Function calculates the magnetisation of the system by summing the spin values over the 
entire matrix and then dividing by the number of elements in the matrix to get the 
average magnetisation, see equation 2 in the report.
------------------------------------------------------------------------------------------
*/
static struct properties * magnetisation(struct system *Sysmem, struct properties *Prop){ 
	
	//initialising the sum to zero
	double sum = 0.0;
	
	for (int x = 0; x < DIM; x++){
		for (int y = 0; y <DIM; y++){
			sum += S(x,y);
		}
	}
	
	Prop->M = sum / (DIM*DIM);
	
	return Prop;	
}

/*
------------------------------------------------------------------------------------------
Function calculates the heat capacity of the system, calculated using equation 3 from
the report. 
------------------------------------------------------------------------------------------
*/
static struct properties * heatcapacity(struct properties *Prop){ 
	
	Prop->C = (Prop->Energy_stddev * Prop->Energy_stddev)/(KB * Prop->Temp * Prop->Temp);
	
	return Prop;	
}


/*
------------------------------------------------------------------------------------------
Function calculates the susceptibility of the system, from equation 4 in the report.
------------------------------------------------------------------------------------------
*/
static struct properties * susceptibility(struct properties *Prop){ 
	
	Prop->Chi = (Prop->Mag_stddev * Prop->Mag_stddev)/(KB * Prop->Temp);
	
	return Prop;	
}


/*
------------------------------------------------------------------------------------------
Function calculates the evolution of the system and whether to flip spins or not. This 
is using equation 1 in the report, and the pseudo-code of the function is detailed in 
algorithm 1 in the report.
------------------------------------------------------------------------------------------
*/
static struct system * evolution(struct system *Sysmem, struct properties *Prop){
	
	double delta_E;
	int x,y;
	

	for(int i = 0; i <DIM*DIM; i++){
		//getting a random integer in range of 0 to DIM-1 from RNG for x and y
		x = gsl_rng_uniform_int(Sysmem->Rangen, DIM); 
		y = gsl_rng_uniform_int(Sysmem->Rangen, DIM);

		//calculating energy required for a flip
		delta_E = S(x,y)*(2*J*(S(x-1,y)+S(x,y-1)+S(x+1,y)+S(x,y+1)) + (MU_B * Prop->B)); 
		
		
		if (delta_E <= 0){
			S(x,y) = -S(x,y);
		}
		else{
			//getting a random number between 0-1 from RNG 
			double U = gsl_rng_uniform(Sysmem->Rangen); 
			
			//flipping the spin based on a Boltzmann factor
			if (U < exp(-delta_E/(KB * Prop->Temp))){ 
				S(x,y) = -S(x,y);
			}
		}
	}

	return Sysmem;
}

/*
------------------------------------------------------------------------------------------
Function calculates the statistics of the system for both the energy and the
magnetisation. The standard deviation of these two quantities is used in calculation of
the heat capacity and the susceptibility of the material, seen in equations 3,4 of 
the report. The first 1000 iterations are just evolved and not used for statistics to 
allow the system to equilibriate first.
------------------------------------------------------------------------------------------
*/
static struct properties * stats(struct system *Sysmem, struct properties *Prop){ 

	//allowing the system to reach equilibrium before calculating statistical values
	for (int n=0; n<1000;n++){
		evolution(Sysmem,Prop);
	}
	
	
	int x=1;
	
	//calculation of initial energy statistics
	energy(Sysmem,Prop);
	Prop->Energy_av = Prop->Energy;
	Prop->Energy_avsquare = Prop->Energy * Prop->Energy;
	
	//calculation of initial magnetisation statistics
	magnetisation(Sysmem,Prop);
	Prop->Mag_av = Prop->M;
	Prop->Mag_avsquare = Prop->M * Prop->M;
	
	//evolution of system 4000 times to allow oscillation about an average
	for(int i=0; i<EVOLVE;i++){
		evolution(Sysmem,Prop);
		
		energy(Sysmem,Prop);
		magnetisation(Sysmem,Prop);
		
		Prop->Energy_av = ((Prop->Energy_av * x) + Prop->Energy) / (x+1);
		//calculating the average of the squared energy <E^2>
		Prop->Energy_avsquare = 
			((Prop->Energy_avsquare * x) + (Prop->Energy*Prop->Energy)) / (x+1);
			
		Prop->Energy_stddev = 
			sqrt((Prop->Energy_avsquare)-(Prop->Energy_av * Prop->Energy_av));
		
		Prop->Mag_av = ((Prop->Mag_av * x) + Prop->M) / (x+1);
		//calculating the average of the squared magnetisation <M^2>
		Prop->Mag_avsquare = ((Prop->Mag_avsquare * x) + (Prop->M*Prop->M)) / (x+1);
		
		Prop->Mag_stddev = sqrt(fabs((Prop->Mag_avsquare)-(Prop->Mag_av * Prop->Mag_av)));
		
		x++;
	}
	return Prop;
}

/*
------------------------------------------------------------------------------------------
Main function of the code: Using getoptlong() to allow user input and choice of functions
to run. The function and options are adapted from Dr. C.D.H. Williams code "mat_gen.c".
------------------------------------------------------------------------------------------
*/
int main(int argc,char ** argv){ 
	
	char command[50];
	srand(time(0));
	static int normal_flg;

	//Allocating Memory based on the lattice size required
	System *Sysmem = malloc(DIM * sizeof(System));
	Properties *Prop = malloc(sizeof(Properties));
	
	//initialising a random number generator of Mersenne Twister
	Sysmem->Rangen = gsl_rng_alloc (gsl_rng_mt19937); 
	//setting the seed randomly
	gsl_rng_set(Sysmem->Rangen, rand()); 
	
	//generating the matrix of spins to be used
	matgen(Sysmem); 
	
	while (1) {
        static struct option long_options[] = {
            // These options set flags. 
            {"verbose", no_argument,      &normal_flg, 1},
            {"normal", no_argument,       &normal_flg, 0},
            // These options don’t set a flag, they are distinguished by their indices. 
            {"B Field", required_argument,  0, 'B'},
			{"Temperature", required_argument,  0, 'T'},
			{"Domains", required_argument, 0, 'D'},
            {0, 0, 0, 0}
        };
        
        // getopt_long needs somewhere to store its option index. 
        int option_index = 0;
        
        int c = getopt_long( argc, argv, "B:T:D:", long_options, &option_index );
        
        // End of options is signalled with '-1' 
        if (c == -1)
            break;
        
        switch (c) {
            case 0:
                // If this option sets a flag we have nothing to do. 
                if (long_options[option_index].flag != 0) {
                    break;
                }
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;
				
			//option for a constant B field and varying temperature	
            case 'B':
				//atof turns the string into a float, and anything not a number into 0.0
				Prop->B = atof(optarg);
				
				FILE *f_en = fopen ("Energy.dat", "w");
				FILE *f_mag = fopen ("Magnetisation.dat", "w");
				FILE *f_hcap = fopen ("HeatCapacity.dat", "w");
				FILE *f_sus = fopen ("Susceptibility.dat", "w");
					
				for (Prop->Temp = 1;Prop->Temp<MAXTEMP;Prop->Temp +=TEMPSTEP){
					//calculating statistical values and evolving the system
					stats(Sysmem,Prop); 
		
		
					energy(Sysmem,Prop);
					fprintf(f_en,"%f\t%e\n", Prop->Temp, Prop->Energy);
		
					magnetisation(Sysmem,Prop);
					fprintf(f_mag,"%f\t%e\n", Prop->Temp, fabs(Prop->M));
	
					heatcapacity(Prop);
					fprintf(f_hcap,"%f\t%e\n", Prop->Temp, Prop->C);
		
					susceptibility(Prop);
					fprintf(f_sus,"%f\t%e\n", Prop->Temp, Prop->Chi);
				}
				
				fclose(f_en);
				fclose(f_mag);
				fclose(f_hcap);
				fclose(f_sus);	
				
				//Command to invoke Gnuplot and plot the output data.
				snprintf(command, sizeof(command), "%s", EPLOT_SCRIPT );
				system(command);
	
				snprintf(command, sizeof(command), "%s", MPLOT_SCRIPT );
				system(command);
	
				snprintf(command, sizeof(command), "%s", HCPLOT_SCRIPT );
				system(command);	
	
				snprintf(command, sizeof(command), "%s", SUSPLOT_SCRIPT );
				system(command);
				break;	
			
			//option for a constant temperature and varying B Field 
			case 'T':
				//atof turns the string into a float, and anything not a number into 0.0
				Prop->Temp = atof(optarg);
				
				//error check for invalid temperature input
				if(Prop->Temp <= 0.0){
					printf("The value input for temperature is not valid");
					return 0;
				}

				FILE *f_brill = fopen ("Brillouin.dat", "w");
				
				//two 'for' loops are run to go from max field to min field and then back
				for (Prop->B = 0;Prop->B<MAXFIELD;Prop->B +=FIELDSTEP){
					//calculating statistical values and evolving the system
					stats(Sysmem,Prop); 
					
					magnetisation(Sysmem,Prop);
					fprintf(f_brill,"%f\t%e\n", Prop->B, Prop->M);
				}
				
				fclose(f_brill);
				
				//Command to invoke Gnuplot and plot the output data.
				snprintf(command, sizeof(command), "%s", BRILLOUIN_SCRIPT );
				system(command);
				break;
			
			//option to vary B field and temperatue and output domain images
			case 'D':;
			//Order of inputs are Temperature, B Field, Number of images
			
				//atof turns the string into a float, and anything not a number into 0.0
				Prop->Temp = atof(argv[optind-1]);
				Prop->B = atof(argv[optind]);
				
				//atoi turns the string into an integer, returning zero for an error 
				//and automatically rounding off any decimal places.
				int output_number = atoi(argv[optind+1]);
				
				optind+=2;
				
				//error check for invalid temperature input
				if(Prop->Temp <= 0.0){
					printf("The value input for temperature is not valid");
					return 0;
				}
				
				//error check for invalid number of images input
				if(output_number <= 0){
					printf("The value input for number of requested images is not valid");
					return 0;
				}
				
				//calculating at what interval an image should be printed
				int image = 5000 / (output_number);
				
				for(int i = 0; i < 5000;i++){
					evolution(Sysmem, Prop);
					
					if(i % image == 0){
						FILE *f_dom = fopen ("Domain.dat", "w");
						for(int x = 0; x<DIM; x++){
							for(int y = 0; y<DIM; y++){
								//outputting the x,y values of the matrix and the spin
								fprintf(f_dom, "%d\t%d\t%d\n", x,y,Sysmem->matrix[x][y]);
							}
						}
						fclose(f_dom);
				
						snprintf(command, sizeof(command), "%s",DOMAINPLOT_SCRIPT);
						system(command);
					}
				}
				break;
        }
    }
	
	// Prints any remaining command line arguments (not options).
    if (optind < argc) {
        fprintf (stderr, "Error: Unrecognised arguments: ");
        while (optind < argc) {
            fprintf (stderr, "%s ", argv[optind++]);
        }
        fprintf (stderr, "\n");
    }

	//freeing the memory of the RNG
	gsl_rng_free(Sysmem->Rangen);
	
	free(Sysmem);
	free(Prop);
}