//
//  mc_sampling.c
//
//  Created by Alessandro Coretti on 11/20/18.
//  Copyright © 2018 Alessandro Coretti. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Maximum number of Monte Carlo iterations
#define _NMAXITER 10000000

//Maximum number of histogram bins
#define _NMAXBINS 1000

//Number of points used to show relaxation
#define _NEQUIL 100

int NITER;
double DX;
double X0;

double L;
int NBINS;
double DBIN;
double MAXhisto=0;

double M;
double W;
double BETA;

double X[_NMAXITER];
double E[_NMAXITER];

void ReadIn(void);
double Energy(double x);
void Sample(void);
void Histogram(void);
void Plot_Histogram(void);
void Plot_Configuration(void);

int main(void) {

    srand48(time(0));

    //Reading input parameters
    ReadIn();

    //Sampling of configurations
    Sample();

    //Producing histogram
    Histogram();

    //Plot histogram
    Plot_Histogram();

    //Plot configurations and energies
    Plot_Configuration();

    return 0;
}

void ReadIn(void) {

    char dummy;

    FILE *fp_input;

    if ((fp_input = fopen("mc_sampling.in", "r+")) == NULL){

        printf("\nReadIn() ERROR: File 'mc_sampling.in' not found.\nExecution aborted.\n\n");
        exit(EXIT_FAILURE);
    }

    fscanf(fp_input, "%c %*[^\n]\n", &dummy);
    fscanf(fp_input, "%d %*[^\n]\n", &NITER); //Number of MC iterations
    if (NITER > _NMAXITER) {
        printf("\nReadIn() ERROR: NITER greater than maximum allowed MC iterations.\nExecution aborted.\n\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fp_input, "%lf %*[^\n]\n", &DX); //Maximum MC displacement (absolute value)
    fscanf(fp_input, "%lf %*[^\n]\n", &X0); //Initial configuration

    fscanf(fp_input, "%c %*[^\n]\n", &dummy);
    fscanf(fp_input, "%lf %*[^\n]\n", &M); //Mass of the particle
    fscanf(fp_input, "%lf %*[^\n]\n", &W); //Characteristic frequency
    fscanf(fp_input, "%lf %*[^\n]\n", &BETA); //Inverse of temperature (in units of kb)

    fscanf(fp_input, "%c %*[^\n]\n", &dummy);
    fscanf(fp_input, "%lf %*[^\n]\n", &L); //Histogram domain
    fscanf(fp_input, "%d %*[^\n]\n", &NBINS); //Number of histogram bins
    if (NBINS > _NMAXBINS) {
        printf("\nReadIn() ERROR: NBINS greater than maximum allowed histogram bins.\nExecution aborted.\n\n");
        exit(EXIT_FAILURE);
    }
    DBIN = (L/(double)NBINS); //Width of the bin

    fclose(fp_input);

    printf("\n** Single-particle harmonic oscillator Monte Carlo sampling **\n");

    printf("\n MC parameters:\n");
    printf("  Number of MC iterations: %d\n", NITER);
    printf("  MC displacement interval: [-%.3lf, +%.3lf]\n", DX, DX);
    printf("  Initial configuration: X = %.2lf\n", X0);

    printf("\n Physical parameters:\n");
    printf("  Mass: %.3lf\n", M);
    printf("  Characteristic frequency: %.3lf\n", W);
    printf("  Inverse of temperature (kb units): %.3lf\n", BETA);

    printf("\n Histogram parameters:\n");
    printf("  Histogram domain: [-%.2lf, +%.2lf]\n", .5*L, .5*L);
    printf("  Number of histogram bins: %d\n", NBINS);

    return;
}

//Potential energy function
double Energy(double x){

    return .5*M*(W*W)*x*x;
}

//Routine for sampling configurations with Metropolis Monte Carlo method
void Sample(void){

    int iter, acc=0;
    double xold, dx, eold, xnew, enew, ranf;

    FILE *fp_conf;

    fp_conf = fopen("configuration.out", "w+");

    //Initialize and saving initial configuration
    X[0] = xold = X0;
    E[0] = Energy(X0);

    fprintf(fp_conf, "%d\t%lf\t%lf\n", 0, X[0], E[0]);

    //Monte Carlo steps
    for (iter=1; iter<NITER; iter++) {

        //Computing energy for "old" configuration
        eold = Energy(xold);

        //Proposing a "new" position
        dx = 2.*DX*((double)lrand48()/(RAND_MAX+1.) - .5);
        xnew = xold + dx;

        //Implementation of periodic boundary conditions (PBC)
//        xnew -= _L*nearbyint(xnew/_L);

        //Computing energy for new proposed configuration
        enew = Energy(xnew);

        //Applying Metropolis algorithm:
        //Extracting a number uniformly distributed between 0 and 1
        ranf = (double)lrand48()/(RAND_MAX+1.);

        //If it is less than the Boltzmann factor of the energy difference between the
        //two configurations than accept the move
        if (ranf < exp(-BETA*(enew-eold))) {

            xold = xnew;
            eold = enew;
            acc++;
        }

        //Saving positions and energies for plotting
        X[iter] = xold;
        E[iter] = enew;

        if (iter < _NEQUIL) {

            fprintf(fp_conf, "%d\t%lf\t%lf\n", iter, X[iter], E[iter]);
        }
    }

    fclose(fp_conf);

    printf("\nAcceptance Ratio: %lf\n\n", (double)acc/(double)NITER);

    return;
}

void Histogram(void){

    int iter, i;
    int count[_NMAXBINS] = {0};

    FILE *fp_histo;

    fp_histo = fopen("histogram.out", "w+");

    //Looping over all Monte Carlo iterations
    for (iter=0; iter<NITER; iter++) {

        //Looping over all bins
        for (i=0; i<NBINS; i++) {

            //Checking which bin the configuration belongs to
            if (X[iter] > -.5*L+i*DBIN && X[iter] < -.5*L+(i+1)*DBIN) {

                //Incrementing counter
                count[i]++;
            }
        }
    }

    //Normalizing and saving data
    for (i=0; i<NBINS; i++) {

        fprintf(fp_histo, "%lf\t%lf\n", -.5*L+.5*DBIN+i*DBIN, (double)count[i]/DBIN/(double)NITER);

        if (MAXhisto < (double)count[i]/DBIN/(double)NITER) MAXhisto = (double)count[i]/DBIN/(double)NITER;
    }

    fclose(fp_histo);

    return;
}

void Plot_Histogram(void){

    FILE *fp_plot;

    fp_plot = fopen("histogram.gp", "w+");

    //Generating script
    fprintf(fp_plot, "reset\n");
    fprintf(fp_plot, "set samples 10000\n");
    fprintf(fp_plot, "set grid\n");
    fprintf(fp_plot, "set xrange[-%lf:%lf]\n", .5*L, .5*L);
    fprintf(fp_plot, "set yrange[0:%lf]\n", MAXhisto);
    fprintf(fp_plot, "set title 'Monte Carlo simulation with %d samples'\n", NITER);
    fprintf(fp_plot, "sigma = %lf\n", sqrt(1./(BETA*M*W*W)));
    fprintf(fp_plot, "k = %lf\n", M*W*W);
    fprintf(fp_plot, "f(x) = 1./(sqrt(2.*pi*sigma**2))*exp(-.5*x**2/sigma**2)\n");
    fprintf(fp_plot, "V(x) = .5*k*x**2\n");
    fprintf(fp_plot, "plot 'histogram.out' u 1:2 w p lc 7 pt 4 t 'Distribution of sampled configurations', f(x) w l lc 8 t 'Analytic result', V(x) w l lc 6 t 'Potential Energy'\n");
    fprintf(fp_plot, "pause -1\n");
        
    fclose(fp_plot);

    //Plotting distribution
    system("gnuplot histogram.gp");

    return;
}

void Plot_Configuration(void) {

    FILE *fp_plot;

    fp_plot = fopen("configuration.gp", "w+");

    //Generating script
    fprintf(fp_plot, "reset\n");
    fprintf(fp_plot, "set multiplot title 'Monte Carlo configurations and energy space'\n");
    fprintf(fp_plot, "set xlabel 'MC iteration'\n");
    fprintf(fp_plot, "set ylabel 'position'\n");
    fprintf(fp_plot, "set origin 0.,0.\n");
    fprintf(fp_plot, "set size 1.,.5\n");
    fprintf(fp_plot, "plot 'configuration.out' u 1:2 w l lc 7 t 'Sampled configurations'\n");
    fprintf(fp_plot, "set origin 0.,.45\n");
    fprintf(fp_plot, "set size 1.,.5\n");
    fprintf(fp_plot, "unset xlabel\n");
    fprintf(fp_plot, "set ylabel 'energy'\n");
    fprintf(fp_plot, "plot 'configuration.out' u 1:3 w l lc 7 t 'Sampled configuration energies'\n");
    fprintf(fp_plot, "unset multiplot\n");
    fprintf(fp_plot, "pause -1\n");

    fclose(fp_plot);

    //Plotting distribution
    system("gnuplot configuration.gp");

    return;
}
