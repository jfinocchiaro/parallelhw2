#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"

using std::vector;

#define _cutoff 0.01 //copied from common.cpp
#define _density 0.0005 //copied from common.cpp
double binSize, gridSize;
int binNum;


//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    //create one vector that will hold all the particle bins
    vector<bin_t> particle_bins;
    bin_t temp;


    set_size( n );
    init_particles( n, particles );
    buildBins(particle_bins, particles, n);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
      	navg = 0;
        davg = 0.0;
      	dmin = 1.0;

        //
        //  compute forces
        //
        //for each bin in grid
        for(int i = 0; i < binNum; ++i)
        {
          for(int j = 0; j < binNum; ++j)
          {
            //current bin
            bin_t& vec = particle_bins[i*binNum+j];
            //for particle in bin, initialize ax and ay to 0
            for(int k = 0; k < vec.size(); ++k)
            {
              vec[k].ax = vec[k].ay = 0;
            }
            //check N/S and E/W neighbors
            for(int dx = -1; dx <= 1; ++dx)
            {
              for(int dy = -1; dy <= 1; ++ dy)
              {
                //if the bin you're checking for is actually in the grid
                if (i+dx >= 0 && i + dx < binNum && j + dy >= 0 && j+dy < binNum)
                {
                  bin_t& vectorholder = particle_bins[(i+dx) *binNum + j + dy];
                  //for every particle in original vector
                  for(int k = 0; k < vec.size(); ++k)
                  {
                    //for every particle in neighboring bin
                    for(int l = 0; l < vectorholder.size(); ++l)
                    {
                      apply_force(vec[k], vectorholder[l], &dmin, &davg, &navg);
                    }
                  }
                }
              }
            }
          }
        }



        //
        //  move particles
        //
        for(int i = 0; i < binNum; ++i)
        {
          for(int j = 0; j < binNum; ++j)
          {
            bin_t& vec = particle_bins[i * binNum + j];
            int tail = vec.size();
            int k = 0;
            for(k; k < tail;)
            {
              move(vec[k]);
              int x = int(vec[k].x / binSize);  //check position of moved particle
              int y = int(vec[k].y / binSize);  //check y-coordinate of moved particle
              if( x == i && y == j)
                ++k;
              else
              {
                temp.push_back(vec[k]);
                vec[k] = vec[--tail];
              }
            }
            vec.resize(k);
          }
        }


        //put moved particles into their new bins
        for(int i = 0; i < temp.size(); ++i)
        {
          int x = int(temp[i].x / binSize);
          int y = int(temp[i].y / binSize);
          particle_bins[x*binNum +y].push_back(temp[i]);
        }
        temp.clear();



        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;

          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    //
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %g\n",n,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
