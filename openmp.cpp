#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

#define _cutoff 0.01
#define _density 0.0005

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0,numthreads;
    double dmin, absmin=1.0,davg,absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    //number of particles
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    vector<bin_t> particle_bins;


    double gridSize = sqrt(n * _density);
    double binSize = _cutoff * 2;
    int binNum = int(gridSize/binSize) + 1;

    set_size( n );
    init_particles( n, particles );
    buildBins(particle_bins, particles, n);


    //  simulate a number of time steps
    double simulation_time = read_timer( );



    #pragma omp parallel private(dmin)
    {
      numthreads = omp_get_num_threads();
      printf("%d threads total", numthreads);
      for( int step = 0; step < 1000; step++ )
      {
          navg = 0;
          davg = 0.0;
  	      dmin = 1.0;
          //
          //  compute all forces
          //  reduction clause lists variables upon which a reduction operation will be done at the end of the parallel region
          #pragma omp for reduction (+:navg) reduction(+:davg)
          {
            //for each bin
            for( int i = 0; i < binNum; i++ )
            {
              for (int j = 0; j < binNum; j ++)
              {
                //get a vector of particles in the current bin
                bin_t& vec = particle_bins[i*binNum + j];
                //for each particle in current bin
                for(int k = 0; k < vec.size(); ++k)
                {
                  //initialize ax and ay to 0
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
          }

          //
          //  move particles
          //  openmp loop
          #pragma omp for
          {
            for( int i = 0; i < n; i++ )
                move( particles[i] );
          }

          if( find_option( argc, argv, "-no" ) == -1 )
          {
            //
            //  compute statistical data
            //  only done by master thread
            #pragma omp master
            if (navg) {
              absavg += davg/navg;
              nabsavg++;
            }

            //only one thread at a time can do it
            #pragma omp critical
            {
    	        if (dmin < absmin)
                {
                  absmin = dmin;
                }
            }
            //
            //  save if necessary
            //
            #pragma omp master
            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
            }//end if -no option is included.
      }//end for loop of steps (1000)
  }//end pragma omp private
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

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
