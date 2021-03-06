#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>
#include <stdio.h>


inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 100;

//
// particle data structure
//
typedef struct
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

typedef std::vector<particle_t> bin_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
double setSize(int n);
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void applyForce(particle_t &particle, particle_t &neighbor);
void move( particle_t &p );


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );
//void MPIsave(FILE *f, int rank, int n, particle_t *p, int * locals, int local_size, MPI_Datatype PARTICLE);

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
