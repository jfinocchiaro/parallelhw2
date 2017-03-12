#ifndef __MPIGRID_H_
#define __MPIGRID_H_

#include "mpicommon.h"
#define cutoff 0.01
#include <math.h>

struct linkedlist 
{
	linkedlist * next;
	int particle_id;
};

typedef struct linkedlist linkedlist_t;

struct grid
{
	int size;
	linkedlist_t ** grid;
	 
};

typedef struct grid grid_t;

//
// grid routines
//

void grid_init(grid_t & grid, int gridsize);
void grid_add(grid_t & grid, particle_t & particle);
bool grid_remove(grid_t & grid, particle_t & p, int gridCoord = -1);
void grid_clear(grid_t & grid);
void grid_clear_row(grid_t & grid, int r);
int  grid_size(grid_t & grid);


//
// Calculate the grid coordinate from a real coordinate
//
inline static int grid_coord(double c)
{
    return (int)floor(c / cutoff);
}
inline static int grid_coord_flat(int size, double x, double y)
{
    return grid_coord(x) * size + grid_coord(y);
}

#endif
