/* Memory allocation and deallocation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to allocate memory to a population */
void allocate_memory_pop (population *pop, int size)
{
    int i;
    pop->ind = new individual[size];
    for (i=0; i<size; i++)
    {
        allocate_memory_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to allocate memory to an individual */
void allocate_memory_ind (individual *ind)
{
    int j;
    if (nreal != 0)
    {
        ind->xreal = new double[nreal];
    }
    if (nbin != 0)
    {
        ind->xbin = new double[nbin];
        ind->gene = new int*[nbin];
        for (j=0; j<nbin; j++)
        {
            ind->gene[j] = new int[nbits[j]];
        }
    }
    ind->obj = new double[nobj];
    if (ncon != 0)
    {
        ind->constr = new double[ncon];
    }
    return;
}

/* Function to deallocate memory to a population */
void deallocate_memory_pop (population *pop, int size)
{
    int i;
    for (i=0; i<size; i++)
    {
        deallocate_memory_ind (&(pop->ind[i]));
    }
    delete [] pop->ind;
    return;
}

/* Function to deallocate memory to an individual */
void deallocate_memory_ind (individual *ind)
{
    int j;
    if (nreal != 0)
    {
        delete [] ind->xreal;
    }
    if (nbin != 0)
    {
        for (j=0; j<nbin; j++)
        {
            delete [] ind->gene[j];
        }
        delete [] ind->xbin;
        delete [] ind->gene;
    }
    delete [] ind->obj;
    if (ncon != 0)
    {
        delete [] ind->constr;
    }
    return;
}
