/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

static double round_to_nearest_integer (double value)
{
    if (value >= 0.0)
    {
        return floor(value + 0.5);
    }
    return ceil(value - 0.5);
}

static double clamp_value (double value, double lower, double upper)
{
    if (value < lower)
    {
        return lower;
    }
    if (value > upper)
    {
        return upper;
    }
    return value;
}

/* Function to print the information of a population in a file */
void report_pop (population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<popsize; i++)
    {
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        if (ncon!=0)
        {
            for (j=0; j<ncon; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
            }
        }
        if (nreal!=0)
        {
            for (j=0; j<nreal; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
            }
        }
        if (nbin!=0)
        {
            for (j=0; j<nbin; j++)
            {
                for (k=0; k<nbits[j]; k++)
                {
                    fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }
        fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
    }
    return;
}

/* Function to print the information of feasible and non-dominated population in a file */
void report_feasible (population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (j=0; j<nobj; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
            if (ncon!=0)
            {
                for (j=0; j<ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nreal!=0)
            {
                for (j=0; j<nreal; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
                }
            }
            if (nbin!=0)
            {
                for (j=0; j<nbin; j++)
                {
                    for (k=0; k<nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }
            fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}

/* Function to print verbose generation level information */
void report_verbose_generation (population *pop, FILE *fpt, int generation, int *member_id)
{
    int i, j;
    if (generation > 1)
    {
        fprintf(fpt,"\n");
    }
    fprintf(fpt,"#Gen %d\n",generation);
    for (i=0; i<popsize; i++)
    {
        fprintf(fpt,"%d, %d",generation,*member_id);
        for (j=0; j<nreal; j++)
        {
            fprintf(fpt,", %e",pop->ind[i].xreal[j]);
        }
        for (j=0; j<nbin; j++)
        {
            fprintf(fpt,", %.0f",clamp_value(round_to_nearest_integer(pop->ind[i].xbin[j]), min_binvar[j], max_binvar[j]));
        }
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,", %e",pop->ind[i].obj[j]);
        }
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,", %e",pop->ind[i].obj_std[j]);
        }
        fprintf(fpt,"\n");
        (*member_id)++;
    }
    return;
}
