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

static int same_reported_design_variables (individual *ind1, individual *ind2)
{
    int j;
    for (j=0; j<nbin; j++)
    {
        double value1;
        double value2;
        value1 = clamp_value(round_to_nearest_integer(ind1->xbin[j]), min_binvar[j], max_binvar[j]);
        value2 = clamp_value(round_to_nearest_integer(ind2->xbin[j]), min_binvar[j], max_binvar[j]);
        if (value1 != value2)
        {
            return 0;
        }
    }
    for (j=0; j<nreal; j++)
    {
        if (fabs(ind1->xreal[j] - ind2->xreal[j]) > 1.0e-9)
        {
            return 0;
        }
    }
    return 1;
}

static int is_feasible_archive_nondominated (population *pop, int size, int index)
{
    int j;
    if (pop->ind[index].constr_violation != 0.0)
    {
        return 0;
    }
    for (j=0; j<size; j++)
    {
        if (index!=j && check_dominance(&(pop->ind[j]), &(pop->ind[index])) == 1)
        {
            return 0;
        }
    }
    return 1;
}

static void report_human_readable_feasible_individual (individual *ind, FILE *fpt)
{
    int j;
    for (j=0; j<nbin; j++)
    {
        fprintf(fpt,"%.0f",clamp_value(round_to_nearest_integer(ind->xbin[j]), min_binvar[j], max_binvar[j]));
        fprintf(fpt,", ");
    }
    fprintf(fpt,OUTPUT_DOUBLE_FORMAT,ind->xreal[0]);
    for (j=0; j<nobj; j++)
    {
        fprintf(fpt,", " OUTPUT_DOUBLE_FORMAT,get_report_objective_value(j, ind->obj[j]));
    }
    for (j=0; j<nobj; j++)
    {
        fprintf(fpt,", " OUTPUT_DOUBLE_FORMAT,ind->obj_std[j]);
    }
    fprintf(fpt,"\n");
}

/* Function to print the information of a population in a file */
void report_pop (population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<popsize; i++)
    {
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].obj[j]);
        }
        if (ncon!=0)
        {
            for (j=0; j<ncon; j++)
            {
                fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].constr[j]);
            }
        }
        if (nreal!=0)
        {
            for (j=0; j<nreal; j++)
            {
                fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].xreal[j]);
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
        fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\n",pop->ind[i].crowd_dist);
    }
    return;
}

/* Function to print the information of feasible and non-dominated population in a file */
void report_feasible (population *pop, FILE *fpt)
{
    int i, j, k;
    if (nobj==2 && nreal==1 && nbin==8)
    {
        fprintf(fpt,"r1, r2, r3, r4, m1, m2, m3, m4, dh, Pt, Q, std_f1, std_f2\n");
        for (i=0; i<popsize; i++)
        {
            if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
            {
                report_human_readable_feasible_individual(&(pop->ind[i]), fpt);
            }
        }
        return;
    }
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (j=0; j<nobj; j++)
            {
                fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].obj[j]);
            }
            if (ncon!=0)
            {
                for (j=0; j<ncon; j++)
                {
                    fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].constr[j]);
                }
            }
            if (nreal!=0)
            {
                for (j=0; j<nreal; j++)
                {
                    fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].xreal[j]);
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
            fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}

void report_archive_feasible (population *pop, int size, FILE *fpt)
{
    int i, j, k;
    int *selected;
    selected = (int *)malloc(size*sizeof(int));
    for (i=0; i<size; i++)
    {
        selected[i] = is_feasible_archive_nondominated(pop, size, i);
    }
    if (nobj==2 && nreal==1 && nbin==8)
    {
        fprintf(fpt,"r1, r2, r3, r4, m1, m2, m3, m4, dh, Pt, Q, std_f1, std_f2\n");
        for (i=0; i<size; i++)
        {
            int duplicate;
            if (!selected[i])
            {
                continue;
            }
            duplicate = 0;
            for (j=0; j<i; j++)
            {
                if (selected[j] && same_reported_design_variables(&(pop->ind[j]), &(pop->ind[i])))
                {
                    duplicate = 1;
                    break;
                }
            }
            if (!duplicate)
            {
                report_human_readable_feasible_individual(&(pop->ind[i]), fpt);
            }
        }
        free(selected);
        return;
    }
    for (i=0; i<size; i++)
    {
        if (!selected[i])
        {
            continue;
        }
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].obj[j]);
        }
        if (ncon!=0)
        {
            for (j=0; j<ncon; j++)
            {
                fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].constr[j]);
            }
        }
        if (nreal!=0)
        {
            for (j=0; j<nreal; j++)
            {
                fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].xreal[j]);
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
        fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,OUTPUT_DOUBLE_FORMAT "\n",pop->ind[i].crowd_dist);
    }
    free(selected);
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
        for (j=0; j<nbin; j++)
        {
            fprintf(fpt,", %.0f",clamp_value(round_to_nearest_integer(pop->ind[i].xbin[j]), min_binvar[j], max_binvar[j]));
        }
        for (j=0; j<nreal; j++)
        {
            fprintf(fpt,", " OUTPUT_DOUBLE_FORMAT,pop->ind[i].xreal[j]);
        }
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,", " OUTPUT_DOUBLE_FORMAT,get_report_objective_value(j, pop->ind[i].obj[j]));
        }
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,", " OUTPUT_DOUBLE_FORMAT,pop->ind[i].obj_std[j]);
        }
        fprintf(fpt,"\n");
        (*member_id)++;
    }
    return;
}
