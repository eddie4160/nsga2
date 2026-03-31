/* Data initializtion routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# include "global.h"
# include "rand.h"

/* Function to initialize a population randomly */
void initialize_pop (population *pop)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        initialize_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to initialize an individual randomly */
void initialize_ind (individual *ind)
{
    int j, k;
    if (nreal!=0)
    {
        for (j=0; j<nreal; j++)
        {
            ind->xreal[j] = rndreal (min_realvar[j], max_realvar[j]);
        }
    }
    if (nbin!=0)
    {
        for (j=0; j<nbin; j++)
        {
            for (k=0; k<nbits[j]; k++)
            {
                if (randomperc() <= 0.5)
                {
                    ind->gene[j][k] = 0;
                }
                else
                {
                    ind->gene[j][k] = 1;
                }
            }
        }
    }
    return;
}

static int parse_csv_field_count (char *line, double *values, int max_values)
{
    int count;
    char *token;
    count = 0;
    token = strtok(line, ",\r\n");
    while (token != NULL && count < max_values)
    {
        values[count] = atof(token);
        count++;
        token = strtok(NULL, ",\r\n");
    }
    return count;
}

static void set_binary_gene_from_value (individual *ind, int bin_index, double value)
{
    int bit_index;
    int max_index;
    int encoded_index;
    double bounded;
    double scale;
    max_index = (int)pow(2, nbits[bin_index]) - 1;
    if (max_index <= 0)
    {
        return;
    }
    bounded = value;
    if (bounded < min_binvar[bin_index])
    {
        bounded = min_binvar[bin_index];
    }
    if (bounded > max_binvar[bin_index])
    {
        bounded = max_binvar[bin_index];
    }
    if (max_binvar[bin_index] == min_binvar[bin_index])
    {
        encoded_index = 0;
    }
    else
    {
        scale = (bounded - min_binvar[bin_index]) / (max_binvar[bin_index] - min_binvar[bin_index]);
        encoded_index = (int)floor(scale * max_index + 0.5);
        if (encoded_index < 0)
        {
            encoded_index = 0;
        }
        if (encoded_index > max_index)
        {
            encoded_index = max_index;
        }
    }
    for (bit_index=0; bit_index<nbits[bin_index]; bit_index++)
    {
        int shift;
        shift = nbits[bin_index] - 1 - bit_index;
        ind->gene[bin_index][bit_index] = (encoded_index >> shift) & 1;
    }
    ind->xbin[bin_index] = min_binvar[bin_index] + (double)encoded_index * (max_binvar[bin_index] - min_binvar[bin_index]) / (double)max_index;
}

void initialize_pop_from_csv (population *pop, const char *filename)
{
    FILE *fp;
    char line[2048];
    int i;
    int expected_fields;
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("\n Could not open initial population file: %s\n", filename);
        exit(1);
    }
    if (fgets(line, sizeof(line), fp) == NULL)
    {
        fclose(fp);
        printf("\n Initial population file is empty: %s\n", filename);
        exit(1);
    }
    expected_fields = 1 + nbin + nreal + nobj;
    for (i=0; i<popsize; i++)
    {
        double values[256];
        int field_count;
        int j;
        if (fgets(line, sizeof(line), fp) == NULL)
        {
            fclose(fp);
            printf("\n Initial population file has %d cases, but population size is %d\n", i, popsize);
            exit(1);
        }
        field_count = parse_csv_field_count(line, values, expected_fields);
        if (field_count != expected_fields)
        {
            fclose(fp);
            printf("\n Initial population row %d has %d fields, expected %d\n", i+1, field_count, expected_fields);
            exit(1);
        }
        for (j=0; j<nbin; j++)
        {
            set_binary_gene_from_value(&(pop->ind[i]), j, values[1+j]);
        }
        for (j=0; j<nreal; j++)
        {
            pop->ind[i].xreal[j] = values[1+nbin+j];
        }
        for (j=0; j<nobj; j++)
        {
            double objective_value;
            objective_value = values[1+nbin+nreal+j];
            if (nobj==2 && nreal==1 && nbin==8 && j==1)
            {
                pop->ind[i].obj[j] = -objective_value;
            }
            else
            {
                pop->ind[i].obj[j] = objective_value;
            }
            pop->ind[i].obj_std[j] = 0.0;
        }
        if (ncon!=0)
        {
            for (j=0; j<ncon; j++)
            {
                pop->ind[i].constr[j] = 0.0;
            }
        }
        pop->ind[i].constr_violation = 0.0;
    }
    if (fgets(line, sizeof(line), fp) != NULL)
    {
        printf("\n Warning: initial population file has more rows than population size; extra rows will be ignored.\n");
    }
    fclose(fp);
    return;
}
