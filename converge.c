/* Convergence analysis utilities */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

# define REF_MARGIN 1.0e-6

typedef struct
{
    double *obj;
}
point;

static int compare_point_obj0 (const void *a, const void *b)
{
    point *pa;
    point *pb;
    pa = (point *)a;
    pb = (point *)b;
    if (pa->obj[0] < pb->obj[0])
    {
        return -1;
    }
    if (pa->obj[0] > pb->obj[0])
    {
        return 1;
    }
    return 0;
}

static int is_feasible_nondominated_index (population *archive_pop, int size, int index)
{
    int j;
    individual *candidate;
    candidate = &(archive_pop->ind[index]);
    if (candidate->constr_violation != 0.0)
    {
        return 0;
    }
    for (j=0; j<size; j++)
    {
        if (j == index)
        {
            continue;
        }
        if (archive_pop->ind[j].constr_violation != 0.0)
        {
            continue;
        }
        if (check_dominance(&(archive_pop->ind[j]), candidate) == 1)
        {
            return 0;
        }
    }
    return 1;
}

static int extract_cumulative_front (population *archive_pop, int size, point **out_front)
{
    int i;
    int count;
    point *front;
    front = (point *)malloc(size*sizeof(point));
    count = 0;
    for (i=0; i<size; i++)
    {
        if (is_feasible_nondominated_index(archive_pop, size, i))
        {
            front[count].obj = archive_pop->ind[i].obj;
            count++;
        }
    }
    *out_front = front;
    return count;
}

static double compute_hypervolume_2d (point *front, int count, double *reference)
{
    int i;
    double hv;
    double prev_f2;
    point *sorted;
    if (count <= 0)
    {
        return 0.0;
    }
    sorted = (point *)malloc(count*sizeof(point));
    for (i=0; i<count; i++)
    {
        sorted[i] = front[i];
    }
    qsort(sorted, count, sizeof(point), compare_point_obj0);
    hv = 0.0;
    prev_f2 = reference[1];
    for (i=0; i<count; i++)
    {
        double width;
        double height;
        width = reference[0] - sorted[i].obj[0];
        height = prev_f2 - sorted[i].obj[1];
        if (width > 0.0 && height > 0.0)
        {
            hv += width * height;
        }
        if (sorted[i].obj[1] < prev_f2)
        {
            prev_f2 = sorted[i].obj[1];
        }
    }
    free(sorted);
    return hv;
}

static double compute_delta_dispersion (point *front, int count)
{
    int i, m;
    int interior_count;
    int *is_boundary;
    double *distances;
    double mean_distance;
    double dispersion;
    if (count < 3)
    {
        return 0.0;
    }
    is_boundary = (int *)calloc(count, sizeof(int));
    distances = (double *)calloc(count, sizeof(double));
    for (m=0; m<nobj; m++)
    {
        int *indices;
        double fmin;
        double fmax;
        indices = (int *)malloc(count*sizeof(int));
        for (i=0; i<count; i++)
        {
            indices[i] = i;
        }
        /* Fallback: simple insertion sort for portability */
        {
            int j;
            for (i=1; i<count; i++)
            {
                int key;
                key = indices[i];
                j = i - 1;
                while (j>=0 && front[indices[j]].obj[m] > front[key].obj[m])
                {
                    indices[j+1] = indices[j];
                    j--;
                }
                indices[j+1] = key;
            }
        }
        fmin = front[indices[0]].obj[m];
        fmax = front[indices[count-1]].obj[m];
        is_boundary[indices[0]] = 1;
        is_boundary[indices[count-1]] = 1;
        if (fmax > fmin)
        {
            for (i=1; i<count-1; i++)
            {
                int idx;
                idx = indices[i];
                distances[idx] += (front[indices[i+1]].obj[m] - front[indices[i-1]].obj[m]) / (fmax - fmin);
            }
        }
        free(indices);
    }
    interior_count = 0;
    mean_distance = 0.0;
    for (i=0; i<count; i++)
    {
        if (!is_boundary[i])
        {
            mean_distance += distances[i];
            interior_count++;
        }
    }
    if (interior_count == 0)
    {
        free(is_boundary);
        free(distances);
        return 0.0;
    }
    mean_distance /= (double)interior_count;
    dispersion = 0.0;
    for (i=0; i<count; i++)
    {
        if (!is_boundary[i])
        {
            dispersion += fabs(distances[i] - mean_distance);
        }
    }
    dispersion /= (double)interior_count;
    free(is_boundary);
    free(distances);
    return dispersion;
}

void report_convergence_metrics (population *archive_pop, int generations, int generation_size, const char *filename)
{
    int g;
    int previous_count;
    double *hv_series;
    double *delta_series;
    FILE *fpt;
    fpt = fopen(filename, "w");
    if (fpt == NULL)
    {
        printf("\n Could not open convergence output file: %s\n", filename);
        return;
    }
    fprintf(fpt,"generation, cumulative_front_size, HV, Delta, hv_change, delta_change\n");
    hv_series = (double *)calloc(generations, sizeof(double));
    delta_series = (double *)calloc(generations, sizeof(double));
    previous_count = 0;
    for (g=0; g<generations; g++)
    {
        point *front;
        int count;
        int cumulative_size;
        int w;
        double *reference;
        double hv_change;
        double delta_change;
        cumulative_size = (g+1) * generation_size;
        count = extract_cumulative_front(archive_pop, cumulative_size, &front);
        reference = (double *)malloc(nobj*sizeof(double));
        for (w=0; w<nobj; w++)
        {
            reference[w] = -INF;
        }
        for (w=0; w<count; w++)
        {
            int m;
            for (m=0; m<nobj; m++)
            {
                if (front[w].obj[m] > reference[m])
                {
                    reference[m] = front[w].obj[m];
                }
            }
        }
        for (w=0; w<nobj; w++)
        {
            if (reference[w] <= -INF/2.0)
            {
                reference[w] = 0.0;
            }
            reference[w] += REF_MARGIN;
        }
        hv_series[g] = (nobj==2) ? compute_hypervolume_2d(front, count, reference) : 0.0;
        delta_series[g] = compute_delta_dispersion(front, count);
        if (previous_count > 0)
        {
            hv_change = ((double)(count - previous_count)) / (double)previous_count;
            delta_change = hv_change;
        }
        else
        {
            hv_change = 0.0;
            delta_change = 0.0;
        }
        fprintf(fpt,"%d, %d, " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT "\n",
                g+1, count, hv_series[g], delta_series[g], hv_change, delta_change);
        previous_count = count;
        free(reference);
        free(front);
    }
    free(hv_series);
    free(delta_series);
    fclose(fpt);
}
