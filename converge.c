/* Convergence analysis utilities */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

# define CONVERGENCE_WINDOW 20
# define REF_MARGIN 1.0e-6
# define STABLE_HV_EPS 1.0e-6
# define STABLE_DELTA_EPS 1.0e-6

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

static int extract_generation_front (population *archive_pop, int generation_size, int generation_index, point **out_front)
{
    int i;
    int count;
    int start;
    point *front;
    start = generation_index * generation_size;
    front = (point *)malloc(generation_size*sizeof(point));
    count = 0;
    for (i=0; i<generation_size; i++)
    {
        individual *ind;
        ind = &(archive_pop->ind[start+i]);
        if (ind->rank == 1 && ind->constr_violation == 0.0)
        {
            front[count].obj = ind->obj;
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
    double *hv_series;
    double *delta_series;
    FILE *fpt;
    fpt = fopen(filename, "w");
    if (fpt == NULL)
    {
        printf("\n Could not open convergence output file: %s\n", filename);
        return;
    }
    fprintf(fpt,"generation, front_size, HV, Delta, window_hv_range, window_delta_range, stable\n");
    hv_series = (double *)calloc(generations, sizeof(double));
    delta_series = (double *)calloc(generations, sizeof(double));
    for (g=0; g<generations; g++)
    {
        point *front;
        int count;
        int ws, we, w;
        double *reference;
        double hv_min;
        double hv_max;
        double d_min;
        double d_max;
        int stable;
        count = extract_generation_front(archive_pop, generation_size, g, &front);
        reference = (double *)malloc(nobj*sizeof(double));
        ws = (g - CONVERGENCE_WINDOW + 1 < 0) ? 0 : (g - CONVERGENCE_WINDOW + 1);
        we = g;
        for (w=0; w<nobj; w++)
        {
            reference[w] = -INF;
        }
        for (w=ws; w<=we; w++)
        {
            point *wfront;
            int wcount;
            int i, m;
            wcount = extract_generation_front(archive_pop, generation_size, w, &wfront);
            for (i=0; i<wcount; i++)
            {
                for (m=0; m<nobj; m++)
                {
                    if (wfront[i].obj[m] > reference[m])
                    {
                        reference[m] = wfront[i].obj[m];
                    }
                }
            }
            free(wfront);
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
        hv_min = hv_max = hv_series[ws];
        d_min = d_max = delta_series[ws];
        for (w=ws; w<=we; w++)
        {
            if (hv_series[w] < hv_min) hv_min = hv_series[w];
            if (hv_series[w] > hv_max) hv_max = hv_series[w];
            if (delta_series[w] < d_min) d_min = delta_series[w];
            if (delta_series[w] > d_max) d_max = delta_series[w];
        }
        stable = (g - ws + 1 >= CONVERGENCE_WINDOW && (hv_max - hv_min) <= STABLE_HV_EPS && (d_max - d_min) <= STABLE_DELTA_EPS) ? 1 : 0;
        fprintf(fpt,"%d, %d, " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", %d\n",
                g+1, count, hv_series[g], delta_series[g], hv_max-hv_min, d_max-d_min, stable);
        free(reference);
        free(front);
    }
    free(hv_series);
    free(delta_series);
    fclose(fpt);
}
