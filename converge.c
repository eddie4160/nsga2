/* Convergence analysis utilities */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

# define CONVERGENCE_WINDOW 20
# define REF_MARGIN 1.0e-6
# define CHANGE_EPS 1.0e-12

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
        if (is_feasible_nondominated_index(archive_pop, start + generation_size, start + i))
        {
            front[count].obj = archive_pop->ind[start+i].obj;
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
    int i;
    int m;
    int interior_count;
    int segment_count;
    int boundary_count;
    int *sorted_idx;
    double *norm_obj;
    double *segment_dist;
    double mean_distance;
    double dispersion;
    double denom;
    if (count < 3)
    {
        return 0.0;
    }
    sorted_idx = (int *)malloc(count*sizeof(int));
    for (i=0; i<count; i++)
    {
        sorted_idx[i] = i;
    }
    for (i=1; i<count; i++)
    {
        int key;
        int j;
        key = sorted_idx[i];
        j = i - 1;
        while (j>=0 && front[sorted_idx[j]].obj[0] > front[key].obj[0])
        {
            sorted_idx[j+1] = sorted_idx[j];
            j--;
        }
        sorted_idx[j+1] = key;
    }
    norm_obj = (double *)malloc(count*nobj*sizeof(double));
    for (m=0; m<nobj; m++)
    {
        double fmin;
        double fmax;
        fmin = front[sorted_idx[0]].obj[m];
        fmax = front[sorted_idx[count-1]].obj[m];
        if (fabs(fmax - fmin) <= CHANGE_EPS)
        {
            for (i=0; i<count; i++)
            {
                norm_obj[i*nobj + m] = 0.0;
            }
        }
        else
        {
            denom = fmax - fmin;
            for (i=0; i<count; i++)
            {
                norm_obj[i*nobj + m] = (front[sorted_idx[i]].obj[m] - fmin) / denom;
            }
        }
    }
    segment_count = count - 1;
    segment_dist = (double *)calloc(segment_count, sizeof(double));
    for (i=0; i<segment_count; i++)
    {
        double d2;
        d2 = 0.0;
        for (m=0; m<nobj; m++)
        {
            double diff;
            diff = norm_obj[(i+1)*nobj + m] - norm_obj[i*nobj + m];
            d2 += diff * diff;
        }
        segment_dist[i] = sqrt(d2);
    }
    boundary_count = 2;
    if (segment_count <= boundary_count)
    {
        free(sorted_idx);
        free(norm_obj);
        free(segment_dist);
        return 0.0;
    }
    interior_count = segment_count - boundary_count;
    mean_distance = 0.0;
    for (i=1; i<segment_count-1; i++)
    {
        mean_distance += segment_dist[i];
    }
    mean_distance /= (double)interior_count;
    if (mean_distance <= CHANGE_EPS)
    {
        free(sorted_idx);
        free(norm_obj);
        free(segment_dist);
        return 0.0;
    }
    dispersion = 0.0;
    for (i=1; i<segment_count-1; i++)
    {
        dispersion += fabs(segment_dist[i] - mean_distance);
    }
    dispersion /= ((double)interior_count * mean_distance);
    free(sorted_idx);
    free(norm_obj);
    free(segment_dist);
    return dispersion;
}

static double compute_relative_change (double current, double previous)
{
    if (fabs(previous) <= CHANGE_EPS)
    {
        return 0.0;
    }
    return (current - previous) / previous;
}

void report_convergence_metrics (population *archive_pop, int generations, int generation_size, const char *filename)
{
    int start_generation;
    int window_size;
    int window_index;
    int g;
    double *hv_series;
    double *delta_series;
    int *front_counts;
    point **fronts;
    double *reference;
    FILE *fpt;
    fpt = fopen(filename, "w");
    if (fpt == NULL)
    {
        printf("\n Could not open convergence output file: %s\n", filename);
        return;
    }
    window_size = (generations < CONVERGENCE_WINDOW) ? generations : CONVERGENCE_WINDOW;
    start_generation = generations - window_size;
    fprintf(fpt,"generation, front_size, HV, Delta, hv_change, delta_change\n");
    hv_series = (double *)calloc(window_size, sizeof(double));
    delta_series = (double *)calloc(window_size, sizeof(double));
    front_counts = (int *)calloc(window_size, sizeof(int));
    fronts = (point **)calloc(window_size, sizeof(point *));
    reference = (double *)malloc(nobj*sizeof(double));
    for (g=0; g<nobj; g++)
    {
        reference[g] = -INF;
    }
    for (g=start_generation; g<generations; g++)
    {
        int count;
        int w;
        window_index = g - start_generation;
        count = extract_generation_front(archive_pop, generation_size, g, &(fronts[window_index]));
        front_counts[window_index] = count;
        for (w=0; w<count; w++)
        {
            int m;
            for (m=0; m<nobj; m++)
            {
                if (fronts[window_index][w].obj[m] > reference[m])
                {
                    reference[m] = fronts[window_index][w].obj[m];
                }
            }
        }
    }
    for (g=0; g<nobj; g++)
    {
        if (reference[g] <= -INF/2.0)
        {
            reference[g] = 0.0;
        }
        reference[g] += REF_MARGIN;
    }
    for (g=start_generation; g<generations; g++)
    {
        double hv_change;
        double delta_change;
        window_index = g - start_generation;
        hv_series[window_index] = (nobj==2) ? compute_hypervolume_2d(fronts[window_index], front_counts[window_index], reference) : 0.0;
        delta_series[window_index] = compute_delta_dispersion(fronts[window_index], front_counts[window_index]);
        if (window_index == 0)
        {
            hv_change = 0.0;
            delta_change = 0.0;
        }
        else
        {
            hv_change = compute_relative_change(hv_series[window_index], hv_series[window_index-1]);
            delta_change = compute_relative_change(delta_series[window_index], delta_series[window_index-1]);
        }
        fprintf(fpt,"%d, %d, " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT "\n",
                g+1, front_counts[window_index], hv_series[window_index], delta_series[window_index], hv_change, delta_change);
    }
    for (g=0; g<window_size; g++)
    {
        free(fronts[g]);
    }
    free(reference);
    free(fronts);
    free(front_counts);
    free(hv_series);
    free(delta_series);
    fclose(fpt);
}
