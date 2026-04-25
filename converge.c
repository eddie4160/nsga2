/* Convergence analysis utilities */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "global.h"
# include "rand.h"

# define REF_MARGIN 1.0e-6
# define CHANGE_EPS 1.0e-12

typedef struct
{
    double *obj;
    int index;
}
point;

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

static int same_reported_design_variables_for_convergence (individual *ind1, individual *ind2)
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
        double value1;
        double value2;
        value1 = floor(ind1->xreal[j] * 1.0e10 + 0.5) / 1.0e10;
        value2 = floor(ind2->xreal[j] * 1.0e10 + 0.5) / 1.0e10;
        if (value1 != value2)
        {
            return 0;
        }
    }
    return 1;
}

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

static int compare_int_asc (const void *a, const void *b)
{
    int ia;
    int ib;
    ia = *((int *)a);
    ib = *((int *)b);
    if (ia < ib)
    {
        return -1;
    }
    if (ia > ib)
    {
        return 1;
    }
    return 0;
}

static void add_candidate_to_nondominated_set (population *archive_pop, int candidate_index, int *nd_indices, int *nd_count)
{
    int j;
    int write_index;
    if (archive_pop->ind[candidate_index].constr_violation != 0.0)
    {
        return;
    }
    for (j=0; j<(*nd_count); j++)
    {
        if (check_dominance(&(archive_pop->ind[nd_indices[j]]), &(archive_pop->ind[candidate_index])) == 1)
        {
            return;
        }
    }
    write_index = 0;
    for (j=0; j<(*nd_count); j++)
    {
        if (check_dominance(&(archive_pop->ind[candidate_index]), &(archive_pop->ind[nd_indices[j]])) == 1)
        {
            continue;
        }
        nd_indices[write_index] = nd_indices[j];
        write_index++;
    }
    *nd_count = write_index;
    nd_indices[*nd_count] = candidate_index;
    (*nd_count)++;
}

static void add_candidate_to_nondominated_set_2d_sorted (population *archive_pop, int candidate_index, int *nd_indices, int *nd_count)
{
    int left;
    int right;
    int pos;
    int remove_end;
    individual *candidate;
    double x0;
    double x1;
    candidate = &(archive_pop->ind[candidate_index]);
    if (candidate->constr_violation != 0.0)
    {
        return;
    }
    x0 = candidate->obj[0];
    x1 = candidate->obj[1];
    left = 0;
    right = *nd_count;
    while (left < right)
    {
        int mid;
        mid = (left + right) / 2;
        if (archive_pop->ind[nd_indices[mid]].obj[0] < x0)
        {
            left = mid + 1;
        }
        else
        {
            right = mid;
        }
    }
    pos = left;
    if (pos > 0)
    {
        individual *prev;
        prev = &(archive_pop->ind[nd_indices[pos-1]]);
        if (prev->obj[1] <= x1)
        {
            return;
        }
    }
    if (pos < *nd_count)
    {
        individual *next;
        next = &(archive_pop->ind[nd_indices[pos]]);
        if (next->obj[0] == x0 && next->obj[1] <= x1)
        {
            return;
        }
    }
    remove_end = pos;
    while (remove_end < *nd_count)
    {
        individual *cur;
        cur = &(archive_pop->ind[nd_indices[remove_end]]);
        if (cur->obj[1] < x1)
        {
            break;
        }
        remove_end++;
    }
    if (remove_end > pos)
    {
        memmove(&(nd_indices[pos]), &(nd_indices[remove_end]), ((*nd_count) - remove_end) * sizeof(int));
        *nd_count -= (remove_end - pos);
    }
    if (pos < *nd_count)
    {
        memmove(&(nd_indices[pos+1]), &(nd_indices[pos]), ((*nd_count) - pos) * sizeof(int));
    }
    nd_indices[pos] = candidate_index;
    (*nd_count)++;
}

static void update_running_pareto_archive_for_generation (population *archive_pop, int generation_size, int generation_index, int *nd_indices, int *nd_count)
{
    int i;
    int start;
    start = generation_index * generation_size;
    for (i=0; i<generation_size; i++)
    {
        int candidate_index;
        candidate_index = start + i;
        if (nobj == 2)
        {
            add_candidate_to_nondominated_set_2d_sorted(archive_pop, candidate_index, nd_indices, nd_count);
        }
        else
        {
            add_candidate_to_nondominated_set(archive_pop, candidate_index, nd_indices, nd_count);
        }
    }
}

static int snapshot_deduplicated_front (population *archive_pop, int *nd_indices, int nd_count, point **out_front)
{
    int i;
    int j;
    int count;
    int duplicate;
    int *sorted_nd;
    point *front;
    if (nd_count <= 0)
    {
        *out_front = (point *)malloc(sizeof(point));
        return 0;
    }
    sorted_nd = (int *)malloc(nd_count*sizeof(int));
    front = (point *)malloc(nd_count*sizeof(point));
    for (i=0; i<nd_count; i++)
    {
        sorted_nd[i] = nd_indices[i];
    }
    qsort(sorted_nd, nd_count, sizeof(int), compare_int_asc);
    count = 0;
    for (i=0; i<nd_count; i++)
    {
        duplicate = 0;
        for (j=0; j<i; j++)
        {
            if (same_reported_design_variables_for_convergence(&(archive_pop->ind[sorted_nd[j]]), &(archive_pop->ind[sorted_nd[i]])))
            {
                duplicate = 1;
                break;
            }
        }
        if (!duplicate)
        {
            front[count].obj = archive_pop->ind[sorted_nd[i]].obj;
            front[count].index = sorted_nd[i];
            count++;
        }
    }
    free(sorted_nd);
    *out_front = front;
    return count;
}

static void report_verbose_style_archive_individual (individual *ind, FILE *fpt, int generation, int member_id)
{
    int j;
    fprintf(fpt,"%d, %d",generation,member_id);
    for (j=0; j<nbin; j++)
    {
        fprintf(fpt,", %.0f",clamp_value(round_to_nearest_integer(ind->xbin[j]), min_binvar[j], max_binvar[j]));
    }
    for (j=0; j<nreal; j++)
    {
        fprintf(fpt,", " OUTPUT_DOUBLE_FORMAT,ind->xreal[j]);
    }
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

static void write_incremental_cumulative_pareto_verbose (population *archive_pop, point **fronts, int *front_counts, int generations, int generation_size, const char *filename)
{
    int g;
    int i;
    int member_id;
    int seen_count;
    int *seen_indices;
    FILE *fpt;
    fpt = fopen(filename, "w");
    if (fpt == NULL)
    {
        printf("\n Could not open cumulative Pareto output file: %s\n", filename);
        return;
    }
    fprintf(fpt,"number of generations, ID, x1(r1), x2(r2), x3(r3), x4(r4), x5(m1), x6(m2), x7(m3), x8(m4), x9(dh), f1(Pt), f2(Q), std_f1, std_f2\n");
    seen_indices = (int *)malloc(generations*generation_size*sizeof(int));
    seen_count = 0;
    member_id = 1;
    for (g=0; g<generations; g++)
    {
        for (i=0; i<front_counts[g]; i++)
        {
            int j;
            int is_new;
            int idx;
            idx = fronts[g][i].index;
            is_new = 1;
            for (j=0; j<seen_count; j++)
            {
                if (same_reported_design_variables_for_convergence(&(archive_pop->ind[seen_indices[j]]), &(archive_pop->ind[idx])))
                {
                    is_new = 0;
                    break;
                }
            }
            if (is_new)
            {
                seen_indices[seen_count] = idx;
                seen_count++;
                report_verbose_style_archive_individual(&(archive_pop->ind[idx]), fpt, g+1, member_id);
                member_id++;
            }
        }
    }
    free(seen_indices);
    fclose(fpt);
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
    int g;
    double *hv_series;
    double *delta_series;
    int *front_counts;
    point **fronts;
    int *nd_indices;
    int nd_count;
    double *reference;
    clock_t convergence_start_clock;
    clock_t phase_clock;
    double elapsed_seconds;
    FILE *fpt;
    convergence_start_clock = clock();
    fpt = fopen(filename, "w");
    if (fpt == NULL)
    {
        printf("\n Could not open convergence output file: %s\n", filename);
        return;
    }
    hv_series = (double *)calloc(generations, sizeof(double));
    delta_series = (double *)calloc(generations, sizeof(double));
    front_counts = (int *)calloc(generations, sizeof(int));
    fronts = (point **)calloc(generations, sizeof(point *));
    nd_indices = (int *)malloc(generations*generation_size*sizeof(int));
    nd_count = 0;
    reference = (double *)malloc(nobj*sizeof(double));
    {
        int total_size;
        int i;
        total_size = generations * generation_size;
        for (g=0; g<nobj; g++)
        {
            reference[g] = -INF;
        }
        for (i=0; i<total_size; i++)
        {
            if (archive_pop->ind[i].constr_violation != 0.0)
            {
                continue;
            }
            for (g=0; g<nobj; g++)
            {
                if (archive_pop->ind[i].obj[g] > reference[g])
                {
                    reference[g] = archive_pop->ind[i].obj[g];
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
    }
    fprintf(fpt,"reference_point");
    for (g=0; g<nobj; g++)
    {
        fprintf(fpt,", " OUTPUT_DOUBLE_FORMAT, reference[g]);
    }
    fprintf(fpt,"\n");
    fprintf(fpt,"generation, front_size, HV, Delta, hv_change, delta_change\n");
    printf("\n[converge] report_convergence_metrics started: generations=%d\n", generations);
    phase_clock = clock();
    for (g=0; g<generations; g++)
    {
        update_running_pareto_archive_for_generation(archive_pop, generation_size, g, nd_indices, &nd_count);
        front_counts[g] = snapshot_deduplicated_front(archive_pop, nd_indices, nd_count, &(fronts[g]));
        if ((g+1) % 10 == 0 || g == generations-1)
        {
            elapsed_seconds = ((double)(clock() - phase_clock)) / (double)CLOCKS_PER_SEC;
            printf("[converge] archive update progress: %d/%d generations processed (phase %.2fs)\n", g+1, generations, elapsed_seconds);
        }
    }
    write_incremental_cumulative_pareto_verbose(archive_pop, fronts, front_counts, generations, generation_size, "cumulative_pareto_verbose.txt");
    phase_clock = clock();
    for (g=0; g<generations; g++)
    {
        double hv_change;
        double delta_change;
        hv_series[g] = (nobj==2) ? compute_hypervolume_2d(fronts[g], front_counts[g], reference) : 0.0;
        delta_series[g] = compute_delta_dispersion(fronts[g], front_counts[g]);
        if (g == 0)
        {
            hv_change = 0.0;
            delta_change = 0.0;
        }
        else
        {
            hv_change = compute_relative_change(hv_series[g], hv_series[g-1]);
            delta_change = compute_relative_change(delta_series[g], delta_series[g-1]);
        }
        fprintf(fpt,"%d, %d, " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT ", " OUTPUT_DOUBLE_FORMAT "\n",
                g+1, front_counts[g], hv_series[g], delta_series[g], hv_change, delta_change);
        if ((g+1) % 10 == 0 || g == generations-1)
        {
            elapsed_seconds = ((double)(clock() - phase_clock)) / (double)CLOCKS_PER_SEC;
            printf("[converge] metric write progress: %d/%d generations processed (phase %.2fs)\n", g+1, generations, elapsed_seconds);
        }
    }
    for (g=0; g<generations; g++)
    {
        free(fronts[g]);
    }
    free(reference);
    free(nd_indices);
    free(fronts);
    free(front_counts);
    free(hv_series);
    free(delta_series);
    fclose(fpt);
    elapsed_seconds = ((double)(clock() - convergence_start_clock)) / (double)CLOCKS_PER_SEC;
    printf("[converge] report_convergence_metrics completed in %.2fs\n", elapsed_seconds);
}
