#ifndef SORT_H
#define SORT_H 1

void sort_ptrackede_(float *ptrackede, int *nparam_ptr, int *ntracke_ptr, int *numberpe_ptr, int *g_iteration_no, int *iteration_no, void *comm, int *err);
void sort_ptrackedi_(float *ptrackede, int *nparam_ptr, int *ntracke_ptr, int *numberpe_ptr, int *g_iteration_no, int *iteration_no, void *comm, int *err);

void report_sorting_timing_();

#endif

