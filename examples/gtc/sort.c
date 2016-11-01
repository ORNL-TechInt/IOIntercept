#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>

#include "mpi.h"

#define MAX_SORTING_COUNT 200

typedef struct _sorting_timing_info {
    double init_start;
    double local_sort_start;
    double shuffle_start;
    double merge_start;
    double merge_end;
    double sort_end;
    uint64_t data_sent;
    uint64_t data_received;
    uint64_t num_sender;
    uint64_t num_receiver;
    int iteration_no;
} sorting_timing_info;

static sorting_timing_info sorting_timing_e[MAX_SORTING_COUNT];
static int sorting_timing_count_e = 0;

static sorting_timing_info sorting_timing_i[MAX_SORTING_COUNT];
static int sorting_timing_count_i = 0;

#define DEFAULT_SEND_BUFFER_SIZE 1024*1024 * 8
#define INSCREMENT_BUFFER_SIZE 1024 * 8

static float **send_buffer = NULL;
static int *send_buffer_size = NULL;
static int *send_buffer_pos = NULL;
static float **sorting_buffer = NULL;
static int *sorting_buffer_size = NULL;
static int *sorting_buffer_pos = NULL;
static int num_pes;
static int mype;


int sort_particles(void *ptrackede_ptr, int *nparam_ptr, int *ntracke_ptr, int *numberpe_ptr, int *g_iteration_no, int *iteration_no, void *comm, sorting_timing_info *sorting_timing, int *sorting_timing_count);

void sort_ptrackede_(void *ptrackede_ptr, int *nparam_ptr, int *ntracke_ptr, int *numberpe_ptr, int *g_iteration_no, int *iteration_no, void *comm, int *err)
{
*err = sort_particles(ptrackede_ptr, nparam_ptr, ntracke_ptr, numberpe_ptr, g_iteration_no, iteration_no, comm, sorting_timing_e, &sorting_timing_count_e);
}

void sort_ptrackedi_(void *ptrackede_ptr, int *nparam_ptr, int *ntracke_ptr, int *numberpe_ptr, int *g_iteration_no, int *iteration_no, void *comm, int *err)
{
*err = sort_particles(ptrackede_ptr, nparam_ptr, ntracke_ptr, numberpe_ptr, g_iteration_no, iteration_no, comm, sorting_timing_i, &sorting_timing_count_i);
}


int sort_particles(void *ptrackede_ptr, int *nparam_ptr, int *ntracke_ptr, int *numberpe_ptr, int *g_iteration_no, int *iteration_no, void *comm, sorting_timing_info *sorting_timing, int *sorting_timing_count)
{
    float *write_buffer;
    double init_start, local_sort_start, shuffle_start, merge_start, merge_end, sort_end;
    uint64_t data_sent = 0, data_received = 0, no_send = 0, no_receive = 0;

    int i, j;
    int nparam = *nparam_ptr;
    int ntracke = *ntracke_ptr;
    int numberpe = *numberpe_ptr;
    float *ptrackede = (float *)ptrackede_ptr;

    int t = *(int *) comm;

    MPI_Comm group_comm = MPI_Comm_f2c(t);

    MPI_Comm_rank(group_comm, &mype);
    num_pes = numberpe;

    float *sorted_data = ptrackede;


    if(!send_buffer) {
        send_buffer = (float **) malloc(num_pes * sizeof(float *));
        if(!send_buffer) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
        memset(send_buffer, 0, num_pes * sizeof(float *));

        send_buffer_size = (int *) malloc(num_pes * sizeof(int));
        if(!send_buffer_size) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

        send_buffer_pos = (int *) malloc(num_pes * sizeof(int));
        if(!send_buffer_pos) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    }

    for(i = 0; i < num_pes; i ++) {
        send_buffer_size[i] = DEFAULT_SEND_BUFFER_SIZE;
        send_buffer[i] = (float *) malloc(send_buffer_size[i] * sizeof(float));
        if(!send_buffer[i]) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
        send_buffer_pos[i] = 0;
    }

    if(!sorting_buffer) {
        sorting_buffer = (float **) malloc(num_pes * sizeof(float *));
        if(!sorting_buffer) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
        memset(sorting_buffer, 0, num_pes * sizeof(float *));

        sorting_buffer_size = (int *) malloc(num_pes * sizeof(int));
        if(!sorting_buffer_size) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
        memset(sorting_buffer_size, 0, num_pes * sizeof(int));

        sorting_buffer_pos = (int *) malloc(num_pes * sizeof(int));
        if(!sorting_buffer_pos) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
        memset(sorting_buffer_pos, 0, num_pes * sizeof(int));
    }



    for(i = 0; i < ntracke; i ++) {
        int pe = (int)ptrackede[i * nparam + 7] - 1; 
        if(send_buffer_pos[pe] + nparam >= send_buffer_size[pe]) {
            float *temp = (float *)realloc(send_buffer[pe], (send_buffer_size[pe] + INSCREMENT_BUFFER_SIZE)*sizeof(float));
            if(!temp) {
                fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
                return -1;
            }
            send_buffer[pe] = temp; 
            send_buffer_size[pe] += INSCREMENT_BUFFER_SIZE;
        } 
        float *destination_pos = send_buffer[pe] + send_buffer_pos[pe]; 

        memcpy(destination_pos, &ptrackede[i * nparam], nparam * sizeof(float));                       
        send_buffer_pos[pe] += nparam;
    }

    MPI_Barrier(group_comm);


    // shuffle 
    MPI_Status status;

    int distance, right_ds, left_ds;
    int num_iterations = num_pes - 1;

    sorting_buffer[mype] = send_buffer[mype];
    sorting_buffer_size[mype] = send_buffer_size[mype];
    sorting_buffer_pos[mype] = send_buffer_pos[mype];
    
    for(i = 0; i < num_iterations; i ++) {
        // advance the distance
        distance = (1 + i) % num_pes;
        left_ds = (mype - distance) % num_pes;
        if(left_ds < 0) {
            left_ds += num_pes;
        }
        right_ds = (mype + distance) % num_pes;
        if(right_ds < 0) {
            right_ds += num_pes;
        }

        // receive from left side and send to right side
        int rc = MPI_Sendrecv(&send_buffer_pos[right_ds], 1, MPI_INTEGER, right_ds, distance,
                     &sorting_buffer_pos[left_ds], 1, MPI_INTEGER, left_ds, distance, group_comm, &status);

        if(send_buffer_pos[right_ds]) {
            sorting_timing[*sorting_timing_count].data_sent += send_buffer_pos[right_ds];
            sorting_timing[*sorting_timing_count].num_sender ++;
        }
        if(sorting_buffer_pos[left_ds]) {
            sorting_timing[*sorting_timing_count].data_received += sorting_buffer_pos[left_ds];
            sorting_timing[*sorting_timing_count].num_receiver ++;
        }

        // allocate receive buffer
        sorting_buffer_size[left_ds] = sorting_buffer_pos[left_ds];
        sorting_buffer[left_ds] = (float *) malloc(sorting_buffer_size[left_ds] * sizeof(float));
        if(!sorting_buffer[left_ds]) {
            fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

        rc = MPI_Sendrecv(send_buffer[right_ds], send_buffer_pos[right_ds], MPI_INTEGER, right_ds, distance,
                     sorting_buffer[left_ds], sorting_buffer_pos[left_ds], MPI_INTEGER, left_ds, distance, group_comm, &status);

        // free send buffer
        free(send_buffer[right_ds]);
        send_buffer[right_ds] = NULL;
        send_buffer_size[right_ds] = 0;
        send_buffer_pos[right_ds] = 0;
    }


    // merge sort
    write_buffer = ptrackede;
    uint64_t ntrack = 0;       
    for(i = 0; i < num_pes; i ++) {
        for(j = 0; j < sorting_buffer_pos[i]; j += nparam) {
            float *par = sorting_buffer[i] + j;
            int index = par[6] - 1;
            float *destination_pos = &(write_buffer[index * nparam]);
            memcpy(destination_pos, par, nparam * sizeof(float));
        }
        ntrack += sorting_buffer_pos[i]/nparam;
    }
    *ntracke_ptr = ntrack;



    if(send_buffer) {
        int i;
        for(i = 0; i < num_pes; i ++) {
            if(send_buffer[i]) {
                free(send_buffer[i]); 
            }
        } 
    }
    free(send_buffer);    
    free(send_buffer_size);    
    free(send_buffer_pos);    
    send_buffer = NULL;
    send_buffer_size = NULL;
    send_buffer_pos = NULL;

    if(sorting_buffer) {
        int i;
        for(i = 0; i < num_pes; i ++) {
            if(sorting_buffer[i] && i != mype) {
                free(sorting_buffer[i]);
            }
        }
    }
    sorting_buffer = NULL;
    sorting_buffer_size = NULL;
    send_buffer_pos = NULL;

    (*sorting_timing_count) ++;

    return 0;
}

void report_sorting_timing_()
{
    char file_name[50];
    sprintf(file_name, "sorting/sorting_timing_e%.5d\0", mype);
    FILE *fptr = fopen(file_name, "w");
    if(!fptr) {
        fprintf(stderr, "cannot open file %s %s:%d\n", file_name, __FILE__, __LINE__);
        return;
    }
    int i;
    for(i = 0; i < sorting_timing_count_e; i ++) {
        fprintf(fptr, "%d\titeration=\t%d\tinit_start=\t%f\tlocal_sort_start=\t%f\tshuffle_start=\t%f\tmerge_start=\t%f\tmerge_end=\t%f\tsort_end=\t%f\tdata_sent=\t%lu\tdata_received=\t%lu\tno_send=\t%lu\tno_receive=\t%lu\n",
        i, sorting_timing_e[i].iteration_no, sorting_timing_e[i].init_start, sorting_timing_e[i].local_sort_start, sorting_timing_e[i].shuffle_start, sorting_timing_e[i].merge_start, sorting_timing_e[i].merge_end, sorting_timing_e[i].sort_end, sorting_timing_e[i].data_sent, sorting_timing_e[i].data_received, sorting_timing_e[i].num_sender, sorting_timing_e[i].num_receiver);
    }
    fclose(fptr);

    sprintf(file_name, "sorting/sorting_timing_i%.5d\0", mype);
    fptr = fopen(file_name, "w");
    if(!fptr) {
        fprintf(stderr, "cannot open file %s %s:%d\n", file_name, __FILE__, __LINE__);
        return;
    }
    for(i = 0; i < sorting_timing_count_i; i ++) {
        fprintf(fptr, "%d\titeration=\t%d\tinit_start=\t%f\tlocal_sort_start=\t%f\tshuffle_start=\t%f\tmerge_start=\t%f\tmerge_end=\t%f\tsort_end=\t%f\tdata_sent=\t%lu\tdata_received=\t%lu\tno_send=\t%lu\tno_receive=\t%lu\n",
        i, sorting_timing_i[i].iteration_no, sorting_timing_i[i].init_start, sorting_timing_i[i].local_sort_start, sorting_timing_i[i].shuffle_start, sorting_timing_i[i].merge_start, sorting_timing_i[i].merge_end, sorting_timing_i[i].sort_end, sorting_timing_i[i].data_sent, sorting_timing_i[i].data_received, sorting_timing_i[i].num_sender, sorting_timing_i[i].num_receiver);
    }
    fclose(fptr);

}

