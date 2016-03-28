#ifndef __BBAPI_H__
#define __BBAPI_H__
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <sched.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>

//New and Fun Data Types
typedef struct {
    char *src;
    char *dest;
    struct transfer_list_t *next;
} transfer_list_t;


typedef struct {
    transfer_list_t *trans_list;
}BBTransferDef_t;

typedef struct {
}BBTransferHandle_t;

typedef struct {
    BBTransferDef_t *work_item;
    BBTransferHandle_t *work_handle;
    //Internal use for the emulate -- used to terminate the emulation layer 
    //Since this layer doesn't make it to production just leave it be. :) 
    bool term_flag;
    struct _internal_work_queue_t *next;
} _internal_work_queue_t;

//Pthread spawned from example code
void *bb_proxy_func(void *arg);
void spawn_bb_proxy();
void term_bb_proxy();
void bb_bind_cpu();

//Emulation of BBAPI 
int BB_CreateTransferDef(BBTransferDef_t **);
int BB_AddFiles(BBTransferDef_t *, char *, char *, int);
int BB_StartTransfer(uint64_t, uint32_t, uint32_t*, BBTransferDef_t *, BBTransferHandle_t *);
int BB_FreeTransferDef(BBTransferDef_t *);

//Internal functions
void transfer_file(char *filename, char *dest);
void move_file(char *src, char *dest);
int enqueue_work(uint64_t tag, uint32_t num_contrib, uint32_t *contrib, BBTransferDef_t *xfer, BBTransferHandle_t *handle, int term_flag);
_internal_work_queue_t* dequeue_work();

//Managment data
extern pthread_t *bb_proxy_thread;
extern pthread_cond_t *cond;
extern pthread_mutex_t *mutex;
extern pthread_mutex_t *wqtex;
extern _internal_work_queue_t *wqhead;
extern _internal_work_queue_t *wqtail;
#endif
