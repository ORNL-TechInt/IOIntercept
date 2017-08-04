#ifndef __EMAPI_H__
#define __EMAPI_H__

//New and Fun Data Types
typedef struct {
    char *src;
    char *dest;
} transfer_list_t;

//Pthread spawned from example code
void *bb_proxy_func(void *arg);
void spawn_bb_proxy();
void term_bb_proxy();
void spawn_bb_proxy_();
void term_bb_proxy_();
void bb_bind_cpu();

//Emulation of BBAPI
int EM_GetTransferHandle(BBTAG, uint64_t, uint32_t [], BBTransferHandle_t *);
int EM_CreateTransferDef(BBTransferDef_t **);
int EM_AddFiles(BBTransferDef_t *, const char *, const char *, BBFILEFLAGS);
int EM_StartTransfer(BBTransferDef_t *, BBTransferHandle_t);
int EM_FreeTransferDef(BBTransferDef_t *);

//Internal functions
void transfer_file(char *filename, char *dest);
void move_file(char *src, char *dest);
int enqueue_work(uint64_t tag, uint32_t num_contrib, uint32_t *contrib, 
        BBTransferDef_t *xfer, BBTransferHandle_t *handle, int term_flag);
void signal_proxy_thread();

//Managment data
extern pthread_t *bb_proxy_thread;
extern pthread_cond_t *cond;
extern pthread_mutex_t *mutex;
extern pthread_mutex_t *wqtex;
#endif
