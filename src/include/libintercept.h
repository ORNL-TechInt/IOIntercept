#ifndef __LIBINTERCEPT_H__
#define __LIBINTERCEPT_H__
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>


#ifdef TITAN
#include <bbapi.h>
#else
#include <bb/include/bbapi.h>
#endif

#ifdef DEBUG
#define PRINTD(MSG, ...) fprintf(stderr, "DEBUG %s:%d: " MSG "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#else
#define PRINTD(MSG, ...)
#endif

#define PFS_DIR "PFS_DIR"
#define PERSIST_DIR "PERSIST_DIR"
#define MANAGE_FILES "MANAGE_FILES"


#ifdef LIBINTERCEPT_PRELOAD
#define __USE_GNU
#include <dlfcn.h>

#define REAL_DECL(func,ret,args) \
    ret (*__real_ ## func)args = NULL;

#define WRAP_DECL(__name) __name

#define MAP_OR_FAIL(func) \
    if (!(__real_ ## func)) \
    { \
       __real_ ## func = dlsym(RTLD_NEXT, #func); \
       if(!(__real_ ## func)) { \
          fprintf(stderr, "lib_io_intercept failed to map symbol: %s\n", #func); \
          exit(1); \
       } \
    }

#else

#define REAL_DECL(func,ret,args) \
    extern ret __real_ ## name args;

#define WRAP_DECL(__name) __wrap_ ## __name

#define MAP_OR_FAIL(func)
#endif

// Handle list struct 
typedef struct {
    BBTransferHandle_t handle;
    BBTransferDef_t *xfer;
    void *next;
} handle_list_t;

//Global variables
char *pfs_dir;
char *persist_dir;
handle_list_t *handle_list;
handle_list_t *tail;

//Function Prototypes
REAL_DECL(close, int, (int fd));
int WRAP_DECL(close)(int fd);
char *readlink_malloc (const char *filename);
void Intercept_ManageFiles();
int Intercept_ExtractFilenames(int fd, char **, char **);
void Intercept_StartTransfer(char *, char *);
int Intercept_StoreHandle(BBTransferHandle_t, BBTransferDef_t *);
int Intercept_RemoveHandle(BBTransferHandle_t);
#endif
