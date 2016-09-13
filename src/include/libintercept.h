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


#define LIBIOINT_PFS_DIR "PFS_DIR"

#define LIBIOINT_PERSIST_DIR "PERSIST_DIR"

typedef struct {
    BBTransferHandle_t handle;
    void *next;
} handle_list_t;

extern char *pfs_dir;
extern char *persist_dir;
extern handle_list_t *head;

REAL_DECL(close, int, (int fd));

int WRAP_DECL(close)(int fd);

extern int Intercept_ExtractFilenames(int fd, char **, char **);
extern void Intercept_StartTransfer(char *, char *);

#endif
