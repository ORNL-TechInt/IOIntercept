#ifndef __LIBINTERCEPT_H__
#define __LIBINTERCEPT_H__

#ifdef LIBINTERCEPT_PRELOAD
#define __USE_GNU
#include <dlfcn.h>
#include <stdlib.h>

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

REAL_DECL(close, int, (int fd));

int WRAP_DECL(close)(int fd);

#endif
