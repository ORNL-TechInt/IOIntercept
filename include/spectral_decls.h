#ifndef __SPECTRAL_DECLS_H__
#define __SPECTRAL_DECLS_H__


#ifdef DEBUG
#define PRINTD(MSG, ...) fprintf(stderr, "DEBUG %s:%d: " MSG "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#else
#define PRINTD(MSG, ...)
#endif

#define PFS_DIR "PFS_DIR"
#define PERSIST_DIR "PERSIST_DIR"
#define MANAGE_FILES "MANAGE_FILES"

#define SUCCESS 0
#define FAILURE -1


#ifdef SPECTRAL_PRELOAD
#include <dlfcn.h>

#define REAL_DECL(func,ret,args) \
    ret (*__real_ ## func)args;

#define WRAP_DECL(__name) __name

#define MAP_OR_FAIL(func) \
    if (!(__real_ ## func)) \
{ \
    __real_ ## func = dlsym(RTLD_NEXT, #func); \
    if(!(__real_ ## func)) { \
        fprintf(stderr, "spectral failed to map symbol: %s\n", #func); \
        exit(1); \
    } \
}

#else

#define REAL_DECL(func,ret,args) \
    extern ret __real_ ## func args;

#define WRAP_DECL(__name) __wrap_ ## __name

#define MAP_OR_FAIL(func)
#endif

enum SPECTRALTRANSFLAGS
{
    SPECEMULATION              = 0x0001, //Emulation based transfer
    SPECBBPROXY                = 0x0002  //BBProxy based transfer
};
typedef enum SPECTRALTRANSFLAGS SPECTRALTRANSFLAGS;

enum SPECTRALFALLBACKFLAGS
{
    NONE                        = 0x0001,
    NEW                         = 0x0002,
    TRANSFERRED                 = 0x0004,
    RSVD                        = 0x0008
};
typedef enum SPECTRALFALLBACKFLAGS SPECTRALFALLBACKFLAGS;




// Handle list struct 
typedef struct {
    SPECTRALFALLBACKFLAGS fallbackstatus;
    BBTransferHandle_t handle;
    BBTransferDef_t *xfer;
    //Emulated transfer handle incase BBProxy fails mid-flight
    BBTransferDef_t *emxfer;
    void *next;
} handle_list_t;

// Function pointers for BB Emulations and BB Proxy
typedef struct {
    /* I may continue to use this transfer def format */
    SPECTRALTRANSFLAGS transflag;
    int (*gettransferhandle) (BBTAG, uint64_t, uint32_t [], BBTransferHandle_t *);
    int (*createtransferdef) (BBTransferDef_t **);
    int (*addfiles) (BBTransferDef_t*, const char *, const char*, BBFILEFLAGS);
    int (*gettransferinfo) (BBTransferHandle_t, BBTransferInfo_t);
    int (*starttransfer) (BBTransferDef_t* , BBTransferHandle_t);
    int (*freetransferdef) (BBTransferDef_t*);
} spectral_functionlist_t;


typedef struct{
    //Global variables
    char *pfs_dir;
    char *persist_dir;
    handle_list_t *handle_list;
    handle_list_t *tail;
    spectral_functionlist_t bbfuncs;
    SPECTRALTRANSFLAGS transfermode;
} spectral_globals_t;

extern spectral_globals_t globals;

#endif
