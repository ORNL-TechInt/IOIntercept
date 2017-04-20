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
    ret (*__real_ ## func)args;

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
    extern ret __real_ ## func args;

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
extern char *pfs_dir;
extern char *persist_dir;
extern handle_list_t *handle_list;
extern handle_list_t *tail;


/* Mapper Macro for close
 *
 * Statically linked: This will point to __real_close() and can 
 * should only be called from __wrap_close()
 *
 * Dynamically linked: This will become a predeclared function pointer
 * that is set in MAP_OR_FAIL. 
 *
 */
REAL_DECL(close, int, (int fd))

/* Wrapper Function for Close
 *
 * \param[in] fd: The file descriptor that is passed to the traditional 
 * close call.
 *
 * \return 0 Success
 * \return -1 Failure
 *
 */
int WRAP_DECL(close)(int fd);

/* Mapper Macro for fclose
 *
 * Statically linked: This will point to __real_fclose() and can 
 * should only be called from __wrap_fclose()
 *
 * Dynamically linked: This will become a predeclared function pointer
 * that is set in MAP_OR_FAIL. 
 *
 */
REAL_DECL(fclose, int, (FILE *fp))



/* Wrapper Function for Close
 *
 * \param[in] fp: The file pointer that is passed to the traditional 
 * fclose call.
 *
 * \return 0 Success
 * \return -1 Failure
 *
 */
int WRAP_DECL(fclose)(FILE *fp);

/* Wrapper Agnostic Close Call
 *
 * \param[in] streaming: indicating if the source call was from
 * fclose or close, dictates certain real callbacks
 *
 * \param[in] fd: file descriptor used to extract the file location
 * from proc
 *
 * \param[in] fp: originally passed in file pointer, only used situationally
 * pass to fclose incase there is additional cleanup surrounding a file-pointer
 *
 * \return 0 Success
 * \return - Failure
 *
 * Most cases don't define errno (YET)
 *
 * This is the primary state driver for handling file persistence.
 * If the environment has not been started or read in this function
 * will initialize BBProxy, if BBProxy init fails, the handler will
 * fall back to libemulate a pthread simulation that will synchronously
 * copy the files out. 
 *
 */
int Intercept_HandleClose(bool streaming, int fd, FILE *fp);


/* File Management 
 *
 * No input parameters
 * No output parameters (TODO)
 *
 * Scans handle list for successful transfers up but not including
 * the most recently enqueued transfer. This cleans up the bbproxy state
 * and unlinks the file avoiding capacity overruns
 *
 * This is only called if the MANAGE_FILES environment variable is set 
 * prior to the calls 
 */
void Intercept_ManageFiles();

/* Creates Source and Destination File locations
 *
 * \param[in] fd: File descriptor from which to extract the source
 * filename
 *
 * \param[out] src: Extracted source filename
 *
 * \param[out] dest: Extracted and created destiation filename
 *
 * \return 1 Success TODO Clean up this for consistency
 * \return -1 Failure
 *
 */
int Intercept_ExtractFilenames(int fd, char **, char **);

/* Starts transfer to parallel filesystem
 *
 * \param[in] source: source file to copy
 *
 * \param[in] destination: destination of transfer
 *
 * \return 0 Success
 * \return BBProxy based failure code
 *
 * This function initiates the actual transfer creating a transfer handle
 * getting transfer info memory allocated
 * creating the transfer definition
 * adding files to the transfer definition
 * and starting the transfer.
 *
 * Any remote bbProxy function call can fail, due to this
 * we can and check rc on all of these calls. If rc indicates failure
 * we immediatly return passing rc to the close handler and determining
 * if we can fall back to the emulation layer.
 *
 */
void Intercept_StartTransfer(char *, char *);


/* Enqueue transfer handle and transfer definition
 *
 * \param[in] handle: Transfer handle 
 *
 * \param[in] xfer: Transfer definition
 *
 * \return 1 Success 
 *
 * \return -1 Failure
 *
 * TODO Find consistency in success codes
 *
 */
int Intercept_StoreHandle(BBTransferHandle_t, BBTransferDef_t *);

/* Dequeue transfer handle and transfer definition
 *
 * \param[in] handle: Transfer handle to remove 
 *
 * \return 0 Success 
 *
 * \return -1 Failure
 *
 * TODO Find consistency in success codes
 *
 */
int Intercept_RemoveHandle(BBTransferHandle_t);



/* Read link malloc
 *
 * \param[in] filename: Symbolic link at /proc/self/fd/<passed in fd number>
 *
 * \return filename of the file pointed to by input link SUCCESS
 * \return NULL when readlink fails
 *
 * This function does not know the lenght of the filename thats pointed to.
 * We read the symbolic link 100 bytes at a time and grow the buffer if we
 * cannot accomodate the string.
 *
 * Once the end is found a string terminator is appended and the final string
 * is returned
 *
 */
char *readlink_malloc (const char *filename);
#endif
