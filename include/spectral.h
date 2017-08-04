#ifndef __SPECTRAL_H__
#define __SPECTRAL_H__

#include <bb/include/bbapi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include "emapi.h"

#define __USE_GNU
#include "spectral_decls.h"

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

/* Sets up transfer to PFS
 *
 * \param[in] fd: file descriptor used to extract the file location
 * from proc
 *
 * \param[out] src: Source file extracted from Fd
 *
 * \param[out] dest: Destination file to be created by transfer agent
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
int spectral_setuptransfer(int, char **, char **);


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
void spectral_managefiles();

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
int spectral_extractfilenames(int fd, char **, char **);

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
void spectral_starttransfer(char *, char *);


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
int spectral_storehandle(BBTransferHandle_t, BBTransferDef_t *, BBTransferDef_t *);

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
int intercept_removehandle(BBTransferHandle_t);



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
