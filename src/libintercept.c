#include "libintercept.h"

char *readlink_malloc (const char *filename);

char *pfs_dir = NULL;
char *persist_dir = NULL;
uint32_t tagctr = 1;
int check(int rc)
{
    if(rc)
    {
	char* errstring;
	BB_GetLastErrorDetails(BBERRORJSON, &errstring);
	printf("Error rc:       %d\n", rc);
	printf("Error details:  %s\n", errstring, errstring);
	free(errstring);

	printf("Aborting due to failures\n");
	exit(-1);
    }
}


int WRAP_DECL(close)(int fd){

    char *dest = NULL;
    char *src = NULL;
    bool early_term = false;
    int rc;

    if (!early_term && !pfs_dir){
	pfs_dir = getenv(LIBIOINT_PFS_DIR);
	persist_dir = getenv(LIBIOINT_PERSIST_DIR);

	/*
	 * On Titan there isn't really an ENV that points to a 
	 * safe place to write.
	 * There is MEMBERWORK but you have to traverse their
	 * ACCT to find a writeable location.
	 *
	 * For now if they don't set it - Just don't do anything
	 *
	 */
	if (!pfs_dir || !persist_dir){
	    printf("Please set PFS_DIR and PERSIST_DIR. You will need to manually extract your files upon application termination\n");
	    early_term = true;
	}else{
	    rc = BB_InitLibrary(getpid(),BBAPI_CLIENTVERSIONSTR);
	    check(rc);
	}
    }


    if (Intercept_ExtractFilenames(fd, &src, &dest) < 0){
	    early_term = true;
    }


    /*
     * If linked dynamically there is no __real_close -- 
     * We will create one and map it to the real close.
     *
     * Statically linked binaries have the __real_close if 
     * linked with option -Wl,-wrap,close
     */
    MAP_OR_FAIL(close);

    /*
     * Execute the real close
     */
    int ret = 0;
    if ((ret = __real_close(fd)) != 0){
	fprintf(stderr,"%s\n",strerror(errno));
	return ret;
    }

    /*
     * If everything else worked cleanly
     * invoke the BBAPI to transfer the files
     * store the handle
     */
    if (!early_term){
	    Intercept_StartTransfer(src, dest);
    }    

    return ret;
}


int Intercept_ExtractFilenames(int fd, char **src, char **dest){
	//TODO If PROCFS ever disapears we'll have to change this.
	if (!access("/proc/self/fd",X_OK)) {
		/* We know path can't be that big */
		char path[255];
		sprintf(path,"/proc/self/fd/%i",fd);
		if (((*src) = readlink_malloc((const char *)path)) != NULL){
			//Extracted Filename now check that path aligns with our persist directory
			/* Your path must write to the persist dir directly
			 * right now won't work if it doesn't start with /ssd/persist/ it won't
			 * get written */
			if (strncmp(persist_dir,*src,strlen(persist_dir)) == 0)
			{
				//Grab portion that's not the ssd persist stuff
				int strdelta = strlen(*src) - strlen(persist_dir);
				char *ptr = *src;
				ptr = ptr + strlen(persist_dir);

				//Gotta account for \0
				//There better be a slash on the pfs_dir
				int destsize = strlen(pfs_dir) + strlen(ptr)+1;
				*dest = (char *)malloc(destsize * sizeof(char));
				strncpy(*dest,pfs_dir,strlen(pfs_dir) + 1);
				strncat(*dest,ptr,strdelta);
			}else{
				return -1;
			}	
		}
	}
	return 0;
}

/*
 * Use BBAPI to transfer the files
 * Failure of any of the bbapi commands cause the application
 * to terminate. In the check(rc) function
 * 
 * No error codes are returned if this returns
 * assume success
 */
void Intercept_StartTransfer(char *src, char *dest){
	int rc;
	BBTransferDef_t *xfer = NULL;
	BBTransferHandle_t handle;

	rc = BB_GetTransferHandle(getpid(), 0, NULL, &handle);
	check(rc);

	rc = BB_CreateTransferDef(&xfer);
	check(rc);        

	rc = BB_AddFiles(xfer, src , dest, 0);
	check(rc);

	rc = BB_StartTransfer(xfer, handle);
	check(rc);
}


char *readlink_malloc (const char *filename)
{
	int size = 100;
	char *buffer = NULL;

	while (1)
	{
		buffer = (char *) realloc (buffer, size);
		int nchars = readlink (filename, buffer, size);
		if (nchars < 0)
		{
			free (buffer);
			return NULL;
		}
		if (nchars < size)
		{
			/* Readlink doesn't null terminate the string it sends back.
			   But a nice property of this is that if nchars == size you can't
			   tell if you read the whole thing so size ALWAYS has to be greater
			   this means you have a slot available to stash your null terminator.

			   Hooray.             
			   */
			buffer[nchars] = '\0';               
			return buffer;
		}
		size *= 2;
	}
}
