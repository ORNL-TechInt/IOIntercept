#include "libintercept.h"


char *readlink_malloc (const char *filename);

char *pfs_dir = NULL;
bool early_term = false;



int WRAP_DECL(close)(int fd){
    PRINTD("Entering Close Wrapper\n");
    
    if (!early_term && !pfs_dir){
        pfs_dir = getenv(LIBIOINT_ENV_VAR);

        /*
         * On Titan there isn't really an ENV that points to a 
         * safe place to write.
         * There is MEMBERWORK but you have to traverse their
         * ACCT to find a writeable location.
         *
         * For now if they don't set it - Just don't do anything
         *
         */
        if (!pfs_dir){
            printf("No PFS Writeout set. You will need to manually extract your files upon application termination\n");
            early_term = true;
        }
    }


    /* Extract filename from FD for generating proxy updates */
    char *filename = NULL;
    //TODO If PROCFS ever disapears we'll have to change this.
    if (!access("/proc/self/fd",X_OK)) {
        /* We know path can't be that big */
        char path[255];
        sprintf(path,"/proc/self/fd/%i",fd);
        filename = readlink_malloc((const char *)path);
        printf("%s---%s\n",path,filename);
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
    int ret = __real_close(fd);

    if (!early_term){
    }

    return ret;
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
