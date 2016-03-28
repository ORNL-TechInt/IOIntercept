#include "libintercept.h"
#include "bbapi.h"

char *readlink_malloc (const char *filename);

char *pfs_dir = NULL;
bool early_term = false;



int WRAP_DECL(close)(int fd){

    char *dest = NULL;

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
            printf("Please set LIBIOINT_ENV_VAR. You will need to manually extract your files upon application termination\n");
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
        if ((filename = readlink_malloc((const char *)path)) != NULL){
            //Extracted Filename now check that path aligns with our persist directory
            /* Your path must write to the persist dir directly
             * right now won't work if it doesn't start with /ssd/persist/ it won't
             * get written */
            if (!strncmp(LIBIOINT_PERSIST_DIR,filename,strlen(LIBIOINT_PERSIST_DIR)))
            {
                //Grab portion that's not the ssd persist stuff
                int strdelta = strlen(filename) - strlen(LIBIOINT_PERSIST_DIR);
                char *ptr = filename;
                ptr = ptr + strlen(LIBIOINT_PERSIST_DIR);

                //Gotta account for \0
                //There better be a slash on the pfs_dir
                int destsize = strlen(pfs_dir) + strlen(ptr)+1;
                dest = (char *)malloc(destsize * sizeof(char));
                strncpy(dest,pfs_dir,strlen(pfs_dir) + 1);
                strncat(dest,ptr,strdelta);
            }else{
                early_term = true;
            }
        }
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

    if (!early_term){
        //OK Real File is close now.
        BBTransferDef_t *xfer = NULL;

        if (dest == NULL){
            perror("Invalid pfs destination\n");
        }

        if (BB_CreateTransferDef(&xfer) < 0 || xfer == NULL){
            //Not checking errno for prototype -- won't until I get a list of errnos from bbproxy
            perror("Failed to create file-transfer definition\n");
        }

        if (BB_AddFiles(xfer, filename, dest, 0) < 0){
            perror("Failed to add files\n");
        }

        /*
         * Contrib passed as null from email line In the latest version of the
         * header (which you don't quite have yet), when numcontrib=1, the
         * contributor list is ignored. I think that should satisfy Scott's
         * file close scenario.  for now with tag, I'm just going to use FD but
         * we probably do need implement a static counter???  
        */

        if (BB_StartTransfer(fd, 1, NULL, xfer, NULL) < 0){
            perror("Failed to start transfer\n");
        }
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
