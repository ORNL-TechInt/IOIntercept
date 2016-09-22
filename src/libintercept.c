#include "libintercept.h"


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

    static bool manage_files = false;
    char *dest = NULL;
    char *src = NULL; 
    char *mf_env = NULL;
    bool early_term = false;
    int rc;

    if (!early_term && !pfs_dir)
    {
        pfs_dir = getenv(PFS_DIR);
        persist_dir = getenv(PERSIST_DIR);
	mf_env = getenv(MANAGE_FILES);

        if (mf_env && strcmp(mf_env,"1") == 0){
            manage_files =  true;
        }

        if (!pfs_dir || !persist_dir)
        {
            printf("Please set PFS_DIR and PERSIST_DIR. You will need to manually extract your files upon application termination\n");
            early_term = true;
        }
        else
        {
            rc = BB_InitLibrary(getpid(),BBAPI_CLIENTVERSIONSTR);
            check(rc);
        }
    }


    if (Intercept_ExtractFilenames(fd, &src, &dest) < 0)
    {
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
    if ((ret = __real_close(fd)) != 0)
    {
        fprintf(stderr,"%s\n",strerror(errno));
        return ret;
    }

    /*
     * Real file has been closed 
     * Clean up previously transmitted files
     */
    if (manage_files)
    {
        Intercept_ManageFiles();
    }


    /*
     * If everything else worked cleanly
     * invoke the BBAPI to transfer the files
     * store the handle
     */
    if (!early_term)
    {
        Intercept_StartTransfer(src, dest);
    } 

    if (dest){
        free(dest);
    }

    if (src){
        free(src);
    }

    return ret;
}

void Intercept_ManageFiles(){
    handle_list_t *itr;
    BBTransferInfo_t info;

    itr = handle_list;
    while (itr)
    {
        BB_GetTransferInfo(itr->handle, &info);
        if (info.status == BBFULLSUCCESS){
            handle_list_t *tmp = itr;
            itr=itr->next;
            BB_FreeTransferDef(tmp->xfer);
            Intercept_RemoveHandle(tmp->handle);
        }
    }
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
                return 0;
            }	
        }
    }
    return 1;
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
    static uint32_t tagctr = 0;
    static uint32_t pid = 0;
    int rc;
    BBTransferDef_t *xfer = NULL;
    BBTransferHandle_t handle;
    BBTransferInfo_t info;

    if (pid == 0){
        pid = getpid();
    }

    int tag = (pid << 16) | (tagctr++);

    rc = BB_GetTransferHandle(tag, 0, NULL, &handle);
    check(rc);

    rc = BB_GetTransferInfo(handle, &info);
    check(rc);

    rc = BB_CreateTransferDef(&xfer);
    check(rc);        

    rc = BB_AddFiles(xfer, src , dest, 0);
    check(rc);

    rc = BB_StartTransfer(xfer, handle);
    check(rc);

    rc = Intercept_StoreHandle(handle,xfer);
}

int Intercept_StoreHandle(BBTransferHandle_t handle, BBTransferDef_t *xfer){
    /* Function invariant tail will always enter and exit being NULL */
    handle_list_t *tmp =  (handle_list_t*)malloc(sizeof(handle_list_t));
    if (!tmp){
        return -1;
    }
    tmp->handle = handle;
    tmp->xfer = xfer;
    tmp->next = NULL;

    /* if tail is null so is head */
    if (tail == NULL){
        tail = tmp;
        handle_list = tail;
    }else
    {
        tail->next = tmp;
        tail = tail->next;
    }

    return 1;
}

int Intercept_RemoveHandle(BBTransferHandle_t handle){
    /* Sanity */
    if (handle_list && handle_list->handle == handle){
        /* We only remove the head of this list */
        handle_list_t *tmp = handle_list;
        if (handle_list == tail){
            /* We're removing the last element in the list */
            handle_list = tail = NULL;
        }else{
            handle_list = handle_list->next;
        }
        free(tmp);
        return 0;
    }    
    return -1;
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
