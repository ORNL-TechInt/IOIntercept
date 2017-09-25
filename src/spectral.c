/*
 *  * Copyright (c) 2017 UT-Battelle, LLC.  All rights reserved.
 *  * Copyright (c) 2017 Oak Ridge National Labs.  All rights reserved.
 *
 */
#include "spectral.h"

/* Global Variables */
spectral_globals_t globals = {false,false, NULL, NULL, NULL, NULL};

/* Init */
int spectral_init(){
    static uint32_t pid = 0;
    int rc = FAILURE;

    //Invariant - either bbProxy or Pthread Model will start
    globals.initialized = true;

    wqtex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
    if (pthread_mutex_init(wqtex,NULL) != 0){
        perror(strerror(errno));
        exit(-1);
    }

    if (pid == 0){
        pid = getpid();
        rc = BB_InitLibrary(pid,BBAPI_CLIENTVERSIONSTR);
        if (check(rc) == SUCCESS){
            globals.transfermode = SPECBBPROXY;
            globals.bbfuncs.gettransferhandle = BB_GetTransferHandle;
            globals.bbfuncs.createtransferdef = BB_CreateTransferDef;
            globals.bbfuncs.addfiles = BB_AddFiles;
            globals.bbfuncs.starttransfer = BB_StartTransfer;
        }else{
            spawn_bb_proxy();
            globals.transfermode = SPECEMULATION;
            globals.bbfuncs.gettransferhandle = EM_GetTransferHandle;
            globals.bbfuncs.createtransferdef = NULL;
            globals.bbfuncs.addfiles = NULL;
            globals.bbfuncs.starttransfer = EM_StartTransfer;
        }
    }

    if (!globals.pfs_dir || !globals.persist_dir)
    {
        printf("Please set PFS_DIR and PERSIST_DIR. You will need to manually\
                extract your files upon application termination\n");
        return FAILURE;
    }

    return SUCCESS;
}

/*  BBProxy File Handlers */
int spectral_setuptransfer(int fd, char **src, char **dest)
{
    char *mf_env = NULL;
    int rc, ret = SUCCESS;

    if (globals.pfs_dir == NULL){
    globals.manage_files = getenv(MANAGE_FILES);
    globals.pfs_dir = getenv(PFS_DIR);
    globals.persist_dir = getenv(PERSIST_DIR);
    }

    /* Extract filenames */
    if (spectral_extractfilenames(fd, src, dest) == FAILURE)
    {
        return FAILURE;
    }

    /* Init */
    if (globals.initialized == false && globals.transfermode == SPECUNINIT)
    {
        rc = spectral_init();
        /* Unable to transmit just close the file and return */
        if (rc == FAILURE){
            return FAILURE;
        }
    }

    return ret;
}

void spectral_starttransfer(char *src, char *dest){
    static uint32_t tagctr = 0;
    int rc;
    int pid = getpid();
    BBTransferDef_t *xfer = NULL;
    BBTransferDef_t *emxfer = NULL;
    BBTransferHandle_t handle = 0;
    BBTransferInfo_t info;

    int tag = (pid << 16) | (tagctr++);

    //Setup Setup Transfer
    rc = globals.bbfuncs.gettransferhandle(tag, 0, NULL, &handle);
    check(rc);

    rc = globals.bbfuncs.createtransferdef(&xfer);
    check(rc);        

    rc = globals.bbfuncs.addfiles(xfer, src , dest, 0);
    check(rc);

    rc = globals.bbfuncs.starttransfer(xfer, handle);
    check(rc);

    //Store copy of the handler for emulation if it's not available.
    EM_CreateTransferDef(&emxfer);
    EM_AddFiles(emxfer,src,dest,0);

    spectral_storehandle(handle,xfer,emxfer);

    /*
     * Real file has been closed 
     * Clean up previously transmitted files
     */
    if (globals.manage_files)
    {
        spectral_managefiles();
    }

    if (dest){
        free(dest);
    }

    if (src){
        free(src);
    }
}


void spectral_managefiles(){
    handle_list_t *itr;
    BBTransferInfo_t info;

    itr = globals.handle_list;
    while (itr)
    {
        bool remove = false;

        /* Base Test */
        if (itr->fallbackstatus == NONE){
            BB_GetTransferInfo(itr->handle, &info);
            if (info.status == BBFULLSUCCESS)
                remove = true;
        }
        else { /* Fall back check */
            if (itr->fallbackstatus == TRANSFERRED)
                remove = true;
        }

        if (remove){
            handle_list_t *tmp = itr;
            BB_FreeTransferDef(tmp->xfer);
            spectral_removehandle(tmp->handle);
        }
        itr=itr->next;
    }
}

int spectral_extractfilenames(int fd, char **src, char **dest){
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
            if (strncmp(globals.persist_dir,*src,strlen(globals.persist_dir)) == 0)
            {
                //Grab portion that's not the ssd persist stuff
                int strdelta = strlen(*src) - strlen(globals.persist_dir);
                char *ptr = *src;
                ptr = ptr + strlen(globals.persist_dir);

                //Gotta account for \0
                //There better be a slash on the pfs_dir
                int destsize = strlen(globals.pfs_dir) + strlen(ptr)+1;
                *dest = (char *)malloc(destsize * sizeof(char));
                strncpy(*dest,globals.pfs_dir,strlen(globals.pfs_dir) + 1);
                strncat(*dest,ptr,strdelta);
                return SUCCESS;
            }
        }
    }
    return FAILURE;
}


int spectral_storehandle(BBTransferHandle_t handle, BBTransferDef_t *xfer, BBTransferDef_t *emxfer){
    /* Function invariant tail will always enter and exit being NULL */
    bool signal_sleeper = false;
    handle_list_t *tmp =  (handle_list_t*)malloc(sizeof(handle_list_t));
    if (!tmp){
        return FAILURE;
    }
    tmp->handle = handle;
    tmp->xfer = xfer;
    tmp->emxfer = emxfer;

    if (globals.transfermode == SPECEMULATION){
        tmp->fallbackstatus = NEW;
        signal_sleeper = true;
    }else{
        tmp->fallbackstatus = NONE;       
    }

    tmp->next = NULL;

    pthread_mutex_lock(wqtex);
    /* if tail is null so is head */
    if (globals.tail == NULL){
        globals.tail = tmp;
        globals.handle_list = globals.tail;
    }else
    {
        globals.tail->next = tmp;
        globals.tail = globals.tail->next;
    }
    pthread_mutex_unlock(wqtex);


    if (signal_sleeper)
        pthread_cond_signal(cond);

    return 1;
}

int spectral_removehandle(BBTransferHandle_t handle){
    /* Sanity */
    pthread_mutex_lock(wqtex);
    if (globals.handle_list && globals.handle_list->handle == handle){
        /* We only remove the head of this list */
        handle_list_t *tmp = globals.handle_list;
        if (globals.handle_list == globals.tail){
            /* We're removing the last element in the list */
            globals.handle_list = globals.tail = NULL;
        }else{
            globals.handle_list = globals.handle_list->next;
        }        
        free(tmp);
        pthread_mutex_unlock(wqtex);
        return 0;
    }

    pthread_mutex_unlock(wqtex);
    return FAILURE;
}
