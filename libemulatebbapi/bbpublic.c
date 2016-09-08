#include "bbapi.h"

int BB_InitLibrary(uint32_t contribID, const char *clientVersion){
	//*Stub do nothing*//
	return 0;
}

int BB_TerminateLibrary(){
	//*Stub do nothing*//
}

int BB_GetTransferHandle(BBTag tag, uint64_t numcontrib, uint32_t contrib[], BBTransferHandle_t* handle){
	handle->tag = tag;
	handle->numcontrib = numcontrib;
	handle->contrib = contrib;
	return 0;
}


int BB_CreateTransferDef(BBTransferDef_t **xfer){
    *xfer = (BBTransferDef_t*)malloc(sizeof(BBTransferDef_t));

    if (*xfer == NULL){
        return -1;
    }

    (*xfer)->trans_list = NULL;
    
    return 0;
}

int BB_AddFiles(BBTransferDef_t *xfer, char *src, char *dest, int flags){
    if (src == NULL || dest == NULL){
        printf("Invalid arguments specified\n");
        return -1;
    }

    if (xfer->trans_list == NULL){
        //Create First Element
        xfer->trans_list = (transfer_list_t*)malloc(sizeof(transfer_list_t));
        xfer->trans_list->src = (char *)malloc(strlen(src)+1);
        xfer->trans_list->dest = (char *)malloc(strlen(dest)+1);
        xfer->trans_list->next = NULL;
    }

    strcpy(xfer->trans_list->src,src);
    strcpy(xfer->trans_list->dest,dest);

    return 0;
}

int BB_StartTransfer(BBTransferDef_t *xfer, BBTransferHandle_t handle){
    /* Does only one rank execute start transfer? Need to clarify this */
    enqueue_work(handle.tag,handle.numcontrib, handle.contrib, xfer, &handle, false);
    return 0;
}


int BB_FreeTransferDef(BBTransferDef_t *xfer){

    return 0;
}
