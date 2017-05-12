#include "spectral.h"

int EM_InitLibrary(uint32_t contribID, const char *clientVersion){
    //*Stub do nothing*//
    return 0;
}

int EM_TerminateLibrary(){
    //*Stub do nothing*//
    return 0;
}

int EM_GetTransferHandle(BBTAG tag, uint64_t numcontrib, uint32_t contrib[],
        BBTransferHandle_t* handle){
    return 0;
}


int EM_CreateTransferDef(BBTransferDef_t **xfer){
    *xfer = (BBTransferDef_t*)malloc(sizeof(transfer_list_t));
    return 0;
}

int EM_AddFiles(BBTransferDef_t *xfer, const char *src, const char *dest, BBFILEFLAGS flags){
    ((transfer_list_t *)xfer)->src = (char *)malloc(strlen(src) + 1);
    ((transfer_list_t *)xfer)->dest = (char *)malloc(strlen(src) + 1);
    strcpy(((transfer_list_t *)xfer)->src,src);
    strcpy(((transfer_list_t *)xfer)->dest,dest);
    return 0;
}

int EM_StartTransfer(BBTransferDef_t *xfer, BBTransferHandle_t handle){
    return 0;
}


int EM_FreeTransferDef(BBTransferDef_t *xfer){
    return 0;
}
