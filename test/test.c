#include <bbapi.h>
#include <stdio.h>
#include <stdint.h>

char src[255];
char dest[255];

int test_bbapi(){
    BBTransferDef_t *xfer = NULL;
    BBTransferDef_t *emxfer = NULL;
    BBTransferHandle_t handle = 0;
    BBTransferInfo_t info;

    BB_InitLibrary(54321, BBAPI_CLIENTVERSIONSTR);

    BB_GetTransferHandle(getpid(), 0, NULL, &handle);

    BB_CreateTransferDef(&xfer);

    BB_AddFiles(xfer, src , dest, 0);

    BB_StartTransfer(xfer, handle);

    //BB_TerminateLibrary();
}


int main(int argc, char **argv)
{
    strcpy(src, argv[1]);
    strcpy(dest, argv[2]);
    test_bbapi();
    return 0;
}
