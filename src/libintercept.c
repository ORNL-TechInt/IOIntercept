#include <stdio.h>
#include "libintercept.h"

int WRAP_DECL(close)(int fd){
    MAP_OR_FAIL(close);
    printf("Got here\n");
    return __real_close(fd);
}
