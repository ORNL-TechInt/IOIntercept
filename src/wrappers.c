#include "spectral.h"

/* Wrappers */
int WRAP_DECL(fclose)(FILE *fp){
    bool cancel = false;
    int fd = fileno(fp);
    int ret = FAILURE;
    char *src = NULL;
    char *dest = NULL;

    if (fd < 0){
        fprintf(stderr, "Invalid file pointer\n");
        return -1;
    }

    if (spectral_setuptransfer(fd, &src, &dest) == FAILURE)
        cancel = true;
   
    
    MAP_OR_FAIL(fclose);
    if ((ret = __real_fclose(fp)) != 0)
    {
        fprintf(stderr,"%s\n",strerror(errno));
        return ret;
    }

    if (!cancel)
        spectral_starttransfer(src,dest);

    return ret;
}

int WRAP_DECL(close)(int fd){
    int ret = FAILURE;
    char *src = NULL;
    char *dest = NULL;
    bool cancel = false;

    if (spectral_setuptransfer(fd, &src, &dest) == FAILURE)
        cancel = true;

    MAP_OR_FAIL(close);
    if ((ret = __real_close(fd)) != 0)
    {
        fprintf(stderr,"%s\n",strerror(errno));
        return ret;
    }

    if (!cancel)
        spectral_starttransfer(src,dest);
    
    return ret;
}

