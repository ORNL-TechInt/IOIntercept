#include <mpi.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <stdbool.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>
#include <openssl/md5.h>

#ifdef TITAN
#include "bbapi.h"
#endif



char *file_md5(char *filename){
    FILE *fptr = fopen(filename,"rb");
    char *temp = (char *)malloc(sizeof(char) * MD5_DIGEST_LENGTH);
    MD5_CTX md5ctx;
    unsigned char data[1024];
    int bytes = 0;

    if (fptr != NULL){
        MD5_Init(&md5ctx);
        while ((bytes = fread(data, 1, 1024, fptr)) != 0){
            MD5_Update(&md5ctx,data,bytes);
        }
        MD5_Final(temp, &md5ctx);
    }
    return temp;
}

/* Persist File
 * Intercept and drain
 * Compare src and dst md5 values
 * 
 * Input parameter file size in KB
 * Writes are performed in 4K chunks
 * size is rounded down
 */
bool simple_file_test(uint32_t filesize){
    char *srcdir = getenv("PERSIST_DIR");
    char *destdir = getenv("PFS_DIR");
    char tfn[256];
    char dfn[256];
    char cmd[256];
    /* 4K contiguous buffer in bytes */
    int buffer4k[1024];
    int testfile = 0;

    /* Determine looping counts*/
    int loop_count = (filesize * 1024) / 4096;

    /* Seed random num gen */
    srand(time(NULL));

    snprintf(tfn,sizeof(tfn), "%s/testfile.%d.%d", srcdir,getpid(),filesize);
    snprintf(dfn,sizeof(tfn), "%s/testfile.%d.%d", destdir,getpid(),filesize);

    if ((testfile = open(tfn, O_WRONLY | O_CREAT | O_TRUNC,
                    S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) == -1)
    {
        perror("Cannot open output file\n"); exit(1);
    }

    for (int lcv = 0;lcv < loop_count;lcv++){
        for (int lcv2 = 0; lcv2 < 1024;lcv2++){
            buffer4k[lcv2] = rand();
        }
        write(testfile, (void *)buffer4k, 4096);
    }
    close(testfile);

    sleep(3);

    char *srcmd5 = file_md5(tfn);
    char *destmd5 = file_md5(dfn);

    for (int lcv = 0;lcv < MD5_DIGEST_LENGTH; lcv++){
        if (srcmd5[lcv] != destmd5[lcv])
        {
            return false;
        }
    }

    if (srcmd5)
        free(srcmd5);

    if (destmd5)
        free(destmd5);


    return true;
}


int main(int argc, char **argv){
#ifdef TITAN
    spawn_bb_proxy();	    
#endif
    MPI_Init(&argc, & argv);

    //128K
    if (!simple_file_test(128)){
        perror("Failed simple_file_test validation");
    }

    //128K
    if (!simple_file_test(256)){
        perror("Failed simple_file_test validation");
    }

    //128K
    if (!simple_file_test(512)){
        perror("Failed simple_file_test validation");
    }

    //128K
    if (!simple_file_test(1024)){
        perror("Failed simple_file_test validation");
    }

#ifdef TITAN
    term_bb_proxy();
#endif
    MPI_Finalize();

    return 0;
}
