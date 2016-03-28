#include <mpi.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <bbapi.h>
#include <errno.h>
int main(int argc, char **argv){

    int testfile = 0;

    spawn_bb_proxy();

    if ((testfile = open("/tmp/persist/a/b/testfile", O_WRONLY | O_CREAT | O_TRUNC,
                        S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) == -1)
    {
            perror("Cannot open output file\n"); exit(1);
    }

    write(testfile,"test",5);

    close(testfile);

    sleep(3);
    
    term_bb_proxy();
    return 0;

}
