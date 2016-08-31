#include <mpi.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
int main(int argc, char **argv){

    int testfile = 0;
    int rank = 0;
    int thread_level;
    char filename[255];


#ifdef TITAN
    MPI_Init_thread(&argc,&argv, MPI_THREAD_FUNNELED, &thread_level);

    if (thread_level != MPI_THREAD_FUNNELED){
        MPI_Finalize();
    }
#else
    MPI_Init(&argc,&argv);
#endif
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    sprintf(filename,"/tmp/persist/a/b/testfile_%d",rank);

    if ((testfile = open(filename, O_WRONLY | O_CREAT | O_TRUNC,
                    S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) == -1)
    {
        fprintf(stderr,"Cannot open output file %s\n",strerror(errno)); exit(1);
    }

    write(testfile,"test",5);

    close(testfile);

    sleep(3);

    MPI_Finalize();
    return 0;

}
