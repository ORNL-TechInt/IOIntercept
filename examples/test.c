#include <mpi.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

#ifdef TITAN
#include "bbapi.h"
#endif

int main(int argc, char **argv){
    int testfile = 0;

#ifdef TITAN
    spawn_bb_proxy();	    
#endif

    char *srcdir = argv[1];
    char tfn[256];

    snprintf(tfn,sizeof(tfn), "%s/testfile.%d", srcdir,getpid());

    if ((testfile = open(tfn, O_WRONLY | O_CREAT | O_TRUNC,
		    S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) == -1)
    {
	perror("Cannot open output file\n"); exit(1);
    }

    write(testfile,"test",5);

    close(testfile);


#ifdef TITAN
    sleep(3); //Don't need the sleep in the non-emulated version bbproxy version
    term_bb_proxy();
#endif

    return 0;
}
