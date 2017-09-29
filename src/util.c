#include "spectral.h"

int check(int rc)
{
    if(rc)
    {
        char* errstring = 0;
        //getLastErrorDetails(BBERRORJSON, &errstring);
        printf("Error rc:       %d\n", rc);
        printf("Error details:  %s\n", errstring);
        free(errstring);
        return FAILURE;
    }
    return SUCCESS;
}


char *readlink_malloc (const char *filename)
{
    int size = 100;
    char *buffer = NULL;

    while (1)
    {
        buffer = (char *) realloc (buffer, size);
        int nchars = readlink (filename, buffer, size);
        if (nchars < 0)
        {
            free (buffer);
            return NULL;
        }
        if (nchars < size)
        {
            /* Readlink doesn't null terminate the string it sends back.
               But a nice property of this is that if nchars == size you can't
               tell if you read the whole thing so size ALWAYS has to be greater
               this means you have a slot available to stash your null terminator.

               Hooray.             
               */
            buffer[nchars] = '\0';               
            return buffer;
        }
        size *= 2;
    }
}
