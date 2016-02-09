# IOIntercept
Repository for LIBIOINTERCEPT 

N:N Burst Buffer solution for applications that require very little modification.

For both dynamic and statically linked applications we will catch the IO close calls and 
automatically invoke the BBAPI to asynchronously transfer the files to the parallel file-system
location that is defined by an environmental variable.

We may also implement neighbor replication extensions handled through CCI

This code will allow the application users to employ the node-local burst buffers without
having to make modifications to their code -- except maybe updating the write to location
