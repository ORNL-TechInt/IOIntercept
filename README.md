# IOIntercept
Repository for LIBIOINTERCEPT 

N:N Burst Buffer solution for applications that require very little modification.

For both dynamic and statically linked applications we will catch the IO close calls and 
automatically invoke the BBAPI to asynchronously transfer the files to the parallel file-system
location that is defined by an environmental variable.

This code will allow the application users to employ the node-local burst buffers without
having to make modifications to their code -- except maybe updating the write to location

To test basic implementation on systems without IBM BBAPI use build option TITAN=ON in cmake.
This builds libbbemulate, which enables the spanwing of an additional pthread to pretends to be
BBPRoxy. Libbbemulate is far less feature rich than the real BBAPI so it should only be used for
basic testing.

CMake Option DEBUG=ON enables diagnostic output and dwarf symbols
