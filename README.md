A C++ implementation of the GSW Fully Homomorphic Encryption (FHE) scheme using the NTL and GMP C++ libraries.

Sources:
NTL C++ library: https://libntl.org/
GMP C++ library: https://gmplib.org/

Articles and references:
GSW scheme: https://eprint.iacr.org/2013/340.pdf

other related articles:
https://web.eecs.umich.edu/~cpeikert/pubs/polyboot.pdf
https://eprint.iacr.org/2021/691.pdf
https://eprint.iacr.org/2020/086.pdf

#################################################################

The program can work both in Microsoft Windows and Linux. In Windows, it is more comfortable to use the MSYS2 unix-like environment (https://www.msys2.org/wiki/Home/).
In linux, the program can be easily compiled and run.

The implementation depends on NTL C++ library. The NTL library, in turn, can be compiled by itself or based on GMP multiprecision C++ library.

Before compiling the program, make sure the flags (inluded in the Makefile) are set based on your system.


#################################################################

To run the program, oopen a terminal and simply run the 'make' and then 'make run' commands in the terminal.


#################################################################

Note: this is an initial implementation of the GSW scheme. More implementations, including two implementations of the RGSW scheme (Ring variant of the GSW scheme) will be provided soon.
