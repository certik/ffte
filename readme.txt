FFTE: A Fast Fourier Transform Package

Description:
    A package to compute Discrete Fourier Transforms of
    1-, 2- and 3- dimensional sequences of length (2^p)*(3^q)*(5^r).

Files:
    fft235.f    : Radix-2,3,4,5 and 8 FFT routine
    kernel.f    : Radix-2,3,4,5 and 8 FFT kernel routine
    mfft235.f   : Radix-2,3,4,5 and 8 multiple FFT routine
    param.h     : Header file for parameters
    readme.txt  : Readme file
    sse2.c      : Radix-2,3,4,5 and 8 FFT kernel routine (for SSE2)
    sse3.c      : Radix-2,3,4,5 and 8 FFT kernel routine (for SSE3)
    vzfft1d.f   : 1-D complex FFT routine (for vector machines)
    vzfft2d.f   : 2-D complex FFT routine (for vector machines)
    vzfft3d.f   : 3-D complex FFT routine (for vector machines)
    zfft1d.f    : 1-D complex FFT routine
    zfft2d.f    : 2-D complex FFT routine
    zfft3d.f    : 3-D complex FFT routine
    tests/      : Test Directory
        Makefile       : Makefile for test programs
        Makefile.sse2  : Makefile for test programs (for SSE2)
        Makefile.sse3  : Makefile for test programs (for SSE3)
        Makefile.vec   : Makefile for test programs (for vector machines)
        speed1d.f      : Speed test program for zfft1d
        speed2d.f      : Speed test program for zfft2d
        speed3d.f      : Speed test program for zfft3d
        test1d.f       : Test program for zfft1d
        test2d.f       : Test program for zfft2d
        test3d.f       : Test program for zfft3d
    mpi/        : MPI version Directory
        pvzfft1d.f : Parallel 1-D complex FFT routine (for vector machines)
        pvzfft2d.f : Parallel 2-D complex FFT routine (for vector machines)
        pvzfft3d.f : Parallel 3-D complex FFT routine (for vector machines)
        pzfft1d.f  : Parallel 1-D complex FFT routine
        pzfft2d.f  : Parallel 2-D complex FFT routine
        pzfft3d.f  : Parallel 3-D complex FFT routine
        pztrans.f  : Global transpose routine
        tests/     : Test Directory
            Makefile       : Makefile for test programs
            Makefile.sse2  : Makefile for test programs (for SSE2)
            Makefile.sse3  : Makefile for test programs (for SSE3)
            Makefile.vec   : Makefile for test programs (for vector machines)
            pspeed1d.f     : Speed test program for pzfft1d
            pspeed2d.f     : Speed test program for pzfft2d
            pspeed3d.f     : Speed test program for pzfft3d
            ptest1d.f      : Test program for pzfft1d
            ptest2d.f      : Test program for pzfft2d
            ptest3d.f      : Test program for pzfft3d

Reference:
    1. Daisuke Takahashi: A Parallel 3-D FFT Algorithm on Clusters of
       Vector SMPs, Proc. 5th International Workshop on Applied
       Parallel Computing (PARA 2000), Lecture Notes in Computer
       Science, No. 1947, Springer-Verlag, pp. 316--323, 2001.

    2. Daisuke Takahashi: A Blocking Algorithm for FFT on Cache-Based
       Processors, Proc. 9th International Conference on High
       Performance Computing and Networking Europe (HPCN Europe 2001),
       Lecture Notes in Computer Science, No. 2110, Springer-Verlag,
       pp. 551--554, 2001.

    3. Daisuke Takahashi: A Blocking Algorithm for Parallel 1-D FFT on
       Shared-Memory Parallel Computers, Proc. 6th International
       Conference on Applied Parallel Computing (PARA 2002),
       Lecture Notes in Computer Science, No. 2367, Springer-Verlag,
       pp. 380--389, 2002.

    4. Daisuke Takahashi: Efficient implementation of parallel
       three-dimensional FFT on clusters of PCs, Computer Physics
       Communications, Vol. 152, pp. 144--150, 2003.

    5. Daisuke Takahashi: A parallel 1-D FFT algorithm for the Hitachi
       SR8000, Parallel Computing, Vol. 29, pp. 679--690, 2003.

Copyright:
    Copyright(C) 2000-2004 Daisuke Takahashi
    Graduate School of Systems and Information Engineering
    University of Tsukuba
    1-1-1 Tennodai, Tsukuba, Ibaraki 305-8573, Japan
    e-mail: daisuke@cs.tsukuba.ac.jp
    You may use, copy, modify this code for any purpose and
    without fee.
    You may distribute this ORIGINAL package.
