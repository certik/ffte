FFTE: A Fast Fourier Transform Package

Description:
    A package to compute Discrete Fourier Transforms of
    1-, 2- and 3- dimensional sequences of length (2^p)*(3^q)*(5^r).

Files:
    zfft1d.f    : 1-dimensional complex FFT routine
    zfft2d.f    : 2-dimensional complex FFT routine
    zfft3d.f    : 3-dimensional complex FFT routine
    kernel.f    : FFT kernel routine (radix-2,3,4,5 and 8)
    param.h     : Header file for parameters
    readme.txt  : Readme file
    tests/      : Test Directory
        Makefile    : Makefile for test programs.
        test1d.f    : Test program for zfft1d
        test2d.f    : Test program for zfft2d
        test3d.f    : Test program for zfft3d
        speed1d.f   : Speed test program for zfft1d
        speed2d.f   : Speed test program for zfft2d
        speed3d.f   : Speed test program for zfft3d
    mpi/        : MPI version Directory
        pzfft1d.f   : 1-dimensional complex FFT routine (MPI version)
        pzfft2d.f   : 2-dimensional complex FFT routine (MPI version)
        pzfft3d.f   : 3-dimensional complex FFT routine (MPI version)
        pztrans.f   : Global transpose routine (MPI version)
        tests/      : Test Directory
            Makefile    : Makefile for test programs.
            ptest1d.f   : Test program for pzfft1d (MPI version)
            ptest2d.f   : Test program for pzfft2d (MPI version)
            ptest3d.f   : Test program for pzfft3d (MPI version)
            pspeed1d.f  : Speed test program for pzfft1d (MPI version)
            pspeed2d.f  : Speed test program for pzfft2d (MPI version)
            pspeed3d.f  : Speed test program for pzfft3d (MPI version)

Reference:
    1. Daisuke Takahashi: A Blocking Algorithm for FFT on Cache-Based
       Processors, Proc. 9th International Conference on High
       Performance Computing and Networking Europe (HPCN Europe 2001),
       Lecture Notes in Computer Science, No. 2110, Springer-Verlag,
       pp. 551--554, 2001.

    2. Daisuke Takahashi: A Blocking Algorithm for Parallel 1-D FFT on
       Shared-Memory Parallel Computers, Proc. 6th International
       Conference on Applied Parallel Computing (PARA 2002),
       Lecture Notes in Computer Science, No. 2367, Springer-Verlag,
       pp. 380--389, 2002.

    3. Daisuke Takahashi: Efficient implementation of parallel
       three-dimensional FFT on clusters of PCs, Computer Physics
       Communications, Vol. 152, pp. 144--150, 2003.

Copyright:
    Copyright(C) 2000-2003 Daisuke Takahashi
    Institute of Information Sciences and Electronics,
    University of Tsukuba
    1-1-1 Tennodai, Tsukuba, Ibaraki 305-8573, Japan
    e-mail: daisuke@is.tsukuba.ac.jp
    You may use, copy, modify this code for any purpose and
    without fee.
    You may distribute this ORIGINAL package.
