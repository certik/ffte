C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000, 2001, 2002, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS,
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@is.tsukuba.ac.jp
C
C
C     HEADER FILE FOR PARAMETERS
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
C The maximum supported 2-D transform length is 32768.
      PARAMETER (NDA2=32768)
C The maximum supported 3-D transform length is 1024.
      PARAMETER (NDA3=1024)
      PARAMETER (NDA4=256)
C The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=8)
C      PARAMETER (NBLK=4)  (for PentiumII and Celeron)
C      PARAMETER (NBLK=8)  (for PentiumIII)
C      PARAMETER (NBLK=16) (for Pentium4, Athlon MP and UltraSPARC-II)
C      PARAMETER (NBLK=32) (for MIPS R14000)
C The parameter NP is a padding parameter to avoid bank and cache
C conflicts in the FFT routines.
C NP=2 appears to work well on most systems.
      PARAMETER (NP=2)
C Size of L2 cache
      PARAMETER (L2SIZE=262144)
