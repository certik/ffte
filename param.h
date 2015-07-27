C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2003, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS
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
      PARAMETER (NBLK=16)
C      PARAMETER (NBLK=8)  (for PentiumIII and Athlon)
C      PARAMETER (NBLK=16) (for Pentium4, Athlon XP, Itanium and
C                           Itanium2)
C The parameter NP is a padding parameter to avoid cache conflicts in
C the FFT routines.
      PARAMETER (NP=8)
C      PARAMETER (NP=2) (for PentiumIII)
C      PARAMETER (NP=4) (for Athlon, Athlon XP, Athlon 64 and Itanium)
C      PARAMETER (NP=8) (for Pentium4 and Itanium2)
C Size of L2 cache
      PARAMETER (L2SIZE=524288)
