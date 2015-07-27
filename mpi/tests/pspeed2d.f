C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     PZFFT2D SPEED TEST PROGRAM
C
C     FORTRAN77 + MPI SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mpif.h'
      PARAMETER (NDA=1048576,LOOP=10)
      COMPLEX*16 A(NDA),B(NDA)
      DIMENSION LNX(3),LNY(3)
      SAVE A,B
C
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPU,IERR)
C
      IF (ME .EQ. 0) THEN
        WRITE(6,*) ' NX,NY ='
        READ(5,*) NX,NY
      END IF
      CALL MPI_BCAST(NX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NY,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
C
      NN=NX*(NY/NPU)
      CALL INIT(A,NN,ME,NPU)
      CALL PZFFT2D(A,B,NX,NY,MPI_COMM_WORLD,NPU,0)
      CALL PZFFT2D(A,B,NX,NY,MPI_COMM_WORLD,NPU,-1)
C
      TIME0=2.0D0**52
      DO I=1,LOOP
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      IF (ME .EQ. 0) THEN
        TIME1=MPI_WTIME()
      END IF
C
      CALL PZFFT2D(A,B,NX,NY,MPI_COMM_WORLD,NPU,-1)
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      IF (ME .EQ. 0) THEN
        TIME2=MPI_WTIME()
        TIME0=DMIN1(TIME0,TIME2-TIME1)
      END IF
      END DO
      IF (ME .EQ. 0) THEN
        FLOPS=(2.5D0*DBLE(LNX(1)+LNY(1))
     1         +4.66666666666666D0*DBLE(LNX(2)+LNY(2))
     2         +6.8D0*DBLE(LNX(3)+LNY(3)))*2.0D0*DBLE(NX)*DBLE(NY)
     3         /TIME0/1.0D6
        WRITE(6,*) ' NPU =',NPU
        WRITE(6,*) ' NX =',NX,' NY =',NY,' TIME =',TIME0,FLOPS,' MFLOPS'
      END IF
C
      CALL MPI_FINALIZE(IERR)
      STOP
      END
      SUBROUTINE INIT(A,NN,ME,NPU)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*)
      INTEGER*8 N
C
      N=INT8(NN)*INT8(NPU)
      DO 10 I=1,NN
C        A(I)=DCMPLX(DBLE(I)+DBLE(NN)*DBLE(ME),
C     1              DBLE(N)-(DBLE(I)+DBLE(NN)*DBLE(ME))+1.0D0)
        A(I)=(0.0D0,0.0D0)
   10 CONTINUE
      RETURN
      END
