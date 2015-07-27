C
C     (C) COPYRIGHT SOFTWARE, 2000, 2001, 2002, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS,
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA-SHI, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@is.tsukuba.ac.jp
C
C
C     PZFFT3D TEST PROGRAM
C
C     FORTRAN77 + MPI SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mpif.h'
      INCLUDE 'param.h'
      PARAMETER (NDA=1048576)
      COMPLEX*16 A(NDA+NP),B(NDA+NP)
C
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPU,IERR)
C
      IF (ME .EQ. 0) THEN
        WRITE(6,*) ' NX,NY,NZ ='
        READ(5,*) NX,NY,NZ
      END IF
      CALL MPI_BCAST(NX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NY,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
C
      NN=NX*NY*NZ/NPU
      CALL INIT(A,NN,ME,NPU)
C
      CALL PZFFT3D(A,B,NX,NY,NZ,ME,NPU,1)
      CALL DUMP(A,NN,ME,NPU)
C
      CALL PZFFT3D(A,B,NX,NY,NZ,ME,NPU,2)
      CALL DUMP(A,NN,ME,NPU)
C
      CALL MPI_FINALIZE(IERR)
      STOP
      END
      SUBROUTINE INIT(A,NN,ME,NPU)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*)
C
      N=NN*NPU
      DO 10 I=1,NN
        A(I)=DCMPLX(DBLE(I-1)*DBLE(NPU)+DBLE(ME)+1.0D0,
     1              DBLE(N)-(DBLE(I-1)*DBLE(NPU)+DBLE(ME)+1.0D0)+1.0D0)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE DUMP(A,NN,ME,NPU)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mpif.h'
      COMPLEX*16 A(*)
C
      DO 20 I=1,NN
        DO 10 J=0,NPU-1
          CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          IF (J .EQ. ME) THEN
            WRITE(6,*) ' ME=',ME,A(I)
          END IF
   10   CONTINUE
   20 CONTINUE
      RETURN
      END