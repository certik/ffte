C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@is.tsukuba.ac.jp
C
C
C     2-DIMENSIONAL COMPLEX FFT ROUTINE (MPI VERSION)
C
C     FORTRAN77 + MPI SOURCE PROGRAM
C
C     CALL PZFFT2D(A,B,NX,NY,ME,NPU,IOPT)
C
C     A(NX/NPU,NY) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(CYCLIC,*)
C     B(NX*NY/NPU) IS WORK VECTOR (COMPLEX*16)
C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
C       ------------------------------------
C         NX = (2**IP) * (3**IQ) * (5**IR)
C         NY = (2**JP) * (3**JQ) * (5**JR)
C       ------------------------------------
C     ME IS PROCESSOR NUMBER (INTEGER*4)
C     NPU IS THE NUMBER OF PROCESSORS (INTEGER*4)
C     IOPT = 1 FOR FORWARD TRANSFORM (INTEGER*4)
C          = 2 FOR INVERSE TRANSFORM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE PZFFT2D(A,B,NX,NY,ME,NPU,IOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(*),B(*)
      COMPLEX*16 C((NDA2+NP)*NBLK+NP),D(NDA2+NP)
      COMPLEX*16 WX(NDA2/2+NP),WY(NDA2/2+NP)
      DATA NX0,NY0/0,0/
      SAVE NX0,NY0,WX,WY
C
      IF (IOPT .EQ. 2) THEN
        NN=NX*NY/NPU
        DO 10 I=1,NN
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
      IF (NX .NE. NX0) THEN
        CALL SETTBL(WX,NX)
        NX0=NX
      END IF
      IF (NY .NE. NY0) THEN
        CALL SETTBL(WY,NY)
        NY0=NY
      END IF
!$OMP PARALLEL PRIVATE(C,D)
      CALL PZFFT2D0(A,A,B,B,C,C,D,WX,WY,NX,NY,NPU)
!$OMP END PARALLEL
C
      IF (IOPT .EQ. 2) THEN
        DN=1.0D0/DBLE(NX*NY)
        DO 20 I=1,NN
          A(I)=DCONJG(A(I))*DN
   20   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE PZFFT2D0(A,AX,B,BX,CX,CY,D,WX,WY,NX,NY,NPU)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(NX/NPU,*),AX(NX/NPU,NY/NPU,*)
      COMPLEX*16 B(NX,*),BX(NX/NPU,NY/NPU,*)
      COMPLEX*16 CX(NX+NP,*),CY(NY+NP,*),D(*)
      COMPLEX*16 WX(*),WY(*)
      DIMENSION LNX(3),LNY(3)
C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
C
      NNX=NX/NPU
      NNY=NY/NPU
      NN=NX*NY/NPU
C
!$OMP DO
      DO 80 II=1,NNX,NBLK
        DO 30 JJ=1,NY,NBLK
          DO 20 I=II,MIN0(II+NBLK-1,NNX)
            DO 10 J=JJ,MIN0(JJ+NBLK-1,NY)
              CY(J,I-II+1)=A(I,J)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
        DO 40 I=II,MIN0(II+NBLK-1,NNX)
          CALL FFT235A(CY(1,I-II+1),D,WY,NY,LNY)
   40   CONTINUE
        DO 70 K=1,NPU
          DO 60 J=1,NNY
            DO 50 I=II,MIN0(II+NBLK-1,NNX)
              BX(I,J,K)=CY(K+(J-1)*NPU,I-II+1)
   50       CONTINUE
   60    CONTINUE
   70   CONTINUE
   80 CONTINUE
!$OMP SINGLE
      CALL PZTRANS(BX,AX,NN,NPU)
!$OMP END SINGLE
!$OMP DO
      DO 160 JJ=1,NNY,NBLK
        DO 120 K=1,NPU
          DO 110 II=1,NNX,NBLK
            DO 100 J=JJ,MIN0(JJ+NBLK-1,NNY)
              DO 90 I=II,MIN0(II+NBLK-1,NNX)
                CX(K+(I-1)*NPU,J-JJ+1)=AX(I,J,K)
   90         CONTINUE
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
        DO 130 J=JJ,MIN0(JJ+NBLK-1,NNY)
          CALL FFT235A(CX(1,J-JJ+1),D,WX,NX,LNX)
  130   CONTINUE
        DO 150 J=JJ,MIN0(JJ+NBLK-1,NNY)
          DO 140 I=1,NX
            B(I,J)=CX(I,J-JJ+1)
  140     CONTINUE
  150   CONTINUE
  160 CONTINUE
!$OMP DO
      DO 210 K=1,NPU
        DO 200 JJ=1,NNY,NBLK
          DO 190 II=1,NNX,NBLK
            DO 180 J=JJ,MIN0(JJ+NBLK-1,NNY)
              DO 170 I=II,MIN0(II+NBLK-1,NNX)
                AX(I,J,K)=B(K+(I-1)*NPU,J)
  170         CONTINUE
  180       CONTINUE
  190     CONTINUE
  200   CONTINUE
  210 CONTINUE
!$OMP SINGLE
      CALL PZTRANS(AX,BX,NN,NPU)
!$OMP END SINGLE
!$OMP DO
      DO 260 K=1,NPU
        DO 250 JJ=1,NNY,NBLK
          DO 240 II=1,NNX,NBLK
            DO 230 J=JJ,MIN0(JJ+NBLK-1,NNY)
              DO 220 I=II,MIN0(II+NBLK-1,NNX)
                A(I,K+(J-1)*NPU)=BX(I,J,K)
  220         CONTINUE
  230       CONTINUE
  240     CONTINUE
  250   CONTINUE
  260 CONTINUE
      RETURN
      END
