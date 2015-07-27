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
C     1-DIMENSIONAL COMPLEX FFT ROUTINE (MPI VERSION)
C
C     FORTRAN77 SOURCE PROGRAM
C
C     CALL PZFFT1D(A,B,WW,WWW,N,ME,NPU,IOPT)
C
C     A(N/NPU) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(CYCLIC)
C     B(N/NPU) IS WORK VECTOR (COMPLEX*16)
C     WW(NDA3*NDA3) IS SINE/COSINE TABLE (COMPLEX*16)
C     WWW(N/NPU) IS SINE/COSINE TABLE (COMPLEX*16)
C     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4)
C       -----------------------------------
C         N = (2**IP) * (3**IQ) * (5**IR)
C       -----------------------------------
C     ME IS PROCESSOR NUMBER (INTEGER*4)
C     NPU IS THE NUMBER OF PROCESSORS (INTEGER*4)
C     IOPT = 1 FOR FORWARD TRANSFORM (INTEGER*4)
C          = 2 FOR INVERSE TRANSFORM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE PZFFT1D(A,B,WW,WWW,N,ME,NPU,IOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(*),B(*),WW(*),WWW(*)
      COMPLEX*16 C((NDA3+NP)*NBLK+NP),D(NDA3+NP)
      COMPLEX*16 WX(NDA3/2+NP),WY(NDA3/2+NP),WZ(NDA3/2+NP)
      DIMENSION IP(3),LNX(3),LNY(3),LNZ(3),LNN(3),LNPU(3)
      DATA N0/0/
      SAVE N0,WX,WY,WZ
C
      NN=N/NPU
      CALL FACTOR(NPU,LNPU)
      CALL FACTOR(NN,LNN)
      DO 10 I=1,3
        IP(I)=LNN(I)+LNPU(I)
        LNZ(I)=MAX0(LNPU(I),(IP(I)+1)/3)
        LNX(I)=MAX0(LNPU(I),(IP(I)-LNZ(I)+1)/2)
        LNY(I)=IP(I)-LNX(I)-LNZ(I)
   10 CONTINUE
      NX=(2**LNX(1))*(3**LNX(2))*(5**LNX(3))
      NY=(2**LNY(1))*(3**LNY(2))*(5**LNY(3))
      NZ=(2**LNZ(1))*(3**LNZ(2))*(5**LNZ(3))
C
      IF (IOPT .EQ. 2) THEN
        DO 20 I=1,NN
          A(I)=DCONJG(A(I))
   20   CONTINUE
      END IF
C
      IF (N .NE. N0) THEN
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL(WZ,NZ)
        CALL SETTBL3(WW,NY,NZ)
        CALL SETTBLP(WWW,NX,NY,NZ,ME,NPU)
        N0=N
      END IF
!$OMP PARALLEL PRIVATE(C,D)
      CALL PZFFT1D0(A,A,A,B,B,C,C,D,WX,WY,WZ,WW,WWW,NX,NY,NZ,NPU)
!$OMP END PARALLEL
C
      IF (IOPT .EQ. 2) THEN
        DN=1.0D0/DBLE(N)
        DO 30 I=1,NN
          A(I)=DCONJG(A(I))*DN
   30   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE PZFFT1D0(A,AX,AZ,B,BX,CY,CZ,D,WX,WY,WZ,WW,WWW,
     1                    NX,NY,NZ,NPU)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(NX/NPU,NY,*),AX(NX/NPU,NY,NZ/NPU,*),AZ(NZ/NPU,NY,*)
      COMPLEX*16 B(NX,NY,*),BX(NX/NPU,NY,NZ/NPU,*)
      COMPLEX*16 CY(NY+NP,*),CZ(NZ+NP,*),D(*)
      COMPLEX*16 WX(*),WY(*),WZ(*),WW(NZ,*),WWW(NY,NX,*)
      DIMENSION LNX(3),LNY(3),LNZ(3)
C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)
C
      NNX=NX/NPU
      NNZ=NZ/NPU
      NN=NX*NY*NZ/NPU
C
!$OMP DO
      DO 90 J=1,NY
        DO 80 II=1,NNX,NBLK
          DO 30 KK=1,NZ,NBLK
            DO 20 I=II,MIN0(II+NBLK-1,NNX)
              DO 10 K=KK,MIN0(KK+NBLK-1,NZ)
                CZ(K,I-II+1)=A(I,J,K)
   10         CONTINUE
   20       CONTINUE
   30     CONTINUE
          DO 40 I=II,MIN0(II+NBLK-1,NNX)
            CALL FFT235A(CZ(1,I-II+1),D,WZ,NZ,LNZ)
   40     CONTINUE
          DO 70 L=1,NPU
            DO 60 K=1,NNZ
              DO 50 I=II,MIN0(II+NBLK-1,NNX)
                BX(I,J,K,L)=CZ(L+(K-1)*NPU,I-II+1)*WW(L+(K-1)*NPU,J)
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE
   80   CONTINUE
   90 CONTINUE
!$OMP SINGLE
      CALL PZTRANS(BX,AX,NN,NPU)
!$OMP END SINGLE
!$OMP DO
      DO 190 K=1,NNZ
        DO 170 L=1,NPU
          DO 160 II=1,NNX,NBLK
            DO 120 JJ=1,NY,NBLK
              DO 110 I=II,MIN0(II+NBLK-1,NNX)
                DO 100 J=JJ,MIN0(JJ+NBLK-1,NY)
                  CY(J,I-II+1)=AX(I,J,K,L)
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
            DO 130 I=II,MIN0(II+NBLK-1,NNX)
              CALL FFT235A(CY(1,I-II+1),D,WY,NY,LNY)
  130       CONTINUE
            DO 150 J=1,NY
              DO 140 I=II,MIN0(II+NBLK-1,NNX)
                B(L+(I-1)*NPU,J,K)=CY(J,I-II+1)*WWW(J,L+(I-1)*NPU,K)
  140         CONTINUE
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
        DO 180 J=1,NY
          CALL FFT235A(B(1,J,K),D,WX,NX,LNX)
  180   CONTINUE
  190 CONTINUE
!$OMP DO
      DO 250 II=1,NX,NBLK
        DO 240 JJ=1,NY,NBLK
          DO 230 KK=1,NNZ,NBLK
            DO 220 I=II,MIN0(II+NBLK-1,NX)
              DO 210 J=JJ,MIN0(JJ+NBLK-1,NY)
                DO 200 K=KK,MIN0(KK+NBLK-1,NNZ)
                  AZ(K,J,I)=B(I,J,K)
  200           CONTINUE
  210         CONTINUE
  220       CONTINUE
  230     CONTINUE
  240   CONTINUE
  250 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL3(W,NY,NZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,NZ,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/DBLE(NY*NZ)
!$OMP PARALLEL DO
      DO 20 J=1,NY
!DIR$ VECTOR ALIGNED
        DO 10 K=1,NZ
          W(1,K,J)=DCOS(PX*DBLE(K-1)*DBLE(J-1))
          W(2,K,J)=DSIN(PX*DBLE(K-1)*DBLE(J-1))
   10   CONTINUE
   20 CONTINUE
!$OMP END PARALLEL DO
      RETURN
      END
      SUBROUTINE SETTBLP(W,NX,NY,NZ,ME,NPU)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,NY,NX,*)
C
      N=NX*NY*NZ
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/DBLE(N)
!$OMP PARALLEL DO
      DO 30 K=1,NZ/NPU
        DO 20 I=1,NX
!DIR$ VECTOR ALIGNED
          DO 10 J=1,NY
            W(1,J,I,K)=DCOS(PX*DBLE(I-1)*(DBLE(ME)+DBLE(K-1)*DBLE(NPU)
     1                                    +DBLE(J-1)*DBLE(NZ)))
            W(2,J,I,K)=DSIN(PX*DBLE(I-1)*(DBLE(ME)+DBLE(K-1)*DBLE(NPU)
     1                                    +DBLE(J-1)*DBLE(NZ)))
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
!$OMP END PARALLEL DO
      RETURN
      END
