C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
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
C     1-D COMPLEX FFT ROUTINE (FOR VECTOR MACHINES)
C
C     FORTRAN77 SOURCE PROGRAM
C
C     CALL ZFFT1D(A,N,IOPT,B)
C
C     A(N) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     B(N*2) IS WORK/COEFFICIENT VECTOR (COMPLEX*16)
C     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4)
C       -----------------------------------
C         N = (2**IP) * (3**IQ) * (5**IR)
C       -----------------------------------
C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
C          = -1 FOR FORWARD TRANSFORM
C          = +1 FOR INVERSE TRANSFORM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT1D(A,N,IOPT,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(*),B(*)
      COMPLEX*16 WX(NDA2),WY(NDA2)
      DIMENSION IP(3),LNX(3),LNY(3)
      SAVE WX,WY
C
      CALL FACTOR(N,IP)
C
      IF (IOPT .EQ. 1) THEN
        DO 10 I=1,N
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
      DO 20 I=1,3
        LNX(I)=(IP(I)+1)/2
        LNY(I)=IP(I)-LNX(I)
   20 CONTINUE
      NX=(2**LNX(1))*(3**LNX(2))*(5**LNX(3))
      NY=(2**LNY(1))*(3**LNY(2))*(5**LNY(3))
C
      IF (IOPT .EQ. 0) THEN
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL2(B(N+1),NY,NX)
        RETURN
      END IF
C
      CALL MFFT235A(A,B,WY,NX,NY,LNY)
      CALL ZTRANSMUL(A,B,B(N+1),NX,NY)
      CALL MFFT235B(B,A,WX,NY,NX,LNX)
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/DBLE(N)
        DO 30 I=1,N
          A(I)=DCONJG(A(I))*DN
   30   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZTRANSMUL(A,B,W,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION LNX(3),LNY(3)
C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
C
      IF (NX .EQ. 1 .OR. NY .EQ. 1) THEN
        DO 10 I=1,NX*NY
          B(I)=A(I)*W(I)
   10   CONTINUE
        RETURN
      END IF
C
      IF (LNX(1)+LNY(1) .LE. 1) THEN
        CALL ZTRANSMULA(A,B,W,NX,NY)
      ELSE
        CALL ZTRANSMULB(A,B,W,NX,NY)
      END IF
      RETURN
      END
      SUBROUTINE ZTRANSMULA(A,B,W,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N1,*),B(N2,*),W(N2,*)
C
      DO 20 I=1,N1
        DO 10 J=1,N2
          B(J,I)=A(I,J)*W(J,I)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE ZTRANSMULB(A,B,W,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N1,*),B(N2,*),W(N2,*)
C
      IF (N2 .GE. N1) THEN
        DO 20 I=0,N1-1
          DO 10 J=1,N1-I
            B(J,I+J)=A(I+J,J)*W(J,I+J)
   10     CONTINUE
   20   CONTINUE
        DO 40 I=1,N2-N1
          DO 30 J=1,N1
            B(I+J,J)=A(J,I+J)*W(I+J,J)
   30     CONTINUE
   40   CONTINUE
        DO 60 I=N2-N1+1,N2-1
          DO 50 J=1,N2-I
            B(I+J,J)=A(J,I+J)*W(I+J,J)
   50     CONTINUE
   60   CONTINUE
      ELSE
        DO 80 I=0,N2-1
          DO 70 J=1,N2-I
            B(I+J,J)=A(J,I+J)*W(I+J,J)
   70     CONTINUE
   80   CONTINUE
        DO 100 I=1,N1-N2
          DO 90 J=1,N2
            B(J,I+J)=A(I+J,J)*W(J,I+J)
   90     CONTINUE
  100   CONTINUE
        DO 120 I=N1-N2+1,N1-1
          DO 110 J=1,N1-I
            B(J,I+J)=A(I+J,J)*W(J,I+J)
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE SETTBL2(W,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,NX,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(NX)*DBLE(NY))
      DO 20 K=1,NY
        DO 10 J=1,NX
          W(1,J,K)=DCOS(PX*DBLE(J-1)*DBLE(K-1))
          W(2,J,K)=DSIN(PX*DBLE(J-1)*DBLE(K-1))
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
