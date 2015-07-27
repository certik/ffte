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
C     1-D COMPLEX FFT ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     CALL ZFFT1D(A,N,IOPT,B)
C
C     A(N) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     B(N) IS WORK VECTOR (COMPLEX*16)
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
      COMPLEX*16 C((NDA2+NP)*NBLK),D(NDA2)
      COMPLEX*16 WX(NDA2),WY(NDA2),W1(NDA2),W2(NDA2),W3(NDA2),W4(NDA2)
      DIMENSION IP(3),LNX(3),LNY(3)
      SAVE WX,WY,W1,W2,W3,W4
C
      CALL FACTOR(N,IP)
C
      IF (IOPT .EQ. 1) THEN
!$OMP PARALLEL DO
        DO 10 I=1,N
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
      IF (N .LE. MIN0(L2SIZE/16/3,NDA2)) THEN
        IF (IOPT .EQ. 0) THEN
          CALL SETTBL(WX,N)
          RETURN
        END IF
        CALL FFT235(A,B,WX,N,IP)
      ELSE
        DO 20 I=1,3
          LNX(I)=(IP(I)+1)/2
          LNY(I)=IP(I)-LNX(I)
   20   CONTINUE
        NX=(2**LNX(1))*(3**LNX(2))*(5**LNX(3))
        NY=(2**LNY(1))*(3**LNY(2))*(5**LNY(3))
        IF (2**LNX(1) .LT. NBLK .OR. 2**LNY(1) .LT. NBLK) THEN
          MX=(2**(LNX(1)/2))*(3**(LNX(2)/2))*(5**(LNX(3)/2))
          MY=(2**(LNY(1)/2))*(3**(LNY(2)/2))*(5**(LNY(3)/2))
        ELSE
          MX=MIN0(NX,MAX0(NBLK,2**(LNX(1)/2)))
          MY=MIN0(NY,MAX0(NBLK,2**(LNY(1)/2)))
        END IF
C
        IF (IOPT .EQ. 0) THEN
          CALL SETTBL(WX,NX)
          CALL SETTBL(WY,NY)
          CALL SETTBL2(W1,W2,W3,W4,NX,NY,MX,MY)
          RETURN
        END IF
C
!$OMP PARALLEL PRIVATE(C,D)
        CALL ZFFT1D0(A,A,B,C,D,WX,WY,W1,W2,W3,W4,NX,NY,MX,MY,LNX,LNY)
!$OMP END PARALLEL
      END IF
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/DBLE(N)
!$OMP PARALLEL DO
        DO 30 I=1,N
          A(I)=DCONJG(A(I))*DN
   30   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZFFT1D0(A,AYX,B,C,D,WX,WY,W1,W2,W3,W4,NX,NY,MX,MY,
     1                   LNX,LNY)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(NX,*),AYX(NY,*),B(NX,*),C(NY+NP,*),D(*)
      COMPLEX*16 WX(*),WY(*)
      COMPLEX*16 W1(MX,*),W2(MX,*),W3(MY,*),W4(NX/MX,*)
      COMPLEX*16 TEMP
      DIMENSION LNX(*),LNY(*)
C
!$OMP DO PRIVATE(IJ,IJ0,IR,J,TEMP)
      DO 110 II=1,NX,NBLK
        DO 30 JJ=1,NY,NBLK
          DO 20 I=II,MIN0(II+NBLK-1,NX)
            DO 10 J=JJ,MIN0(JJ+NBLK-1,NY)
              C(J,I-II+1)=A(I,J)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
        DO 40 I=II,MIN0(II+NBLK-1,NX)
          CALL FFT235(C(1,I-II+1),D,WY,NY,LNY)
   40   CONTINUE
        IF (2**LNX(1) .LT. NBLK .OR. 2**LNY(1) .LT. NBLK) THEN
          DO 70 IS=1,NY/MY
            DO 60 IK=1,MY
              J=IK+(IS-1)*MY
              DO 50 I=II,MIN0(II+NBLK-1,NX)
                IR=(I-1)/MX+1
                IJ=MOD(I-1,MX)+1
                B(I,J)=C(J,I-II+1)*(W1(IJ,IK)*W2(IJ,IS)
     1                *W3(IK,IR)*W4(IR,IS))
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE
        ELSE
          IR=(II-1)/MX+1
          IJ0=MOD(II-1,MX)+1
          DO 100 IS=1,NY/MY
            DO 90 IK=1,MY
              TEMP=W3(IK,IR)*W4(IR,IS)
              J=IK+(IS-1)*MY
              IJ=IJ0
              DO 80 I=II,MIN0(II+NBLK-1,NX)
                B(I,J)=C(J,I-II+1)*(W1(IJ,IK)*W2(IJ,IS)*TEMP)
                IJ=IJ+1
   80         CONTINUE
   90       CONTINUE
  100     CONTINUE
        END IF
  110 CONTINUE
!$OMP DO
      DO 150 JJ=1,NY,NBLK
        DO 120 J=JJ,MIN0(JJ+NBLK-1,NY)
          CALL FFT235(B(1,J),C,WX,NX,LNX)
  120   CONTINUE
        DO 140 I=1,NX
          DO 130 J=JJ,MIN0(JJ+NBLK-1,NY)
            AYX(J,I)=B(I,J)
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL2(W1,W2,W3,W4,NX,NY,MX,MY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W1(2,MX,*),W2(2,MX,*),W3(2,MY,*),W4(2,NX/MX,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(NX)*DBLE(NY))
!$OMP PARALLEL DO
      DO 30 K=1,MY
!DIR$ VECTOR ALIGNED
        DO 10 J=1,MX
          W1(1,J,K)=DCOS(PX*DBLE(J-1)*DBLE(K-1))
          W1(2,J,K)=DSIN(PX*DBLE(J-1)*DBLE(K-1))
   10   CONTINUE
!DIR$ VECTOR ALIGNED
        DO 20 IR=1,NX/MX
          W3(1,K,IR)=DCOS(PX*DBLE(K-1)*DBLE(IR-1)*DBLE(MX))
          W3(2,K,IR)=DSIN(PX*DBLE(K-1)*DBLE(IR-1)*DBLE(MX))
   20   CONTINUE
   30 CONTINUE
!$OMP PARALLEL DO
      DO 60 IS=1,NY/MY
!DIR$ VECTOR ALIGNED
        DO 40 J=1,MX
          W2(1,J,IS)=DCOS(PX*DBLE(J-1)*DBLE(IS-1)*DBLE(MY))
          W2(2,J,IS)=DSIN(PX*DBLE(J-1)*DBLE(IS-1)*DBLE(MY))
   40   CONTINUE
!DIR$ VECTOR ALIGNED
        DO 50 IR=1,NX/MX
          W4(1,IR,IS)=DCOS(PX*DBLE(IR-1)*DBLE(MX)*DBLE(IS-1)*DBLE(MY))
          W4(2,IR,IS)=DSIN(PX*DBLE(IR-1)*DBLE(MX)*DBLE(IS-1)*DBLE(MY))
   50   CONTINUE
   60 CONTINUE
      RETURN
      END
