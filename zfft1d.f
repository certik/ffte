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
C     1-DIMENSIONAL COMPLEX FFT ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     CALL ZFFT1D(A,B,N,IOPT)
C
C     A(N) IS COMPLEX INPUT VECTOR (COMPLEX*16)
C     B(N) IS COMPLEX OUTPUT VECTOR (COMPLEX*16)
C     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4)
C       -----------------------------------
C         N = (2**IP) * (3**IQ) * (5**IR)
C       -----------------------------------
C     IOPT = 1 FOR FORWARD TRANSFORM (INTEGER*4)
C          = 2 FOR INVERSE TRANSFORM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT1D(A,B,N,IOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(*),B(*)
      COMPLEX*16 C((NDA2+NP)*NBLK+NP),D(NDA2+NP)
      COMPLEX*16 W1(NDA2/2+NP),W2(NDA2/2+NP)
      COMPLEX*16 WW1(NDA2+NDA4*NP+NP),WW2(NDA2+NDA4*NP+NP)
      COMPLEX*16 WW3(NDA2+NDA4*NP+NP),WW4(NDA2+NDA4*NP+NP)
      DIMENSION IP(3),IP1(3),IP2(3)
      DATA N0/0/
      SAVE N0,W1,W2,WW1,WW2,WW3,WW4
C
      CALL FACTOR(N,IP)
C
      IF (N .LE. MIN0(L2SIZE/16/3,NDA2)) THEN
        IF (N .NE. N0) THEN
          CALL SETTBL(W1,N)
          N0=N
        END IF
        IF (IOPT .EQ. 2) THEN
          DO 10 I=1,N
            A(I)=DCONJG(A(I))
   10     CONTINUE
        END IF
        CALL FFT235B(A,B,W1,N,IP)
        IF (IOPT .EQ. 2) THEN
          DN=1.0D0/DBLE(N)
          DO 20 I=1,N
            B(I)=DCONJG(B(I))*DN
   20     CONTINUE
        END IF
      ELSE
        IF (IOPT .EQ. 2) THEN
          DO 30 I=1,N
            A(I)=DCONJG(A(I))
   30     CONTINUE
        END IF
        DO 40 I=1,3
          IP1(I)=(IP(I)+1)/2
          IP2(I)=IP(I)-IP1(I)
   40   CONTINUE
        N1=(2**IP1(1))*(3**IP1(2))*(5**IP1(3))
        N2=(2**IP2(1))*(3**IP2(2))*(5**IP2(3))
        IF (2**IP1(1) .LT. NBLK .OR. 2**IP2(1) .LT. NBLK) THEN
          M1=MIN0(N1,(2**(IP1(1)/2))*(3**(IP1(2)/2))*(5**(IP1(3)/2)))
          M2=MIN0(N2,(2**(IP2(1)/2))*(3**(IP2(2)/2))*(5**(IP2(3)/2)))
        ELSE
          M1=MIN0(N1,MAX0(NBLK,2**(IP1(1)/2)))
          M2=MIN0(N2,MAX0(NBLK,2**(IP2(1)/2)))
        END IF
C
        IF (N .NE. N0) THEN
          CALL SETTBL(W1,N1)
          CALL SETTBL(W2,N2)
          CALL SETTBLS(WW1,WW2,WW3,WW4,N1,N2,M1,M2)
          N0=N
        END IF
!$OMP PARALLEL PRIVATE(C,D)
        CALL ZFFT1D0(A,B,C,D,W1,W2,WW1,WW2,WW3,WW4,N1,N2,M1,M2,IP1,IP2)
!$OMP END PARALLEL
        IF (IOPT .EQ. 2) THEN
          DN=1.0D0/DBLE(N)
          DO 50 I=1,N
            B(I)=DCONJG(B(I))*DN
   50     CONTINUE
        END IF
      END IF
      RETURN
      END
      SUBROUTINE ZFFT1D0(A,B,C,D,W1,W2,WW1,WW2,WW3,WW4,N1,N2,M1,M2,
     1                   IP1,IP2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A(N1,*),B(N2,*),C(N2+NP,*),D(*)
      COMPLEX*16 W1(*),W2(*)
      COMPLEX*16 WW1(M2+NP,*),WW2(N2/M2+NP,*)
      COMPLEX*16 WW3(M2+NP,*),WW4(N2/M2+NP,*)
      COMPLEX*16 TEMP
      DIMENSION IP1(*),IP2(*)
C
!$OMP DO PRIVATE(IJ,IJ0,IR,J,TEMP)
      DO 110 II=1,N1,NBLK
        DO 30 JJ=1,N2,NBLK
          DO 20 I=II,MIN0(II+NBLK-1,N1)
            DO 10 J=JJ,MIN0(JJ+NBLK-1,N2)
              C(J,I-II+1)=A(I,J)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
        DO 40 I=II,MIN0(II+NBLK-1,N1)
          CALL FFT235A(C(1,I-II+1),D,W2,N2,IP2)
   40   CONTINUE
        IF (2**IP1(1) .LT. NBLK .OR. 2**IP2(1) .LT. NBLK) THEN
          DO 70 IS=1,N2/M2
            DO 60 IK=1,M2
              J=IK+(IS-1)*M2
              DO 50 I=II,MIN0(II+NBLK-1,N1)
                IR=(I-1)/M1+1
                IJ=MOD(I-1,M1)+1
                A(I,J)=C(J,I-II+1)*(WW1(IK,IJ)*WW2(IS,IJ)
     1                *WW3(IK,IR)*WW4(IS,IR))
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE
        ELSE
          IR=(II-1)/M1+1
          IJ0=MOD(II-1,M1)+1
          DO 100 IS=1,N2/M2
            DO 90 IK=1,M2
              TEMP=WW3(IK,IR)*WW4(IS,IR)
              J=IK+(IS-1)*M2
              IJ=IJ0
              DO 80 I=II,MIN0(II+NBLK-1,N1)
                A(I,J)=C(J,I-II+1)*(WW1(IK,IJ)*WW2(IS,IJ)*TEMP)
                IJ=IJ+1
   80         CONTINUE
   90       CONTINUE
  100     CONTINUE
        END IF
  110 CONTINUE
!$OMP DO
      DO 150 JJ=1,N2,NBLK
        DO 120 J=JJ,MIN0(JJ+NBLK-1,N2)
          CALL FFT235A(A(1,J),C,W1,N1,IP1)
  120   CONTINUE
        DO 140 I=1,N1
          DO 130 J=JJ,MIN0(JJ+NBLK-1,N2)
            B(J,I)=A(I,J)
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBLS(W1,W2,W3,W4,N1,N2,M1,M2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      DIMENSION W1(2,M2+NP,*),W2(2,N2/M2+NP,*)
      DIMENSION W3(2,M2+NP,*),W4(2,N2/M2+NP,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/DBLE(N1*N2)
!$OMP PARALLEL
!$OMP DO
      DO 30 J=1,M1
!DIR$ VECTOR ALIGNED
        DO 10 K=1,M2
          W1(1,K,J)=DCOS(PX*DBLE(J-1)*DBLE(K-1))
          W1(2,K,J)=DSIN(PX*DBLE(J-1)*DBLE(K-1))
   10   CONTINUE
!DIR$ VECTOR ALIGNED
        DO 20 IS=1,N2/M2
          W2(1,IS,J)=DCOS(PX*DBLE(J-1)*DBLE(IS-1)*DBLE(M2))
          W2(2,IS,J)=DSIN(PX*DBLE(J-1)*DBLE(IS-1)*DBLE(M2))
   20   CONTINUE
   30 CONTINUE
!$OMP DO
      DO 60 IR=1,N1/M1
!DIR$ VECTOR ALIGNED
        DO 40 K=1,M2
          W3(1,K,IR)=DCOS(PX*DBLE(IR-1)*DBLE(M1)*DBLE(K-1))
          W3(2,K,IR)=DSIN(PX*DBLE(IR-1)*DBLE(M1)*DBLE(K-1))
   40   CONTINUE
!DIR$ VECTOR ALIGNED
        DO 50 IS=1,N2/M2
          W4(1,IS,IR)=DCOS(PX*DBLE(IR-1)*DBLE(M1)*DBLE(IS-1)*DBLE(M2))
          W4(2,IS,IR)=DSIN(PX*DBLE(IR-1)*DBLE(M1)*DBLE(IS-1)*DBLE(M2))
   50   CONTINUE
   60 CONTINUE
!$OMP END PARALLEL
      RETURN
      END
