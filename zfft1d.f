C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
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
      COMPLEX*16 C((NDA2+NP)*(NBLK+1)+NP)
      COMPLEX*16 W1(NDA2/2+NP),W2(NDA2/2+NP)
      COMPLEX*16 WW((NDA2+NP)*4+NP)
      DIMENSION IP(3),IP1(3),IP2(3)
      SAVE W1,W2,WW
C
      CALL FACTOR(N,IP)
C
      IF (IOPT .EQ. 1) THEN
        DO 10 I=1,N
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
      IF (N .LE. MIN0(L2SIZE/16/3,NDA2)) THEN
        IF (IOPT .EQ. 0) THEN
          CALL SETTBL(W1,N)
          RETURN
        END IF
        CALL FFT235(A,B,W1,N,IP)
      ELSE
        DO 20 I=1,3
          IP1(I)=(IP(I)+1)/2
          IP2(I)=IP(I)-IP1(I)
   20   CONTINUE
        N1=(2**IP1(1))*(3**IP1(2))*(5**IP1(3))
        N2=(2**IP2(1))*(3**IP2(2))*(5**IP2(3))
        IF (2**IP1(1) .LT. NBLK .OR. 2**IP2(1) .LT. NBLK) THEN
          M1=MIN0(N1,(2**(IP1(1)/2))*(3**(IP1(2)/2))*(5**(IP1(3)/2)))
          M2=MIN0(N2,(2**(IP2(1)/2))*(3**(IP2(2)/2))*(5**(IP2(3)/2)))
        ELSE
          M1=MIN0(N1,MAX0(NBLK,2**(IP1(1)/2)))
          M2=MIN0(N2,MAX0(NBLK,2**(IP2(1)/2)))
        END IF
        NW2=M1*M2+NP
        NW3=NW2+M1*(N2/M2)+NP
        NW4=NW3+M2*(N1/M1)+NP
C
        IF (IOPT .EQ. 0) THEN
          CALL SETTBL(W1,N1)
          CALL SETTBL(W2,N2)
          CALL SETTBLS(WW,WW(NW2+1),WW(NW3+1),WW(NW4+1),N1,N2,M1,M2)
          RETURN
        END IF
C
        ND=(N2+NP)*NBLK+NP
!$OMP PARALLEL PRIVATE(C)
        CALL ZFFT1D0(A,A,B,C,C(ND+1),W1,W2,WW,WW(NW2+1),WW(NW3+1),
     1               WW(NW4+1),N1,N2,M1,M2,IP1,IP2)
!$OMP END PARALLEL
      END IF
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/DBLE(N)
        DO 30 I=1,N
          A(I)=DCONJG(A(I))*DN
   30   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZFFT1D0(A1,A2,B,C,D,W1,W2,WW1,WW2,WW3,WW4,N1,N2,M1,M2,
     1                   IP1,IP2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      COMPLEX*16 A1(N1,*),A2(N2,*),B(N1,*),C(N2+NP,*),D(*)
      COMPLEX*16 W1(*),W2(*)
      COMPLEX*16 WW1(M1,*),WW2(M1,*),WW3(M2,*),WW4(N1/M1,*)
      COMPLEX*16 TEMP
      DIMENSION IP1(*),IP2(*)
C
!$OMP DO PRIVATE(IJ,IJ0,IR,J,TEMP)
      DO 110 II=1,N1,NBLK
        DO 30 JJ=1,N2,NBLK
          DO 20 I=II,MIN0(II+NBLK-1,N1)
            DO 10 J=JJ,MIN0(JJ+NBLK-1,N2)
              C(J,I-II+1)=A1(I,J)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
        DO 40 I=II,MIN0(II+NBLK-1,N1)
          CALL FFT235(C(1,I-II+1),D,W2,N2,IP2)
   40   CONTINUE
        IF (2**IP1(1) .LT. NBLK .OR. 2**IP2(1) .LT. NBLK) THEN
          DO 70 IS=1,N2/M2
            DO 60 IK=1,M2
              J=IK+(IS-1)*M2
              DO 50 I=II,MIN0(II+NBLK-1,N1)
                IR=(I-1)/M1+1
                IJ=MOD(I-1,M1)+1
                B(I,J)=C(J,I-II+1)*(WW1(IJ,IK)*WW2(IJ,IS)
     1                *WW3(IK,IR)*WW4(IR,IS))
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE
        ELSE
          IR=(II-1)/M1+1
          IJ0=MOD(II-1,M1)+1
          DO 100 IS=1,N2/M2
            DO 90 IK=1,M2
              TEMP=WW3(IK,IR)*WW4(IR,IS)
              J=IK+(IS-1)*M2
              IJ=IJ0
              DO 80 I=II,MIN0(II+NBLK-1,N1)
                B(I,J)=C(J,I-II+1)*(WW1(IJ,IK)*WW2(IJ,IS)*TEMP)
                IJ=IJ+1
   80         CONTINUE
   90       CONTINUE
  100     CONTINUE
        END IF
  110 CONTINUE
!$OMP DO
      DO 150 JJ=1,N2,NBLK
        DO 120 J=JJ,MIN0(JJ+NBLK-1,N2)
          CALL FFT235(B(1,J),C,W1,N1,IP1)
  120   CONTINUE
        DO 140 I=1,N1
          DO 130 J=JJ,MIN0(JJ+NBLK-1,N2)
            A2(J,I)=B(I,J)
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBLS(W1,W2,W3,W4,N1,N2,M1,M2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      DIMENSION W1(2,M1,*),W2(2,M1,*),W3(2,M2,*),W4(2,N1/M1,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(N1)*DBLE(N2))
!$OMP PARALLEL
!$OMP DO
      DO 30 K=1,M2
!DIR$ VECTOR ALIGNED
        DO 10 J=1,M1
          W1(1,J,K)=DCOS(PX*DBLE(J-1)*DBLE(K-1))
          W1(2,J,K)=DSIN(PX*DBLE(J-1)*DBLE(K-1))
   10   CONTINUE
!DIR$ VECTOR ALIGNED
        DO 20 IR=1,N1/M1
          W3(1,K,IR)=DCOS(PX*DBLE(K-1)*DBLE(IR-1)*DBLE(M1))
          W3(2,K,IR)=DSIN(PX*DBLE(K-1)*DBLE(IR-1)*DBLE(M1))
   20   CONTINUE
   30 CONTINUE
      DO 60 IS=1,N2/M2
!DIR$ VECTOR ALIGNED
        DO 40 J=1,M1
          W2(1,J,IS)=DCOS(PX*DBLE(J-1)*DBLE(IS-1)*DBLE(M2))
          W2(2,J,IS)=DSIN(PX*DBLE(J-1)*DBLE(IS-1)*DBLE(M2))
   40   CONTINUE
!DIR$ VECTOR ALIGNED
        DO 50 IR=1,N1/M1
          W4(1,IR,IS)=DCOS(PX*DBLE(IR-1)*DBLE(M1)*DBLE(IS-1)*DBLE(M2))
          W4(2,IR,IS)=DSIN(PX*DBLE(IR-1)*DBLE(M1)*DBLE(IS-1)*DBLE(M2))
   50   CONTINUE
   60 CONTINUE
!$OMP END PARALLEL
      RETURN
      END
