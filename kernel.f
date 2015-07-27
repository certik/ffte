C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2003, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS,
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@is.tsukuba.ac.jp
C
C
C     RADIX-2, 3, 4, 5 AND 8 FFT KERNEL ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE FFT23458(A,B,W,N,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION IP(*)
C
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
C
      KEY=1
      J=1
      L=N
      M=1
      DO 10 K=1,KP8
        L=L/8
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,B,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,A,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
        END IF
        M=M*8
        J=J+L
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,B,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,A,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
        END IF
        M=M*5
        J=J+L
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,B,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,A,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
        END IF
        M=M*4
        J=J+L
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,B,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,A,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
        END IF
        M=M*3
        J=J+L
   40 CONTINUE
      IF (IP(1) .EQ. 1) THEN
        IF (KEY .GE. 0) THEN
          CALL FFT2(A,A,M)
        ELSE
          CALL FFT2(B,A,M)
        END IF
      END IF
      RETURN
      END
      SUBROUTINE FFT2(A,B,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,*),B(2,M,*)
C
      DO 10 I=1,M
        X0=A(1,I,1)
        Y0=A(2,I,1)
        X1=A(1,I,2)
        Y1=A(2,I,2)
        B(1,I,1)=X0+X1
        B(2,I,1)=Y0+Y1
        B(1,I,2)=X0-X1
        B(2,I,2)=Y0-Y1
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT3A(A,B,W,L)
      ELSE
        CALL FFT3B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT3A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,3,*),W(2,*)
C
C      PI=4.0D0*DATAN(1.0D0)
C      S60=DSIN(PI/3.0D0)
      S60=0.86602540378443865D0
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        X0=A(1,J,2)+A(1,J,3)
        Y0=A(2,J,2)+A(2,J,3)
        X1=A(1,J,1)-0.5D0*X0
        Y1=A(2,J,1)-0.5D0*Y0
        X2=S60*(A(2,J,2)-A(2,J,3))
        Y2=S60*(A(1,J,3)-A(1,J,2))
        B(1,1,J)=A(1,J,1)+X0
        B(2,1,J)=A(2,J,1)+Y0
        B(1,2,J)=WR1*(X1+X2)-WI1*(Y1+Y2)
        B(2,2,J)=WR1*(Y1+Y2)+WI1*(X1+X2)
        B(1,3,J)=WR2*(X1-X2)-WI2*(Y1-Y2)
        B(2,3,J)=WR2*(Y1-Y2)+WI2*(X1-X2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,3,*),W(2,*)
C
C      PI=4.0D0*DATAN(1.0D0)
C      S60=DSIN(PI/3.0D0)
      S60=0.86602540378443865D0
      DO 10 I=1,M
        X0=A(1,I,1,2)+A(1,I,1,3)
        Y0=A(2,I,1,2)+A(2,I,1,3)
        X1=A(1,I,1,1)-0.5D0*X0
        Y1=A(2,I,1,1)-0.5D0*Y0
        X2=S60*(A(2,I,1,2)-A(2,I,1,3))
        Y2=S60*(A(1,I,1,3)-A(1,I,1,2))
        B(1,I,1,1)=A(1,I,1,1)+X0
        B(2,I,1,1)=A(2,I,1,1)+Y0
        B(1,I,2,1)=X1+X2
        B(2,I,2,1)=Y1+Y2
        B(1,I,3,1)=X1-X2
        B(2,I,3,1)=Y1-Y2
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        DO 20 I=1,M
          X0=A(1,I,J,2)+A(1,I,J,3)
          Y0=A(2,I,J,2)+A(2,I,J,3)
          X1=A(1,I,J,1)-0.5D0*X0
          Y1=A(2,I,J,1)-0.5D0*Y0
          X2=S60*(A(2,I,J,2)-A(2,I,J,3))
          Y2=S60*(A(1,I,J,3)-A(1,I,J,2))
          B(1,I,1,J)=A(1,I,J,1)+X0
          B(2,I,1,J)=A(2,I,J,1)+Y0
          B(1,I,2,J)=WR1*(X1+X2)-WI1*(Y1+Y2)
          B(2,I,2,J)=WR1*(Y1+Y2)+WI1*(X1+X2)
          B(1,I,3,J)=WR2*(X1-X2)-WI2*(Y1-Y2)
          B(2,I,3,J)=WR2*(Y1-Y2)+WI2*(X1-X2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT4A(A,B,W,L)
      ELSE
        CALL FFT4B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT4A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,4,*),W(2,*)
C
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        X0=A(1,J,1)+A(1,J,3)
        Y0=A(2,J,1)+A(2,J,3)
        X1=A(1,J,1)-A(1,J,3)
        Y1=A(2,J,1)-A(2,J,3)
        X2=A(1,J,2)+A(1,J,4)
        Y2=A(2,J,2)+A(2,J,4)
        X3=A(2,J,2)-A(2,J,4)
        Y3=A(1,J,4)-A(1,J,2)
        B(1,1,J)=X0+X2
        B(2,1,J)=Y0+Y2
        B(1,3,J)=WR2*(X0-X2)-WI2*(Y0-Y2)
        B(2,3,J)=WR2*(Y0-Y2)+WI2*(X0-X2)
        B(1,2,J)=WR1*(X1+X3)-WI1*(Y1+Y3)
        B(2,2,J)=WR1*(Y1+Y3)+WI1*(X1+X3)
        B(1,4,J)=WR3*(X1-X3)-WI3*(Y1-Y3)
        B(2,4,J)=WR3*(Y1-Y3)+WI3*(X1-X3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,4,*),W(2,*)
C
      DO 10 I=1,M
        X0=A(1,I,1,1)+A(1,I,1,3)
        Y0=A(2,I,1,1)+A(2,I,1,3)
        X1=A(1,I,1,1)-A(1,I,1,3)
        Y1=A(2,I,1,1)-A(2,I,1,3)
        X2=A(1,I,1,2)+A(1,I,1,4)
        Y2=A(2,I,1,2)+A(2,I,1,4)
        X3=A(2,I,1,2)-A(2,I,1,4)
        Y3=A(1,I,1,4)-A(1,I,1,2)
        B(1,I,1,1)=X0+X2
        B(2,I,1,1)=Y0+Y2
        B(1,I,3,1)=X0-X2
        B(2,I,3,1)=Y0-Y2
        B(1,I,2,1)=X1+X3
        B(2,I,2,1)=Y1+Y3
        B(1,I,4,1)=X1-X3
        B(2,I,4,1)=Y1-Y3
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        DO 20 I=1,M
          X0=A(1,I,J,1)+A(1,I,J,3)
          Y0=A(2,I,J,1)+A(2,I,J,3)
          X1=A(1,I,J,1)-A(1,I,J,3)
          Y1=A(2,I,J,1)-A(2,I,J,3)
          X2=A(1,I,J,2)+A(1,I,J,4)
          Y2=A(2,I,J,2)+A(2,I,J,4)
          X3=A(2,I,J,2)-A(2,I,J,4)
          Y3=A(1,I,J,4)-A(1,I,J,2)
          B(1,I,1,J)=X0+X2
          B(2,I,1,J)=Y0+Y2
          B(1,I,3,J)=WR2*(X0-X2)-WI2*(Y0-Y2)
          B(2,I,3,J)=WR2*(Y0-Y2)+WI2*(X0-X2)
          B(1,I,2,J)=WR1*(X1+X3)-WI1*(Y1+Y3)
          B(2,I,2,J)=WR1*(Y1+Y3)+WI1*(X1+X3)
          B(1,I,4,J)=WR3*(X1-X3)-WI3*(Y1-Y3)
          B(2,I,4,J)=WR3*(Y1-Y3)+WI3*(X1-X3)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT5A(A,B,W,L)
      ELSE
        CALL FFT5B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT5A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,5,*),W(2,*)
C
C      PI=4.0D0*DATAN(1.0D0)
C      S72=DSIN(2.0D0*PI/5.0D0)
C      S3672=DSIN(PI/5.0D0)/S72
C      SQR5Q=DSQRT(5.0D0)*0.25D0
      S72=0.95105651629515357D0
      S3672=0.61803398874989485D0
      SQR5Q=0.55901699437494742D0
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        X0=A(1,J,2)+A(1,J,5)
        Y0=A(2,J,2)+A(2,J,5)
        X1=A(1,J,3)+A(1,J,4)
        Y1=A(2,J,3)+A(2,J,4)
        X2=S72*(A(1,J,2)-A(1,J,5))
        Y2=S72*(A(2,J,2)-A(2,J,5))
        X3=S72*(A(1,J,3)-A(1,J,4))
        Y3=S72*(A(2,J,3)-A(2,J,4))
        X4=X0+X1
        Y4=Y0+Y1
        X5=SQR5Q*(X0-X1)
        Y5=SQR5Q*(Y0-Y1)
        X6=A(1,J,1)-0.25D0*X4
        Y6=A(2,J,1)-0.25D0*Y4
        X7=X6+X5
        Y7=Y6+Y5
        X8=X6-X5
        Y8=Y6-Y5
        X9=Y2+S3672*Y3
        Y9=-X2-S3672*X3
        X10=S3672*Y2-Y3
        Y10=X3-S3672*X2
        B(1,1,J)=A(1,J,1)+X4
        B(2,1,J)=A(2,J,1)+Y4
        B(1,2,J)=WR1*(X7+X9)-WI1*(Y7+Y9)
        B(2,2,J)=WR1*(Y7+Y9)+WI1*(X7+X9)
        B(1,3,J)=WR2*(X8+X10)-WI2*(Y8+Y10)
        B(2,3,J)=WR2*(Y8+Y10)+WI2*(X8+X10)
        B(1,4,J)=WR3*(X8-X10)-WI3*(Y8-Y10)
        B(2,4,J)=WR3*(Y8-Y10)+WI3*(X8-X10)
        B(1,5,J)=WR4*(X7-X9)-WI4*(Y7-Y9)
        B(2,5,J)=WR4*(Y7-Y9)+WI4*(X7-X9)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,5,*),W(2,*)
C
C      PI=4.0D0*DATAN(1.0D0)
C      S72=DSIN(2.0D0*PI/5.0D0)
C      S3672=DSIN(PI/5.0D0)/S72
C      SQR5Q=DSQRT(5.0D0)*0.25D0
      S72=0.95105651629515357D0
      S3672=0.61803398874989485D0
      SQR5Q=0.55901699437494742D0
      DO 10 I=1,M
        X0=A(1,I,1,2)+A(1,I,1,5)
        Y0=A(2,I,1,2)+A(2,I,1,5)
        X1=A(1,I,1,3)+A(1,I,1,4)
        Y1=A(2,I,1,3)+A(2,I,1,4)
        X2=S72*(A(1,I,1,2)-A(1,I,1,5))
        Y2=S72*(A(2,I,1,2)-A(2,I,1,5))
        X3=S72*(A(1,I,1,3)-A(1,I,1,4))
        Y3=S72*(A(2,I,1,3)-A(2,I,1,4))
        X4=X0+X1
        Y4=Y0+Y1
        X5=SQR5Q*(X0-X1)
        Y5=SQR5Q*(Y0-Y1)
        X6=A(1,I,1,1)-0.25D0*X4
        Y6=A(2,I,1,1)-0.25D0*Y4
        X7=X6+X5
        Y7=Y6+Y5
        X8=X6-X5
        Y8=Y6-Y5
        X9=Y2+S3672*Y3
        Y9=-X2-S3672*X3
        X10=S3672*Y2-Y3
        Y10=X3-S3672*X2
        B(1,I,1,1)=A(1,I,1,1)+X4
        B(2,I,1,1)=A(2,I,1,1)+Y4
        B(1,I,2,1)=X7+X9
        B(2,I,2,1)=Y7+Y9
        B(1,I,3,1)=X8+X10
        B(2,I,3,1)=Y8+Y10
        B(1,I,4,1)=X8-X10
        B(2,I,4,1)=Y8-Y10
        B(1,I,5,1)=X7-X9
        B(2,I,5,1)=Y7-Y9
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        DO 20 I=1,M
          X0=A(1,I,J,2)+A(1,I,J,5)
          Y0=A(2,I,J,2)+A(2,I,J,5)
          X1=A(1,I,J,3)+A(1,I,J,4)
          Y1=A(2,I,J,3)+A(2,I,J,4)
          X2=S72*(A(1,I,J,2)-A(1,I,J,5))
          Y2=S72*(A(2,I,J,2)-A(2,I,J,5))
          X3=S72*(A(1,I,J,3)-A(1,I,J,4))
          Y3=S72*(A(2,I,J,3)-A(2,I,J,4))
          X4=X0+X1
          Y4=Y0+Y1
          X5=SQR5Q*(X0-X1)
          Y5=SQR5Q*(Y0-Y1)
          X6=A(1,I,J,1)-0.25D0*X4
          Y6=A(2,I,J,1)-0.25D0*Y4
          X7=X6+X5
          Y7=Y6+Y5
          X8=X6-X5
          Y8=Y6-Y5
          X9=Y2+S3672*Y3
          Y9=-X2-S3672*X3
          X10=S3672*Y2-Y3
          Y10=X3-S3672*X2
          B(1,I,1,J)=A(1,I,J,1)+X4
          B(2,I,1,J)=A(2,I,J,1)+Y4
          B(1,I,2,J)=WR1*(X7+X9)-WI1*(Y7+Y9)
          B(2,I,2,J)=WR1*(Y7+Y9)+WI1*(X7+X9)
          B(1,I,3,J)=WR2*(X8+X10)-WI2*(Y8+Y10)
          B(2,I,3,J)=WR2*(Y8+Y10)+WI2*(X8+X10)
          B(1,I,4,J)=WR3*(X8-X10)-WI3*(Y8-Y10)
          B(2,I,4,J)=WR3*(Y8-Y10)+WI3*(X8-X10)
          B(1,I,5,J)=WR4*(X7-X9)-WI4*(Y7-Y9)
          B(2,I,5,J)=WR4*(Y7-Y9)+WI4*(X7-X9)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT8A(A,B,W,L)
      ELSE
        CALL FFT8B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT8A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,8,*),W(2,*)
C
C      C21=DSQRT(0.5D0)
      C21=0.70710678118654752D0
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        WR5=WR2*WR3-WI2*WI3
        WI5=WR2*WI3+WI2*WR3
        WR6=WR3*WR3-WI3*WI3
        WI6=WR3*WI3+WR3*WI3
        WR7=WR3*WR4-WI3*WI4
        WI7=WR3*WI4+WI3*WR4
        X0=A(1,J,1)+A(1,J,5)
        Y0=A(2,J,1)+A(2,J,5)
        X1=A(1,J,1)-A(1,J,5)
        Y1=A(2,J,1)-A(2,J,5)
        X2=A(1,J,3)+A(1,J,7)
        Y2=A(2,J,3)+A(2,J,7)
        X3=A(2,J,3)-A(2,J,7)
        Y3=A(1,J,7)-A(1,J,3)
        U0=X0+X2
        V0=Y0+Y2
        U1=X0-X2
        V1=Y0-Y2
        X4=A(1,J,2)+A(1,J,6)
        Y4=A(2,J,2)+A(2,J,6)
        X5=A(1,J,2)-A(1,J,6)
        Y5=A(2,J,2)-A(2,J,6)
        X6=A(1,J,4)+A(1,J,8)
        Y6=A(2,J,4)+A(2,J,8)
        X7=A(1,J,4)-A(1,J,8)
        Y7=A(2,J,4)-A(2,J,8)
        U2=X4+X6
        V2=Y4+Y6
        U3=Y4-Y6
        V3=X6-X4
        B(1,1,J)=U0+U2
        B(2,1,J)=V0+V2
        B(1,5,J)=WR4*(U0-U2)-WI4*(V0-V2)
        B(2,5,J)=WR4*(V0-V2)+WI4*(U0-U2)
        B(1,3,J)=WR2*(U1+U3)-WI2*(V1+V3)
        B(2,3,J)=WR2*(V1+V3)+WI2*(U1+U3)
        B(1,7,J)=WR6*(U1-U3)-WI6*(V1-V3)
        B(2,7,J)=WR6*(V1-V3)+WI6*(U1-U3)
        U0=X1+C21*(X5-X7)
        V0=Y1+C21*(Y5-Y7)
        U1=X1-C21*(X5-X7)
        V1=Y1-C21*(Y5-Y7)
        U2=X3+C21*(Y5+Y7)
        V2=Y3-C21*(X5+X7)
        U3=X3-C21*(Y5+Y7)
        V3=Y3+C21*(X5+X7)
        B(1,2,J)=WR1*(U0+U2)-WI1*(V0+V2)
        B(2,2,J)=WR1*(V0+V2)+WI1*(U0+U2)
        B(1,6,J)=WR5*(U1+U3)-WI5*(V1+V3)
        B(2,6,J)=WR5*(V1+V3)+WI5*(U1+U3)
        B(1,4,J)=WR3*(U1-U3)-WI3*(V1-V3)
        B(2,4,J)=WR3*(V1-V3)+WI3*(U1-U3)
        B(1,8,J)=WR7*(U0-U2)-WI7*(V0-V2)
        B(2,8,J)=WR7*(V0-V2)+WI7*(U0-U2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,8,*),W(2,*)
C
C      C21=DSQRT(0.5D0)
      C21=0.70710678118654752D0
      DO 10 I=1,M
        X0=A(1,I,1,1)+A(1,I,1,5)
        Y0=A(2,I,1,1)+A(2,I,1,5)
        X1=A(1,I,1,1)-A(1,I,1,5)
        Y1=A(2,I,1,1)-A(2,I,1,5)
        X2=A(1,I,1,3)+A(1,I,1,7)
        Y2=A(2,I,1,3)+A(2,I,1,7)
        X3=A(2,I,1,3)-A(2,I,1,7)
        Y3=A(1,I,1,7)-A(1,I,1,3)
        U0=X0+X2
        V0=Y0+Y2
        U1=X0-X2
        V1=Y0-Y2
        X4=A(1,I,1,2)+A(1,I,1,6)
        Y4=A(2,I,1,2)+A(2,I,1,6)
        X5=A(1,I,1,2)-A(1,I,1,6)
        Y5=A(2,I,1,2)-A(2,I,1,6)
        X6=A(1,I,1,4)+A(1,I,1,8)
        Y6=A(2,I,1,4)+A(2,I,1,8)
        X7=A(1,I,1,4)-A(1,I,1,8)
        Y7=A(2,I,1,4)-A(2,I,1,8)
        U2=X4+X6
        V2=Y4+Y6
        U3=Y4-Y6
        V3=X6-X4
        B(1,I,1,1)=U0+U2
        B(2,I,1,1)=V0+V2
        B(1,I,5,1)=U0-U2
        B(2,I,5,1)=V0-V2
        B(1,I,3,1)=U1+U3
        B(2,I,3,1)=V1+V3
        B(1,I,7,1)=U1-U3
        B(2,I,7,1)=V1-V3
        U0=X1+C21*(X5-X7)
        V0=Y1+C21*(Y5-Y7)
        U1=X1-C21*(X5-X7)
        V1=Y1-C21*(Y5-Y7)
        U2=X3+C21*(Y5+Y7)
        V2=Y3-C21*(X5+X7)
        U3=X3-C21*(Y5+Y7)
        V3=Y3+C21*(X5+X7)
        B(1,I,2,1)=U0+U2
        B(2,I,2,1)=V0+V2
        B(1,I,6,1)=U1+U3
        B(2,I,6,1)=V1+V3
        B(1,I,4,1)=U1-U3
        B(2,I,4,1)=V1-V3
        B(1,I,8,1)=U0-U2
        B(2,I,8,1)=V0-V2
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        WR5=WR2*WR3-WI2*WI3
        WI5=WR2*WI3+WI2*WR3
        WR6=WR3*WR3-WI3*WI3
        WI6=WR3*WI3+WR3*WI3
        WR7=WR3*WR4-WI3*WI4
        WI7=WR3*WI4+WI3*WR4
        DO 20 I=1,M
          X0=A(1,I,J,1)+A(1,I,J,5)
          Y0=A(2,I,J,1)+A(2,I,J,5)
          X1=A(1,I,J,1)-A(1,I,J,5)
          Y1=A(2,I,J,1)-A(2,I,J,5)
          X2=A(1,I,J,3)+A(1,I,J,7)
          Y2=A(2,I,J,3)+A(2,I,J,7)
          X3=A(2,I,J,3)-A(2,I,J,7)
          Y3=A(1,I,J,7)-A(1,I,J,3)
          U0=X0+X2
          V0=Y0+Y2
          U1=X0-X2
          V1=Y0-Y2
          X4=A(1,I,J,2)+A(1,I,J,6)
          Y4=A(2,I,J,2)+A(2,I,J,6)
          X5=A(1,I,J,2)-A(1,I,J,6)
          Y5=A(2,I,J,2)-A(2,I,J,6)
          X6=A(1,I,J,4)+A(1,I,J,8)
          Y6=A(2,I,J,4)+A(2,I,J,8)
          X7=A(1,I,J,4)-A(1,I,J,8)
          Y7=A(2,I,J,4)-A(2,I,J,8)
          U2=X4+X6
          V2=Y4+Y6
          U3=Y4-Y6
          V3=X6-X4
          B(1,I,1,J)=U0+U2
          B(2,I,1,J)=V0+V2
          B(1,I,5,J)=WR4*(U0-U2)-WI4*(V0-V2)
          B(2,I,5,J)=WR4*(V0-V2)+WI4*(U0-U2)
          B(1,I,3,J)=WR2*(U1+U3)-WI2*(V1+V3)
          B(2,I,3,J)=WR2*(V1+V3)+WI2*(U1+U3)
          B(1,I,7,J)=WR6*(U1-U3)-WI6*(V1-V3)
          B(2,I,7,J)=WR6*(V1-V3)+WI6*(U1-U3)
          U0=X1+C21*(X5-X7)
          V0=Y1+C21*(Y5-Y7)
          U1=X1-C21*(X5-X7)
          V1=Y1-C21*(Y5-Y7)
          U2=X3+C21*(Y5+Y7)
          V2=Y3-C21*(X5+X7)
          U3=X3-C21*(Y5+Y7)
          V3=Y3+C21*(X5+X7)
          B(1,I,2,J)=WR1*(U0+U2)-WI1*(V0+V2)
          B(2,I,2,J)=WR1*(V0+V2)+WI1*(U0+U2)
          B(1,I,6,J)=WR5*(U1+U3)-WI5*(V1+V3)
          B(2,I,6,J)=WR5*(V1+V3)+WI5*(U1+U3)
          B(1,I,4,J)=WR3*(U1-U3)-WI3*(V1-V3)
          B(2,I,4,J)=WR3*(V1-V3)+WI3*(U1-U3)
          B(1,I,8,J)=WR7*(U0-U2)-WI7*(V0-V2)
          B(2,I,8,J)=WR7*(V0-V2)+WI7*(U0-U2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL(W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(*)
      DIMENSION IP(3)
C
      CALL FACTOR(N,IP)
C
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
C
      J=1
      L=N
      DO 10 K=1,KP8
        L=L/8
        CALL SETTBL2(W(J),L,8)
        J=J+L
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        CALL SETTBL2(W(J),L,5)
        J=J+L
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        CALL SETTBL2(W(J),L,4)
        J=J+L
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        CALL SETTBL2(W(J),L,3)
        J=J+L
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL2(W,L,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/DBLE(L*M)
      DO 10 I=1,L
        W(1,I)=DCOS(PX*DBLE(I-1))
        W(2,I)=DSIN(PX*DBLE(I-1))
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FACTOR(N,IP)
      DIMENSION IP(*)
C
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND.
     1    MOD(N,5) .NE. 0) RETURN
   10 IF (N2 .LE. 1) RETURN
      IF (MOD(N2,2) .EQ. 0) THEN
        IP(1)=IP(1)+1
        N2=N2/2
        GO TO 10
      ELSE IF (MOD(N2,3) .EQ. 0) THEN
        IP(2)=IP(2)+1
        N2=N2/3
        GO TO 10
      ELSE IF (MOD(N2,5) .EQ. 0) THEN
        IP(3)=IP(3)+1
        N2=N2/5
        GO TO 10
      END IF
      RETURN
      END
