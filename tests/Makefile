F77 = g77
FFLAGS = -O3 -ffast-math -fomit-frame-pointer -malign-double -I..

all: test1d test2d test3d speed1d speed2d speed3d

test1d : test1d.o zfft1d.o kernel.o
	$(F77) test1d.o zfft1d.o kernel.o -o test1d

test2d : test2d.o zfft2d.o kernel.o
	$(F77) test2d.o zfft2d.o kernel.o -o test2d

test3d : test3d.o zfft3d.o kernel.o
	$(F77) test3d.o zfft3d.o kernel.o -o test3d

test1d.o : test1d.f ../param.h
	$(F77) $(FFLAGS) -c test1d.f -o test1d.o

test2d.o : test2d.f ../param.h
	$(F77) $(FFLAGS) -c test2d.f -o test2d.o

test3d.o : test3d.f ../param.h
	$(F77) $(FFLAGS) -c test3d.f -o test3d.o

speed1d : speed1d.o zfft1d.o kernel.o
	$(F77) speed1d.o zfft1d.o kernel.o -o speed1d

speed2d : speed2d.o zfft2d.o kernel.o
	$(F77) speed2d.o zfft2d.o kernel.o -o speed2d

speed3d : speed3d.o zfft3d.o kernel.o
	$(F77) speed3d.o zfft3d.o kernel.o -o speed3d

speed1d.o : speed1d.f ../param.h
	$(F77) $(FFLAGS) -c speed1d.f -o speed1d.o

speed2d.o : speed2d.f ../param.h
	$(F77) $(FFLAGS) -c speed2d.f -o speed2d.o

speed3d.o : speed3d.f ../param.h
	$(F77) $(FFLAGS) -c speed3d.f -o speed3d.o

zfft1d.o : ../zfft1d.f ../param.h
	$(F77) $(FFLAGS) -c ../zfft1d.f -o zfft1d.o

zfft2d.o : ../zfft2d.f ../param.h
	$(F77) $(FFLAGS) -c ../zfft2d.f -o zfft2d.o

zfft3d.o : ../zfft3d.f ../param.h
	$(F77) $(FFLAGS) -c ../zfft3d.f -o zfft3d.o

kernel.o : ../kernel.f
	$(F77) $(FFLAGS) -c ../kernel.f -o kernel.o

clean:
	rm -f *.o
