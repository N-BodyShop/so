#
# Makefile for SO (Spherical Overdensity calculator).
#

#UNCOMMENT APPROPRIATE CFLAGS HERE
#Generic CFLAGS:
#CFLAGS	= -O2
#gcc CFLAGS:
CFLAGS = -O3 -funroll-loops
#CFLAGS = -ggdb
#R10k CFLAGS:
#CFLAGS = -O3 -mips4 -64 -r10000

CC=cc
LIBS	=   -lm

default: so

clean:
	rm -f *.o
	rm so

dist:
	cd ..; tar cvf so.tar so/*.c so/*.h so/Makefile; gzip so.tar	

so: so.o kd2.o smooth2.o cosmo.o romberg.o nr.o
	$(CC) $(CFLAGS) -o so so.o kd2.o smooth2.o cosmo.o romberg.o nr.o $(LIBS)

so.o: so.c kd2.h smooth2.h cosmo.h

kd2.o: kd2.c kd2.h tipsydefs.h cosmo.h

smooth2.o: smooth2.c kd2.h smooth2.h cosmo.h

nr.o: nr.c kd2.h

#
#	May need to specify -lrpc on some systems
#

