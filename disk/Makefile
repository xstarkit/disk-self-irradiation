# $Id: Makefile.def,v 1.3 2005/07/18 12:58:34 bursa Exp $



CC = gcc
SDLIB = /home/bursa/projects/slimdisk
SIM5LIB = /home/bursa/sim5/lib

CFLAGS = -I$(SDLIB) -I$(SIM5LIB) -L$(SDLIB) -Wall -O3 
LFLAGS = $(CFLAGS) -lm -ldl -lgsl -lgslcblas


default: sd-radial-model


# compile object files as a "position independent code" 
%.o: %.c
	$(CC) $(CFLAGS) -fpic -c -o $@ $<

clean:
	rm -f *.o


sd-radial-model-src = \
    $(SDLIB)/fluxmodel.c \
    $(SIM5LIB)/sim5kerr.c \
    sd-radial-model.c

sd-radial-model-obj = $(sd-radial-model-src:.c=.o)

sd-radial-model: $(sd-radial-model-obj)
	$(CC) $(LFLAGS) -o $@ $(sd-radial-model-obj)



sd-opacity-src = \
    $(SDLIB)/fluxmodel.c \
    $(SIM5LIB)/sim5kerr.c \
    sd-opacity.c

sd-opacity-obj = $(sd-opacity-src:.c=.o)

sd-opacity: $(sd-opacity-obj)
	$(CC) $(LFLAGS) -o $@ $(sd-opacity-obj)


