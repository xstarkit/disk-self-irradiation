# $Id: Makefile.def,v 1.3 2005/07/18 12:58:34 bursa Exp $



CC = gcc
SIM5LIB = sim5lib

#SDLIB = /home/bursa/projects/slimdisk
#SIM5LIB = /home/bursa/projects/sim5/lib
#FLAGS = -I$(SDLIB) -I$(SIM5LIB) -L$(SDLIB) -L$(SIM5LIB) -Wall -O3  
CFLAGS = -I$(SIM5LIB) -L$(SIM5LIB) -Wall -O3 -w
LFLAGS = $(CFLAGS) -lm -ldl -lgsl -lgslcblas


default: disk-model


# compile object files as a "position independent code" 
%.o: %.c
	$(CC) $< $(CFLAGS) -fpic -c -o $@

clean:
	rm -f *.o


sd-radial-model-src = \
    $(SDLIB)/fluxmodel.c \
    $(SIM5LIB)/sim5kerr.c \
    sd-radial-model.c

sd-radial-model-obj = $(sd-radial-model-src:.c=.o)

sd-radial-model: $(sd-radial-model-obj)
	$(CC) $(sd-radial-model-obj) $(LFLAGS) -o $@ 





sd-opacity-src = \
    $(SDLIB)/fluxmodel.c \
    $(SIM5LIB)/sim5kerr.c \
    sd-opacity.c

sd-opacity-obj = $(sd-opacity-src:.c=.o)

sd-opacity: $(sd-opacity-obj)
	$(CC) $(LFLAGS) -o $@ $(sd-opacity-obj)






disk-model-src = \
    $(SIM5LIB)/sim5lib.c \
    disk-model.c

disk-model-obj = $(disk-model-src:.c=.o)

disk-model: $(disk-model-obj)
	$(CC) $(disk-model-obj) $(LFLAGS) -o $@ 





