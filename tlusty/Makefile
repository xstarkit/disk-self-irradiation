default: t205


# compiler
CC = gcc

#-O3 ... optmization level
#-ggdb .. include debugging symbols
CFLAGS = -O3 -ggdb

# linker flags
LFLAGS = -lm $(CFLAGS)

.SUFFIXES:
.SUFFIXES: .c .o

t205:  clean-bin
	gfortran -w --no-automatic t205f.f -o t205 $(LFLAGS)

#run: compile
#	@rm -f fort.*
#	@echo 1 > fort.1
#	./t205 < t1.5 > t1.6
#	cp fort.141 t1.141
#
#run-ir: compile
#	@rm -f fort.*
#	@echo 1 > fort.1
#	cp t1.142 fort.2
#	./t205 < t1.5i > t1.6i
#	cp fort.141 t1.141i

clean-bin:
	@rm -f t205 a.out


