#----!!!! Operating system specific variables.  Override these by making a
#         file named makeincl with definitions of the variables
#         listed below. 

VERSION=1.5
# C compiler options
CC= gcc
COPT= -g -O3 -ffast-math # -pg
INCLUDE= -I./CVODE/include -I./Mathutils
INCLUDE= -I${SUNPREFIX}/include -I./Mathutils
DEFINE= -DAPPEND_UNDERSCORE -DNOBLAS_E -DRUN_NETWORK_VERSION=\"$(VERSION)\" 
CFLAGS =  $(COPT) $(INCLUDE) $(DEFINE)

# Sundial library directory
SUNPREFIX = ./sundials

# Linker options
LD= $(CC)
SUNDIALS_LIB_DIR = ${SUNPREFIX}/lib
SUNDIALS_LIBS    = -lsundials_cvode -lsundials_nvecserial -lsundials_shared

LIBS= -L. -lmathutils -L${SUNDIALS_LIB_DIR} ${SUNDIALS_LIBS} -lm
# Some code modifications may require additional blas and lapack libraries
# listed below.
#LIBS= -L. -lmathutils -lcvode -llapack -lblas -lm -lc

#----!!!! Below this point should be independent of operating system. 

# Output directory for executables
BINDIR=../bin

# Optional include file to override default variables
-include makeincl

OBJ= network.o 

EXEC=run_network

$(EXEC): %: %.o $(OBJ) 
	mkdir -p $(BINDIR)
	$(COMPILER_PREFIX) $(CC) $(CFLAGS) -c $@.c
	$(LINKER_PREFIX) $(CC) -o $(BINDIR)/$@ $@.o $(OBJ) $(LIBS) 

$(OBJ): libmathutils.a sundials

libmathutils.a:
	cd Mathutils; make CC="$(CC)" CFLAGS="$(CFLAGS)"

sundials:
	gzip -dc cvode_ser.tar.gz | tar xvf -
	cd $(SUNPREFIX); ./configure; make install

clean:
	rm -f *.o *.a
	rm -rf sundials

#----!!!! For making distributions of source code

FILES=LICENSE INSTALL Makefile network.h network.c Mathutils/ run_network.c HOW_TO_INSTALL README cvode_ser.tar.gz
SOURCES=LICENSE INSTALL Makefile network.h network.c run_network.c HOW_TO_INSTALL README

dist: $(FILES)
	tar czf dist.tgz $(FILES)

ci: $(FILES)
	ci -zLT -s- -u $(SOURCES)
	rcs -U $(SOURCES)
	chmod u+w $(SOURCES)

co: $(FILES)
	co -zLT -s- -u $(SOURCES)
