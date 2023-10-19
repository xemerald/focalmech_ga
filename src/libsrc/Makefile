#
#    Revision history:
#
#     Revision 1.0  2018/03/19 17:38:40  Benjamin Yang
#     Initial version.
#
CFLAGS = $(GLOBALFLAGS) -O3 -g -I../../include -flto

LIB_DIR = ../../lib

LOCALSRCS = thrd_pool.c raytracing.c
LOCALOBJS = $(LOCALSRCS:%.c=%.o)

all: $(LOCALOBJS) pgplot
	@mv *.o $(LIB_DIR)

#
# PGPlot library
#
pgplot: PHONY
	@echo Checking the environment...
ifdef PGPLOT_DIR
	@echo PGPLOT_DIR is setted at $(PGPLOT_DIR)
	-@cd ./pgplot && \
		echo Making pgplot library... ; \
		mkdir -p build ; cd build ; \
		cp ../drivers.list . ; \
		../makemake .. linux gfortran_gcc ; \
		$(MAKE) ; \
		$(MAKE) cpg ; \
		mv *.a *.so ../../$(LIB_DIR) ; \
		mv grfont.dat rgb.txt $(PGPLOT_DIR) ; \
		mv *.h ../../../../include ; \
		cd ..
else
	$(error Need to setup the PGPLOT_DIR!!!)
endif

#
# Compile rule for Object
#
.c.o:
	@echo Compiling $<...
	@$(CC) $(CFLAGS) -c $<

#
# Clean-up rules
#
clean:
	@echo Cleaning build objects...
	@rm -f a.out core *.o *.obj *% *~
	@cd ./pgplot && rm -rf build\

clean_lib:
	@echo Removing libraries...
	@cd $(LIB_DIR) && \
		rm -f *.o *.a *.so

PHONY:
