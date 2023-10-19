#
#
#
CFLAGS = $(GLOBALFLAGS) -O3 -g -I../include -flto
LIBS = -lm -lpthread $(MT_FLAG) -lpng -lz -lrt
FCOMPL = gfortran

EW_BINDIR = $(EW_HOME)/$(EW_VERSION)/bin
EW_LIBDIR = $(EW_HOME)/$(EW_VERSION)/lib
LC_LIBDIR = ../lib

#
#
#
EWLIBS = $(EW_LIBDIR)/lockfile_ew.o $(EW_LIBDIR)/lockfile.o $(EW_LIBDIR)/libew_mt.a

LCLIBS = $(LC_LIBDIR)/thrd_pool.o $(LC_LIBDIR)/raytracing.o

CPGPLOT_LIB = $(LC_LIBDIR)/libcpgplot.a $(LC_LIBDIR)/libpgplot.a

#
#
#
OBJS = \
	focalmech_ga.o \
	plot_focal.o \
	fpl_func.o \
	tko_azi_cal.o \
	ga_core.o

#
# Compile rule for focalmech_ga module
#
focalmech_ga: $(OBJS)
	@echo Creating $@...
	@$(FCOMPL) -o $(EW_BINDIR)/$@ $(OBJS) $(CPGPLOT_LIB) $(EWLIBS) $(LCLIBS) $(LIBS)

pfile2fevent: pfile2fevent.o
	@echo Creating $@...
	@$(CC) -o $(EW_BINDIR)/$@ pfile2fevent.o $(EWLIBS) $(LIBS)

#
#
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

clean_bin:
	@echo Removing binary execution file...
	@rm -f $(EW_BINDIR)/focalmech_ga*
