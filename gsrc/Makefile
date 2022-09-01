FC = gfortran
OBJ = detsirh.o stosirh.o ignbin.o mcmc.o # util.o
FFLAGS = -c -O3 
all: $(OBJ)


util.o:  util.f
	$(FC) $(FFLAGS) $(PKG_LIBS) util.f

mcmc.o:  mcmc.f
	$(FC) $(FFLAGS) $(PKG_LIBS) mcmc.f

detsirh.o:  detsirh.f
	$(FC) $(FFLAGS) $(PKG_LIBS) detsirh.f

stosirh.o:  stosirh.f
	$(FC) $(FFLAGS) $(PKG_LIBS) stosirh.f

ignbin.o: ignbin.f
	$(FC) $(FFLAGS) $(PKG_LIBS) ignbin.f






