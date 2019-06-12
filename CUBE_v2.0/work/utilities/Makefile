OPTIONS=-DFFTFINE
OPTIONS+=-Dsigma8
#OPTIONS+=-Danalysis
#OPTIONS+=-Ddiv3 -Dhackdsp

#OPTIONS+=-w2
OPTIONS+=-DPID
#OPTIONS+=-Dpositive
#OPTIONS+=-DEmode -DdspE

all: ic.x cicpower.x #cicrsd.x cicrsd1.x halofinder.x dsp.x lpt.x lpt1.x ang_mom_corr.x #vel.x qspace.x corr.x #convert.x
	@echo "done"

ic.x:    parameters.o variables.o pencil_fft.o initial_conditions.o powerspectrum.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
ic_nu.x: parameters.o variables.o pencil_fft.o initial_conditions_nu.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)

lpt.x: parameters.o variables.o pencil_fft.o lpt.o halofind.o powerspectrum.o indexedsort.o
	$(FC) $(XFLAG_NO_OMP) $(OPTIONS) $^ -o $@ $(FFTFLAG)
lpt1.x: parameters.o lpt1.o
	$(FC) $(XFLAG_NO_OMP) $(OPTIONS) $^ -o $@ $(FFTFLAG)
dsp.x: parameters.o pencil_fft.o powerspectrum.o displacement.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
cicpower.x: parameters.o variables.o pencil_fft.o powerspectrum.o cicpower.o hydrodf.o hydrogr.o hydro3d.o hydro1d.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
cicrsd.x: parameters.o cicrsd.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
cicrsd1.x: parameters.o cicrsd1.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
halofinder.x: parameters.o variables.o halofind.o indexedsort.o buffer_grid.o buffer_particle.o particle_initialization.o halofinder.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
ang_mom_corr.x: parameters.o variables.o ang_mom_corr.o indexedsort.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
vel.x:  parameters.o pencil_fft.o powerspectrum.o velocity_field.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
corr.x: parameters.o pencil_fft.o powerspectrum.o correlation_qspace.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@ $(FFTFLAG)
convert.x: parameters.o convert_zip_to_xv.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
qspace.x: parameters.o xva_qspace.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@

variables.o: parameters.o

parameters.o: ../main/parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
hydro1d.o: ../main/hydro1d.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
hydro3d.o: ../main/hydro3d.f90 hydro1d.o
	$(FC) $(OFLAG) $(OPTIONS) $<
hydrogr.o: ../main/hydrogr.f90 hydro3d.o
	$(FC) $(OFLAG) $(OPTIONS) $<
hydrodf.o: ../main/hydrodf.f90 hydrogr.o
	$(FC) $(OFLAG) $(OPTIONS) $<
variables.o: ../main/variables.f90 parameters.o hydrodf.o
	$(FC) $(OFLAG) $(OPTIONS) $<
neutrinos.o: ../main/neutrinos.f90 variables.o
	$(FC) $(OFLAG) $(OPTIONS) $<
particle_initialization.o: ../main/particle_initialization.f90 variables.o neutrinos.o
	$(FC) $(OFLAG) $(OPTIONS) $<
pencil_fft.o: ../main/pencil_fft.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
powerspectrum.o: powerspectrum.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)

initial_conditions.o:    initial_conditions.f90    parameters.o variables.o pencil_fft.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
initial_conditions_nu.o: initial_conditions_nu.f90 parameters.o variables.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)

lpt.o: lpt.f90 parameters.o variables.o pencil_fft.o powerspectrum.o
	$(FC) $(OFLAG_NO_OMP) $(OPTIONS) $< $(FFTFLAG)
lpt1.o: lpt1.f90 parameters.o
	$(FC) $(OFLAG_NO_OMP) $(OPTIONS) $< $(FFTFLAG)
displacement.o: displacement.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
xva_qspace.o: xva_qspace.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $<
correlation_qspace.o: correlation_qspace.f90 parameters.o pencil_fft.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
cicpower.o: cicpower.f90 parameters.o variables.o pencil_fft.o powerspectrum.o hydrodf.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
cicrsd.o: cicrsd.f90 parameters.o 
	$(FC) $(OFLAG) $(OPTIONS) $<
cicrsd1.o: cicrsd1.f90 parameters.o 
	$(FC) $(OFLAG) $(OPTIONS) $<
ang_mom_corr.o: ang_mom_corr.f90 parameters.o indexedsort.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
halofinder.o: halofinder.f90 parameters.o variables.o neutrinos.o halofind.o buffer_grid.o buffer_particle.o particle_initialization.o
	$(FC) $(OFLAG) $(OPTIONS) $<
halofind.o: ../main/halofind.f90 variables.o neutrinos.o indexedsort.o
	$(FC) $(OFLAG) $(OPTIONS) $<
indexedsort.o: ../main/indexedsort.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
buffer_grid.o: ../main/buffer_grid.f90 variables.o
	$(FC) $(OFLAG) $(OPTIONS) $<
buffer_particle.o: ../main/buffer_particle.f90 variables.o
	$(FC) $(OFLAG) $(OPTIONS) $<
velocity_field.o: velocity_field.f90 parameters.o pencil_fft.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
convert_zip_to_xv.o: convert_zip_to_xv.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $<
*.o: Makefile

clean:
	rm -f *.mod *.o *.out *.x *~
