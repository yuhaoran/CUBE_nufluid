OPTIONS+=-DPID
#OPTIONS+=-DNEUTRINOS
OPTIONS+=-DHALOFIND

#OPTIONS+=-Danalysis
#OPTIONS+=-DFORCETEST
#OPTIONS+=-Wall

EXTRA=
OPTIONS+=$(EXTRA)

MODFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o)))

all: main.x
	@echo "done"
main.x: $(OBJFILE)
	@echo "Link files:"
	$(FC) $(XFLAG) $(OPTIONS) $(OBJFILE) -o $@ $(FFTFLAG)

## For addtional subroutines, add its dependency on variables.o, and dependency of main.o on it.

hydro3d.o: hydro1d.o
hydrogr.o: hydro3d.o parameters.o
hydrodf.o: hydrogr.o

variables.o: parameters.o hydrodf.o
neutrinos.o: parameters.o
pencil_fft.o: parameters.o

buffer_grid.o: variables.o neutrinos.o
buffer_particle.o: variables.o neutrinos.o
update_particle.o: variables.o neutrinos.o
checkpoint.o: variables.o neutrinos.o
cubefft.o: variables.o
ext_pp_force.o: variables.o neutrinos.o
kernel_f.o: variables.o
kernel_c.o: variables.o
particle_initialization.o: variables.o neutrinos.o
timestep.o: variables.o neutrinos.o
halofind.o: variables.o neutrinos.o indexedsort.o

projection.o: variables.o

initialize.o: variables.o cubefft.o pencil_fft.o z_checkpoint.txt z_halofind.txt
finalize.o: variables.o cubefft.o pencil_fft.o
pm.o: variables.o neutrinos.o cubefft.o pencil_fft.o

main.o: parameters.o variables.o initialize.o ext_pp_force.o cubefft.o pencil_fft.o kernel_f.o kernel_c.o timestep.o particle_initialization.o buffer_grid.o buffer_particle.o pm.o update_particle.o checkpoint.o projection.o halofind.o finalize.o indexedsort.o hydrogr.o hydrodf.o

%.o: %.f90 Makefile
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)

clean:
	rm -f *.mod *.o *.out *.x *~
