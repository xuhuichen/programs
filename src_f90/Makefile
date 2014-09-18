#FC = mpif90 -mcmodel=medium -i-dynamic
# -check bound -check uninit
# FC = mpif90
FC = mpif90 -traceback

# FFLAGS = -O2
# FFLAGS = -O2 -axW -132

F1 = compton2d.f90 \
	reader.f90 \
	setup2d.f90 \
	xec2d.f90 \
	imcgen2d.f90 \
	volume2d.f90 \
	gamma1_2d.f90 \
	nontherm2d.f90 \
	imcsurf2d_para.f90 \
	planck2d.f90 \
	pp2d.f90 \
	imctrk2d.f90 \
	census2d.f90 \
	compb_2d.f90 \
	comtot2d.f90 \
	imcdate2d.f90 \
	ref_matrix.f90 \
	imcleak2d.f90 \
	graphics2d.f90 \
	imcvol2d_para.f90 \
	imcfield2d.f90 \
        imcredist.f90 \
	icloss2d.f90 \
	update2d.f90 \
	rand.f90 \
	fp_mpi.f90 \
	surf_mpi.f90 \
	vol_mpi.f90 \
        write_record.f90 \
        read_record.f90 

O1 = compton2d.o \
	reader.o \
	setup2d.o \
	xec2d.o \
	imcgen2d.o \
	volume2d.o \
	gamma1_2d.o \
	nontherm2d.o \
	imcsurf2d_para.o \
	planck2d.o \
	pp2d.o \
	imctrk2d.o \
	census2d.o \
	compb_2d.o \
	comtot2d.o \
	imcdate2d.o \
	ref_matrix.o \
	imcleak2d.o \
	graphics2d.o \
	imcvol2d_para.o \
	imcfield2d.o \
        imcredist.o \
	icloss2d.o \
	update2d.o \
	rand.o \
	fp_mpi.o \
	surf_mpi.o \
	vol_mpi.o \
        write_record.o \
        read_record.o 

compton2d : $(F1) $(O1)
	$(FC) -o compton2d $(O1)

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c -o $*.o $<

clean: 
	rm -f $(O1)

