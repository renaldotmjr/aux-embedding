#FC = ifort
FC = gfortran

#FFLAGS = -O2 -r8
FFLAGS = 
#FFLAGS2 = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FFLAGS2 = 

OBJ_DIR = object.std

#INCLUDE = -I./include

OBJ =$(OBJ_DIR)/embparam.o \
	$(OBJ_DIR)/embtol.o \
	$(OBJ_DIR)/global_variables.o \
	$(OBJ_DIR)/inpdefvars.o \
	$(OBJ_DIR)/readinput.o \
	$(OBJ_DIR)/dealloc_emb.o \
	$(OBJ_DIR)/defaults.o \
	$(OBJ_DIR)/mathfun.o \
	$(OBJ_DIR)/unitconv.o \
	$(OBJ_DIR)/physcon.o \
	$(OBJ_DIR)/stringsman.o \
	$(OBJ_DIR)/findzatom.o \
	$(OBJ_DIR)/latconv.o \
	$(OBJ_DIR)/rho_pwscf.o \
	$(OBJ_DIR)/embgrid.o \
	$(OBJ_DIR)/mountzet.o \
	$(OBJ_DIR)/getbasis.o \
	$(OBJ_DIR)/genshl.o \
	$(OBJ_DIR)/embnorm.o \
	$(OBJ_DIR)/volume.o \
	$(OBJ_DIR)/cart2crys.o \
	$(OBJ_DIR)/crys2cart.o \
	$(OBJ_DIR)/crys2index.o \
	$(OBJ_DIR)/index2crys.o \
	$(OBJ_DIR)/putingrid.o \
	$(OBJ_DIR)/funchergauss.o \
	$(OBJ_DIR)/limit_s.o \
	$(OBJ_DIR)/limit_p.o \
	$(OBJ_DIR)/limit_d.o \
	$(OBJ_DIR)/calcnstep.o \
	$(OBJ_DIR)/findatomgrid.o \
	$(OBJ_DIR)/simpson.o \
	$(OBJ_DIR)/simpson_5.o \
	$(OBJ_DIR)/trapez.o \
	$(OBJ_DIR)/simplesum.o \
	$(OBJ_DIR)/intnumcall.o \
	$(OBJ_DIR)/hergauss_s.o \
	$(OBJ_DIR)/hergauss_p.o \
	$(OBJ_DIR)/hergauss_d.o \
	$(OBJ_DIR)/intauxrho.o \
	$(OBJ_DIR)/overlap_ss.o \
	$(OBJ_DIR)/overlap_ps.o \
	$(OBJ_DIR)/overlap_ds.o \
	$(OBJ_DIR)/overlap_fs.o \
	$(OBJ_DIR)/overlap_gs.o \
	$(OBJ_DIR)/overlap_sp.o \
	$(OBJ_DIR)/overlap_sd.o \
	$(OBJ_DIR)/overlap_pp.o \
	$(OBJ_DIR)/overlap_dp.o \
	$(OBJ_DIR)/overlap_dd.o \
	$(OBJ_DIR)/overlap_pd.o \
	$(OBJ_DIR)/normalization_s.o \
	$(OBJ_DIR)/drvoverlap.o \
	$(OBJ_DIR)/latticevec.o \
	$(OBJ_DIR)/calcbastrv.o \
	$(OBJ_DIR)/transatom.o \
	$(OBJ_DIR)/intoverlap.o \
	$(OBJ_DIR)/pythag.o \
	$(OBJ_DIR)/jacobi.o \
	$(OBJ_DIR)/reslinjac.o \
	$(OBJ_DIR)/reslinjac_NC.o \
	$(OBJ_DIR)/checksameatm.o \
	$(OBJ_DIR)/drvhgf.o \
	$(OBJ_DIR)/rhofit.o \
	$(OBJ_DIR)/fiterror.o \
	$(OBJ_DIR)/nelecfit.o \
	$(OBJ_DIR)/nelecfitatm.o \
	$(OBJ_DIR)/gencube.o \
	$(OBJ_DIR)/mounttvecs.o \
	$(OBJ_DIR)/printvec.o \
	$(OBJ_DIR)/printextembch.o \
	$(OBJ_DIR)/printdemoncube.o \
	$(OBJ_DIR)/printauxis.o \
	$(OBJ_DIR)/embcubreadauxis.o

SUB = modules/embparam.f90 \
	modules/embtol.f90 \
	modules/global_variables.f90 \
	modules/inpdefvars.f90 \
	util/readinput.f90 \
	util/dealloc_emb.f90 \
	util/defaults.f90 \
	util/mathfun.f90 \
	util/unitconv.f90 \
	util/physcon.f90 \
	util/stringsman.f90 \
	util/findzatom.f90 \
	lattice/latconv.f90 \
	util/rho_pwscf.f90 \
	util/embgrid.f90 \
	util/mountzet.f90 \
	util/getbasis.f90 \
	util/genshl.f90 \
	util/embnorm.f90 \
	util/volume.f90 \
	lattice/cart2crys.f90 \
	lattice/crys2cart.f90 \
	lattice/crys2index.f90 \
	lattice/index2crys.f90 \
	lattice/putingrid.f90 \
	hergauss/funchergauss.f90 \
	hergauss/limit_s.f90 \
	hergauss/limit_p.f90 \
	hergauss/limit_d.f90 \
	hergauss/calcnstep.f90 \
	hergauss/findatomgrid.f90 \
	hergauss/simpson.f90 \
	hergauss/simpson_5.f90 \
	hergauss/trapez.f90 \
	hergauss/simplesum.f90 \
	hergauss/intnumcall.f90 \
	hergauss/hergauss_s.f90 \
	hergauss/hergauss_p.f90 \
	hergauss/hergauss_d.f90 \
	hergauss/intauxrho.f90 \
	overlap/overlap_ss.f90 \
	overlap/overlap_ps.f90 \
	overlap/overlap_ds.f90 \
	overlap/overlap_fs.f90 \
	overlap/overlap_gs.f90 \
	overlap/overlap_sp.f90 \
	overlap/overlap_sd.f90 \
	overlap/overlap_pp.f90 \
	overlap/overlap_dp.f90 \
	overlap/overlap_dd.f90 \
	overlap/overlap_pd.f90 \
	overlap/normalization_s.f90 \
	overlap/drvoverlap.f90 \
	overlap/latticevec.f90 \
	overlap/calcbastrv.f90 \
	overlap/transatom.f90 \
	overlap/intoverlap.f90 \
	util/pythag.f90 \
	util/jacobi.f90 \
	util/reslinjac.f90 \
	util/reslinjac_NC.f90 \
	util/checksameatm.f90 \
	analysis/drvhgf.f90 \
	analysis/rhofit.f90 \
	analysis/fiterror.f90 \
	analysis/nelecfit.f90 \
	analysis/nelecfitatm.f90 \
	analysis/gencube.f90 \
	util/mounttvecs.f90 \
	prints/printvec.f90 \
	prints/printextembch.f90 \
	prints/printdemoncube.f90 \
	prints/printauxis.f90 \
	EMBdeMon/embcubreadauxis.f90

default:
	mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $(SUB) 
	mv *.o $(OBJ_DIR)
	$(FC) $(FFLAGS) $(OBJ) PWDE.f90 -o PWDE.x $(FFLAGS2)
	mv *.mod $(OBJ_DIR)
	mv PWDE.x bin/
