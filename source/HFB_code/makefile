# compiler 
#XLF= ifort -O3 -parallel 
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= Hf
EXE_PATH=../${EMPM_BIN_DIR}
#
all: technical.o math.o geom.o int_arr.o inter_stor_dd.o hp_field.o hn_field.o fermi_energy_n.o fermi_energy_p.o Dp_field.o Dn_field.o make_sp_levels.o make_qsp_levels.o dimm.o read_input.o HFB_energy.o check_densities.o initialize_basis.o square_root.o deallocate_all.o gauss_int.o integral_R.o occupations.o interaction_unformatted_dd.o hfb_iteration_dd.o diagonalization.o diag_subspace.o kinetic.o power2.o ln_matrix.o R_val.o exp_matrix.o make_densities.o transf_interaction_new2_dd.o TDA_new2_DD.o phonon_density_calc.o transit_calc.o transit_calc_spur.o spur_vec.o transition.o NatOrb.o  NatOrb_transf_unf.o HFB_energy_NAT.o transit_calc_NO.o transit_calc_spur_NO.o radial.o TDA_NO.o HFB.o    
	${XLF} technical.o math.o geom.o int_arr.o inter_stor_dd.o hp_field.o hn_field.o fermi_energy_n.o fermi_energy_p.o Dp_field.o Dn_field.o make_sp_levels.o make_qsp_levels.o dimm.o read_input.o HFB_energy.o check_densities.o initialize_basis.o square_root.o deallocate_all.o gauss_int.o integral_R.o occupations.o interaction_unformatted_dd.o hfb_iteration_dd.o diagonalization.o diag_subspace.o kinetic.o power2.o ln_matrix.o R_val.o exp_matrix.o make_densities.o transf_interaction_new2_dd.o TDA_new2_DD.o phonon_density_calc.o transit_calc.o transit_calc_spur.o spur_vec.o transition.o NatOrb.o  NatOrb_transf_unf.o HFB_energy_NAT.o transit_calc_NO.o transit_calc_spur_NO.o radial.o TDA_NO.o HFB.o ${OPT} -o ${EXE_PATH}${EXF}

technical.o: technical.f
	${XLF} -c technical.f

math.o: math.f
	${XLF} -c math.f

geom.o: geom.f
	${XLF} -c geom.f

int_arr.o: int_arr.f90
	${XLF} -c int_arr.f90

inter_stor_dd.o: inter_stor_dd.f90  
	${XLF} -c inter_stor_dd.f90

hp_field.o: hp_field.f
	${XLF} -c hp_field.f

hn_field.o: hn_field.f
	${XLF} -c hn_field.f

fermi_energy_n.o: fermi_energy_n.f
	${XLF} -c fermi_energy_n.f

fermi_energy_p.o: fermi_energy_p.f
	${XLF} -c fermi_energy_p.f

Dp_field.o: Dp_field.f
	${XLF} -c Dp_field.f

Dn_field.o: Dn_field.f
	${XLF} -c Dn_field.f

make_sp_levels.o: make_sp_levels.f
	${XLF} -c make_sp_levels.f

make_qsp_levels.o: make_qsp_levels.f
	${XLF} -c make_qsp_levels.f

dimm.o: dimm.f
	${XLF} -c dimm.f

read_input.o: read_input.f
	${XLF} -c read_input.f

HFB_energy.o: HFB_energy.f
	${XLF} -c HFB_energy.f

check_densities.o: check_densities.f
	${XLF} -c check_densities.f

initialize_basis.o: initialize_basis.f
	${XLF} -c initialize_basis.f

square_root.o: square_root.f
	${XLF} -c square_root.f

deallocate_all.o: deallocate_all.f
	${XLF} -c deallocate_all.f

gauss_int.o: gauss_int.f
	${XLF} -c gauss_int.f

integral_R.o: integral_R.f
	${XLF} -c integral_R.f

occupations.o: occupations.f
	${XLF} -c occupations.f

interaction_unformatted_dd.o: interaction_unformatted_dd.f
	${XLF} -c interaction_unformatted_dd.f

hfb_iteration_dd.o: hfb_iteration_dd.f
	${XLF} -c hfb_iteration_dd.f

diagonalization.o: diagonalization.f
	${XLF} -c diagonalization.f

diag_subspace.o: diag_subspace.f
	${XLF} -c diag_subspace.f

kinetic.o: kinetic.f
	${XLF} -c kinetic.f

power2.o: power2.f
	${XLF} -c power2.f 

ln_matrix.o: ln_matrix.f
	${XLF} -c ln_matrix.f

R_val.o: R_val.f
	${XLF} -c R_val.f

exp_matrix.o: exp_matrix.f
	${XLF} -c exp_matrix.f 

make_densities.o: make_densities.f
	${XLF} -c make_densities.f

transf_interaction_new2_dd.o: transf_interaction_new2_dd.f90
	${XLF} -c transf_interaction_new2_dd.f90

TDA_new2_DD.o: TDA_new2_DD.f
	${XLF} -c TDA_new2_DD.f

phonon_density_calc.o: phonon_density_calc.f
	${XLF} -c phonon_density_calc.f
  
transit_calc.o: transit_calc.f
	${XLF} -c transit_calc.f

transit_calc_spur.o: transit_calc_spur.f
	${XLF} -c transit_calc_spur.f

spur_vec.o: spur_vec.f
	${XLF} -c spur_vec.f

transition.o: transition.f
	${XLF} -c transition.f

NatOrb.o: NatOrb.f
	${XLF} -c NatOrb.f

NatOrb_transf_unf.o:  NatOrb_transf_unf.f90
	${XLF} -c  NatOrb_transf_unf.f90

HFB_energy_NAT.o: HFB_energy_NAT.f
	${XLF} -c  HFB_energy_NAT.f 

transit_calc_NO.o: transit_calc_NO.f
	${XLF} -c transit_calc_NO.f

transit_calc_spur_NO.o: transit_calc_spur_NO.f
	${XLF} -c transit_calc_spur_NO.f

radial.o : radial.f90
	${XLF} -c radial.f90

TDA_NO.o: TDA_NO.f
	${XLF} -c TDA_NO.f

HFB.o: HFB.f
	${XLF} -c HFB.f

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}    




	
