all:
	gfortran -c m_rdmd.F90
	gfortran -w -ffixed-form -fd-lines-as-comments -c m_lbfgs.F
	gfortran -c m_gammatodm2.F90
	gfortran -c m_E_grad_occ.F90
	gfortran -c m_optocc.F90
	gfortran -c m_noft_driver.F90
clean:
	/bin/rm -rf *.o *.mod
tar:
	tar -cvf module_noft.tar *F90 *F
