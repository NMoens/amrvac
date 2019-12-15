SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
TESTS := rhd_3d_2step_tvdlf_mm.log

#  rhd_3d_2step_tvdmu_al.log		\
# rhd_3d_3step_hll_cada.log rhd_3d_4step_hll_mc.log rhd_3d_4step_hllc_ko.log	\
# rhd_3d_rk4_tvdlf_cada.log

rhd_3d_2step_tvdlf_mm.log: rhd_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
# rhd_3d_2step_tvdmu_al.log: rhd_3d.par $(SCHEME_DIR)/2step_tvdmu_al.par
# rhd_3d_3step_hll_cada.log: rhd_3d.par $(SCHEME_DIR)/3step_hll_cada.par
# rhd_3d_4step_hll_mc.log: rhd_3d.par $(SCHEME_DIR)/4step_hll_mc.par
# rhd_3d_4step_hllc_ko.log: rhd_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
# rhd_3d_rk4_tvdlf_cada.log: rhd_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
