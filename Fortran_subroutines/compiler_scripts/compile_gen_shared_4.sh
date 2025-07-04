#!/bin/bash

# Exit on any error
set -e

# Define the source files as a space-separated list
SOURCE_FILES=(
  commondat1_mod.f90
  commondat_mod.f90
  final_str_mod.f90
  orb_mod.f90
  rum_rad_mod.f90
  check_bd_mod.f90
  quality_mod.f90
  loops_mod.f90
  infosymm_mod.f90
  factorial.f90
  comb.f90
  wigner.f90
  mod_vector_rep.f90
  mod_rumer_structures.f90
  mod_intra_bond_factor.f90
  mod_symm_break_factor.f90
  mod_nnat_bond_cal.f90
  mod_nnat_bond_cal_2.f90
  mod_prio_rad_str.f90
  mod_main_bond_cal.f90
  mod_quality_factor.f90
  mod_Rumer_set_id.f90
  mod_check_str_bond.f90
  mod_symm_loops.f90
  mod_coordination_val.f90
  mod_nnat_bond_sig.f90
  geocal.f90
  Invmat.f90
  mod_mat_ind.f90
  mod_eq_dst_check.f90
  mod_eq_dstr_set.f90
  mod_print_rumer.f90
  mod_All_Rumer_set.f90
  mod_write_rumer_xmi.f90
  mod_write_symm_xmi_1.f90
  mod_write_symm_xmi_new.f90
  mod_qult_str_arrange.f90
  mod_symmetry_cal_sig.f90
  mod_symmetry_cal_pi.f90
  mod_str_selection.f90
  mod_cov_struc.f90
  get_ctrl_inputs.f90
)

# Define module name for Python
MODULE_NAME="symm_str"

# Define the output directory
OUTPUT_DIR="/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface"

# Ensure MinGW-w64 cross-compiler is installed
if ! command -v x86_64-w64-mingw32-gfortran &> /dev/null; then
    echo "Error: x86_64-w64-mingw32-gfortran not found. Please install MinGW-w64."
    exit 1
fi


source /home/sourav/Desktop/environment/myenv/bin/activate
echo "Compiling Fortran source files into a Windows DLL..."

# Compile using MinGW-w64 Fortran compiler
x86_64-w64-mingw32-gfortran -shared -o "${MODULE_NAME}.dll" "${SOURCE_FILES[@]}"

# Check if the DLL was created
if [[ ! -f "${MODULE_NAME}.dll" ]]; then
  echo "Error: Compilation failed. DLL file not created."
  exit 1
fi

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Move the generated DLL file to the target directory
mv "${MODULE_NAME}.dll" "$OUTPUT_DIR"

echo "Compilation complete. Windows-compatible DLL file created: ${OUTPUT_DIR}/${MODULE_NAME}.dll"

