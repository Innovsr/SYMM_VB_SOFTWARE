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

# Ensure f2py is installed
if ! command -v f2py &> /dev/null; then
    echo "Error: f2py not found. Please install NumPy."
    exit 1
fi

# Activate virtual environment if it exists
VENV_PATH="/home/sourav/Desktop/environment/myenv/bin/activate"
if [ -f "$VENV_PATH" ]; then
    source "$VENV_PATH"
fi

echo "Compiling Fortran source files into a Python module using f2py..."

# Compile using f2py
f2py -c --lower -m "$MODULE_NAME" "${SOURCE_FILES[@]}"
#f2py -c -m "$MODULE_NAME" "${SOURCE_FILES[@]}"

# Check if the shared object file was created
SHARED_OBJECT=$(ls "${MODULE_NAME}".*.so 2>/dev/null || true)
if [[ -z "$SHARED_OBJECT" ]]; then
  echo "Error: Compilation failed. Shared object file not created."
  exit 1
fi

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Move the generated shared object file (.so) to the target directory
mv "$SHARED_OBJECT" "$OUTPUT_DIR"

echo "Compilation complete. Python-compatible shared object file created: ${OUTPUT_DIR}/${SHARED_OBJECT}"

# Deactivate virtual environment if it was activated
if [[ -n "$VIRTUAL_ENV" ]]; then
    deactivate
fi
