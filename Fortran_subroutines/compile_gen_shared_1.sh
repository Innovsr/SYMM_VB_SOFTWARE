#!/bin/bash

# Exit on any error
set -e

# Define the source files as a space-separated list
SOURCE_FILES=(
  commondat1.f90
  commondat.f90
  orb.f90
  mod_rum_rad.f90
  check_mod.f90
  quality.f90
  final_str.f90
  loops.f90
  infosymm.f90
  factorial.f90
  comb.f90
  close_file.f90
  wigner.f90
  sym_check.f90
  TwoSideJacobi.f90
  MatTran.f90
  mod_MatLDR.f90
  vector_rep.f90
  rumer_structures.f90
  main_bond.f90
  check_str_bond.f90
  mod_symm_loops.f90
  mod_coordination_val.f90
  nnat_bond_sig.f90
  geocal.f90
  lowercase.f90
  main_bond_str.f90
  Invmat.f90
  ind_matrix.f90
  eq_dst_check.f90
  eq_dstr_set.f90
  write_symm_xmi_1.f90
  write_symm_xmi_new.f90
  qult_str_arrange.f90
  intra_bond_fac.f90
  symm_break_fac.f90
  nnat_bond.f90
  nnat_bond_2.f90
  prio_rad.f90
  mod_quality_factor.f90
  Rum_set_id.f90
  mod_print_rumer.f90
  All_Rumer_set.f90
  write_rumer_xmi.f90
  symm_cal_sig.f90
  symm_cal_pi.f90
  str_select.f90
  cov_struc.f90
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
f2py -c -m "$MODULE_NAME" "${SOURCE_FILES[@]}"

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
