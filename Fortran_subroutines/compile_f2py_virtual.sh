#!/bin/bash

# Exit on any error
set -e

# Define the source files
SOURCE_FILES="
commondat1.f90
commondat.f90
quality.f90
factorial.f90
new_row.f90
close_file.f90
wigner.f90
sym_check.f90
TwoSideJacobi.f90
MatTran.f90
Keycheck.f90
MatLDR.f90
vector_rep.f90
rumer_structures.f90
Rumer_set_id.f90
main_bond.f90
check_str_bond.f90
symm_loops.f90
coordination_val.f90
nnat_bond_sig.f90
geocal.f90
lowercase.f90
main_bond_str.f90

Invmat.f90
mat_ind.f90
eq_dst_check.f90
eq_dstr_set.f90
rumer.f90
All_Rumer_set.f90
write_rumer_xmi.f90
write_symm_xmi_1.f90
write_symm_xmi_new.f90
qult_str_arrange.f90
intra_bond_fac.f90
symm_break_fac.f90
nnat_bond.f90
nnat_bond_2.f90
prio_rad.f90
quality_fac.f90
symm_cal_sig.f90
symm_cal_pi.f90
str_select.f90
read_spl_key.f90
ol.f90
cov_struc.f90
read_info.f90
get_ctrl_inputs.f90
get_geometry_info.f90
get_orbs_info.f90
get_splkeywds.f90
get_input_data.f90
"

# Define the module name (used in Python imports)
MODULE_NAME="symm_str"

# Define the output directory
OUTPUT_DIR="/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface"

# Create the concatenated string of all source files for f2py
SOURCE_FILES_LIST=$(echo $SOURCE_FILES | tr '\n' ' ')

# Activate the virtual environment if required
# Uncomment the following line if you're using a virtual environment
source /home/sourav/Desktop/environment/myenv/bin/activate

echo "Compiling Fortran source files into a Python module using f2py..."

# Use f2py to compile all source files into a Python module
f2py -c -m $MODULE_NAME $SOURCE_FILES_LIST

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
# Uncomment the following line if a virtual environment was activated
deactivate

