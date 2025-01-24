#!/bin/bash

# Define the list of targets
# targets=("C12" "Ca40" "Be9" "B10" "B11" "Ca48" "Ti48" "Fe54" "Ni58" "Ni64" "Ag108" "Th232" "Li6" "Li7" "Al27" "Cu63" "Au197" "He4")
targets=("C12" "Ca40")

# Loop through each target
for ntg in "${targets[@]}"; do
    # Run the notebook with the current target name as the parameter
    papermill sesaps_xsec.ipynb output_notebook_$ntg.ipynb -p ntg $ntg
done