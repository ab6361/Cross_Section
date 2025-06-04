#!/bin/bash

# This line initiates an if statement that checks the number of arguments passed to the script.
# The special variable $# holds the count of arguments provided. -ne operator stands for "not equal"
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <filename_without_extension>"
  exit 1
fi
# The $0 variable represents the name of the script itself. fi ends the if block
FILENAME="$1"

# The || operator is a logical OR operator in Bash. It is used here to handle the case where the cd command fails
cd src/ || { echo "Directory src/ not found. Ensure this script is located in mc-single-arm/"; exit 1; }

# Run the simulation with the desired input file
echo "$FILENAME" | ./mc_single_arm

# The exit 1 command will terminate the script with an exit status of 1, which is a common convention to indicate that an error has occurred
cd ../util/ntuple/ || { echo "Directory util/ntuple/ not found."; exit 1; }
# |: Takes the output of the command on the left and uses it as the input for the command on the right
echo "$FILENAME" | ./make_ntuple

# Navigate to worksim for final conversion to ROOT
cd ../../worksim || { echo "Directory worksim not found."; exit 1; }
# $ works like f in python
h2root "${FILENAME}.rzdat"
