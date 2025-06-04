#!/bin/bash

# This line initiates an if statement that checks the number of arguments passed to the script.
# The special variable $# holds the count of arguments provided. -ne operator stands for "not equal"
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <filename_without_extension>"
  exit 1
fi
# The $0 variable represents the name of the script itself. fi ends the if block
FILENAME="$1"
# |: Takes the output of the command on the left and uses it as the input for the command on the right
# Runs the executable simc and provides the infile name
echo "$FILENAME" | ./simc

# The || operator is a logical OR operator in Bash. It is used here to handle the case where the cd command fails
# The exit 1 command will terminate the script with an exit status of 1, which is a common convention to indicate that an error has occurred
cd util/ntuple/ || { echo "Directory util/ntuple/ not found."; exit 1; }

echo "$FILENAME" | ./make_ntuple
cd ../../worksim || { echo "Directory worksim not found."; exit 1; }

# $ works like f in python
h2root "${FILENAME}.rzdat"
