#!/bin/bash
# script to clean files from a nautilus simulation
# version 1.0

function clean {
    rm output_1D.*
    rm rates1D.*
    
    # To delete the stderr and stdout of a bash scheduler of the server. 
    # The last "." is very important, in order to avoid suppression of the submission script itself.
    rm *.sh.* 
}





echo "deleting files..."
clean
echo "done"


