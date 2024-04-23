#!/bin/bash


# define In- and output locations through mounting
InputFileLocation="/mnt/c/Users/riber/Programming/masters/algorithms/project5/unique_distance_matrices"
OutputFileLocation="/mnt/c/Users/riber/Programming/masters/algorithms/project5/output"

###
### quicktree
###

cd ../quicktree-2.0/

# loop over the files in the input directory
for file in "$InputFileLocation"/*
do
    # check if the file is a regular file and not a ubunto id file
    if [ -f "$file" ] && [[ ! "$file" == *:Zone.identifier ]]; then
        # get the filename without the path
        filename=$(basename "$file")

        # run quicktree on the file
        quicktree -in m "$file" > "$OutputFileLocation/quicktree_$filename"
        
        # print status message
        echo "Ran quicktree on $filename"
    
    else
        # print message indicating that the file is skipped
        echo "Skipping file: $file"
    fi
done

echo "All files processed for quicktree."

###
### rapidNJ
###

# Move into rapidnj folder
cd ../rapidNJ-latest


# Loop over the files in the Windows directory
for file in "$InputFileLocation"/*
do
    # Check if the file is a regular file
    if [ -f "$file" ] && [[ ! "$file" == *:Zone.identifier ]]; then
        # Get the filename without the path
        filename=$(basename "$file")

        # Run quicktree on the file
        ./bin/rapidnj -i pd "$file" > "$OutputFileLocation/rapidnj_$filename"
        
        # Print status message
        echo "Ran rapidnj on $filename"
    
    else
        # Print message indicating that the file is skipped
        echo "Skipping file: $file"
    fi
done


echo "All files processed for rapidnj."