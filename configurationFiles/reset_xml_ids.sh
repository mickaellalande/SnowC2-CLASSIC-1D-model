#!/bin/bash

# Input and output file paths
input_file="outputVariableDescriptors.xml"
output_file="outputVariableDescriptors_clean.xml"

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input XML file '$input_file' not found."
    exit 1
fi

# Create a temporary file for storing modified XML content
tmp_file=$(mktemp)

# Initialize a counter for the new id values
new_id=0

# Process the input XML file with awk and sed
awk -v new_id="$new_id" '
    {
        # Match and replace id="X" within <variable> elements
        if ($0 ~ /<variable id="[0-9]+"/) {
            # Extract the current id value
            match($0, /id="[0-9]+"/)
            current_id = substr($0, RSTART + 4, RLENGTH - 5)

            # Replace the current id with the new id
            gsub("id=\"" current_id "\"", "id=\"" new_id "\"")

            # Increment the new id for the next occurrence
            new_id++
        }
        # Print the modified or unmodified line
        print
    }
' "$input_file" > "$tmp_file"

# Move the temporary file to the output file
mv "$tmp_file" "$output_file"

echo "Updated XML content has been saved to '$output_file'."





