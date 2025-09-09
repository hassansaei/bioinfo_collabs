#!/bin/bash
# Define the flags and arguments
while getopts "hd:i:x:" flag; do
 case $flag in
   h) # Handle the -h flag
   # Display script help information
   ;;
   d) # Handle the -d flag with an argument
   START_DIR=$OPTARG
   # Process the specified file
   ;;
   \?)
   # Handle invalid options
   ;;
   i) # Handel -i flag with an argument
   salmon_index=$OPTARG
   # Salmon index file
   ;;
   x) # Handel -x flag with an argument
   # Handel Hisat2 index file
   HISAT_IND=$OPTARG
   ;;
 esac
done

# Create new Project directory
current_date=$(date +"%Y-%m-%d")
main_dir="RNASeqProject_$current_date"

mkdir "$START_DIR/$main_dir"

# Move into the main directory
cd "$START_DIR/$main_dir" || exit
# Create subdirectories
dir1="0_Rawdata_Filtered"
dir2="1_Salmon_out"
dir3="2_DESeq2_out"

mkdir "$dir1"
mkdir "$dir2"
mkdir "$dir3"
cd "../../"
echo $PWD

# Find all R1 FASTQ files and assume each has an R2 pair
find $PWD/$START_DIR -name '*_1.fq.gz' | while read -r R1_FILE; do

  # Determine the corresponding R2 file based on the R1 file name
  R2_FILE="${R1_FILE/_1/_2}"
  # Check if the R2 file actually exists
  if [ -f "$R2_FILE" ]; then

    # Extract the sample name to be used for output
    OUTPUT_NAME=$(basename "${R1_FILE}")

    # Run fastp
    /root/tools/./fastp \
      -i "$R1_FILE" \
      -I "$R2_FILE" \
      -o "$PWD/$START_DIR"$main_dir"/"0_Rawdata_Filtered"/${OUTPUT_NAME}_R1_filtered.fastq.gz" \
      -O "$PWD/$START_DIR"$main_dir"/"0_Rawdata_Filtered"/${OUTPUT_NAME}_R2_filtered.fastq.gz" \
      -q 20 \
      -l 50 \
      -h "$PWD/$START_DIR/0_Rawdata_Filtered/${OUTPUT_NAME}_report.html" \
      -R "${OUTPUT_NAME}"

  else
    echo "Corresponding R2 file for $R1_FILE not found. Skipping."
  fi
done


# Salmon mapping and transcript abundance

DIR="$PWD/$START_DIR"$main_dir"/"0_Rawdata_Filtered""

# Find all R1 FASTQ files and assume each has an R2 pair
find $DIR -name '*_1.fq.gz' | while read -r R1_FILE; do

  # Determine the corresponding R2 file based on the R1 file name
  R2_FILE="${R1_FILE/_1/_2}"

  # Check if the R2 file actually exists
  if [ -f "$R2_FILE" ]; then

    # Extract the sample name to be used for output
    OUTPUT_NAME=$(basename "${R1_FILE}")

    # Run salmon
    /root/tools/salmon/bin/./salmon quant -i "$salmon_index" -l A -1 "$R1_FILE" -2 "$R2_FILE" -p 12 --validateMappings -o "$PWD/$START_DIR"$main_dir"/"1_Salmon_out"/"$OUTPUT_NAME""

  else
    echo "Corresponding R2 file not found. Skipping."
  fi
done

