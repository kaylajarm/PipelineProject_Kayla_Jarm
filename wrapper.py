import os
from Bio import SeqIO
import subprocess
import glob
import statistics


# Define your first and last name
first_name = "Kayla"
last_name = "Jarm"

# Create directory
output_dir = f"PipelineProject_{first_name}_{last_name}"
os.system(f"mkdir {output_dir}")
os.chdir(output_dir)

# Function to extract CDS features from GenBank file and generate FASTA file
def extract_cds(genbank_file, output_fasta):
    with open(genbank_file, "r") as gb_fh, open(output_fasta, "w") as fasta_fh:
        for record in SeqIO.parse(gb_fh, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    protein_id = feature.qualifiers.get("protein_id", ["Unknown"])[0]
                    fasta_fh.write(f">{protein_id}\n{feature.location.extract(record).seq}\n")

# Paths
genbank_file = "/Users/kaylajarm/Desktop/Python-Pipeline-Kayla-Jarm/NC_006273.2.gb"
output_fasta = "/Users/kaylajarm/Desktop/Python-Pipeline-Kayla-Jarm/HCMV_CDS.fasta"
log_file = "/Users/kaylajarm/Desktop/Python-Pipeline-Kayla-Jarm/PipelineProject.log"

# Step 2: Extract CDS features and generate FASTA file
extract_cds(genbank_file, output_fasta)

# Step 3: Count the number of CDS in the HCMV genome
num_cds = sum(1 for record in SeqIO.parse(genbank_file, "genbank") for feature in record.features if feature.type == "CDS")

# Step 4: Write number of CDS to log file
with open(log_file, "w") as log_fh:
    log_fh.write(f"The HCMV genome (NC_006273.2) has {num_cds} CDS.\n")

# Step 5: Build Kallisto index
kallisto_index_output = "HCMV_index.idx"
subprocess.run(["kallisto", "index", "-i", kallisto_index_output, output_fasta])

# Function to parse abundance.tsv file and calculate TPM statistics
def calculate_tpm_stats(abundance_file):
    tpm_values = []
    with open(abundance_file, "r") as file:
        next(file)  # Skip header
        for line in file:
            fields = line.strip().split("\t")
            tpm = float(fields[4])  # TPM value is at index 4
            tpm_values.append(tpm)
    min_tpm = min(tpm_values)
    med_tpm = statistics.median(tpm_values)
    mean_tpm = statistics.mean(tpm_values)
    max_tpm = max(tpm_values)
    return min_tpm, med_tpm, mean_tpm, max_tpm

# List of donor FASTQ files
donor_fastq_files = glob.glob("Donor*_*.fastq")  # Assumes fastq files are in the current directory

# Path to the log file
log_file_path = "/Users/kaylajarm/Desktop/Python-Pipeline-Kayla-Jarm/PipelineProject.log"

# Write header to log file
with open(log_file_path, "a") as log_file:  # Open in append mode ("a")
    log_file.write("sample condition min_tpm med_tpm mean_tpm max_tpm\n")

# Iterate over each donor FASTQ file
for fastq_file in donor_fastq_files:
    sample_condition = os.path.splitext(os.path.basename(fastq_file))[0]  # Extract sample condition from file name
    
    # Assuming the output directory for each sample is named after the sample condition
    abundance_file = "/Users/kaylajarm/Desktop/Python-Pipeline-Kayla-Jarm/abundance.tsv"
    
    if not os.path.exists(abundance_file):
        print(f"Abundance file not found for sample {sample_condition}. Skipping...")
        continue
    
    # Calculate TPM statistics
    min_tpm, med_tpm, mean_tpm, max_tpm = calculate_tpm_stats(abundance_file)
    
    # Append details to log file
    with open(log_file_path, "a") as log_file:  # Open in append mode ("a")
        log_file.write(f"{sample_condition} {min_tpm} {med_tpm} {mean_tpm} {max_tpm}\n")

print("TPM statistics calculation and logging complete.")
