from Bio import Entrez, SeqIO
import os
import subprocess
import numpy as np
import glob


def create_project_directory(directory_name):
    # Create project directory if it doesn't exist
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)

def retrieve_genbank(accession):
    Entrez.email = "kayla.jarm@example.com"  # Replace with your email address
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    return handle.read()

def extract_cds_features(genbank_data):
    cds_features = []
    for record in SeqIO.parse(genbank_data, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                cds_features.append(feature)
    return cds_features

def write_fasta_file(cds_features, output_file, record):
    with open(output_file, "w") as fasta_out:
        for feature in cds_features:
            protein_id = feature.qualifiers.get("protein_id", ["unknown_protein_id"])[0]
            sequence = feature.location.extract(record.seq)
            fasta_out.write(f">{protein_id}\n{sequence}\n")

def build_kallisto_index(fasta_file, index_name):
    os.system(f"kallisto index -i {index_name} {fasta_file}")

def write_to_log(log_file, num_cds):
    with open(log_file, "a") as log:
        log.write(f"The HCMV genome (NC_006273.2) has {num_cds} CDS.\n")

def process_fastq_files(fastq_files, srr_numbers, output_directory, index_name, log_file):
    os.makedirs(output_directory, exist_ok=True)  # Ensure output directory exists
    
    # Write header to log file
    # Write header to log file with the desired label for condition
    with open(log_file, "a") as log:
        log.write("sample\texperimental_condition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")


    for fastq_file, srr_number in zip(fastq_files, srr_numbers):
        # Run kallisto quantification
        output_file = run_kallisto_quant(fastq_file, output_directory, index_name)
        
        # Check if kallisto quantification was successful
        if output_file is None:
            print(f"Skipping processing of {fastq_file} due to kallisto quantification failure")
            continue
        
        # Extract abundance information
        abundances = extract_abundance(output_file)
        
        # Calculate TPM statistics
        tpm_values = [float(tpm) for tpm in abundances.values()]
        min_tpm = np.min(tpm_values)
        med_tpm = np.median(tpm_values)
        mean_tpm = np.mean(tpm_values)
        max_tpm = np.max(tpm_values)
        
        # Write TPM statistics to log file
        with open(log_file, "a") as log:
            log.write(f"{srr_number}\tcondition\t{min_tpm}\t{med_tpm}\t{mean_tpm}\t{max_tpm}\n")
        
        print(f"Processed {fastq_file} (SRR: {srr_number}) and calculated TPM statistics")



def run_kallisto_quant(fastq_file, output_directory, index_name):
    # Construct output file path
    output_file = os.path.join(output_directory, "abundance.tsv")

    # Run kallisto quant with single-end flag
    result = subprocess.run(["kallisto", "quant", "-i", index_name, "-o", output_directory, "--single", "-l", "200", "-s", "20", fastq_file])

    
    # Check if kallisto quantification was successful
    if result.returncode != 0:
        print(f"Error: kallisto quant failed for {fastq_file}")
        return None

    # Return path to abundance file
    return output_file


def extract_abundance(output_file):
    abundances = {}
    with open(output_file, "r") as f:
        next(f)  # Skip header line
        for line in f:
            fields = line.strip().split("\t")
            transcript_id = fields[0]
            abundance = fields[4]  # Adjust this index according to kallisto output format
            abundances[transcript_id] = abundance
    return abundances

def write_abundance_tsv(abundances, output_file):
    with open(output_file, "w") as f:
        f.write("Transcript ID\tAbundance\n")
        for transcript_id, abundance in abundances.items():
            f.write(f"{transcript_id}\t{abundance}\n")


import glob

def main():
    # Define project directory and filenames
    project_directory = "PipelineProject_Kayla_Jarm"
    log_file = os.path.join(project_directory, "PipelineProject.log")
    accession = "NC_006273.2"
    genbank_filename = os.path.join(project_directory, f"{accession}.gb")
    fasta_filename = os.path.join(project_directory, f"{accession}.fasta")
    index_name = os.path.join(project_directory, f"{accession}_index.idx")

    # Define pattern to match fastq files
    fastq_pattern = "/Users/kaylajarm/Desktop/Python-Pipeline-Kayla-Jarm/SRR*.fastq"

    # Retrieve list of fastq files using glob
    fastq_files = sorted(glob.glob(fastq_pattern))

    srr_numbers = []
    for fastq_file in fastq_files:
        srr_number = os.path.basename(fastq_file).split("_")[0]
        srr_numbers.append(srr_number)

    fastq_output_directory = project_directory  # You can change this if you want to save output files in a different directory

    # Create project directory
    create_project_directory(project_directory)

    # Retrieve GenBank file
    genbank_data = retrieve_genbank(accession)
    with open(genbank_filename, "w") as gb_out:
        gb_out.write(genbank_data)

    # Extract CDS features
    cds_features = extract_cds_features(genbank_filename)

    # Get the first record for extracting sequences
    first_record = SeqIO.read(genbank_filename, "genbank")

    # Write CDS features to FASTA file
    write_fasta_file(cds_features, fasta_filename, first_record)

    # Build kallisto index
    build_kallisto_index(fasta_filename, index_name)

    # Write to log file
    num_cds = len(cds_features)
    write_to_log(log_file, num_cds)

    # Process FASTQ files
    process_fastq_files(fastq_files, srr_numbers, fastq_output_directory, index_name, log_file)

    print("Transcriptome index built successfully.")

if __name__ == "__main__":
    main()
