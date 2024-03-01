#installed homebrew and kallisto from terminal
from Bio import Entrez, SeqIO  # BioPython library for handling biological data
import os  # Operating System module for file system operations
import subprocess  # Subprocess module for running shell commands
import numpy as np  # NumPy library for numerical operations
import glob  #glob module for file path manipulation

#function to create a project directory
def create_project_directory(directory_name):
    #create project directory if it doesn't exist
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)

#retrieve GenBank data
def retrieve_genbank(accession):
    Entrez.email = "kayla.jarm@example.com" 
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    return handle.read()

#extract CDS features from GenBank data
def extract_cds_features(genbank_data):
    cds_features = []
    for record in SeqIO.parse(genbank_data, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":  # Check if feature type is CDS (Coding DNA Sequence)
                cds_features.append(feature)
    return cds_features

#write CDS features to a FASTA file
def write_fasta_file(cds_features, output_file, record):
    with open(output_file, "w") as fasta_out:
        for feature in cds_features:
            protein_id = feature.qualifiers.get("protein_id", ["unknown_protein_id"])[0]
            sequence = feature.location.extract(record.seq)
            fasta_out.write(f">{protein_id}\n{sequence}\n")

#build a Kallisto index from a FASTA file
def build_kallisto_index(fasta_file, index_name):
    os.system(f"kallisto index -i {index_name} {fasta_file}")

#write information to a log file
def write_to_log(log_file, num_cds):
    with open(log_file, "a") as log:
        log.write(f"The HCMV genome (NC_006273.2) has {num_cds} CDS.\n")

#process FASTQ files using kallisto quantification
def process_fastq_files(fastq_files, srr_numbers, output_directory, index_name, log_file):
    os.makedirs(output_directory, exist_ok=True)  # Ensure output directory exists
    
    # Write header to log file with the desired label for condition
    with open(log_file, "a") as log:
        log.write("sample\texperimental_condition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")

    #loop through each FASTQ file and corresponding SRR number
    for fastq_file, srr_number in zip(fastq_files, srr_numbers):
        #run kallisto quantification
        output_file = run_kallisto_quant(fastq_file, output_directory, index_name)
        
        #check if kallisto quantification was successful
        if output_file is None:
            print(f"Skipping processing of {fastq_file} due to kallisto quantification failure")
            continue
        
        #extract abundance information
        abundances = extract_abundance(output_file)
        
        #calculate TPM statistics
        tpm_values = [float(tpm) for tpm in abundances.values()]
        min_tpm = np.min(tpm_values)
        med_tpm = np.median(tpm_values)
        mean_tpm = np.mean(tpm_values)
        max_tpm = np.max(tpm_values)
        
        # Write TPM statistics to log file
        with open(log_file, "a") as log:
            log.write(f"{srr_number}\tcondition\t{min_tpm}\t{med_tpm}\t{mean_tpm}\t{max_tpm}\n")
        
        print(f"Processed {fastq_file} (SRR: {srr_number}) and calculated TPM statistics")

# Function to run kallisto quantification
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

# Function to extract abundance information from kallisto output
def extract_abundance(output_file):
    abundances = {}
    with open(output_file, "r") as f:
        next(f)  # Skip header line
        for line in f:
            fields = line.strip().split("\t")
            transcript_id = fields[0]
            abundance = fields[4] 
            abundances[transcript_id] = abundance
    return abundances

#write abundance information to a TSV file
def write_abundance_tsv(abundances, output_file):
    with open(output_file, "w") as f:
        f.write("Transcript ID\tAbundance\n")
        for transcript_id, abundance in abundances.items():
            f.write(f"{transcript_id}\t{abundance}\n")

def main():
    #define project directory and filenames
    project_directory = "PipelineProject_Kayla_Jarm"
    log_file = os.path.join(project_directory, "PipelineProject.log")  #path to log file
    accession = "NC_006273.2"  #accession number for the genome
    genbank_filename = os.path.join(project_directory, f"{accession}.gb")  #genBank file name
    fasta_filename = os.path.join(project_directory, f"{accession}.fasta")  #FASTA file name
    index_name = os.path.join(project_directory, f"{accession}_index.idx")  #index file name for kallisto

    #define pattern to match fastq files
    fastq_pattern = "/Users/kaylajarm/Desktop/Python-Pipeline-Kayla-Jarm/SRR*.fastq"  # Pattern for FASTQ files

    #retrieve list of fastq files using glob
    fastq_files = sorted(glob.glob(fastq_pattern))  # List of FASTQ files

    #extract SRR numbers from fastq file names
    srr_numbers = []
    for fastq_file in fastq_files:
        srr_number = os.path.basename(fastq_file).split("_")[0]  # Extract SRR number from file name
        srr_numbers.append(srr_number)

    fastq_output_directory = project_directory  #output directory for processed FASTQ files

    #create project directory if it doesn't exist
    create_project_directory(project_directory)

    #retrieve GenBank file
    genbank_data = retrieve_genbank(accession)  # Retrieve GenBank data
    with open(genbank_filename, "w") as gb_out:
        gb_out.write(genbank_data)  # Write GenBank data to file

    # Extract CDS features from GenBank file
    cds_features = extract_cds_features(genbank_filename)  # Extract CDS features

    # Get the first record from GenBank file for extracting sequences
    first_record = SeqIO.read(genbank_filename, "genbank")

    # Write CDS features to FASTA file
    write_fasta_file(cds_features, fasta_filename, first_record)  # Write CDS features to FASTA file

    # Build kallisto index from FASTA file
    build_kallisto_index(fasta_filename, index_name)  # Build kallisto index

    # Write information about the number of CDS features to log file
    num_cds = len(cds_features)  # Number of CDS features
    write_to_log(log_file, num_cds)  # Write to log file

    # Process FASTQ files using kallisto quantification
    process_fastq_files(fastq_files, srr_numbers, fastq_output_directory, index_name, log_file)  # Process FASTQ files

    print("Transcriptome index built successfully.")  # Print success message

if __name__ == "__main__":
    main()  # Call main function if script is executed directly

