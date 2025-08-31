#Update and Install the Necessary Tools
!apt-get update
!apt-get install -y sra-toolkit

!pip install biopython pandas matplotlib seaborn

!pip install Bio

from collections import Counter
from Bio import SeqIO
import pandas as pd

sra_accessions = ["SRR34219405"] # Replace with your SRA accessions

for accession in sra_accessions:
    # Download SRA file
    !prefetch {accession}
    # Convert SRA to FASTQ
    !fasterq-dump {accession} --split-files --outdir /content/fastq_files # Adjust output directory as needed


from Bio import SeqIO
import os

fastq_directory = "/content/fastq_files"
fastq_files = [f for f in os.listdir(fastq_directory) if f.endswith(".fastq")]

for fastq_file in fastq_files:
    file_path = os.path.join(fastq_directory, fastq_file)
    print(f"Loading {fastq_file}...")
    try:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # Process each FASTQ record (e.g., access record.id, record.seq, record.letter_annotations["phred_quality"])
                # Example: print(f"ID: {record.id}, Length: {len(record.seq)}")
                pass # Replace with your desired processing
        print(f"Successfully loaded {fastq_file}")
        except Exception as e:
            print(f"Error loading {fastq_file}: {e}")

#Download FASTQ Files Using fasterq-dump from SRA
#Sequence Read Archive

SRA_ID_1 = "SRR34219405"  # Example SRA ID 1
!fasterq-dump {SRA_ID_1} --split-files --threads 4

from Bio import SeqIO
import os

SRA_ID_1 = "SRR34219405"  # Example SRA ID 1
output_dir = "." # Assuming the fastq files are in the current directory, adjust if necessary

# List all fastq files for the given SRA ID
fastq_files = [f for f in os.listdir(output_dir) if f.startswith(SRA_ID_1) and f.endswith(".fastq")]

for fastq_file in fastq_files:
    fastq_path = os.path.join(output_dir, fastq_file)
    fasta_file = fastq_file.replace(".fastq", ".fasta")
    fasta_path = os.path.join(output_dir, fasta_file)

    print(f"Converting {fastq_file} to {fasta_file}...")

    try:
        with open(fastq_path, "r") as fastq_handle, open(fasta_path, "w") as fasta_handle:
            sequences = SeqIO.parse(fastq_handle, "fastq")
            SeqIO.write(sequences, fasta_handle, "fasta")
        print(f"Successfully converted {fastq_file} to {fasta_file}")
    except FileNotFoundError:
        print(f"Error: The file '{fastq_path}' was not found.")
    except Exception as e:
        print(f"An error occurred during conversion: {e}")

print(f"Conversion process for {SRA_ID_1} completed.")
