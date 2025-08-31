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
