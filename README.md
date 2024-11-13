# BioInf_Mega_Tools

The repo consists of tools for processing of bioinformatic objects.

It includes 2 modules:

1. The module **bioinf_mega_tools.py** for working with DNA and RNA sequences with 2 functions: *run_dna_rna_tools* and *filter_fastq*
2. The module **bio_files_processor.py** for working with bioinformatic files: .fasta, .gbk, .txt (with BLAST results). Functions *convert_multiline_fasta_to_oneline*, *parse_blast_output* and *select_genes_from_gbk_to_fasta* are included in the module

## Module bioinf_mega_tools.py
### Function - *run_dna_rna_tools*

The function returns DNA/RNA sequence or sequences after one of the following operations made for input DNA/RNA sequence or sequences:
 - transcription (works only for DNA)
 - reversion
 - complementation
 - reversion and complementation together

### Function - *filter_fastq*

The function filters DNA and RNA reads from fastq file by parameters:
 - GC content in sequence
 - sequence length
 - average quality score encoding

## Module bio_files_processor.py
### Function - *convert_multiline_fasta_to_oneline*

The function rewrites the sequences in fasta file from multiline format to oneline.

### Function - *parse_blast_output*

The function writes to txt file the descriptions of sequences with the highest identity percent for every query from txt file with BLAST results.

### Function - *select_genes_from_gbk_to_fasta*

The function selects neighboring genes of gene/genes of interest from gbk file and writes them to fasta file.

## Running instructions

1. Clone repository to your work directory.
2. Import the module/modules.
3. Check out the information about functions.
4. Call the required function with the required arguments.

## Example of use 

```python
import bioinf_mega_tools
bioinf_mega_tools.run_dna_rna_tools('ATG', 'transcribe') # Returns 'AUG'
bioinf_mega_tools.run_dna_rna_tools('ATG', 'reverse') # Returns 'GTA'
bioinf_mega_tools.run_dna_rna_tools('AtG', 'complement') # Returns 'TaC'
bioinf_mega_tools.run_dna_rna_tools('ATg', 'reverse_complement') # Returns 'cAT'
bioinf_mega_tools.run_dna_rna_tools('ATG', 'aT', 'reverse') # Returns ['GTA', 'Ta']
```

## System requirements
- Python 3.x
