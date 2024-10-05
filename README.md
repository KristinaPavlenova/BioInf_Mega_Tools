# BioInf_Mega_Tools

The module for for working with DNA and RNA sequences. It includes 2 functions: *run_dna_rna_tools* and *filter_fastq*

## Function - *run_dna_rna_tools*

The function returns DNA/RNA sequence or sequences after one of the following operations made for input DNA/RNA sequence or sequences:
 - transcription (works only for DNA)
 - reversion
 - complementation
 - reversion and complementation together

## Function - *filter_fastq*

The function filters input DNA and RNA reads by parameters:
 - GC content in sequence
 - sequence length
 - average quality score encoding

## Running instructions

1. Clone repository to your work directory.
2. Import the module.
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
