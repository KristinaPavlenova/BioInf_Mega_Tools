# BioInf_Mega_Tools

The repo consists of tools for processing of bioinformatic objects.

It includes 2 modules:

1. The module **bioinf_mega_tools.py** for working with biological sequences:
classes *DNASequence*, *RNASequence*, *AminoAcidSequence* and function *filter_fastq*
2. The module **bio_files_processor.py** for working with bioinformatic files: .fasta, .gbk, .txt (with BLAST results). Functions *convert_multiline_fasta_to_oneline*, *parse_blast_output* and *select_genes_from_gbk_to_fasta* are included in the module

## Module bioinf_mega_tools.py
### Classes *DNASequence*, *RNASequence*, *AminoAcidSequence*

*DNASequence*, *RNASequence* implement the interface of parent class *NucleicAcidSequence* and abstract class *BiologicalSequence*. *DNASequence*, *RNASequence* have methods, which return instances of *DNASequence* and *RNASequence* classes after one of the following operations:
 - transcription (works only for *DNASequence*)
 - reversion
 - complementation
 - reversion and complementation together

 *AminoAcidSequence* implement the interface of abstract class *BiologicalSequence* and has method *molecular_weight*, which returns molecular weight of amino acid sequence.

 *BiologicalSequence* class has attribute *sequence* and method *check_alphabet* to check sequence alphabet for correctness.

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
2. Install Biopython (https://biopython.org/wiki/Download) for *filter_fastq* use.
3. Import the module/modules or classes.
4. Check out the information about functions/classes.

## Example of use 

```python
from bioinf_mega_tools import BiologicalSequence, DNASequence
my_dna = DNASequence('ATG')
my_dna.transcribe() # Returns 'AUG'
my_dna.reverse() # Returns 'GTA'
my_dna.complement() # Returns 'TAC'
my_dna.reverse_complement() # Returns 'CAT'
BiologicalSequence.check_alphabet(my_dna.sequence, my_dna.alphabet) # Returns True
```

## System requirements
- Python 3.x
