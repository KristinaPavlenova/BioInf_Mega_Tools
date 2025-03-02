from scripts.filter_params import is_filter_gc
from scripts.filter_params import is_filter_length, is_filter_quality
import scripts.rw_fastq

from abc import ABC

class BiologicalSequence(ABC):
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        
    def __len__(self):
        return len(self.sequence)
    
    def __getitem__(self, index):
         return self.sequence[index]

    def __str__(self):
        return self.sequence
    
    def __repr__(self):
        return f"{self.__class__.__name__}: {self.sequence}"
    
    @staticmethod
    def check_alphabet(sequence, alphabet):
        return set(sequence).issubset(alphabet)
    

class NucleicAcidSequence(BiologicalSequence):
    complement_nucl = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}

    def complement(self):
        return self.__class__("".join([self.complement_nucl[i] for i in self.sequence]))
    
    def reverse(self):
        return self.__class__(self.sequence[::-1])
    
    def reverse_complement(self):
        return self.__class__(self.complement().reverse().sequence)
    
    
class DNASequence(NucleicAcidSequence):
    alphabet = {'A', 'T', 'G', 'C'}

    def transcribe(self):
        return RNASequence(self.sequence.replace('T', 'U'))
    

class RNASequence(NucleicAcidSequence):
    alphabet = {'A', 'G', 'C', 'U'}
    pass


class AminoAcidSequence(BiologicalSequence):
    alphabet = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}
    amino_acid_weights = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
        'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.18,
        'L': 131.18, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
        'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
    }
    
    def molecular_weight(self):
        return sum(self.amino_acid_weights.get(amino_acid, 0) for amino_acid in self.sequence)



def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | float | int = (0, 100),
    length_bounds: tuple | int = (0, 2**32),
    quality_threshold: float | int = 0,
) -> None:
    """
    The function filters DNA and RNA reads from fastq file by parameters:
    GC content in sequence, sequence length, average quality score encoding.
    All boundaries are within the parameters ranges.

    Args:

    input_fastq - name of fastq file to read
    output_fastq - name of fastq file to write (path: ./filtered/output_fastq)

    gc_bounds: percentage range for GC content in 'dna_rna_sequence';
    in the case gc_bounds = (a: float | int, b: float | int)
    percentage range equels (a, b), default = (0, 100);
    in the case gc_bounds = a: float | int percentage range equels (0, a)

    length_bounds: range for 'dna_rna_sequence' length;
    in the case length_bounds = (n0: int, n1: int)
    length range equels (n0, n1), default = (0, 2**32)
    in the case length_bounds = n1: int length range equels (0, n1)

    quality_threshold: lower limit of average read quality threshold
    for filter in the scale phred33,
    default = 0
    """
    filtered_seqs = dict()
    seqs = scripts.rw_fastq.read_fastq(input_fastq)

    for seq_id, (seq, qual) in seqs.items():
        if (
            is_filter_gc(seq, gc_bounds)
            and is_filter_length(seq, length_bounds)
            and is_filter_quality(qual, quality_threshold)
        ):
            filtered_seqs[seq_id] = (seq, qual)
    scripts.rw_fastq.write_fastq(output_fastq, filtered_seqs)
