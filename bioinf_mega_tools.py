from abc import ABC
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os, os.path
import click
import logging


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

    def check_alphabet(self):
        return set(self.sequence).issubset(self.alphabet)


class NucleicAcidSequence(BiologicalSequence):
    complement_nucl = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}

    def complement(self):
        return self.__class__("".join([self.complement_nucl[i] for i in self.sequence]))

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    alphabet = {"A", "T", "G", "C"}

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    alphabet = {"A", "G", "C", "U"}
    pass


class AminoAcidSequence(BiologicalSequence):
    alphabet = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}
    amino_acid_weights = {
        "A": 89.09,
        "R": 174.20,
        "N": 132.12,
        "D": 133.10,
        "C": 121.15,
        "E": 147.13,
        "Q": 146.15,
        "G": 75.07,
        "H": 155.16,
        "I": 131.18,
        "L": 131.18,
        "K": 146.19,
        "M": 149.21,
        "F": 165.19,
        "P": 115.13,
        "S": 105.09,
        "T": 119.12,
        "W": 204.23,
        "Y": 181.19,
        "V": 117.15,
    }

    def molecular_weight(self):
        return sum(
            self.amino_acid_weights.get(amino_acid, 0) for amino_acid in self.sequence
        )


logging.basicConfig(
    filename='filter_fastq.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


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
    os.makedirs("./filtered/", exist_ok=True)

    logging.info(f"Starting filtering {input_fastq}.")

    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)
    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)
    try:
        with open(input_fastq) as input_handle:
            with open(os.path.join("./filtered/", output_fastq), "w") as output_handle:
                for record in SeqIO.parse(input_handle, "fastq"):
                    qualities = record.letter_annotations["phred_quality"]
                    if (
                        gc_bounds[0] <= gc_fraction(record.seq) * 100 <= gc_bounds[1]
                        and length_bounds[0] <= len(record.seq) <= length_bounds[1]
                        and sum(qualities) / len(qualities) >= quality_threshold
                    ):
                        SeqIO.write(record, output_handle, "fastq")
        logging.info(f"Filtering complete. Output saved to {os.path.join('./filtered/', output_fastq)}.")
    except Exception as e:
        logging.error(f"Error during filtering: {e}")
        raise



def parse_range(_, __, value):
    if value is None:
        return None
    if "," in value:
        a, b = map(float, value.split(","))
        return (a, b)
    return float(value)


@click.command()
@click.argument('input_fastq', type=click.Path(exists=True))
@click.argument('output_fastq')
@click.option('--gc_bounds', callback=parse_range, default="0,100")
@click.option('--length_bounds', callback=parse_range, default="0,4294967296")
@click.option('--quality_threshold', type=float, default=0)
def cli(input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold):
    """CLI wrapper for filter_fastq function"""
    filter_fastq(
        input_fastq=input_fastq,
        output_fastq=output_fastq,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold
    )


if __name__ == "__main__":
    cli()
