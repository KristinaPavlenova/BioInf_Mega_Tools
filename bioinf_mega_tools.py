from scripts.filter_params import filter_gc, filter_length, filter_quality
import scripts.dna_rna_tools

operations = {
    "transcribe": scripts.dna_rna_tools.transcribe,
    "reverse": scripts.dna_rna_tools.reverse,
    "complement": scripts.dna_rna_tools.complement,
    "reverse_complement": scripts.dna_rna_tools.reverse_complement,
}


def run_dna_rna_tools(*args: str):
    """
    The function returns DNA or RNA sequence
    after one of the following operations:
    - transcribe (works only for DNA)
    - reverse
    - complement
    - reverse_complement (reverse and complement together)

    Args:
    last one - operation name
    others - DNA/RNA sequence or sequences to operate
    (if the value is not DNA/RNA returns a warning, not a sequence)
    """
    seqs = args[:-1]
    operation = args[-1]

    for i in range(len(seqs)):
        is_dna = scripts.dna_rna_tools.is_dna
        is_rna = scripts.dna_rna_tools.is_rna
        if not is_dna(seqs[i]) and not is_rna(seqs[i]):
            seqs[i] = "Warning, input is not DNA/RNA"

    result = operations[operation](seqs)

    if len(result) == 1:
        return result[0]
    return result


def filter_fastq(
    seqs: dict,
    gc_bounds: tuple | float | int = (0, 100),
    length_bounds: tuple | int = (0, 2**32),
    quality_threshold: float | int = 0,
) -> dict:
    """
    The function filters DNA and RNA reads by parameters:
    GC content in sequence, sequence length, average quality score encoding.
    All boundaries are within the parameters ranges.

    Args:

    seqs: {'name': ('dna_rna_sequence', 'quality_score_encoding')}
    is a collection with reads

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
    result = dict()

    for key, value in seqs.items():
        if (
            filter_gc(value[0], gc_bounds)
            and filter_length(value[0], length_bounds)
            and filter_quality(value[1], quality_threshold)
        ):
            result[key] = value
    return result
