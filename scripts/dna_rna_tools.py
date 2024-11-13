TRANSCRIBE_NUCL = {"T": "U", "t": "u"}
COMPLEMENT_NUCL = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "U": "A",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "u": "a",
}


def is_dna(seq: str) -> bool:
    return set(seq.upper()).issubset({"A", "T", "G", "C"})


def is_rna(seq: str) -> bool:
    seq = seq.upper()
    return "U" in seq and set(seq).issubset({"A", "G", "C", "U"})


def transcribe(seq: str) -> str:
    if "T" in seq.upper() and is_dna(seq):
        for nucl in TRANSCRIBE_NUCL:
            seq = seq.replace(nucl, TRANSCRIBE_NUCL[nucl])
    return seq


def reverse(seq: str) -> str:
    if is_dna(seq) or is_rna(seq):
        result = seq[::-1]
    return result


def complement(seq: str) -> str:
    if is_dna(seq) or is_rna(seq):
        result = "".join([COMPLEMENT_NUCL[i] for i in seq])
    return result


def reverse_complement(seq: str) -> str:
    return reverse(complement(seq))
