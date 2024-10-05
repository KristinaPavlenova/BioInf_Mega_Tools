transcribe_nucl = {"T": "U", "t": "u"}
complement_nucl = {
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


def transcribe(seqs: tuple) -> list:
    result = []
    for seq in seqs:
        if "T" in seq.upper() and is_dna(seq):
            for key in transcribe_nucl:
                seq = seq.replace(key, transcribe_nucl[key])
            result.append(seq)
        else:
            result.append(seq)
    return result


def reverse(seqs: tuple) -> list:
    result = []
    for seq in seqs:
        if is_dna(seq) or is_rna(seq):
            result.append(seq[::-1])
        else:
            result.append(seq)
    return result


def complement(seqs: tuple) -> list:
    result = []
    for seq in seqs:
        if is_dna(seq) or is_rna(seq):
            result.append("".join([complement_nucl[i] for i in seq]))
        else:
            result.append(seq)
    return result


def reverse_complement(seqs: tuple) -> list:
    return reverse(complement(seqs))
