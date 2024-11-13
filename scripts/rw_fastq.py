import os, os.path


def read_fastq(input_fastq: str) -> dict:
    seqs = dict()
    with open(input_fastq) as fastq:
        for line in fastq:
            if line.startswith("@") and " " in line:
                seq = fastq.readline()
                fastq.readline()
                scores = fastq.readline().strip()
                seqs[line] = (seq, qual)
    return seqs


def write_fastq(output_fastq: str, filtered_seqs: dict) -> None:
    if not os.path.exists("./filtered/"):
        os.mkdir("./filtered/")
    with open(os.path.join("./filtered/", output_fastq), "w") as fastq:
        for seq_id, (seq, qual) in filtered_seqs.items():
            fastq.write(seq_id)
            fastq.write(seq)
            fastq.write("+" + seq_id[1:])
            fastq.write(qual + "\n")
