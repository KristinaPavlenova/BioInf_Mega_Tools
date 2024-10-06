import os, os.path


def read_fastq(input_fastq: str) -> dict:
    seqs = dict()
    with open(input_fastq) as fastq:
        for line in fastq:
            if line.startswith("@") and " " in line:
                seq = fastq.readline()
                fastq.readline()
                scores = fastq.readline().strip()
                seqs[line] = (seq, scores)
    return seqs


def write_fastq(output_fastq: str, filtered_seqs: dict):
    if not os.path.exists("./filtered/"):
        os.mkdir("./filtered/")
    with open(os.path.join("./filtered/", output_fastq), "w") as fastq:
        for key, value in filtered_seqs.items():
            fastq.write(key)
            fastq.write(value[0])
            fastq.write("+" + key[1:])
            fastq.write(value[1] + "\n")
    pass
