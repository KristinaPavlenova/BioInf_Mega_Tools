def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = "x_oneline.fasta"
) -> None:
    """
    The function rewrites the sequences in fasta file from multiline format to oneline.

    Args:
    input_fasta - path of the input fasta file
    output_fasta - path of the output fasta file
    """
    if output_fasta == "x_oneline.fasta":
        output_fasta = input_fasta.split(".")[0] + output_fasta[1:]
    with open(input_fasta, "r") as in_fasta, open(output_fasta, "w") as out_fasta:
        for line in in_fasta:
            if not line.startswith(">"):
                seq = []
                while line.strip().isalpha():
                    seq.append(line.strip())
                    line = in_fasta.readline()
                out_fasta.write("".join(seq) + "\n")
            out_fasta.write(line)


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    The function writes to txt file the descriptions of sequences
    with the highest identity percent for every query
    from txt file with BLAST results.

    Args:
    input_file - path of the input txt file
    output_file - path of the output txt file
    """
    with open(input_file, "r") as blast_results, open(
        output_file, "w"
    ) as best_ident_discriptions:
        best_ident_list = []
        for line in blast_results:
            if "Sequences producing significant alignments:" in line:
                next_line = blast_results.readline()
                len_discription_col = len(next_line.split("Scientific")[0])
                blast_results.readline()
                line = blast_results.readline()
                best_ident_list.append(line[:len_discription_col].strip())
        best_ident_list.sort()
        for discription in best_ident_list:
            print(discription, file=best_ident_discriptions)


def select_genes_from_gbk_to_fasta(
    input_gbk: str, output_fasta: str, *args: str, n_before: int = 1, n_after: int = 1
) -> None:
    """
    The function selects neighboring genes of gene/genes of interest from gbk file
    and writes then to fasta file.

    Args:
    input_gbk - path of the input gbk file
    output_fasta - path of the output fasta file
    *args - gene or genes of interest
    n_before - number of genes before gene of interest to select (default = 1)
    n_after - number of genes after gene of interest to select (default = 1)
    """
    with open(input_gbk, "r") as input_gbk, open(output_fasta, "w") as output_fasta:
        genes = args
        all_genes = []
        for line in input_gbk:
            if '/gene="' in line:
                gene = line.split('"')[1]
                while '/translation="' not in line:
                    line = input_gbk.readline()
                seq = []
                seq.append(line.strip().split('"')[1])
                line = input_gbk.readline()
                while line.split('"')[0].strip().isalpha():
                    seq.append(line.split('"')[0].strip())
                    line = input_gbk.readline()
                all_genes.append((gene, "".join(seq)))
        for gene_of_interest in genes:
            for gene, seq in all_genes:
                if gene == gene_of_interest:
                    gene_index = all_genes.index((gene, seq))
            for i in range(n_before, 0, -1):
                if gene_index - i >= 0:
                    output_fasta.write(">" + all_genes[gene_index - i][0] + "\n")
                    output_fasta.write(all_genes[gene_index - i][1] + "\n")
            for i in range(1, n_after + 1):
                if gene_index + i <= (len(all_genes) - 1):
                    output_fasta.write(">" + all_genes[gene_index + i][0] + "\n")
                    output_fasta.write(all_genes[gene_index + i][1] + "\n")
