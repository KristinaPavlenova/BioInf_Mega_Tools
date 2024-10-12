def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = '_oneline.fasta'):
    if output_fasta == '_oneline.fasta':
        output_fasta = input_fasta.split('.')[0] + output_fasta
    with open(input_fasta, 'r') as in_fasta, open(output_fasta, 'w') as out_fasta:
        for line in in_fasta:
            if not line.startswith('>'):
                seq = []
                while line.strip().isalpha():
                    seq.append(line.strip())
                    line = in_fasta.readline()
                out_fasta.write(''.join(seq) + '\n')
            out_fasta.write(line)


def parse_blast_output(input_file: str, output_file: str):
    with open(input_file, 'r') as blast_results, open(output_file, 'w') as best_ident_discriptions:
        best_ident_list = []
        for line in blast_results:
            if 'Sequences producing significant alignments:' in line:
                blast_results.readline()
                blast_results.readline()
                line = blast_results.readline()
                best_ident_list.append(line.split('  ')[0])
        best_ident_list.sort()
        for discription in best_ident_list:
            print(discription, file=best_ident_discriptions)                
