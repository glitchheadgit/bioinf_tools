def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str = 'shifted') -> None:
    """
    Changes fasta start position.

    Arguments
    ---------
    input_fasta: str
        Path to the input fasta file
    shift: int
        Number of shift
    output_fasta:
        Path to output file
    """
    with open(input_fasta, 'r') as file_input:
        with open(output_fasta + '.fasta', 'w') as file_output:
            file_output.write(file_input.readline())
            seq = file_input.readline().strip()
            seq = seq[shift:] + seq[:shift]
            file_output.write(seq)


def parse_blast_output(input_file: str, output_file: str = 'blast_proteins') -> None:
    """
    Parses blast output. Saves list of the most accurate proteins on each query.

    Arguments
    ---------
    input_file: str
        Path to blast output
    output_file: str
        Path to output file
    """
    proteins = set()
    with open(input_file) as file_input:
        with open(output_file + '.txt', 'w') as file_output:
            for line in file_input:
                if line.startswith('Alignments'):
                    file_input.readline()
                    protein = file_input.readline()[1:]
                    if 'MULTISPECIES' in protein:
                        protein = protein.strip('MULTISPECIES: ')
                    proteins.add(protein)
            for protein in sorted(proteins):
                file_output.write(protein)
