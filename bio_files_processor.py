import os
from typing import List, Tuple


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Converts multiline fasta to oneline fasta.

    Arguments
    ---------
    input_fasta: str
        Path to the multiline fasta
    output_fasta: str (optional)
        Name of the output oneline fasta file. If not specified name will be the same as the original file.
    """
    if output_fasta is None:
        output_fasta = os.path.basename(input_fasta)
    else:
        output_fasta += '.fasta'
    with open(input_fasta) as file_input:
        with open(output_fasta, 'w') as file_output:
            seq = []
            for line in file_input:
                if line.startswith('>'):
                    if seq:
                        file_output.write(''.join(seq) + '\n')
                    seq = []
                    file_output.write(line)
                else:
                    seq.append(line.strip())
            else:
                file_output.write(''.join(seq))


def select_gene_from_gbk(input_gbk: str, target_gene: str, n_before: int, n_after: int) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """
    Selects genes around gene of interest in gbk file.

    Arguments
    ---------
    input_gbk: str
        Path to the gbk file
    target_gene: str
        Name of the target gene
    n_before: int
        Number of genes we want to get before target
    n_after: int
        Number of genes we want to get after target

    Returns
    -------
    Tuple of two lists: genes before and after target, which consists of tuples with name of a gene and it's aminoacid sequence
    """
    with open(input_gbk) as file_input:
        genes_before = n_before * [None]
        genes_after = n_after * [None]
        while True:
            line = file_input.readline().strip()
            if line == '':
                return (None, None)
            if line.startswith('/gene'):
                if target_gene in line:
                    break
                gene = line[7:-1] + f'_before_{target_gene}'
                seq = None
                while seq is None:
                    line = file_input.readline().strip()
                    if line.startswith('/translation'):
                        subseqs = []
                        subseqs.append(line[1:].strip('translation="'))
                        while True:
                            line = file_input.readline().strip()
                            if line.endswith('"'):
                                subseqs.append(line[:-1])
                                break
                            subseqs.append(line)
                        seq = ''.join(subseqs)
                genes_before.append((gene, seq))
                del genes_before[0]
        while True:
            line = file_input.readline().strip()
            if line == '':
                return [gene for gene in genes_before if not gene is None], [gene for gene in genes_after if not gene is None]
            if line.startswith('/gene'):
                gene = line[7:-1] + f'_after_{target_gene}'
                seq = None
                while seq is None:
                    line = file_input.readline().strip()
                    if line.startswith('/translation'):
                        subseqs = []
                        subseqs.append(line[1:].strip('translation="'))
                        while True:
                            line = file_input.readline().strip()
                            if line.endswith('"'):
                                subseqs.append(line[:-1])
                                break
                            subseqs.append(line)
                        seq = ''.join(subseqs)
                genes_after.append((gene, seq))
                del genes_after[0]
                if not None in genes_after:
                    return genes_before, genes_after


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: List[str], n_before: int, n_after: int, output_fasta: str = None) -> None:
    """
    Writes output of select_gene_from_gbk to fasta file

    Arguments
    ---------
    input_gbk: str
        Path to gbk file
    genes: List[str]
        List of genes of interest
    n_before: int
        Number of genes to get before target
    n_after: int
        Number of genes to get after target
    output_fasta: str (optional)
        Name of output file. If not specified name will be the same as input file.
    """
    if output_fasta is None:
        output_fasta = os.path.basename(input_gbk).split('.')[0] + '_target_proteins.fasta'
    else:
        output_fasta += '.fasta'
    for gene in genes:
        with open(output_fasta, 'w') as file:
            for target in genes:
                genes_before, genes_after = select_gene_from_gbk(input_gbk, target, n_before, n_after)
                for gene, seq in genes_before:
                    file.write('>' + gene + '\n')
                    file.write(seq + '\n')
                for gene, seq in genes_after:
                    file.write('>' + gene + '\n')
                    file.write(seq + '\n')
