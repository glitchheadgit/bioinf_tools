from modules import fastq_module
from modules import protein_module
from modules import dna_module
from typing import Union, Dict, List, Tuple
import os
import sys


NUCLEOTIDES = {'a', 't', 'g', 'c', 'u'}


def filter_fastq(input_path: str, gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]]= (0, 100), length_bounds: Union[int, Tuple[int, int]] = (0, 2**32), quality_threshold: int = 0, output_filename: str = None) -> None:
    """
    Filters appropriate sequences.

    Arguments
    ----------
    input_path: str
        Path to fastq file
    output_filename: str
        Name of filtered fastq file
    gc_bounds: Union[int, float, Tuple[Union[int,float], Union[int,float]]
        Filter parameter (in percents) of gc composition. Tuple associated with lowest and highest levels of gc content in  sequences, also it can get on input int or float associated with highest level of gc content, lowest level will be set to 0.
    length_bounds: Union[int, Tuple[int, int]]
        Filter parameter of sequences length. Logic the same as in gc_bounds.
    quality_threshold: int
        Filter parameter of mean sequence quality. Sets lowest mean quality of sequences.

    Saves filtered sequences in a fastq file named input_path/output_filename.
    """
    seqs = fastq_module.read_fastq(input_path)
    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int) or isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)
    seqs_filtered = {}
    for name, seq in seqs.items():
        seq_length = len(seq[0])
        if seq_length < length_bounds[0] or seq_length > length_bounds[1]:
            continue
        gc_count = fastq_module.count_gc(seq[0])
        if gc_count < gc_bounds[0] or gc_count > gc_bounds[1]:
            continue
        mean_quality = fastq_module.calculate_mean_quality(seq[1])
        if mean_quality < quality_threshold:
            continue
        seqs_filtered[name] = seq
    if output_filename is None:
        output_filename = os.path.basename(input_path)
    else:
        output_filename += '.fastq'
    fastq_module.write_fastq(output_filename, seqs_filtered)


def protein_tool(*args: str) -> Union[str, List[Union[Dict[str, int], str]]]:
    """
    Receives a request from the user and runs the desired function.

    Arguments
    ----------
    seq : str
        Amino acid sequences.
    operation : str
        Type of user's request.

    Returns
    -------
    str
        If a single sequence is supplied, outputs the result as a string or or identify a problem with a specific sequence.
    list
        If several sequences are supplied, outputs the result as a list.
    """
    *seqs, operation = args
    operations = {'one letter': protein_module.change_abbreviation, 'RNA': protein_module.to_rna, 'DNA': protein_module.to_dna, 'charge': protein_module.define_charge, 'polarity': protein_module.define_polarity}
    output = []
    for seq in seqs:
        seq_check = protein_module.is_correct_seq(seq.upper())
        if seq_check:
            function_output = operations[operation](seq.upper())
            output.append(function_output)
        else:
            print(f'Something wrong with {seq}', file=sys.stderr)
            continue
    if len(output) == 1 and (operation == 'RNA' or operation == 'DNA' or operation == 'one letter'):
        return ''.join(output)
    else:
        return output


def rna_dna_tool(*args: List[str]) -> List[Union[bool, int, str]]:
    """
    Processes incoming sequences with the specified command.

    Arguments
    ----------
    seq : List[str]
        Nucleotide sequences.
    operation : str
        Command name.

    Returns
    -------
    str
        If a single sequence is supplied, outputs the result as a string or or identify a problem with a specific sequence.
    list
        If several sequences are supplied, outputs the result as a list.
    """
    *sequences_list, command = args
    sequences = []
    for sequence in sequences_list:
        sequence_low = sequence.lower()
        if ('u' in sequence_low and 't' in sequence_low):
            print(f'{sequence} is not NA. Deleted.')
            continue
        if not all([x in NUCLEOTIDES for x in sequence_low]):
            print(f'{sequence} is not NA. Deleted.')
            continue
        sequences.append(sequence)
    match command:
        case 'transcribe':
            result = dna_module.transcribe(sequences)
        case 'reverse':
            result = dna_module.reverse(sequences)
        case 'complement':
            result = dna_module.complement(sequences)
        case 'reverse_complement':
            result = dna_module.reverse_complement(sequences)
        case 'is_palindrome':
            result = dna_module.is_palindrome(sequences)
        case 'length':
            result = dna_module.length(sequences)
        case 'is_dna':
            result = dna_module.is_dna(sequences)
        case _:
            raise ValueError("The specified command is invalid!")
    return result[0] if len(result) == 1 else result
