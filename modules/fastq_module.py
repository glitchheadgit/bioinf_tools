from typing import Union, Dict, Tuple


def read_fastq(input_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Reads a fastq file and makes from data a dict of sequences.

    Arguments
    ---------
    input_path: str
        Path to a fastq file

    Returns
    -------
    Dictionary of seq names as keys and a tuple of sequences and it's quality as values.
    """
    seqs = {}
    with open(input_path, 'r') as file:
        for line in file:
            if line.startswith("@"):
                counter = 0
                name = line[1:]
            if counter == 1:
                seq = line
            if counter == 3:
                quality = line
                seqs[name] = (seq, quality)
            counter += 1
    return seqs


def write_fastq(output_filename: str, seqs: Dict[str, Tuple[str, str]]):
    """
    Writes dict of sequence names as keys and tuple of sequence and it's quality as value to a file in a fastq format.

    Arguments
    ---------
    output_filename: str
        Name of an output fastq file
    seqs: Dict[str, Tuple[str, str]]
        Dict of sequences
    """
    with open(output_filename, 'w') as file:
        for name, seq in seqs.items():
            file.write(f'@{name}')
            file.write(f'{seq[0]}')
            file.write(f'+{name}')
            file.write(f'{seq[1]}')


def count_gc(seq: str) -> Union[float, int]:
    """
    Counts GC content of nucleotide sequence in percentages.

    Arguments
    ----------
    seq: str
        Nucleotide sequence
    Returns
    -------
        int or float
    """
    return (seq.count('G') + seq.count('C'))/len(seq)*100


def calculate_mean_quality(seq_quality: str) -> Union[float, int]:
    """
    Counts mean quality of nucleotide sequence.

    Arguments
    ----------
    seq: str
        Nucleotide sequence
    Returns
    -------
        int or float
    """
    mean_quality = 0
    for nucleotide_quality in seq_quality:
        mean_quality += ord(nucleotide_quality) - 33
    return mean_quality/len(seq_quality)



