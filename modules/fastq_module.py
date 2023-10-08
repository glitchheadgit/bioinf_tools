from typing import Union


def count_gc(seq: str) -> Union[float, int]:
    """
    Counts GC content of nucleotide sequence in percentages.

    Arguments
    ----------
    seq: str
        Nucleotide sequence
    Returns
    int or float
    -------
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
    int or float
    -------
    """
    mean_quality = 0
    for nucleotide_quality in seq_quality:
        mean_quality += ord(nucleotide_quality) - 33
    return mean_quality/len(seq_quality)



