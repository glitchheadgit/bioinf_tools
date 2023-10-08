from typing import List, Union


DNA = {
    "T": "A",
    "A": "T",
    "C": "G",
    "G": "C",
    "t": "a",
    "a": "t",
    "c": "g",
    "g": "c"
}
RNA = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",
    "a": "u",
    "u": "a",
    "c": "g",
    "g": "c"
}
TRANSCRIBE_LIST = {
    "C": "C",
    "G": "G",
    "c": "c",
    "g": "g",
    "A": "A",
    "a": "a",
    "T": "U",
    "t": "u",
}


def transcribe(sequences: List[str]) -> List[str]:
    """
    Gets transcribed sequences

    Arguments
    ---------
    List of sequences

    Returns
    -------
    List of transcribed sequences
    """
    transcribed = []
    for num, sequence in enumerate(sequences):
        if is_dna((sequence,))[0]:
            transcribed.append("")
            for nucleotide in sequence:
                transcribed[num] += TRANSCRIBE_LIST[nucleotide]
        else:
            raise ValueError(f"{sequence} is not dna!")
    return transcribed


def reverse(sequences: List[str]) -> List[str]:
    """
    Gets reverse sequences

    Arguments
    ---------
    List of sequences

    Returns
    -------
    List of reverse sequences
    """
    reversed_seq = [sequence[::-1] for sequence in sequences]
    return reversed_seq


def complement(sequences: List[str]) -> List[str]:
    """
    Gets complement sequences

    Arguments
    ---------
    List of sequences

    Returns
    -------
    List of complement sequences
    """
    complemented = []
    for num, sequence in enumerate(sequences):
        complemented.append("")
        if "u" in sequence.lower():
            for nucleotide in sequence:
                complemented[num] += RNA[nucleotide]
        else:
            for nucleotide in sequence:
                complemented[num] += DNA[nucleotide]
    return complemented


def reverse_complement(sequences: List[str]) -> List[str]:
    """
    Gets reverse complement sequences

    Arguments
    ---------
    List of sequences

    Returns
    -------
    List of reverse complement sequences
    """
    return reverse(complement(sequences))


def is_dna(sequences: List[str]) -> List[Union[bool, str]]:
    """
    Checks if sequences are DNA.

    Arguments
    ---------
    List of sequences

    Returns
    -------
    List of True/False/'Unclear' values for according sequence
    """
    bool_list = [
        (True if "t" in sequence.lower() else "Unclear")
        if "u" not in sequence.lower()
        else False
        for sequence in sequences
    ]
    return bool_list


def length(sequences: List[str]) -> List[int]:
    """
    Measures sequences length

    Arguments
    ---------
    List of sequences

    Returns
    ---------
    List of integers
    """
    length_list = [len(sequence) for sequence in sequences]
    return length_list


def is_palindrome(sequences: List[str]) -> List[bool]:
    """
    Checks if sequences are palindromes.

    Arguments
    ---------
    List of sequences

    Returns
    -------
    List of True/False values for according sequence
    """
    bool_list = [
        complement((sequence[::-1].lower(),))[0] == sequence.lower()
        for sequence in sequences
    ]
    return bool_list
