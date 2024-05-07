import os
import datetime
import tempfile
import shutil
import re

from typing import List, Tuple
from dataclasses import dataclass
from typing import List


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> None:
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
        output_fasta += ".fasta"
    with open(input_fasta) as file_input:
        with open(output_fasta, "w") as file_output:
            seq = []
            for line in file_input:
                if line.startswith(">"):
                    if seq:
                        file_output.write("".join(seq) + "\n")
                    seq = []
                    file_output.write(line)
                else:
                    seq.append(line.strip())
            else:
                file_output.write("".join(seq))


def select_gene_from_gbk(
    input_gbk: str, target_gene: str, n_before: int, n_after: int
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
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
            if line == "":
                return (None, None)
            if line.startswith("/gene"):
                if target_gene in line:
                    break
                gene = line[7:-1] + f"_before_{target_gene}"
                seq = None
                while seq is None:
                    line = file_input.readline().strip()
                    if line.startswith("/translation"):
                        subseqs = []
                        subseqs.append(line[1:].strip('translation="'))
                        while True:
                            line = file_input.readline().strip()
                            if line.endswith('"'):
                                subseqs.append(line[:-1])
                                break
                            subseqs.append(line)
                        seq = "".join(subseqs)
                genes_before.append((gene, seq))
                del genes_before[0]
        while True:
            line = file_input.readline().strip()
            if line == "":
                return [gene for gene in genes_before if not gene is None], [
                    gene for gene in genes_after if not gene is None
                ]
            if line.startswith("/gene"):
                gene = line[7:-1] + f"_after_{target_gene}"
                seq = None
                while seq is None:
                    line = file_input.readline().strip()
                    if line.startswith("/translation"):
                        subseqs = []
                        subseqs.append(line[1:].strip('translation="'))
                        while True:
                            line = file_input.readline().strip()
                            if line.endswith('"'):
                                subseqs.append(line[:-1])
                                break
                            subseqs.append(line)
                        seq = "".join(subseqs)
                genes_after.append((gene, seq))
                del genes_after[0]
                if not None in genes_after:
                    return genes_before, genes_after


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: List[str],
    n_before: int,
    n_after: int,
    output_fasta: str = None,
) -> None:
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
        output_fasta = (
            os.path.basename(input_gbk).split(".")[0] + "_target_proteins.fasta"
        )
    else:
        output_fasta += ".fasta"
    for gene in genes:
        with open(output_fasta, "w") as file:
            for target in genes:
                genes_before, genes_after = select_gene_from_gbk(
                    input_gbk, target, n_before, n_after
                )
                for gene, seq in genes_before:
                    file.write(">" + gene + "\n")
                    file.write(seq + "\n")
                for gene, seq in genes_after:
                    file.write(">" + gene + "\n")
                    file.write(seq + "\n")


def change_fasta_start_pos(
    input_fasta: str, shift: int, output_fasta: str = "shifted"
) -> None:
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
    with open(input_fasta, "r") as file_input:
        with open(output_fasta + ".fasta", "w") as file_output:
            file_output.write(file_input.readline())
            seq = file_input.readline().strip()
            seq = seq[shift:] + seq[:shift]
            file_output.write(seq)


def parse_blast_output(input_file: str, output_file: str = "blast_proteins") -> None:
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
        with open(output_file + ".txt", "w") as file_output:
            for line in file_input:
                if line.startswith("Alignments"):
                    file_input.readline()
                    protein = file_input.readline()[1:]
                    if "MULTISPECIES" in protein:
                        protein = protein.strip("MULTISPECIES: ")
                    proteins.add(protein)
            for protein in sorted(proteins):
                file_output.write(protein)


class MeasureTime:
    """Context manager that prints its process duration in seconds and hh:mm:ss format"""

    def __init__(self):
        self.starttime = None
        self.endtime = None
        self.elapsedtime = None

    def __enter__(self):
        self.starttime = datetime.datetime.now()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.endtime = datetime.datetime.now()
        self.elapsedtime = (self.endtime - self.starttime).total_seconds()
        print(
            "Process took in hh:mm:ss:",
            datetime.timedelta(seconds=self.elapsedtime // 1),
        )
        print(f"Process took in seconds: {self.elapsedtime:.2f}")


class TmpDir:
    """Context manager that creates temporary dir for its subprocess files

    Arguments
    ---------
    - working_dir - str, sets dir in which temporary dir will be created
    - to_delete - bool, choose if you want to delete temporary dir after subprocesses are over
    """

    def __init__(self, working_dir: str = ".", to_delete: bool = True):
        self.working_dir = working_dir
        self.to_delete = to_delete

    def __enter__(self):
        self.temp_dir = tempfile.mkdtemp(dir=self.working_dir)
        return self.temp_dir

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.to_delete:
            shutil.rmtree(self.temp_dir)


class MyDict:
    """My dictionary class that duplicates standard dictionary functionality, except that it iterates over tuples of key and value"""

    def __init__(self):
        self._data = {}

    def __setitem__(self, key, value):
        self._data[key] = value

    def __getitem__(self, key):
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        self._keys = iter(self._data)
        return self

    def __next__(self):
        key = next(self._keys)
        return key, self._data[key]


class FastaFormatError(ValueError):
    pass


@dataclass
class FastaRecord:
    """Data class for fasta record returned with OpenFasta context manager.

    Arguments
    ---------
    - id - str, ID of fasta record (example: "GTD326487.1")
    - description - str, additional information about this record (example: "Species anonymous 24 chromosome")
    - sequence - biological sequence (example: "ATCGACTACGACTAGCATCACGATCACGATACGATGCATCAGTAGCACTAGATCA")
    """

    id: str
    description: str
    sequence: str

    def __repr__(self):
        return f"FastaRecord(id={self.id}, description={self.description}, sequence={self.sequence})"


class OpenFasta:
    def __init__(self, file_path):
        self.file_path = file_path
        self.handler = None

    def __enter__(self):
        if is_singleline_fasta(self.file_path):
            self.handler = open(self.file_path)
            return self
        raise FastaFormatError(
            'Check fasta format, it should be single line.\nYou can change multiline to single line fasta with "translate_multiline_to_singleline" function!'
        )

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handler.close()

    def __iter__(self):
        return self

    def __next__(self):
        id, description, sequence = None, None, None
        while True:
            line = self.handler.readline()
            if not line:
                raise StopIteration
            line = line.strip()
            if line.startswith(">"):
                id, description = line[1:].split(None, 1)
            else:
                sequence = line
                return FastaRecord(id, description, sequence)

    def read_record(self) -> "FastaRecord":
        """readline method analogue from open context manager"""
        return next(self)

    def read_records(self) -> List["FastaRecord"]:
        """readlines method analogue from open context manager"""
        records = []
        for record in self:
            records.append(record)
        return records


def is_singleline_fasta(file: str) -> bool:
    """Function that checks if fasta file is single line

    Arguments:
    ---------
    - file - str, path to the fasta file

    Returns
    -------
    - bool
    """
    count = 0
    with open(file) as fasta:
        for line in fasta:
            line = line.strip()
            if not line:
                return False
            if line.startswith(">"):
                count += 1
                continue
            else:
                count += 1
                if count % 2 == 0:
                    continue
                else:
                    return False
        return True


def translate_multiline_to_singleline(input: str) -> None:
    """Function that creates single line fasta file from multiline fasta

    Arguments:
    ---------
    - input - str, path to the multiline fasta
    """
    with open(input) as mfasta:
        output = re.sub(r"(.*)\..*$", "\g<1>", input) + ".singleline.fasta"
        with open(output, "w") as sfasta:
            header = ""
            sequence = ""

            for line in mfasta:
                line = line.strip()
                if line.startswith(">"):
                    if sequence:
                        sfasta.write(f"{header}\n{sequence}\n")
                        sequence = ""

                    header = line
                else:
                    sequence += line
            if header:
                sfasta.write(f"{header}\n{sequence}")
