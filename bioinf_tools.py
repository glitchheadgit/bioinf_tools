import numpy as np
import os
import sys
import time
import requests
import re

from abc import ABC, abstractmethod
from typing import Union, Dict, Tuple, Set
from typing_extensions import Self
from io import StringIO
from typing import Dict, List, Callable, Union, Optional
from dataclasses import dataclass
from bs4 import BeautifulSoup
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class ExpressedSequences(BiologicalSequence):
    def __init__(self, seq: str, alphabet: Set[str]) -> None:
        self.letters = set(seq)
        self.seq = seq
        self._alphabet = alphabet

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, s):
        if not s:
            raise ValueError("Sequence can't be empty!")
        if not self.check_alphabet():
            raise ValueError("Sequence contains incorrect symbols!")
        self._seq = s

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, index: int) -> "Self":
        return self.__class__(self.seq[index])

    def __str__(self) -> str:
        return self.seq

    def check_alphabet(self) -> bool:
        return all([x in self._alphabet for x in self.letters]) and not (
            "U" in self.letters and "T" in self.letters
        )


class NucleicAcidSequence(ExpressedSequences):
    def __init__(self, seq: str) -> None:
        self._alphabet = {"A", "T", "G", "C", "U"}
        super().__init__(seq, self._alphabet)
        self._complement_table = None

    def complement(self) -> "Self":
        if self._complement_table is None:
            raise NotImplementedError
        complemented = []
        for nucleotide in self.seq:
            complemented.append(self._complement_table[nucleotide])
        return self.__class__("".join(complemented))

    def gc_content(self) -> int:
        return (self.seq.count("G") + self.seq.count("C")) / len(self.seq) * 100


class DNASequence(NucleicAcidSequence):
    def __init__(self, seq: str) -> None:
        super().__init__(seq)
        self._complement_table = {"T": "A", "A": "T", "C": "G", "G": "C"}
        self._transcribe_table = {
            "C": "C",
            "G": "G",
            "A": "A",
            "T": "U",
        }

    def transcribe(self) -> "RNASequence":
        transcribed = []
        for nucleotide in self.seq:
            transcribed.append(self._transcribe_table[nucleotide])
        return RNASequence("".join(transcribed))


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq: str):
        super().__init__(seq)
        self._complement_table = {"A": "U", "U": "A", "C": "G", "G": "C"}


class AminoAcidSequence(ExpressedSequences):
    def __init__(self, seq: str) -> None:
        self._alphabet = {
            "A",
            "G",
            "D",
            "L",
            "N",
            "P",
            "C",
            "Y",
            "S",
            "I",
            "H",
            "W",
            "E",
            "F",
            "R",
            "T",
            "V",
            "K",
            "M",
            "Q",
        }
        super().__init__(seq, self._alphabet)
        self._positive_aa = {"R", "K", "H"}
        self._negative_aa = {"D", "E"}

    def define_charge(self) -> Dict[str, int]:
        positive_count = 0
        negative_count = 0
        neutral_count = 0
        for aa in self.seq:
            if aa in self._positive_aa:
                positive_count += 1
            elif aa in self._negative_aa:
                negative_count += 1
            else:
                neutral_count += 1
        result = {
            "Positive": positive_count,
            "Negative": negative_count,
            "Neutral": neutral_count,
        }
        return result


def filter_fastq(
    input_path: str,
    output_filename: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (
        0,
        100,
    ),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
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
    records = SeqIO.parse(input_path, "fastq")
    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int) or isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)
    records_filtered = []
    for rec in records:
        seq_length = len(rec)
        if seq_length < length_bounds[0] or seq_length > length_bounds[1]:
            continue
        gc_count = gc_fraction(rec.seq) * 100
        if gc_count < gc_bounds[0] or gc_count > gc_bounds[1]:
            continue
        mean_quality = np.mean(rec.letter_annotations["phred_quality"])
        if mean_quality < quality_threshold:
            continue
        records_filtered.append(rec)
    SeqIO.write(records_filtered, output_filename, "fasta")


def telegram_logger(chat_id: str) -> Callable:
    """Decorator to log function execution details to a Telegram chat.

    Args:
        chat_id (str): The ID of the Telegram chat to which the logs will be sent.

    Returns:
        function with additional logging functionality

    Example:
        @telegram_logger("123456789")
        def my_function():
            # Function implementation
    """
    token = os.getenv("TG_API_TOKEN")
    url = f"https://api.telegram.org/bot{token}/"

    def decorator(func: Callable) -> Callable:
        def wrapper(*args, **kwargs) -> Union[object, None]:
            buffer = StringIO()
            sys.stdout = buffer
            sys.stderr = buffer
            start_time = time.time()
            func_name = func.__name__
            log_file = f"{func_name}.log"
            try:
                result = func(*args, **kwargs)
                end_time = time.time()
                elapsed_time = end_time - start_time
                if elapsed_time < 86400:  # Время выполнения менее 1 дня
                    time_str = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
                else:
                    days = int(elapsed_time // 86400)
                    time_str = time.strftime(
                        f"{days} days, %H:%M:%S", time.gmtime(elapsed_time)
                    )
                message = f"&#x2728; Function <code>{func_name}</code> executed successfully in <code>{time_str}</code>"
                buffer.seek(0)
                r = requests.post(
                    url + "sendDocument",
                    data={"chat_id": chat_id, "caption": message, "parse_mode": "html"},
                    files={"document": (log_file, buffer)},
                )
                if r.status_code // 100 == 2:
                    return result
                else:
                    r = requests.post(
                        url + "sendMessage",
                        data={
                            "chat_id": chat_id,
                            "text": message + "\nFunction doesn't have logs!",
                            "parse_mode": "html",
                        },
                    )
                    return result
            except Exception as e:
                end_time = time.time()
                elapsed_time = end_time - start_time

                if elapsed_time < 86400:
                    time_str = time.strftime("%H:%M:%S.%f", time.gmtime(elapsed_time))
                else:
                    days = int(elapsed_time // 86400)
                    time_str = time.strftime(
                        f"{days} days, %H:%M:%S", time.gmtime(elapsed_time)
                    )

                error_type = type(e).__name__
                error_text = str(e)
                error_message = f"&#128548; Function <code>{func_name}</code> encountered an error:\n\n<code>{error_type}: {error_text}</code>"
                buffer.seek(0)
                r = requests.post(
                    url + "sendDocument",
                    data={
                        "chat_id": chat_id,
                        "caption": error_message,
                        "parse_mode": "html",
                    },
                    files={"document": (log_file, buffer)},
                )
                if r.status_code // 100 == 2:
                    return result
                else:
                    r = requests.post(
                        url + "sendMessage",
                        data={
                            "chat_id": chat_id,
                            "text": error_message + "\nFunction doesn't have logs!",
                            "parse_mode": "html",
                        },
                    )
                    return result
            finally:
                buffer.close()

        return wrapper

    return decorator


@dataclass
class GenscanOutput:
    """Data class to hold output of the Genscan prediction.

    Attributes:
        status (str): The status code of the Genscan prediction.
        cds_list (Dict[str, str]): Predicted coding sequences.
        intron_list (List[Dict[str, int]]): Predicted introns.
        exon_list (List[Dict[str, int]]): Predicted exons.

    Example:
        output = GenscanOutput(status=200, cds_list={}, intron_list=[], exon_list=[])
    """

    status: str
    cds_list: Dict[str, str]
    intron_list: List[Dict[str, int]]
    exon_list: List[Dict[str, int]]


def run_genscan(
    sequence: Optional[str] = None,
    sequence_file: Optional[str] = None,
    organism: str = "Vertebrate",
    exon_cutoff: float = 1.00,
    sequence_name: str = "",
) -> GenscanOutput:
    """Run Genscan prediction for a given sequence.

    Args:
        sequence (str): The DNA sequence to predict genes from.
        sequence_file (str): Path to a file containing the DNA sequence.
        organism (str): The organism for which to run the prediction.
        exon_cutoff (float): The exon cutoff value for the prediction.
        sequence_name (str): Name of the sequence.

    Returns:
        GenscanOutput: Data Class instance with the Genscan prediction.

    Raises:
        ValueError: If incorrect organism or exon_cutoff is provided.
        ValueError: If neither sequence nor sequence_file is provided.

    Example:
        output = run_genscan(sequence="ATGCATGCATGC")
    """
    ORGANISMS = {"Vertebrate", "Arabidopsis", "Maize"}
    EXON_CUTOFFS = {1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01}
    url = "http://hollywood.mit.edu/GENSCAN.html"
    job_url = "http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi/"
    with requests.Session() as s:
        req = s.get(url)
        payload = {}
        file = {}
        soup = BeautifulSoup(req.text, "lxml")
        if organism not in ORGANISMS:
            raise ValueError(f"Wrong organism provided. Choose from {ORGANISMS}")
        if exon_cutoff not in EXON_CUTOFFS:
            raise ValueError(f"Wrong exon cutoff provided. Choose from {EXON_CUTOFFS}")
        if sequence:
            seq = soup.find("textarea", attrs={"name": "-s"})
            payload[seq["name"]] = sequence
        elif sequence_file:
            seq = soup.find("input", attrs={"name": "-u"})
            file = {
                seq["name"]: ("input.fasta", open(sequence_file, "rb"), "text/plain")
            }
        else:
            raise ValueError("You must specify a sequence or a file with it!")
        org = soup.find("select", attrs={"name": "-o"})
        payload[org["name"]] = organism
        cut = soup.find("select", attrs={"name": "-e"})
        payload[cut["name"]] = exon_cutoff
        name = soup.find("input", attrs={"name": "-n"})
        payload[name["name"]] = sequence_name
        payload["-p"] = "Predicted peptides only"
        r = s.post(job_url, data=payload, files=file)
        if file:
            file[seq["name"]][1].close()

    soup = BeautifulSoup(r.text, "lxml")
    results = soup.find("pre").text
    exons_list = []
    exons = re.search(r"Predicted genes/exons.*\n([^<]*)Suboptimal", results).group(1)
    exons = re.sub(r" ?- ?-", "", exons.strip())
    exons = [s.strip() for s in exons.strip().split("\n") if s.strip()]
    exons = map(lambda x: x.split(), exons)
    i = 0
    for exon in exons:
        if i == 0:
            sequence_number = "1"
            i += 1
            exons_list.append([])
            continue
        if sequence_number != exon[0].split(".")[0]:
            sequence_number = exon[0].split(".")[0]
            i = 1
            exons_list.append([])
        if exon[0] == "NO":
            print("No exons found")
            return None
        exons_list[-1].append({"index": i, "start": int(exon[3]), "end": int(exon[4])})
        i += 1
    introns_list = []
    exon1 = None
    for sequence_exons in exons_list:
        i = 1
        exon1 = None
        introns_list.append([])
        for exon in sequence_exons:
            if exon1 is None:
                exon1 = exon
            else:
                introns_list[-1].append(
                    {"index": i, "start": exon1["end"] + 1, "end": exon["start"] - 1}
                )
                exon1 = exon
                i += 1
    cds_list = []
    pps = (
        re.search(r"Predicted peptide sequence.*\n([^<]*)", results)
        .group(1)
        .strip()
        .split("\n\n\n\n")
    )
    for i, pp in enumerate(pps):
        name, seq = pp.split("\n\n")[0][1:], "".join(pp.split("\n\n")[1:])
        if sequence and sequence_name:
            name = sequence_name
        cds_list.append({"index": i + 1, "name": name, "sequence": seq})
    return GenscanOutput(r.status_code, cds_list, introns_list, exons_list)
