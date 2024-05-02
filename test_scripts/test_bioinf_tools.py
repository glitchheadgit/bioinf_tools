import pytest

from bioinf_tools import DNASequence, RNASequence, AminoAcidSequence, run_genscan


def are_files_identical(file1, file2):
    with open(file1, "rb") as f1, open(file2, "rb") as f2:
        return f1.read() == f2.read()


def test_wrong_sequence():
    with pytest.raises(ValueError) as e_info:
        DNASequence("ATUTCG")


def test_complementary():
    assert isinstance(DNASequence("ATGC").complement(), DNASequence)


def test_transcription():
    assert isinstance(DNASequence("ATGC").transcribe(), RNASequence)


def test_genscan():
    assert run_genscan(sequence_file="data/NUDT9_gene_multiline.fasta").status == 200


def test_aminoacid_sequence():
    assert (AminoAcidSequence('ADKEH').check_alphabet() and AminoAcidSequence('ADKEH').define_charge()['Positive'] == 2)
