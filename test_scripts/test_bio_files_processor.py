import pytest

from bio_files_processor import OpenFasta, FastaRecord, FastaFormatError, is_singleline_fasta


def test_open_fasta():
    with pytest.raises(FastaFormatError):
        with OpenFasta('data/NUDT9_gene_multiline.fasta'):
            pass


def test_singleline_check():
    assert is_singleline_fasta('data/NUDT9_gene_oneline.fasta')


def test_fasta_record_return():
    with OpenFasta('data/NUDT9_gene_oneline.fasta') as f:
        assert isinstance(f.read_record(), FastaRecord)