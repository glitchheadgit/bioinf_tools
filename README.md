# HW18: API, OOP, ML, and File Handling practices

This repository serves as a comprehensive Bioinformatics Institute homework refactoring focusing on Object-Oriented Programming (OOP) in general, Machine Learning (ML), API and file handling.

## Overview
In this project, you will find implementations for the following functionalities:

* **RandomForestClassifierCustom**: A custom implementation of Random Forest classifier that supports parallel fitting and prediction.
* **Classes for RNA, DNA, Amino Acids**: These classes offer functionality for essential molecular biology tasks including obtaining complementary sequences, transcriptions (DNA to RNA), and calculating the count of charge groups in amino acids.
* **Fastq Filter**: Filters FASTQ files based on quality, length, and GC content.
* **Telegram Logger Decorator**: A decorator for logging information to Telegram.
* **Genscan Custom API**: An API wrapper for running CDS, exon, and intron prediction using Genscan on given nucleotide sequences or FASTA files.
* **Fasta File Context Manager**: A context manager for working with FASTA files, capable of reading records and returning custom objects with metadata and sequences.
* **Fasta Format Checker**: Checks whether a FASTA file is in one-line or multi-line format.
* **Fasta Multiline to Oneline Converter**: Converts multi-line FASTA files to one-line format.
* **Context manager with tmpdir functionality**: Creates temporary dir for subprocesses
* **Test Suite**: Includes a suite of tests covering key functionalities of the modules.
* **Showcases Notebook**: A Jupyter notebook (showcases.ipynb) containing examples demonstrating the usage of the implemented functions and modules.



## Usage
To utilize the functionalities provided in this repository, simply clone the repository and import the required modules into your Python environment. Refer to the showcases.ipynb notebook for detailed examples and usage scenarios.


## Dependencies
You can install them via pip using the provided requirements.txt file:
```bash
pip install -r requirements.txt
```