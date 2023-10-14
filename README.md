# BioInf Tools

This Python utility allows you to work with protein, DNA, RNA sequences and fastq format dictionaries. You can perform various operations on sequences, such as translating them into RNA, DNA, aminoacid sequences, counting charged, uncharged, hydrophobic, hydrophilic amino acids and many more... 

## Table of Contents

- [Installation](#installation)
- [Bio files processor](#bio_files_processor)
- [Additional scripts](#change_fasta_start_pos)
- [Fastq tool](#fastq_tool)
- [RNA DNA tool](#rna_dna_tool)
- [Protein tool](#protein_tool)
- [Troubleshooting](#troubleshooting)

## Installation

You can clone this repository or download the source code. 

##### Requirements:

Python3

# fastq_tool

## filter_fastq

Function gets on input a fastq file. Function filters sequences according to defined parameters: GC content, length, mean quality - and saves filtered sequences in a fastq file named input_path or output_filename.

Required arguments:

- input_path - str of a path to the fastq file. 

Optional arguments:

- output_filename - str of a filename in which filtered sequences will be writed.

- gc_bounds - int, float of max value or tuple with min and max values. Default = (0, 100)

- length_bounds - int of max value or tuple with min and max values. Defaul = (0, 2**32)

- quality_threshold - int, minimum mean quality of a sequence

```python
filter_fastq({'BEST_SEQ_EVER':('GG', '!!')}, gc_bounds = 50)
## Blank fasta
filter_fastq({'BEST_SEQ_EVER':('GG', '!!')}, quality_threshold = 0)
## fasta file
# >BEST_SEQ_EVER
# GG
# +BEST_SEQ_EVER
# !!
```

# bio_files_processor

## convert_multiline_fasta_to_oneline

Converts multiline fasta to oneline fasta.

**Arguments**:

- input_fasta: str
    Path to the multiline fasta
- output_fasta: str (optional)
    Name of the output oneline fasta file. If not specified name will be the same as the original file.

```python
    convert_multiline_fasta_to_oneline('bad.fasta', 'good')
## bad.fasta
# >my_genome
# AGCTA
# CG
## good.fasta
# >my_genome
# AGCTACG
```

## select_gene_from_gbk

Selects genes around gene of interest in gbk file.

**Arguments**:

- input_gbk: str
    Path to the gbk file
- target_gene: str
    Name of the target gene
- n_before: int
    Number of genes we want to get before target
- n_after: int
    Number of genes we want to get after target

**Returns**:
Tuple of two lists: genes before and after target, which consists of tuples with name of a gene and it's aminoacid sequence

```python
select_gene_from_gbk('input.gbk', 'yeeW', 3, 3)
# Returns:
# ([('galT_before_modE', 'MTQFNPVDHPHRRYNPLTGQWILVSPHRAKRPWQGAQETPAKQVLPAHDPDCFLCAGNVRVTGDKNPDYTRTYVFTNDFAALMSDTPDAPESNDPLMRCQSARGTSRVICFSPDHSKTLPELSVAALTEIVKTWQEQTAELGKTYPWVQVFENKGAAMGCSNPHPHGQIWANSFLPNEAEREDRLQKEYFAEQKSPMLVDYVQRELADGSRTVVETEHWLAVVPYWAAWPFETLLLPKAHVLRITDLTDAQRSDLALALKKLTSRYDNLFQCSFPYSMGWHGAPFNGEENQHWQLHAHFYPPLLRSATVRKFMVGYEMLAETQRDLTAEQAAERLRAVSDIHFRESGV'), ('galE_before_modE', 'MRVLVTGGSGYIGSHTCVQLLQNGHDVIILDNLCNSKRSVLPVIERLGGKHPTFVEGDIRNEALMTEILHDHAIDTVIHFAGLKAVGESVQKPLEYYDNNVNGTLRLISAMRAANVKNFIFSSSATVYGDQPKIPYVESFPTGTPQSPYGKSKLMVEQILTDLQKAQPDWSIALLRYFNPVGAHPSGDMGEDPQGIPNNLMPYIAQVAVGRRDSLAIFGNDYPTEDGTGVRDYIHVMDLADGHVVAMEKLANKPGVHIYNLGAGVGNSVLDVVNAFSKACGKPVNYHFAPRREGDLPAYWADASKADRELNWRVTRTLDEMAQDTWHWQSRHPQGYPD'), ('modF_before_modE', 'MSSLQILQGTFRLSDTKTLQLPQLTLNAGDSWAFVGSNGSGKSALARALAGELPLLKGERQSQFSHITRLSFEQLQKLVSDEWQRNNTDMLGPGEDDTGRTTAEIIQDEVKDAPRCMQLAQQFGITALLDRRFKYLSTGETRKTLLCQALMSEPDLLILDEPFDGLDVASRQQLAERLASLHQSGITLVLVLNRFDEIPEFVQFAGVLADCTLAETGAKEELLQQALVAQLAHSEQLEGVQLPEPDEPSARHALPTNEPRIVLNNGVVSYNDRPILNNLSWQVNPGEHWQIVGPNGAGKSTLLSLITGDHPQGYSNDLTLFGRRRGSGETIWDIKKHIGYVSSSLHLDYRVSTTVRNVILSGYFDSIGIYQAVSDRQQKLVQQWLDILGIDKRTADAPFHSLSWGQQRLALIVRALVKHPTLLILDEPLQGLDPLNRQLIRRFVDVLISEGETQLLFVSHHAEDAPACITHRLEFVPDGGLYRYVLTKIY')], [('acrZ_after_modE', 'MLELLKSLVFAVIMVPVVMAIILGLIYGLGEVFNIFSGVGKKDQPGQNH'), ('modA_after_modE', 'MARKWLNLFAGAALSFAVAGNALADEGKITVFAAASLTNAMQDIATQYKKEKGVDVVSSFASSSTLARQIEAGAPADLFISADQKWMDYAVDKKAIDTASRQTLLGNSLVVVAPKASEQKDFTIDSKTNWTSLLNGGRLAVGDPEHVPAGIYAKEALQKLGAWDTLSPKLAPAEDVRGALALVERNEAPLGIVYGSDAVASKGVKVVATFPEDSHKKVEYPVAVVEGHNNATVKAFYDYLKGPQAAEIFKRYGFTTK'), ('modB_after_modE', 'MILTDPEWQAVLLSLKVSSLAVLFSLPFGIFFAWLLVRCTFPGKALLDSVLHLPLVLPPVVVGYLLLVSMGRRGFIGERLYDWFGITFAFSWRGAVLAAAVMSFPLMVRAIRLALEGVDVKLEQAARTLGAGRWRVFFTITIPLTLPGIIVGTVLAFARSLGEFGATITFVSNIPGETRTIPSAMYTLIQTPGGESGAARLCIISIALAMISLLISEWLARISRERAGR')])
```

## select_genes_from_gbk_to_fasta

Writes output of select_gene_from_gbk to fasta file

**Arguments**:

- input_gbk: str
    Path to gbk file
- genes: List[str]
    List of genes of interest
- n_before: int
    Number of genes to get before target
- n_after: int
    Number of genes to get after target
- output_fasta: str (optional)
    Name of output file. If not specified name will be the same as input file.

```python
select_genes_from_gbk_to_fasta('input.gbk', ['modE'], 3, 3, 'output.gbk']
# output.gbk
##>galT_before_modE
##MTQFNPVDHPHRRYNPLTGQWILVSPHRAKRPWQGAQETPAKQVLPAHDPDCFLCAGNVRVTGDKNPDYTRTYVFTNDFAALMSDTPDAPESNDPLMRCQSARGTSRVICFSPDHSKTLPELSVAALTEIVKTWQEQTAELGKTYPWVQVFENKGAAMGCSNPHPHGQIWANSFLPNEAEREDRLQKEYFAEQKSPMLVDYVQRELADGSRTVVETEHWLAVVPYWAAWPFETLLLPKAHVLRITDLTDAQRSDLALALKKLTSRYDNLFQCSFPYSMGWHGAPFNGEENQHWQLHAHFYPPLLRSATVRKFMVGYEMLAETQRDLTAEQAAERLRAVSDIHFRESGV
##>galE_before_modE
##MRVLVTGGSGYIGSHTCVQLLQNGHDVIILDNLCNSKRSVLPVIERLGGKHPTFVEGDIRNEALMTEILHDHAIDTVIHFAGLKAVGESVQKPLEYYDNNVNGTLRLISAMRAANVKNFIFSSSATVYGDQPKIPYVESFPTGTPQSPYGKSKLMVEQILTDLQKAQPDWSIALLRYFNPVGAHPSGDMGEDPQGIPNNLMPYIAQVAVGRRDSLAIFGNDYPTEDGTGVRDYIHVMDLADGHVVAMEKLANKPGVHIYNLGAGVGNSVLDVVNAFSKACGKPVNYHFAPRREGDLPAYWADASKADRELNWRVTRTLDEMAQDTWHWQSRHPQGYPD
##>modF_before_modE
##MSSLQILQGTFRLSDTKTLQLPQLTLNAGDSWAFVGSNGSGKSALARALAGELPLLKGERQSQFSHITRLSFEQLQKLVSDEWQRNNTDMLGPGEDDTGRTTAEIIQDEVKDAPRCMQLAQQFGITALLDRRFKYLSTGETRKTLLCQALMSEPDLLILDEPFDGLDVASRQQLAERLASLHQSGITLVLVLNRFDEIPEFVQFAGVLADCTLAETGAKEELLQQALVAQLAHSEQLEGVQLPEPDEPSARHALPTNEPRIVLNNGVVSYNDRPILNNLSWQVNPGEHWQIVGPNGAGKSTLLSLITGDHPQGYSNDLTLFGRRRGSGETIWDIKKHIGYVSSSLHLDYRVSTTVRNVILSGYFDSIGIYQAVSDRQQKLVQQWLDILGIDKRTADAPFHSLSWGQQRLALIVRALVKHPTLLILDEPLQGLDPLNRQLIRRFVDVLISEGETQLLFVSHHAEDAPACITHRLEFVPDGGLYRYVLTKIY
##>acrZ_after_modE
##MLELLKSLVFAVIMVPVVMAIILGLIYGLGEVFNIFSGVGKKDQPGQNH
##>modA_after_modE
##MARKWLNLFAGAALSFAVAGNALADEGKITVFAAASLTNAMQDIATQYKKEKGVDVVSSFASSSTLARQIEAGAPADLFISADQKWMDYAVDKKAIDTASRQTLLGNSLVVVAPKASEQKDFTIDSKTNWTSLLNGGRLAVGDPEHVPAGIYAKEALQKLGAWDTLSPKLAPAEDVRGALALVERNEAPLGIVYGSDAVASKGVKVVATFPEDSHKKVEYPVAVVEGHNNATVKAFYDYLKGPQAAEIFKRYGFTTK
##>modB_after_modE
##MILTDPEWQAVLLSLKVSSLAVLFSLPFGIFFAWLLVRCTFPGKALLDSVLHLPLVLPPVVVGYLLLVSMGRRGFIGERLYDWFGITFAFSWRGAVLAAAVMSFPLMVRAIRLALEGVDVKLEQAARTLGAGRWRVFFTITIPLTLPGIIVGTVLAFARSLGEFGATITFVSNIPGETRTIPSAMYTLIQTPGGESGAARLCIISIALAMISLLISEWLARISRERAGR
```

# Additional scripts

## change_fasta_start_pos

Changes fasta start position.

**Arguments**:

- input_fasta: str
    Path to the input fasta file
- shift: int
    Number of shift
- output_fasta:
    Path to output file

```python
change_fasta_start_pos('input.fasta', 2, 'output)
## input.fasta
# >genome
# ATGCA
## output.fasta
# >genome
# GCAAT
```

## parse_blast_output

Parses blast output. Saves list of the most accurate proteins findings on each query.

**Arguments**:

- input_file: str
    Path to blast output
- output_file: str
    Path to output file

```python
parse_blast_output('input.txt', 'output')
#output.fasta
##DNA methylase [Enterobacteriaceae]
##DUF1380 domain-containing protein [Escherichia coli]
##DUF1380 family protein [Enterobacteriaceae]
##DUF4158 domain-containing protein [Klebsiella pneumoniae]
##DUF905 domain-containing protein, partial [Shigella sonnei]
##DinI-like family protein [Escherichia coli]
##KlcA [Escherichia coli]
##PilI type IV pilus biogenesis protein, partial [Salmonella enterica]
##PilK [Escherichia coli]
##PsiB [Escherichia coli]
```

# rna_dna_tool

Function gets on input a command name and a sequence or an array of sequences. Function delets non nucleotide sequences and sends sequences to subfunction(command). Prints returned value from subfunction.
**Commands**:

## transcribe

```python
rna_dna_tool('AcTg', 'transcribe') # Returns 'AcUg'
```

Function gets on input a sequence or an array of sequences and returns transcribed sequence or an array of sequences

## reverse

```python
 rna_dna_tool('AttG', 'reverse') # Returns GttA
```

Function gets on input a sequence or an array of sequences and returns reverse sequence or an array of sequences

## complement

```python
 rna_dna_tool('AttG', 'AcAc', 'complement') # Returns ['TaaC', 'TgTg']
```

Function gets on input a sequence or an array of sequences and returns complement sequence or an array of sequences

## reverse_complement

```python
rna_dna_tools('AttG', 'AcAc', 'complement')# Returns ['CaaT', 'gTgT']
```

Function gets on input a sequence or an array of sequences and returns reverse complement sequence or an array of sequences

## is_dna

```python
dna_rna_tools('AtAT','AGCA','Augc','is_dna') # returns [True, 'Unclear', False]
```

Function gets on input a sequence or an array of sequences and check if a sequence/ sequences are DNA. Returns True, False, 'Unclear' value or an array of this values.

## length

```python
length(['AtAT','ATGcTG','AugcA']) # returns [4, 6, 5]
```

Function gets on input a sequence or an array of sequences and returns sequence length value or an array of values.

## is_palindrome

```python
is_palindrome(['AtAT','ATGc','Augc']) # returns [True, False, False]
```

Function gets on input a sequence or an array of sequences checks if a sequence/sequences are palindrome and returns True, False value or an array of values.

# protein_tool

The program can process one or more amino acid sequences written in a one-letter format and also does not take into account the size of the input and output letters. Tool can work with amino acids which are mentioned in table below.

| One letter amino acid | Three letter amino acid |
| --------------------- | ----------------------- |
| A                     | Ala                     |
| C                     | Cys                     |
| D                     | Asp                     |
| E                     | Glu                     |
| F                     | Phe                     |
| G                     | Gly                     |
| H                     | His                     |
| I                     | Ile                     |
| K                     | Lys                     |
| L                     | Leu                     |
| M                     | Met                     |
| N                     | Asn                     |
| P                     | Pro                     |
| Q                     | Gln                     |
| R                     | Arg                     |
| S                     | Ser                     |
| T                     | Tre                     |
| V                     | Val                     |
| W                     | Trp                     |
| Y                     | Tyr                     |

## change_abbreviation(seq)

Name's operation: "one letter".
Sequences are written in three-letter format, can be converted to one-letter format. Amino acids must be separeted by "-". 

- seq: Amino acid sequence (str).

- Returns: string of one-letter format of sequence or list of sequences.
  
  ### Example:
  
  ```python
  protein_tool('aLa-CyS', 'one letter') #input ignore letter's size
  # Returns 'AC'
  protein_tool('Ala-Cys', 'Ala', 'one letter')
  # Returns ['AC', 'A']
  ```
  
  ## to_dna(seq)
  
  Name's operation: "DNA".
  Transforms aminoacid sequence to according DNA sequence. 

Arguments:

- sequence: aminoacid sequence to transform into DNA.

Returns:

- String of according DNA sequence.
  
  ### Example:
  
  ```python
  protein_tool('AsDr', 'DNA')
  # Returns 'GCN(TCN or AGY)GAYAGY'
  protein_tool('YWNGAS', 'DNA')
  # Returns'TAY(CGN or AGR)AAYGGNGCN(TCN or AGY)'
  ```
  
  ## to_rna(seq, rna_dict)
  
  Name's operation: "RNA".
  Translates an amino acid sequence into an RNA sequence. 

Arguments:

- seq: Amino acid sequence (str).
- rna_dict: Dictionary defining the correspondence of amino acids to RNA triplets (default, standard code).

Returns: 

- String or list of RNA sequences.
  
  ### Example:
  
  ```python
  protein_tool('FM', 'RNA')
  # Returns'UUYAUG'
  ```
  
  ## define_charge(seq, positive_charge, negative_charge)
  
  Name's operation: "charge".
  Counts the number of amino acids with positive charge, negative charge, and neutral amino acids in the sequence. 

Arguments:

- seq: Amino acid sequence (string).
- positive_charge: List of amino acids with positive charge (default is ['R', 'K', 'H']).
- negative_charge: List of amino acids with negative charge (default is ['D', 'E']).

Returns: 

- Dictionary (or list of dictionaries) containing the counts of amino acids and their labels.

### Example:

```python
protein_tool('ASDRKHDE', 'charge')
# Returns {'Positive': 3, 'Negative': 3, 'Neutral': 2}
```

## define_polarity(seq)

Name's operation: "polarity".
Counts polar and nonpolar aminoacids in sequence. 

Arguments:

- sequence: sequence in which we count polar and nonpolar aminoacids. \newline

Returns:

- Dictionary with keys 'Polar', 'Nonpolar' and appropriate aminoacid counters as values.
  
  ### Example:
  
  ```python
  protein_tool('ASDR', 'polarity')
  # Returns[{'Polar': 3, 'Nonpolar': 1}]
  ```
  
  ## Troubleshooting
  
  Sequences containing characters that do not code for amino acids will be removed from the analysis. The program will write an error and display the sequence with which the problem occurred.
  
  ### Example:
  
  ```python
  protein_tool('ASDR', 'ala1', 'polarity')
  # Prints Something wrong with ala1
  # Returns [{'Polar': 3, 'Nonpolar': 1}]
  ```
  
  Sequences specified in a three-letter format are accepted only by the function for changing the recording format. In other cases, the program will produce either an error or incorrect calculations.
  
  ### Example:
  
  ```python
  protein_tool('AGFHGF', 'Ala-ala', 'DNA') # key error "-"
  protein_tool('AGFHGF', 'Ala-ala', 'polarity')
  # Returns [{'Polar': 1, 'Nonpolar': 5}, {'Polar': 0, 'Nonpolar': 7}] # "-" is counted as non-polar
  ```
