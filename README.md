 anarci_info_cdrs

Simple pipeline to run ANARCI on antibody or TCR sequences and extract gene and CDR information.

------------------------------------------------------------
OVERVIEW
------------------------------------------------------------

This project runs ANARCI on input sequences and extracts:

- V gene
- J gene
- CDR1
- CDR2
- CDR3
- Numbered sequence

------------------------------------------------------------
INPUT
------------------------------------------------------------

Place your input CSV file inside:

data/ab_chains/

The CSV file must contain the following columns:

chain      : chain type (H, L, A, B, etc.)
sequence   : amino acid sequence
species    : species name (human, mouse, etc.)

Example:

chain,sequence,species
A,CAVRDSNYQLIW,human

------------------------------------------------------------
OUTPUT
------------------------------------------------------------

The output will be saved automatically in:

output/

The output file contains:

- original sequence
- V gene
- J gene
- CDR1
- CDR2
- CDR3
- numbering information

------------------------------------------------------------
INSTALLATION
------------------------------------------------------------

Install dependencies using conda (recommended):

conda install -c bioconda hmmer -y
conda install -c conda-forge biopython -y
pip install anarci

Or install ANARCI from source:

git clone https://github.com/oxpig/ANARCI.git
cd ANARCI
python setup.py install

------------------------------------------------------------
USAGE
------------------------------------------------------------

1. Place your CSV file in:

data/ab_chains/

2. Run:

python src/main.py

3. Output will appear in:

output/

------------------------------------------------------------
PROJECT STRUCTURE
------------------------------------------------------------

project/

data/
    ab_chains/
        input.csv

output/

src/
    main.py

------------------------------------------------------------
REQUIREMENTS
------------------------------------------------------------

- Python 3.8+
- ANARCI
- hmmer
- biopython
- numpy
- pandas

------------------------------------------------------------
NOTES
------------------------------------------------------------

Supports antibody and TCR sequences.
Uses ANARCI IMGT numbering scheme.
