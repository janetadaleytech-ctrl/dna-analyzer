# DNA Sequence Analyzer & ORF Finder

A Python bioinformatics tool that analyzes DNA sequences and identifies Open
Reading Frames (ORFs) are regions of DNA that could code for proteins. Built with 
an interactive Streamlit web interface.

## Overview

This project simulates the core biological process of gene expression from raw DNA sequence all the way to a predicted protein and searches for regions of the
sequence that could realistically encode a gene.

## Features

- **Sequence validation** — checks that input only contains valid DNA bases (A, T, C, G)
- **Base composition analysis** — counts occurrences of each nucleotide
- **GC content calculation** — calculates the percentage of G and C bases in the 
  sequence. This matters biologically because G-C base pairs form 3 hydrogen bonds 
  compared to only 2 for A-T pairs, making GC-rich DNA more thermally stable.
- **Reverse complement generation** — DNA is double-stranded and antiparallel, 
  meaning the second strand runs in the opposite direction and pairs base-for-base 
  with the first (A-T, C-G). This function simulates reading the opposite strand, 
  which matters because genes can be encoded on either strand of DNA.
- **Transcription** — converts a DNA sequence into RNA by replacing every T with U, 
  simulating the first step of gene expression.
- **Translation** — converts an RNA sequence into a chain of amino acids by reading 
  it in groups of three bases (codons), using the genetic code. Translation stops 
  as soon as a stop codon is reached, since that's the signal that tells a ribosome 
  the protein is complete.
- **ORF Finder** — scans a sequence for regions that start with a start codon (ATG) 
  and end with a stop codon (TAA, TAG, or TGA), representing candidate protein-coding 
  regions.
- **Six-frame search** — a single sequence can be read in 3 different reading frames 
  depending on where you start counting groups of 3, and genes can appear on either 
  the forward strand or the reverse complement strand. This tool checks all 6 
  possible frames (3 forward + 3 reverse) to make sure no potential gene is missed.
- **Interactive visualizations** — bar chart of nucleotide composition, pie chart 
  of GC vs AT content, and a chart showing where each ORF sits along the sequence.

## Biology Background

DNA carries genetic information through a process called the **central dogma**: 
DNA → RNA → Protein. First, DNA is transcribed into RNA (a working copy that can 
leave the nucleus). Then RNA is translated into a protein by reading it three bases 
at a time, each 3-base "codon" corresponds to one amino acid, the building blocks 
of proteins.

An **ORF** is a stretch of sequence that starts with a start 
codon (ATG) and ends with a stop codon, with a clean, uninterrupted run of codons 
in between. ORFs are candidate regions that could code for a real protein, which 
makes ORF-finding a basic but genuinely useful gene-finding technique.

Because a sequence can be read starting from three different reading 
frames, and because genes can sit on either DNA strand, a thorough search checks 
all **6 reading frames**, in order to avoid missing real genes.

## How to Run

1. Clone this repository
2. Install dependencies:
pip install -r requirements.txt
3. Run the app:
```
    streamlit run app.py
```

## Project Structure

```
DNA_Analyzer/
├── app.py # Streamlit web interface
├── analyzer.py # Core DNA analysis functions
├── requirements.txt
└── README.md
```
## What I Learned

Building this project taught me core Python fundamentals: functions, loops, dictionaries, string manipulation, and working with external libraries like matplotlib and Streamlit, while also applying them to real bioinformatics concepts like the genetic code, reading frames, and the central dogma of molecular biology. 
It also gave me practice structuring a project properly, separating core logic 
from the user interface, and debugging real errors along the way.

## Example

Input: `ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG`

- GC Content: 56.41%
- Reverse Complement: `CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT`
- Transcribed RNA: `AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG`
- Protein: `Met-Ala-Ile-Val-Met-Gly-Arg`