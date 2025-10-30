# Protein Sequence Analysis

A simple Python project that analyzes protein sequences from a FASTA file using Pandas.  
It calculates basic statistics about sequence lengths and generates outputs for data exploration.

## Features

- Reads sequences from a `.fasta` file
- Calculates the length of each protein sequence
- Exports a `.csv` file with IDs and sequence lengths
- Generates a text report summarizing:
  - Minimum sequence length
  - Maximum sequence length
  - Mean sequence length
  - Median sequence length

## Tech Stack

- Python 3
- Pandas

## How to Run

1. Install dependencies:
   ```bash
   pip install pandas
   ```

2. Run the script:
   ```bash
   python Protein-Sequence-Analysis.py
   ```

4. The program will create:
   - `protein_lengths.csv` — containing each sequence ID and its length
   - `summary_report.txt` — with basic statistical summaries

## Example Output

**CSV:**
| ID | Length |
|----|---------|
| sp|P12345|  | 350 |
| sp|Q67890|  | 289 |

**Report:**
```
Total sequences: 250
Min length: 45
Max length: 1012
Mean length: 276.8
Median length: 251
```