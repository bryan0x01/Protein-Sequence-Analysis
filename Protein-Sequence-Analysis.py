#!/usr/bin/env python3
# protein FASTA analysis with pandas
# Usage: python protein_pipeline.py proteins.fasta

import argparse
import pandas as pd

def read_fasta(path):
    """Yield (header, seq) for each record in a FASTA file."""
    header = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)

# bio notes ---
HYDROPHOBIC = set(list("AILMFWV"))
CHARGED_POS = set(list("KRH"))
CHARGED_NEG = set(list("DE"))

def aa_fraction(seq, aa_set):
    if not seq:
        return 0.0
    hits = sum(1 for x in seq if x in aa_set)
    return hits / len(seq)

# --- main ---
def main():
    p = argparse.ArgumentParser(description="Protein Sequence Analysis (lengths + simple notes)")
    p.add_argument("--fasta", required=True, help="Protein FASTA file")
    p.add_argument("--out-csv", default="protein_lengths.csv", help="Output CSV with per-sequence info")
    p.add_argument("--out-report", default="protein_report.txt", help="Short text report with summary & notes")
    args = p.parse_args()

    # collect rows
    rows = []
    for header, seq in read_fasta(args.fasta):
        # id is everything up to first space if present
        prot_id = header.split()[0]
        length = len(seq)
        frac_hydrophobic = aa_fraction(seq, HYDROPHOBIC)
        frac_pos = aa_fraction(seq, CHARGED_POS)
        frac_neg = aa_fraction(seq, CHARGED_NEG)
        rows.append({
            "id": prot_id,
            "length": length,
            "frac_hydrophobic": round(frac_hydrophobic, 4),
            "frac_pos_charged": round(frac_pos, 4),
            "frac_neg_charged": round(frac_neg, 4),
        })

    if not rows:
        pd.DataFrame(columns=["id","length","frac_hydrophobic","frac_pos_charged","frac_neg_charged"]).to_csv(args.out_csv, index=False)
        with open(args.out_report, "w") as out:
            out.write("No sequences found.\n")
        return

    # make a small DataFrame
    df = pd.DataFrame(rows)

    # basic descriptive stats (lengths only)
    n = len(df)
    min_len = int(df["length"].min())
    max_len = int(df["length"].max())
    mean_len = float(df["length"].mean())
    median_len = float(df["length"].median())

    # ids for shortest / longest (first match if ties)
    shortest_id = df.loc[df["length"].idxmin(), "id"]
    longest_id  = df.loc[df["length"].idxmax(), "id"]

    #biological notes
    # - if proteins skew short, maybe many peptides/signal-like; if longer, maybe multi-domain
    if mean_len < 150:
        length_note = "Proteins skew short on average (<150 aa). Could indicate many small/peptide-like proteins."
    elif mean_len > 400:
        length_note = "Proteins skew long on average (>400 aa). Could indicate multi-domain or complex proteins."
    else:
        length_note = "Average protein length is in a typical range."

    # hydrophobicity
    avg_hyd = float(df["frac_hydrophobic"].mean())
    if   avg_hyd >= 0.35:
        hyd_note = "Relatively hydrophobic on average; may include many membrane-associated proteins."
    elif avg_hyd <= 0.22:
        hyd_note = "Relatively hydrophilic on average; may include many soluble proteins."
    else:
        hyd_note = "Hydrophobicity looks moderate overall."

    # save CSV
    df.to_csv(args.out_csv, index=False)

    # write text report
    with open(args.out_report, "w") as out:
        out.write("Protein Sequence Analysis (simple)\n")
        out.write("----------------------------------\n")
        out.write(f"Total proteins: {n}\n")
        out.write(f"Min length: {min_len} (id: {shortest_id})\n")
        out.write(f"Max length: {max_len} (id: {longest_id})\n")
        out.write(f"Mean length: {mean_len:.1f}\n")
        out.write(f"Median length: {median_len:.1f}\n")
        out.write("\n")
        out.write("Quick biological notes:\n")
        out.write(f"- {length_note}\n")
        out.write(f"- {hyd_note}\n")
        out.write("\n")
        out.write("Columns in CSV:\n")
        out.write("id,length,frac_hydrophobic,frac_pos_charged,frac_neg_charged\n")

if __name__ == "__main__":
    main()