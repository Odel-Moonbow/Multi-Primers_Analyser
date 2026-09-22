# Multiple Primer Analysis (multi_pa) software: Self-dimer and heterodimer analysis for DNA and RNA oligonucleotides, including degenerate primers.
# This script uses the CS50 library to interact with a SQLite database to store primer information and perform analyses. 

from Bio.SeqUtils import gc_fraction
from colorama import Fore, Back, Style
from cs50 import SQL
from degenerate_primer import degenerate_primer_combinations, calculate_degenerate_tm, deg_self_dimer, self_dimer_check, hetero_dimer_check, deg_hetero_dimer2
from itertools import product
import re
from tabulate import tabulate

# Define a function to prompt the user for a positive integer
def get_pos_int(prompt):
    while True:
        n = input(prompt)
        if n.isdigit() and int(n) > 0:
            return int(n)
        print("Please enter a positive integer.")


# Define IUPAC degenerate nucleotide dictionary
IUPAC_CODES = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
    'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

# Configure CS50 Library to use SQLite database
db = SQL("sqlite:///primers_table.db")

# Define the main function for multiple primer analysis
def multi_pa():
    # Prompt a valid number of primers
    n = get_pos_int("Number of primers: ")


    for _ in range(n):
        # Re-prompt user for a valid name
        while True: 
            name = input("Primer's name: ").strip()             
            if not name:
                print("Please enter a valid name.")             
            elif len(db.execute("SELECT * FROM primers_table WHERE name = ?", name)) > 0:
                print("Please choose a different name. This name already exists.") 
            else:
                break

        # Re-prompt user for a valid sequence
        while True:
            prompt_primer_seq = input(f"{name} sequence: ").strip().replace("-", "")
            if not prompt_primer_seq:
                print("Please enter a sequence")

            #Analyse the valid sequence and store info into a table
            elif re.search("^[ACGTUWSMKRYBDHVN]+$", prompt_primer_seq.upper()):
                primer_seq = prompt_primer_seq.upper()
                length = len(primer_seq)
                
                # Generate all possible combinations (in case a degenerate primer is submitted)
                sequences = degenerate_primer_combinations(primer_seq, IUPAC_CODES)
                # Calculates min, max, and average Tm for a degenerate or non-degenerate primer pool
                results = calculate_degenerate_tm(sequences)
                tm = round(results['avg_tm'], 2)

                # Calculate the GC content for a degenerate or non-degenerate primer pool
                gc_perc = round(gc_fraction(primer_seq, ambiguous='weighted') * 100, 2)
                
                # Populate the table
                db.execute("INSERT INTO primers_table (name, sequence, length, Tm, GC_percentage) VALUES (?, ?, ?, ?, ?)",
                        name, prompt_primer_seq, length, tm, gc_perc)
                break
            else:
                print(f"{name} is invalid. Please enter a valid IUPAC nucleotide code-based sequence")

    # Run the analysis on each primer
    for row in db.execute("SELECT name, sequence, length, Tm, GC_percentage, self_dimer FROM primers_table"):
        print()
        print(Fore.BLUE + Style.BRIGHT + f"{row['name']} primer", end="")
        print(Style.RESET_ALL)
        print(f"Sequence: {row['sequence']}")
        print(f"Length: {row['length']}")
        print(f"Tm: {row['Tm']}°C")
        print(f"GC: {row['GC_percentage']}%")
        print()
        print(Fore.BLUE + Style.BRIGHT + f"{row['name']} self-dimer analysis", end="")
        print(Style.RESET_ALL)
        seq1 = row["sequence"]
        seq1_upper = seq1.upper()
        # For each alignment, generate the lists of total number of matched bases and store them in a variable
        score1 = deg_self_dimer(seq1_upper)
        score2 = deg_self_dimer(seq1_upper[::-1])

        # Display dimers with 3 or more annealing bases in a row only
        if self_dimer_check(score1, score2, seq1) == 0:
            print(f"No self-dimers found in {row['name']}")
            # Populate the table
            db.execute("UPDATE primers_table SET self_dimer = ? WHERE name = ?", "Pass", row['name'])
        else:
            # Populate the table
            db.execute("UPDATE primers_table SET self_dimer = ? WHERE name = ?", "Fail", row['name'])


    # Split primers' names and sequences in two separate lists
    total_seq = []
    total_seq_name = []
    for row in db.execute("SELECT name, sequence FROM primers_table"):
        total_seq.append(row["sequence"])
        total_seq_name.append(row["name"])


    # Analyse hetero-dimers on each possible primers pair in the database
    # Keep track of sequences that form hetero-dimers
    hetero_dimer_forming_seq = []
    for j in range(len(total_seq)):
        for z in range(j + 1, len(total_seq)):
            print()
            print(Fore.BLUE + Back.YELLOW + Style.BRIGHT + f"Cross-dimer analysis between {total_seq_name[j]} and {total_seq_name[z]}", end="")
            print(Style.RESET_ALL)
            score3, score4 = deg_hetero_dimer2(total_seq[j], total_seq[z])
            if hetero_dimer_check(score3, score4, total_seq[j], total_seq[z]) == 0:
                print(f"No hetero-dimers found in {total_seq_name[j]} and {total_seq_name[z]} sequences")
            else:
                hetero_dimer_forming_seq.append(total_seq_name[j])
                hetero_dimer_forming_seq.append(total_seq_name[z])
        # Populate the table
        if total_seq_name[j] not in hetero_dimer_forming_seq:
            db.execute("UPDATE primers_table SET hetero_dimer = ? WHERE name = ?", "Pass", total_seq_name[j])
        else:
            db.execute("UPDATE primers_table SET hetero_dimer = ? WHERE name = ?", "Fail", total_seq_name[j])

    # Display summary table of all primers
    print()    
    print(Style.BRIGHT + "Table of submitted primers.", end="")
    print(Style.RESET_ALL)
    print(tabulate(db.execute("SELECT name, sequence, length, Tm, GC_percentage, self_dimer, hetero_dimer FROM primers_table"), headers='keys', tablefmt='fancy_grid', numalign='right'))
    # Prompt user for more primers
    if input("Would you like to analyse more oligonucleotides? (y/n): ").lower() == "y":
        multi_pa()       
    else:
        print("Analysis completed. Thanks for using the Multiple Primer Analysis (multi_pa) software! Good luck with your research!")
        db.execute("DELETE FROM primers_table")
multi_pa()