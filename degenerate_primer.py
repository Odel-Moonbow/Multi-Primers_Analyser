from itertools import product
from Bio.SeqUtils import MeltingTemp as mt
import re


def main():
# Define IUPAC degenerate nucleotide dictionary
    IUPAC_CODES = {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
        'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
        'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }


# Generates all possible sequence combinations for a degenerate primer.
def degenerate_primer_combinations(primer, d):
    options = [d[base.upper()] for base in primer]
    return [''.join(comb) for comb in product(*options)]

# Calculates min, max, and average Tm for a degenerate primer pool.
def calculate_degenerate_tm(sequences, na_conc=50.0, mg_conc=0.0):
    tm_values = [
        mt.Tm_NN(seq, Na=na_conc, Mg=mg_conc) 
        for seq in sequences
    ]
    
    return {
        "min_tm": min(tm_values),
        "max_tm": max(tm_values),
        "avg_tm": sum(tm_values) / len(tm_values),
        "total_variants": len(sequences)
    }


#return a reverse of the sequence
def rev_st(seq):
    # reverse strand
    seq_r = seq[::-1]
    return seq_r


# Find all possible self-dimers in the degenerate primer and save their score (= matches found) in a list
def deg_self_dimer(seq):
    scores = []
    x = 0
    for z in range(len(seq)):
        score = 0
        for i in range(x, len(seq)):
            x = 0
            if seq[i] in ("A", "W", "M", "R", "D", "H", "V", "N"):
                if seq[-(i + 1 - z)] in ("T", "W", "K", "Y", "B", "D", "H", "N", "U"):
                    score += 1
            elif seq[i] in ("T", "U", "W", "K", "Y", "B", "D", "H", "N"):
                if seq[-(i + 1 - z)] in ("A", "W", "M", "R", "D", "H", "V", "N"):
                    score += 1
            elif seq[i] in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                if seq[-(i + 1 - z)] in ("G", "S", "K", "R", "B", "D", "V", "N"):
                    score += 1
            elif seq[i] in ("G", "S", "K", "R", "B", "D", "V", "N"):
                if seq[-(i + 1 - z)] in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                    score += 1
            elif seq[i] == "U":
                if seq[-(i + 1 - z)] in ("A", "W", "M", "R", "D", "H", "V", "N"):
                    score += 1
            x += 1 + z
        scores.append(score)
    return scores

# Define a function to check self dimers and print only those with matches =>3 in a row
def self_dimer_check(score1, score2, seq1):
    # Group scores > 0 in 2 new lists
    spaces1, spaces2 = get_index2(score1), get_index2(score2)
    # Define a counting variable
    count = 0
    for i in range(len(spaces1)):        
        # Check for consecutive ||| string
        matches_pattern = self_matches2(seq1, spaces1[i])   
        if re.search("\|\|\|\|", matches_pattern) or re.search("\|\|\|", matches_pattern):
            count += 1
            print()
            print(f"Self-dimer number {count}")
            print()
            seq2 = " " * spaces1[i] + rev_st(seq1)
            print(f"5' {seq1} 3'")
            match_assign_ch = print_matches2(seq1, spaces1[i])        
            print(match_assign_ch)
            print(f"3' {seq2} 5'")
        else:
            pass

    # Start from the second element to avoid double counting the first alignment           
    for i in range(1, len(spaces2)):
        # Check for consecutive ||| string
        matches_pattern = self_matches2(rev_st(seq1), spaces2[i])
        if re.search("\|\|\|\|", matches_pattern) or re.search("\|\|\|", matches_pattern):
            count += 1
            print()
            print(f"Self-dimer number {count}")
            print()
            seq2 = " " * spaces2[i] + seq1
            print(f"5' {seq2} 3'")
            match_assign_ch = print_matches2(rev_st(seq1), spaces2[i])
            print(match_assign_ch)
            print(f"3' {rev_st(seq2)} 5'")
        else:
            pass
    return count
    

# Define a function to check hetero dimers and print only those with matches =>3 in a row
def hetero_dimer_check(score1, score2, seq1, seq3):  
    # Group scores > 0 in 2 new lists
    spaces1, spaces2 = get_index2(score1), get_index2(score2)

    # Define a counting variable
    count = 0
    for i in range(len(spaces1)):         
        # Check for consecutive ||| string
        matches_pattern = hetero_matches2(seq1, rev_st(seq3), spaces1[i])
        # Display matches with three or more | in a row      
        if re.search("\|\|\|\|", matches_pattern) or re.search("\|\|\|", matches_pattern):
            count += 1
            print()
            print(f"Hetero-dimer number {count}")
            print()
            seq2 = " " * spaces1[i] + rev_st(seq3)
            print(f"5' {seq1} 3'")
            match_assign_ch = print_hetero_matches2(seq1, rev_st(seq3), spaces1[i])        
            print(match_assign_ch)
            print(f"3' {seq2} 5'")
        else:
            pass
                            
    for i in range(len(spaces2)):
        # Check for consecutive ||| string
        matches_pattern = hetero_matches2(rev_st(seq3), seq1, spaces2[i])
        if re.search("\|\|\|\|", matches_pattern) or re.search("\|\|\|", matches_pattern):
            count += 1
            print()
            print(f"Hetero-dimer number {count}")
            print()
            seq2 = " " * spaces2[i] + seq1
            print(f"5' {seq2} 3'")
            match_assign_ch = print_hetero_matches2(rev_st(seq3), seq1, spaces2[i])
            print(match_assign_ch)
            print(f"3' {rev_st(seq3)} 5'")
        else:
            pass       
    return count

# Define a function that prints "|" characters to visualise matches in all self-dimers (degenerate and non)
def print_matches2(seq, x):
    seq2 = " " * x + seq
    print(" " * x, end="   ")
    match_assign_ch = ""
    for i in range(x, len(seq)):
        if seq[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N") and seq2[-(i + 1 - x)].upper() in ("T", "W", "K", "Y", "B", "D", "H", "N", "U"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() in ("T", "U", "W", "K", "Y", "B", "D", "H", "N") and seq2[-(i + 1 - x)].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N") and seq2[-(i + 1 - x)].upper() in ("G", "S", "K", "R","U","W","K","Y","B","D","H","N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() in ("G", "S", "K", "R", "B", "D", "V", "N") and seq2[-(i + 1 - x)].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() == "U" and seq2[-(i + 1 - x)].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        else:
            match_assign_ch = "".join([match_assign_ch, " "])
    return match_assign_ch

# Define a function that assigns "|" characters to matches found in all self-dimers (degenerate and non)
def self_matches2(seq, x):
    seq2 = " " * x + seq
    match_assign_ch = ""
    for i in range(x, len(seq)):
        if seq[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N") and seq2[-(i + 1 - x)].upper() in ("T", "W", "K", "Y", "B", "D", "H", "N", "U"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() in ("T", "U", "W", "K", "Y", "B", "D", "H", "N") and seq2[-(i + 1 - x)].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N") and seq2[-(i + 1 - x)].upper() in ("G", "S", "K", "R","U","W","K","Y","B","D","H","N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() in ("G", "S", "K", "R", "B", "D", "V", "N") and seq2[-(i + 1 - x)].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        elif seq[i].upper() == "U" and seq2[-(i + 1 - x)].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
            match_assign_ch = "".join([match_assign_ch, "|"])
        else:
            match_assign_ch = "".join([match_assign_ch, " "])
    return match_assign_ch

# Define a function that prints "|" characters to visualise matches in all hetero-dimers (degenerate and non)
def print_hetero_matches2(seq, seq2, x):
    seq2 = " " * x + seq2
    print(" " * x, end="   ")
    match_assign_ch = ""
    for i in range(x, len(seq)):
        try:
            if seq[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N") and seq2[i].upper() in ("T", "W", "K", "Y", "B", "D", "H", "N", "U"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() in ("T", "U", "W", "K", "Y", "B", "D", "H", "N") and seq2[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N") and seq2[i].upper() in ("G", "S", "K", "R","U","W","K","Y","B","D","H","N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() in ("G", "S", "K", "R", "B", "D", "V", "N") and seq2[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() == "U" and seq2[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            else:
                match_assign_ch = "".join([match_assign_ch, " "])
        except IndexError:
                pass
    return match_assign_ch

# Define a function that assigns "|" characters to matches found in hetero-dimers (degenerate and non)
def hetero_matches2(seq, seq2, x):
    seq2 = " " * x + seq2
    match_assign_ch = ""
    for i in range(x, len(seq)):
        try:
            if seq[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N") and seq2[i].upper() in ("T", "W", "K", "Y", "B", "D", "H", "N", "U"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() in ("T", "U", "W", "K", "Y", "B", "D", "H", "N") and seq2[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N") and seq2[i].upper() in ("G", "S", "K", "R","U","W","K","Y","B","D","H","N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() in ("G", "S", "K", "R", "B", "D", "V", "N") and seq2[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            elif seq[i].upper() == "U" and seq2[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                match_assign_ch = "".join([match_assign_ch, "|"])
            else:
                match_assign_ch = "".join([match_assign_ch, " "])
        except IndexError:
                pass
    return match_assign_ch




# Group nucleotide base positions (0-based index) where at least one match is found
def get_index2(scores):
    list = []
    for i in range(len(scores)):
	    if scores[i] > 0:
                list.append(i)
    return list

# Given two sequences, return two lists with scores for all possible hetero-dimers.
def deg_hetero_dimer2(seq, seq2):
    scores1, scores2 = [], []
    x = 0
    for z in range(len(seq)):        
        score1 = 0
        for i in range(x, len(seq)):
            x = 0
            try:
                if seq[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                    if seq2[-(i + 1 - z)].upper() in ("T", "W", "K", "Y", "B", "D", "H", "N", "U"):
                        score1 += 1
                elif seq[i].upper() in ("T", "U", "W", "K", "Y", "B", "D", "H", "N"):
                    if seq2[-(i + 1 - z)].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                        score1 += 1
                elif seq[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                    if seq2[-(i + 1 - z)].upper() in ("G", "S", "K", "R", "B", "D", "V", "N"):
                        score1 += 1
                elif seq[i].upper() in ("G", "S", "K", "R", "B", "D", "V", "N"):
                    if seq2[-(i + 1 - z)].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                        score1 += 1
                elif seq[i].upper() == "U":
                    if seq2[-(i + 1 - z)].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                        score1 += 1
            except IndexError:
                pass
            x = x + 1 + z
        scores1.append(score1)

    # Reset count
    x = 0
    #Compare seq2 vs seq1
    seq2 = seq2[::-1]
    for z in range(len(seq2)):        
        score = 0
        for i in range(x, len(seq2)):            
            x = 0
            try:
                if seq2[i].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                    if seq[i - z].upper() in ("T", "W", "K", "Y", "B", "D", "H", "N", "U"):
                        score += 1
                elif seq2[i].upper() in ("T", "U", "W", "K", "Y", "B", "D", "H", "N"):
                    if seq[i - z].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                        score += 1
                elif seq2[i].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                    if seq[i - z].upper() in ("G", "S", "K", "R", "B", "D", "V", "N"):
                        score += 1
                elif seq2[i].upper() in ("G", "S", "K", "R", "B", "D", "V", "N"):
                    if seq[i - z].upper() in ("C", "S", "M", "Y", "B", "H", "V", "N"):
                        score += 1
                elif seq2[i].upper() == "U":
                    if seq[i - z].upper() in ("A", "W", "M", "R", "D", "H", "V", "N"):
                        score += 1
            except IndexError:
                pass
            x = x + 1 + z
        scores2.append(score)
    return scores1, scores2


if __name__ == "__main__":
    main()