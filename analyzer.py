
def validate_sequence(sequence):
    for base in sequence:
        if base not in "ATCG":
            return False
    return True


def count_bases(sequence):

    base_count ={"A": 0, "T": 0, "C": 0, "G": 0}
    for base in sequence:
        if base in base_count:
            base_count[base] += 1
    return base_count


def gc_content(sequence):

    if not validate_sequence(sequence):
        return "Invalid DNA sequence"
    gc_count = sequence.count("G") + sequence.count("C")
    total_count  = len(sequence)
    gc_percentage = (gc_count/ total_count) * 100
    return gc_percentage


def reverse_complement(sequence):
    if not validate_sequence(sequence):
        return "Invalid DNA sequence"
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_comp = "".join(complement[base] for base in reversed(sequence))
    return reverse_comp


def transcribe(sequence):
    if not validate_sequence(sequence):
        return "Invalid DNA sequence"
    rna_sequence = sequence.replace("T", "U")
    return rna_sequence


def translate(sequence):
    if not validate_sequence(sequence):
        return "Invalid DNA sequence"
    codon_table = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}
    rna_sequence = transcribe(sequence)
    protein = []
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        amino_acid = codon_table.get(codon, "Unknown")
        if amino_acid == "STOP":
            break
        protein.append(amino_acid)
    return "-".join(protein)


def find_start_codons(sequence):
    if not validate_sequence(sequence):
        return "Invalid DNA sequence"
    start_codons = []
    for i in range (len(sequence) - 2):
        codon = sequence[i:i+3]
        if codon == "ATG":
            start_codons.append(i)
    return start_codons


def find_orf_from_start(sequence, start_pos):
    if not validate_sequence(sequence):
        return "Invalid DNA Sequence"
    for i in range(start_pos, len(sequence) -2, 3):
        codon = sequence[i:i+3]
        if codon in ["TAA", "TAG", "TGA"]:
            return sequence[start_pos:i+3]
    return sequence[start_pos:]


def find_all_orfs(sequence):
    if not validate_sequence(sequence):
        return "Invalid DNA Sequence"
    orfs = []
    start_positions = find_start_codons(sequence)
    for start in start_positions:
        orf = find_orf_from_start(sequence, start)
        if orf[-3:] in ["TAA", "TAG", "TGA"]:
            orfs.append(orf)
    return orfs

# Finds all ORFs within a single reading frame
def find_orfs_in_frame(sequence, frame):
    if not validate_sequence(sequence):
        return "Invalid DNA Sequence"
    orfs = []
    for i in range(frame, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon == "ATG":
            orf = find_orf_from_start(sequence, i)
            if orf[-3:] in ["TAA", "TAG", "TGA"]:
                orfs.append(orf)
    return orfs


def find_all_orfs_six_frames(sequence):
    if not validate_sequence(sequence):
        return "Invalid DNA Sequence"
    
    results = {}
    rev_comp = reverse_complement(sequence)

    for frame in range(3):
        label = f"Frame +{frame + 1}"
        results[label] = find_orfs_in_frame(sequence, frame)
    
    for frame in range(3):
        label = f"Frame -{frame + 1}"
        results[label] = find_orfs_in_frame(rev_comp, frame)
    
    return results


