import sys
import subprocess
import os
import re


class Translation:
    dna = ""
    genetic_code = {}
    abbreviation = {}

    def __init__(self):
        global genetic_code
        genetic_code = {
            # Phe
            'UUU': 'Phe', 'UUC': 'Phe',
            # Leu
            'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
            # Ser
            'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
            # Tyr
            'UAU': 'Tyr', 'UAC': 'Tyr',
            # STOP
            'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
            # Cys
            'UGU': 'Cys', 'UGC': 'Cys',
            # Trp
            'UGG': 'Trp',
            # Pro
            'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
            # His
            'CAU': 'His', 'CAC': 'His',
            # Gln
            'CAA': 'Gln', 'CAG': 'Gln',
            # Arg
            'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
            # Ile
            'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
            # Met
            'AUG': 'Met',
            # Thr
            'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
            # Asn
            'AAU': 'Asn', 'AAC': 'Asn',
            # Lys
            'AAA': 'Lys', 'AAG': 'Lys',
            # Val
            'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
            # Ala
            'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
            # Asp
            'GAU': 'Asp', 'GAC': 'Asp',
            # Glu
            'GAA': 'Glu', 'GAG': 'Glu',
            # Gly
            'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
        }
        global abbreviation
        abbreviation = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q',
            'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
            'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'STOP': 'STOP'
        }

    @staticmethod
    # Converts DNA to Protein
    def dna_to_protein(in_dna, title):
        # Transcription (DNA to RNA)
        print("Transcribing frame number ", title, " ...")
        rna = in_dna.replace('T', 'U')

        # Translation (RNA to Protein)
        my_start = my_end = 0
        protein = ""
        print("Translating frame number ", title, " ...")
        while (len(rna) - my_end) / 3 >= 1:
            my_end = my_start + 3
            codon = rna[my_start:my_end]
            amino = genetic_code[str(codon)]
            protein += abbreviation[amino]
            my_start = my_end
        return protein

    @staticmethod
    # TODO: This method might be redundant
    # Converts files into a format that is compatible with HMMER
    def fasta_format(title, sequence):
        counter = 1
        fasta = ""
        for char in str(sequence):
            if (counter % 50) == 0:
                fasta += char + '\n'
            else:
                fasta += char
            counter += 1
        return title + "\n" + fasta

    @staticmethod
    # Computes a score for a set of hits considering:
    # overlap as penalty, and match as positive score
    def scoring_function(current):
        current = list(map(list, current))
        overlap = 0
        if len(current) == 1:
            return current[0][1] - current[0][0]

        for i in range(0, len(current)):
            for j in range(i + 1, len(current)):
                if current[i][1] <= current[j][0] or current[i][0] >= current[j][1]:
                    if current[i][1] <= current[j][0]:
                        overlap += current[j][0] - current[i][1] + 1
                    else:
                        overlap += current[i][0] - current[j][1] + 1
                elif current[i][0] <= current[j][0] and current[i][1] <= current[j][1]:
                    overlap += current[i][1] - current[j][0] + 1
                elif current[i][0] >= current[j][0] and current[i][1] <= current[j][1]:
                    overlap += current[i][1] - current[i][0]
                elif current[i][0] >= current[j][0] and current[i][1] >= current[j][1]:
                    overlap += current[j][1] - current[i][0] + 1
                elif current[i][0] <= current[j][0] and current[i][1] >= current[j][1]:
                    overlap += current[j][1] - current[j][0]
                else:
                    print("UNKNOWN EXCEPTION!")

        coverage = 0
        for i in range(0, len(current)):
            coverage += current[i][1] - current[i][0]

        return coverage - (2*overlap)

    # Creates a binary tree in which a hit is selected or not: 2^hits
    def naive(self, hits, current, n):
        # Base Case
        if n == 0:
            return self.scoring_function(current), current

        # return the maximum of two cases:
        # (1) n-1 th item not included
        # (2) n-1 th item included
        old_current = current.copy()
        current.add(hits[n - 1])
        left_score, left_set = self.naive(hits, current, n - 1)
        right_score, right_set = self.naive(hits, old_current, n - 1)

        if left_score > right_score:
            return left_score, left_set
        else:
            return right_score, right_set

    def hmmsearch(self):
        # Reading HMM protein profile and the DNA sequence as arguments
        print("Reading the protein file ...")
        protein_profile_file = sys.argv[1]

        print("Reading the dna file ...")
        dna_sequence_file = open(sys.argv[2], 'r')

        # Reading DNA sequence line by line
        self.dna = dna_sequence_file.read().split('\n')
        title = self.dna[0]
        self.dna = self.dna[1:]
        new_dna = ""
        for line in self.dna:
            new_dna += line

        # Creating the protein_file with 6 different sequences (ORFs) with 6 different recognizable titles
        protein_fasta = str()
        reverse = 1
        for frame_number in range(3):
            protein = self.dna_to_protein(new_dna[frame_number::reverse], frame_number)
            protein_fasta = protein_fasta + self.fasta_format('>' + str(frame_number) + '_' + title.replace('>', '') + " [reading frame number: " + str(frame_number) + ", reverse: No" + "]", protein) + '\n'
        # The reverse 3 ORFs is not yet needed.
        # reverse = -1
        # for frame_number in range(3):
        #     protein = self.dna_to_protein(new_dna[reverse*(frame_number+1)::reverse])
        #     protein_fasta = protein_fasta + self.fasta_format('>' + str(frame_number) + 'r' + '_' + title.replace('>', '')
        #                                                       + " [reading frame number: " + str(frame_number) + ", reverse: Yes" + "]", protein) + '\n'

        # Writing the protein_file to a file for HMMER use
        protein_file = open("protein.fa", 'w+')
        protein_file.write(protein_fasta)
        seq_path = os.path.abspath("protein.fa")
        protein_file.seek(0, 0)

        print("Running HMMER ...")
        process = subprocess.run(["hmmsearch", protein_profile_file, seq_path], stdout=subprocess.PIPE,
                                 universal_newlines=True)
        # process = subprocess.run(["hmmsearch", protein_profile_file, "/home/ahhajianpour/Downloads/hmmer-3.2/tutorial/some1.fa"], stdout=subprocess.PIPE,
        #                          universal_newlines=True)

        # Reading HMMER output and extracting necessary information
        output = process.stdout
        line_counter = 0
        current_hit = 0
        number_of_hits = 0
        hmm = list()
        alignment = list()
        line_by_line = output.split('\n')
        for line in line_by_line:
            reformed_line = re.sub(' +', ' ', line).strip().split(' ')
            # Finding the number of hits
            if reformed_line[0] == 'Scores' and reformed_line[1] == 'for' and reformed_line[2] == 'complete' and reformed_line[3] == 'sequences':
                read = line_counter + 4
                if len(re.sub(' +', ' ', line_by_line[read]).strip().split(' ')) < 9:
                    return "No hit found!", "No hit found!"
                while len(re.sub(' +', ' ', line_by_line[read]).strip().split(' ')) >= 9:
                    number_of_hits += int(re.sub(' +', ' ', line_by_line[read]).strip().split(' ')[7])
                    read += 1
                hmm = list(range(0, number_of_hits))
                alignment = list(range(0, number_of_hits))
            # Finding hits in each sequence
            if reformed_line[0] == '>>':
                read = line_counter + 3
                # Reading the start and the end index of each hit
                while len(re.sub(' +', ' ', line_by_line[read]).strip().split(' ')) == 16:
                    hmm[current_hit] = list(map(int, re.sub(' +', ' ', line_by_line[read]).strip().split(' ')[6:8]))
                    alignment[current_hit] = list(map(int, re.sub(' +', ' ', line_by_line[read]).strip().split(' ')[9:11]))
                    read += 1
                    current_hit += 1
            line_counter += 1
        return hmm, alignment


some = Translation()
hmm , alignment = some.hmmsearch()
if hmm != "No hit found!":
    my_out = some.naive(list(map(tuple, hmm)), set(), len(hmm))
    print("Highest score is: ", my_out[0], ", and the best set is: ", my_out[1])
else:
    print(hmm)