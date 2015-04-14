'''
Created on Apr 14, 2015
@author: Jessica Mullins
'''
import random

charList = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}

# ML = positive integer for motif length
# NM = non-negative integer for number of variable positions
# SL = positive integer for sequence length
# SC = positive integer for sequence count
def benchmark(ML, NM, SL, SC):
    # Step 1: Generate SC random sequences with uniform nucleotide frequencies. 
    # Each random sequence has length SL
    sequnces = createSequence(SC, SL);
    
    # Step 2: Generate a random DNA string (motif) of length ML.
    # Mark a random subset of the NM positions in this motif as 'variable'
    motif = createSequence(1, ML);
    markPositions(); # TODO

def createSequence(SC, SL):
    sequences = [];
    currSequence = [];
    # Create SC number of sequences
    for i in range(0, SC):
        # Create a random DNA sequence of length SL
        for j in range(SL):
            randomNumber = random.randint(0,3); # Is this uniformly distributed?
            currSequence.append(charList[randomNumber]);
        sequences.append(currSequence);
        currSequence = [];
    #  print(sequences);
    return sequences;

def markPositions():
    return;       
        
if __name__ == '__main__':
    benchmark(1, 2, 3, 4);
    pass