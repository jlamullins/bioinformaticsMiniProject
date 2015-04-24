'''
Created on Apr 14, 2015
@author: Jessica Mullins
'''
import random
import os

charList = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}
nmPositions = [];
dirName = "dataset";
dirNum = 1;

# ML = positive integer for motif length
# NM = non-negative integer for number of variable positions
# SL = positive integer for sequence length
# SC = positive integer for sequence count
def benchmark(ML, NM, SL, SC):
    # Step 2: Generate SC random sequences with uniform nucleotide frequencies. 
    # Each random sequence has length SL
    sequences = createSequences(SC, SL);
    
    # Step 3: Generate a random DNA string (motif) of length ML.
    # Mark a random subset of the NM positions in this motif as 'variable'
    motif = createMotif(ML, NM);

    # Step 4: Generate SC 'binding sites', each of which matches the motif 
    #exactly, except possibly at the positions marked variable. At these variable positions, each binding site will have a randomly chosen nucleotide
    bindingSites = addBindingSites(motif, SC, NM);
    
    # Step 5: "Plant" one sampled site at a random location in each random sequence generated in step 2
    # Planting means overwriting the substring at that location with the site
    sequences = plantSampledSites(sequences, bindingSites, SC, SL, ML);
    
    # Step 6: Write out the SC sequences into a FASTA format file called
    # "sequences.fa"
    directory = "/bioData/" + dirName + str(dirNum) + "/";
    if not os.path.exists(directory):
        os.makedirs(directory);
    writeSequencesToFiles(sequences, directory + "sequences.fa", SC, SL);
    
    # Step 7: In a separate text file called "sites.txt" write down the location of 
    # the planted site in each sequence
    writeSitesToFile(nmPositions, directory + "sites.txt", NM, True);
    
    # Step 8: In a separate file called "motif.txt" write the motif that was
    # generated in step 3.
    writeMotifToFile(motif, directory + "motif.txt", ML);
    
    # Step 9: In a separate file called "motiflength.txt" write down the motif length
    writeToFile(directory + "motiflength.txt", ML);
    
def createSequences(SC, SL):
    sequences = [];
    # Create SC number of sequences
    for i in range(0, SC):
        # Create a random DNA sequence of length SL
        sequences.append(createSequence(SL));
    #  print(sequences);
    return sequences;

def createSequence (len):
    result = [];
    for j in range(0, len):
        randomNumber = random.randint(0,3); # Is this uniformly distributed?
        result.append(charList[randomNumber]);
    return result;

def createMotif(ML, NM):
    motif = createSequence(ML);
    i = 0;
    while i < NM:
        # Mark a random subset of the NM positions as variable
        randomNumber = random.randint(0, ML-1); # Is this uniformly distributed?
        if (motif[randomNumber] != '*'):
            motif[randomNumber] = '*';
            nmPositions.append(randomNumber);
            i = i + 1;
    return motif;       

def addBindingSites(motif, SC, NM):
    scBindingSites = [];
    temp = [];
    for i in range(0, SC):
        temp = list(motif);
        for j in range(0, NM):
            randomNumber = random.randint(0,3);
            temp[nmPositions[j]]  = charList[randomNumber];
        scBindingSites.append(temp)
    return scBindingSites;

def plantSampledSites(sequences, bindingSites, SC, SL, ML):
    for i in range(0, SC):
        j = random.randint(0, SL - ML - 1);
        k = 0;
        for k in range(0, ML -1):
            sequences[i][j+k] = bindingSites[i][k];
            k = k + 1;
        i = i + 1;
    return sequences;

def writeSequencesToFiles(sequences,name, SC, SL ):
    f = open(name, "w"); # opens the file with the given name
    f.write(">" + name + "\n");
    for i in range(0, SC):
        for j in range(0, SL):
            f.write(sequences[i][j]);
        f.write("\n");
    f.close();
    
def writeSitesToFile(nmPositions, name, NM, ifComma):
    f = open(name, "w");
    for i in range(0,NM):
        f.write(str(nmPositions[i]));
        if (ifComma):
            f.write(",");
    f.close();

def writeMotifToFile(motif, name, ML):
    f = open(name, "w");
    f.write("MOTIF1");
    f.write("\t");
    f.write(str(ML));
    f.write("\t");
    for i in range(0, ML):
        f.write(motif[i]);
        
def writeToFile(name, ML):
    f = open(name, "w");
    f.write(str(ML));

if __name__ == '__main__':
    # ML = 8, NM = 1, SL = 500, SC = 10 are the default parameters
    ML = 8;
    NM = 1;
    SL = 500;
    SC = 10;
    i = 0;
    
    # 10 with default
    while i < 10:
        benchmark(ML, NM, SL, SC);
        dirNum = dirNum + 1;
        i = i + 1;
    
    paramList = [0, 2];
    j = 0;
    while j < 2:
        NM = paramList[j];
        i = 0;
        while i < 10:
            benchmark(ML, NM, SL, SC);
            dirNum = dirNum + 1;
            i = i + 1;
        j = j+ 1; 
    NM = 1;      
    paramList = [6, 7];
    j = 0;
    while j < 2:
        ML = paramList[j];
        i = 0;
        while i < 10:
            benchmark(ML, NM, SL, SC);
            dirNum = dirNum + 1;
            i = i + 1;
        j = j+ 1;  
    ML = 8;
    paramList = [5, 20];
    j = 0;
    while j < 2:
        SC = paramList[j];
        i = 0;
        while i < 10:
            benchmark(ML, NM, SL, SC);
            dirNum = dirNum + 1;
            i = i + 1;
        j = j+ 1;  
    
    pass