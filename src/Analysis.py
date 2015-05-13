'''
Created on May 11, 2015

@author: SethT
'''
from __future__ import division
import os
import numpy as np


charList = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}

def readMotif(filename, isPredicted):
    f = open (filename);
    motifdata = [];
    result = [];
    for line in f.readlines():
        if (not "<" in line and not ">" in line): 
            motifdata = line.strip().split("\t");
            result.append(motifdata);
    if isPredicted:
    ##    print(result)
        return result;
    else:
        return motifdata;

def convertToPWM(motifData):
    len = motifData[1];
    col = 0;
    row = 0;
    mPWM = [[0 for col in range(int(len))] for row in range(4)];
 ##   print (mPWM)
    motif = motifData[2];
##    print (motif)
    pos = 0;
    for letter in motif:
        if letter == '*':
            mPWM[0][pos] = .25;
            mPWM[1][pos] = .25;
            mPWM[2][pos] = .25;
            mPWM[3][pos] = .25;
        else:
            let = charList[letter]
            mPWM[let][pos] = 1;
        pos = pos + 1;
 ##   print(mPWM)
    return mPWM;

def computeRelativeEntropy(mPWM, pPWM, len):
    ## Build distribution
    summationP = np.sum(pPWM, 0);
    normP = normalize(pPWM, summationP, len);
    ## Probabilities of each nucleotide 
    p = np.asarray(normP, dtype=np.float);
##    print p
    q = np.asarray(mPWM, dtype=np.float);
##    print q
    ##result = np.sum(np.where(p != 0, p * np.log(p / q), 0)) ## should this be log base 2??
    result = kl(p, q, len);   
    return result;
    
def kl(p, q, len):
    result = [];
    for i in range(len):
        val = 0;
        if p[0][i] != 0 and q[0][i] != 0: ## A probability
            val = val + (p[0][i] * np.log(p[0][i]/q[0][i]));
        if p[1][i] != 0 and q[1][i] != 0: ## C probability
            val = val + (p[1][i]*  np.log(p[1][i]/q[1][i]));
        if p[2][i] != 0 and q[2][i] != 0: ## G probability
            val = val + (p[2][i] * np.log(p[2][i]/q[2][i]));
        if p[3][i] != 0 and q[3][i] != 0: ## T probability
            val = val + (p[3][i] * np.log(p[3][i]/q[3][i]));
        result.append(val);
    return result;
 
 
def getResults(mFileName, pFileName):
    ## Read motif data
    motifData = readMotif(mFileName + "/motif.txt", False);
    len = int(motifData[1]);
    ## Read predicted motif
    pPWM = readMotif(pFileName + "/predictedmotif.txt", True)
    
    ## Convert the motif to a PWM format
    mPWM = convertToPWM(motifData);
    pPWM = strToNum(pPWM, len); 

    ## Compute relative entropy
    return computeRelativeEntropy(mPWM, pPWM, len); 


def readSites(filename, tab):
    f = open (filename);
    for line in f.readlines(): 
        if not tab:
            return  line.strip().split(",");
        else:
            return  line.strip().split()

def getOverlap(mFileName, pFileName, SC):
    ## Read motif data
    sitesData = readSites(mFileName + "/sites.txt", False);
    ## Read predicted motif
    pPWM = readSites(pFileName + "/predictedsites.txt", True)
  #  print pPWM
  #  print sitesData
    count = 0.0;
    total = 0.0;
    for j in range (0, SC):
    #    print j
        if int(sitesData[j]) == int(pPWM[j]):
            count = count + 1;
        total = total + 1;
    
 #   print count
 #   print total
    x = count / total;
    return x ## returns overlap percentage

def strToNum(pPWM, len):
    for i in range(4):
        for j in range(len):
            pPWM[i][j] = float(pPWM[i][j])
    return pPWM;

def normalize(pPWM, summationP, len):
    for i in range(4):
        for j in range(len):
            pPWM[i][j] = (pPWM[i][j])/summationP[0]
    return pPWM;

if __name__ == '__main__':
    ## Entropy is - summation of qklog2qk, where qk is the probablilty of k
    ## Using Shannon entropy
    
    motifFileName = "../bioData/dataset";
    predictedMotifFileName = "../results/dataset";
    
    ## Read the default data
    dataSet = 1;
    defaultResults = [];
    defaultOverlappingResults = [];
    for dataSet in range(1,11):
        defaultResults.append(getResults(motifFileName + str(dataSet), predictedMotifFileName + str(dataSet)));
        defaultOverlappingResults.append(getOverlap(motifFileName+str(dataSet), predictedMotifFileName + str(dataSet), 10));
    print ("DEFAULT")
    print ("Relative entropy")
    print defaultResults;
    print ("% of overlapping Sites")
    print defaultOverlappingResults
     
    ## Check when NM is = 0
    dataSet = 11;
    nm0Results = [];
    nm0OverLappingResults = []
    for dataSet in range(11,21):
        nm0Results.append(getResults(motifFileName + str(dataSet), predictedMotifFileName + str(dataSet)));
        nm0OverLappingResults.append(getOverlap(motifFileName+str(dataSet), predictedMotifFileName + str(dataSet), 10));
    print("NM = 0")
    print ("Relative Entropy")
    print nm0Results;
    print ("% of overlapping sites")
    print nm0OverLappingResults;
         
    ## Check when NM is = 2
    dataSet = 21;
    nm2Results = [];
    nm2OverLappingResults = []
    for dataSet in range(21,31):
        nm2Results.append(getResults(motifFileName + str(dataSet), predictedMotifFileName + str(dataSet)));
        nm2OverLappingResults.append(getOverlap(motifFileName+str(dataSet), predictedMotifFileName + str(dataSet), 10));
    print("NM = 2")
    print ("Relative entropy")
    print nm2Results;
    print ("% of overlapping sites")
    print nm2OverLappingResults;
 
    ## Check when ML is = 6
    dataSet = 31;
    ml6Results = [];
    ml6OverLappingResults = []
    for dataSet in range(31,41):
        ml6Results.append(getResults(motifFileName + str(dataSet), predictedMotifFileName + str(dataSet)));
        ml6OverLappingResults.append(getOverlap(motifFileName+str(dataSet), predictedMotifFileName + str(dataSet), 10));
    print("ML = 6")
    print ("Relative Entropy")
    print ml6Results;
    print ("% of overlapping sites")
    print ml6OverLappingResults;
 
    ## Check when ML is = 7
    dataSet = 41;
    ml7Results = [];
    ml7OverLappingResults = []
    for dataSet in range(41,51):
        ml7Results.append(getResults(motifFileName + str(dataSet), predictedMotifFileName + str(dataSet)));
        ml7OverLappingResults.append(getOverlap(motifFileName+str(dataSet), predictedMotifFileName + str(dataSet), 10));
    print("ML = 7")
    print ("Relative Entropy")
    print ml7Results;
    print ("% of overlapping sites")
    print ml7OverLappingResults;
    
    ## Check when SC is = 5
    dataSet = 51;
    sc5Results = [];
    sc5OverLappingResults = []
    for dataSet in range(51,61):
        sc5Results.append(getResults(motifFileName + str(dataSet), predictedMotifFileName + str(dataSet)));
        sc5OverLappingResults.append(getOverlap(motifFileName+str(dataSet), predictedMotifFileName + str(dataSet), 5));
    print("SC = 5")
    print ("Relative Entropy")
    print sc5Results;
    print ("% of overlapping Sites")
    print sc5OverLappingResults;

    ## Check when SC is = 20
    dataSet = 61;
    sc20Results = [];
    sc20OverLappingResults = []
    for dataSet in range(61,71):
        sc20Results.append(getResults(motifFileName + str(dataSet), predictedMotifFileName + str(dataSet)));
        sc20OverLappingResults.append(getOverlap(motifFileName+str(dataSet), predictedMotifFileName + str(dataSet), 20));
    print("SC = 20")
    print ("Relative Entropy")
    print sc20Results;
    print ("% of overlapping Sites")
    print sc20OverLappingResults;

    
    pass