import sys
import time
import numpy as np
import pandas as pd


# (bterry7) Read in and parse the data into dictionary
# Ensure that the expected number of arguments were input

# Analysis:
# No loops, just reading and storing data
# Timing: O(1)
if len(sys.argv) != 3:
    print('Usage: ./projectc.py <URL> <output file or ->')

else:
    # Set flags and information for writing the data to text file or console later
    if str(sys.argv[2]) == "-":
        outputToFileFlag = 0
    else:
        outputToFileFlag = 1
        fileName = sys.argv[2] + '.txt'

    # Try to read in data and report exception as needed
    # Analysis:
    # While there is a loop, it's not for iteration, just validation
    # Timing: O(1)
    while True:
        try:
            sourceData = pd.read_json(sys.argv[1])
            break
        except:
            print('Invalid URL provided \nPlease provide valid URL or quit (type "quit")')
            a = input('Input URL ')
            if a == "quit":
                sys.exit()
            else:
                sys.argv[1] = a
        
    # Reformat Pandas into dictionary where
    # Analysis:
    # 1 loop over the indexes (patient IDs). Loop is O(n)
    # Timing: O(n)
    pids = list(sourceData.index) # patient IDs
    emr = {}
    dna = {}

    for (p,pid) in enumerate(pids):
        emr[pid] = list(sourceData.values[p])[0]
        dna[pid] = list(sourceData.values[p])[1]

    # (bterry7 & rali1)
    # Create a dictionary of diseases, for reference
    # Second number will be total number of people with that disease to be used in correlation calcualtions
    # Analysis:
    # Nested for loop--every disease for every patient.
    # First for loop is O(n)
    # Second for loop iterates O(n) per n in the first; this n is likely to be smaller though
    # Timing: O(n^2). 
    diseaseCodes = "abcdeABCD"
    diseaseNames = ["Pancreatic cancer", "Breast cancer", "Lung cancer","Lymphoma", "Leukemia", "Gastro-reflux", \
                    "Hyperlipidemia", "High blood pressure", "Macular degeneration (any degree)"] 
    
    diseases = {}
    for i in range(len(diseaseCodes)):
        diseases[diseaseCodes[i]] = [diseaseNames[i], 0]

    for pid in pids:
        curDisease = emr[pid]
        for code in curDisease:
            diseases[code][1] += 1
   
    ######################################################################################
    # (rali1) ############################################################################
    ######################################################################################
    # Removes empty subdicts within a dict
    def dictClean(dict):
        keys = list(dict.keys())
        for key in keys:
            if dict[key] == {}:
                dict.pop(key)

    # Adds a patient ID to permissable categories for a sequence
    def addPID(pid, emr, matches, size, combo):
        matches[size][combo]['all'].append(pid)
        for diseaseCode in emr:
            if diseaseCode not in matches[size][combo].keys():
                matches[size][combo][diseaseCode] = []
            matches[size][combo][diseaseCode].append(pid)

    ndna = len(dna[pids[0]])
    # Every possible match size is initialized
    matches = {}
    for size in range(3,ndna+1):
        matches[size] = {}

    # Keep track of possible size combinations left to match
    # and how many matches there are for each size
    sizes2Match  = list(range(3,ndna+1))
    matchesCount = list(np.zeros((ndna-2)))

    start = time.time()
    # All IDs except last one will be forward looking for matches
    for (p,pid1) in enumerate(pids[:-1]):
        sizesLeft = sizes2Match.copy()
        # Check sequences at indexes that can fit remaining sizes to check
        for i in range(ndna-sizes2Match[0]+1):
            # Don't check sequences that exist beyond the final character
            if i + sizesLeft[-1] > ndna:
                sizesLeft = sizesLeft[:-1]
            # All subsequent IDs are potential matches for a sequence at this ID's i
            pidsLeft = pids[p+1:]
            # Check the sizes, w, that are left and can fit at index i
            for size in sizesLeft:
                combo = dna[pid1][i:i+size] # sequence to match
                # If we've already checked that combination for matches in a list of sub IDs
                if combo in matches[size].keys():
                    continue
                # For the IDs that have not been removed from potential match list
                for pid2 in pidsLeft.copy():
                    # Check for a match for combo in pid2
                    if combo in dna[pid2]:
                        # Add sequence to matches if not already there
                        if combo not in matches[size].keys():
                            # Introduce the combo to matches
                            matchesCount[size-3]       += 1
                            matches[size][combo]        = {}
                            matches[size][combo]['all'] = []
                            # Record ID where combo was pulled from
                            addPID(pid1, emr[pid1], matches, size, combo)
                            # Remove size from search list if all combos of size are found
                            if matchesCount[size-3] == 2**size:
                                sizes2Match.remove(size)                      
                        # Add matched ID to reported matches
                        addPID(pid2, emr[pid2], matches, size, combo)
                    # Otherwise remove patient from comparison at that index because
                    # if a sequence at index i is not matched, then no larger sequences
                    # at index i will match because they contain the smaller sub-sequence
                    else:
                        pidsLeft.remove(pid2)
    # Remove empty dicts (unmatched sizes)
    dictClean(matches)
    end = time.time()
    print("# of Patients: %i" % (len(pids)))
    print("Runtime: %5.3fs" % (end-start))    
    ######################################################################################
    #(rali1)#(end)########################################################################
    ######################################################################################

    # (bterry7) (set up to be written, did not create text)
    # Print Report, to file or console based on original input
    # Analysis:
    # Write to file. No looping or recursion. Will create file in current directory, so no searching
    # Timing: O(1)
    outputText = 'test'
    if outputToFileFlag:
        print(fileName)
        f = open(fileName,'w')
        f.write(outputText)
        f.close()
    else:
        print(outputText)
