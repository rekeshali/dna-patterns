# COCS 505 - Project C
# Group 10: rali1, bterry7, twall4, cbrice1
# Correlate disease frequency with matches in DNA sequences.
import sys
import time
import numpy as np
import pandas as pd

# (Cbrice 1): Big O analysis for Rekesh's section 
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
    # However, reading in the data will depend on the size of the data, increasing linearly with the size (number of samples) of the .json file
    # Validating the input and inputting/outputting is a constant time operation
    # Timing: O(n)
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
    # Converting to list depends on the number of patient IDs
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
    diseaseCodes = "abcdeABCD"
    diseaseNames = ["Pancreatic cancer", "Breast cancer", "Lung cancer","Lymphoma", "Leukemia", "Gastro-reflux", \
                    "Hyperlipidemia", "High blood pressure", "Macular degeneration (any degree)"] 
    
    # Analysis:
    # Loop over the number of possible diseases. If this value is fixed, then it is a fixed time.
    # This framework allows for easy addition of diseases; if the number of disease is a variable, the runtime becomes O(n). That is not the case for this project
    # Timing: O(1)
    diseases = {}
    for i in range(len(diseaseCodes)):
        diseases[diseaseCodes[i]] = [diseaseNames[i], 0]
    
    # Analysis:
    # Nested for loop--every disease for every patient.
    # First for loop is O(n)
    # Second for loop iterates O(n) per n in the first; this n is likely to be smaller though
    # Timing: O(n^2). 
    for pid in pids:
        curDisease = emr[pid]
        for code in curDisease:
            diseases[code][1] += 1
    
    ######################################################################################
    # (rali1) ############################################################################
    ######################################################################################
    # Removes empty subdicts within a dict
    # my_func()
    # Time Analysis:(Cbrice1)
    # For loop is O(n)
    #Timing: O(n) end (Cbrice1)
    def dictClean(dict):
        keys = list(dict.keys())
        for key in keys:
            if dict[key] == {}:
                dict.pop(key)

    # Adds a patient ID to permissable categories for a sequence
    # my_func()
    # Analysis:(Cbrice1)
    # For loop is O(n) where n is num diseases of patient
    #Timing: O(n) end (Cbrice1)
    def addPID(pid, emr, matches, size, combo):
        matches[size][combo]['all'].append(pid)
        for diseaseCode in emr:
            if diseaseCode not in matches[size][combo].keys():
                matches[size][combo][diseaseCode] = []
            matches[size][combo][diseaseCode].append(pid)

    ndna = len(dna[pids[0]])
    # Every possible match size is initialized
    # Analysis:(Cbrice1)
    # For loop is O(n) where n is possible sizes
    #Timing: O(n) (Cbrice1)
    matches = {}
    for size in range(3,ndna+1):
        matches[size] = {}

    # Keep track of possible size combinations left to match
    # and how many matches there are for each size
    sizes2Match  = list(range(3,ndna+1))
    matchesCount = list(np.zeros((ndna-2)))

    start = time.time()
    # All IDs except last one will be forward looking for matches
    # Analysis:(Cbrice1)
    # Nested for loop that scales with n patients.
    # First loop iterates first (n-1) pids is O(n)
    # Second for loop iterates (n-1) pids after the first which is O(n) per n in the first.
    # Other loops (index and size) do not scale with # of patients so are O(1)
    #Timing: O(n^2) end (Cbrice1)
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
    ######################################################################################
    #(rali1)#(end)########################################################################
    ######################################################################################
    
    
    #(twall4)#(start)#####################################################################
    #This portion of code  calculates the percantages for correlations and formats data 
    #for use in the report or to be printed  
    
    #intializes the outputText variable to be populated with the "report" formatted information
    outputText ='' 
    
    #Nested for loops --
    # 3rd (most inner) loop: operates on n-#matches --> time = O(n)
    # 2nd (middle) loop: operates on m-#diseases and 3rd loop --> time = O(m)
    # 1st (outer) loop: operates on k-#DNAsequences and inner loops --> time = O(k)
    # total time = O(n*m*k)
    
    #Nested for loops sort down through the nested dictionaries in "matches"
    for key, seq in sorted(matches.items()):
        for dna, sick in sorted(seq.items()):
            #places the DNA sequnces into outputText variable
            outputText = outputText + ("{}:\n".format(dna))            
            for condition, afflictID in sorted(sick.items()):        
                #sets up collection of all patients with DNA sequence
                allID = sick["all"]    
                #creates a dictionary to switch from disease code to full name
                switch = { "a":"Pancreatic cancer", "b":"Breast cancer", "c":"Lung cancer",
                           "d":"Lymphoma", "e":"Leukemia", "A":"Gastro-reflux", "B":"Hyperlipidemia",
                           "C":"High blood pressure", "D":"Macular degeneration (any degree)",
                         "all":"All IDs with sequence" }
                                                           
                #calculates the correlation percentage, makes use of the diseases dictionary defined earlier in the code
                # Correlation = (ppl with disease & DNA)/(ppl with DNA)
                if condition != 'all':
                    per = len(afflictID)/len(allID)
                    if per >= .40 and per < .60 :
                        cor = "Slightly Correlated"
                    elif per >= .60 and per < .80 :
                        cor = "Moderately Correlated"
                    elif per >= .80:
                        cor = "Significantly Correlated"
                    else:
                        cor = ''
                else:
                    cor = ''
                
                #produced formatted output depending on wether a correlation exists or if the entry specifies the "All IDs with seq" option
                if cor == '':
                    outputText = outputText
                else:
                    stringID = ''
                    for pid in afflictID:
                        stringID += str(pid) + ", "
                    outputText = outputText + ('\t{}: {}\n\t\t{}\n'.format(switch[condition], cor, stringID[:-2]))
            stringID = ''
            for pid in allID:
                stringID += str(pid) + ", "
            outputText = outputText + ('\t{}: \n\t\t{}\n'.format(switch["all"], stringID[:-2]))  
    
    #(twall4)#(end)############################################################################
    ###########################################################################################
    
    
    # (bterry7) (set up to be written, did not create text)
    # Print Report, to file or console based on original input
    # Analysis:
    # Write to file. No looping or recursion. Will create file in current directory, so no searching
    # Time to write will depend on size of output string. However, this prints based on matches, not number of samples.
    # If the length of the DNA sequence is known and constant regardless of sample size, there is a set upper limit for this timing that does not change with sample size
    # So, while it can run faster, there is a definitive set maximum time that does not scale with the number of samples (assuming the DNA sequence lenght is set)
    # Timing: O(1)
    if outputToFileFlag:
        print(fileName)
        f = open(fileName,'w')
        f.write(outputText)
        f.close()
    else:
        print(outputText)

    print("# of Patients: %i" % (len(pids)))
    print("Runtime: %5.3fs" % (end-start))    

