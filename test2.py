import sys 
import pandas as pd


# (bterry7) Read in and parse the data into dictionary
# Ensure that the expected number of arguments were input

# Analysis:
# No loops
# Taking user input is a constant time operation
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
    # 1 loop over the indexes (patient IDs). Loop is O(n)
    # Timing: O(n)
    indexes = list(sourceData.index)
    data = {}

    for i in range(len(indexes)):
        data[indexes[i]] = list(sourceData.values[i])



    # (bterry7)
    # Create a dictionary of diseases, for reference
    # Second number will be total number of people with that disease to be used in correlation calcualtions
    # Analysis:
    # Nested for loop--every disease for every patient.
    # First for loop is O(n)
    # Second for loop iterates O(n) per n in the first; this n is likely to be smaller though
    # Timing: O(n^2). 
    diseases = {'a': ['Pancreatic cancer', 0],
                'b': ['Breast cancer', 0],
                'c': ['Lung cancer', 0],
                'd': ['Lymphoma', 0],
                'e': ['Leukemia', 0],
                'A': ['Gastro-reflux', 0],
                'B': ['Hyperlipidemia', 0],
                'C': ['High blood Pressure', 0],
                'D': ['Macular dengeneration (any degree)', 0]}

    for i in data.keys():
        curDisease = data[i][0]
        for j in curDisease:
            diseases[j][1] += 1


    outputText = 'test'
    
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
