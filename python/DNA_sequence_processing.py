# Question 1: Prompt user to select certain analysis methods performed on a DNA sequence --------------------------------------------------
import os.path
import re
print("Question 1:")

# A. DNA_Comp
def DNA_Comp(sequence):
    countA = sequence.count('A')
    countT = sequence.count('T')
    countG = sequence.count('G')
    countC = sequence.count('C')
    countN = sequence.count('N')

    print("A Count: " + str(countA) + "\n" + 
          "T Count: " + str(countT) + "\n" +
          "G Count: " + str(countG) + "\n" +
          "C Count: " + str(countC) + "\n" +
          "N (unknown) Count: " + str(countN))

# B. AT_Content
def AT_Content(sequence):
    AT_Number = sequence.count('A') + sequence.count('T')
    AT_Percent = (AT_Number / len(sequence)) * 100
    #print(str(AT_Number))
    #print(str(len(sequence)))
    #print(str(AT_Percent))

    formatAT = ("{:.2f}".format(AT_Percent))
    print("AT Content: " + str(formatAT) + "%")

# C. GC_Content
def GC_Content(sequence):
    GC_Number = sequence.count('G') + sequence.count('C')
    GC_Percent = (GC_Number / len(sequence)) * 100
    #print(str(GC_Number))
    #print(str(len(sequence)))
    #print(str(GC_Percent))

    formatGC = ("{:.2f}".format(GC_Percent))     # format to two decimal places 
    print("GC Content: " + str(formatGC) + "%")

# D. Complement of DNA sequence
def Complement(sequence):
    flip = str.maketrans('ATGC','TACG')
    complement = sequence.translate(flip)
    print("The complement of your DNA sequence is: " + "\n" + complement)
    

# E. Reverse Complement of DNA sequence
def rev_Complement(sequence):
    flip = str.maketrans('ATGC','TACG')
    revcomplement = sequence.translate(flip)[::-1] # list slicing to get reverse of string
    print("The reverse complement of your DNA sequence is: " + "\n" + revcomplement) 
    
#----------------------------------------------------------------------------------------------------


# check if file exists 
while True:
    filename = input("Please enter your chosen FASTA file: ")
    fileExists = os.path.exists(filename)

    if fileExists == True:
        print ("Your file exists.")
        break
    else:
        print("Cannot find the file you specified. Please enter a new filename.")

# store the sequence (without the header) as a string 
sequence = ""


# now open the file and check if it is in FASTA format 
while True:
    with open(filename) as myfile:
        if myfile.readline()[0] != '>':
            print("The file you entered is not in FASTA format.")
        else:
            for line in myfile:
                sequence = sequence + line.strip() # removes whitespace(s) from the beginning and end of string; getting rid of header 
            break

    while True:
        filename = input("Please enter your chosen FASTA file: ") # as a backup, check again if the file exists in case user enters a different file name 
        fileExists = os.path.exists(filename)

        if fileExists == True:
            print ("Your file exists.")
            break
        else:
            print("Cannot find the file you specified. Please enter a new filename.")

print("Here is the sequence from your file: " + "\n" + sequence + "\n\n")
        
# loop a user through each menu item, allowing them to select an item
# exit loop when user selects item corresponding to a 'break' statement

print("A. Calculate DNA composition: This will print to the screen the numbers of A, G, C and T nucleotides, and any unknowns (N’s)." + "\n" + 
"B. Calculate AT content: Prints to the screen the percentage of AT in the sequence." + "\n" +  
"C. Calculate GC content: Prints to the screen the percentage of GC in the sequence." + "\n" +  
"D. Compliment: Prints to the screen the compliment of the DNA sequence." + "\n" + 
"E. Reverse compliment: Prints to the screen the reverse compliment.")
print("\n\n")


while True:
    item = input("Enter 'A', 'B', 'C', 'D', or 'E' from the above list. Enter 'F' to exit the module: ").upper()
    print("\n\n")

    if item == 'A':
        DNA_Comp(sequence)
    elif item == 'B':
        AT_Content(sequence)
    elif item == 'C':
        GC_Content(sequence)
    elif item == 'D':
        Complement(sequence)
    elif item == 'E':
        rev_Complement(sequence)
    elif item == 'F':
        print("Exiting module.")
        break
    else:
        print("Enter letters from the list above.")
        
'''
# Function calls
DNA_Comp(sequence)
AT_Content(sequence)
GC_Content(sequence)
Complement(sequence)
rev_Complement(sequence)
'''
print("\n\n")
###########################################################################################################################
###########################################################################################################################


# Question 2: read in sequences and accession numbers from two files and create separate output files with accession:sequence associations ----------------------------

print("Question 2:")

def formatSeq(seq):         # format each line from the sequences file: get rid of newlines, whitespaces, special characters 
    formattedSeq = ""
    for x in seq:
        x = x.upper()
        formattedSeq = formattedSeq + x
    form = formattedSeq.replace('-','')
    form = form.strip()
    return form

def formatAcc(acc):      # format each line from the accession number file 
    acc = acc.strip()
    return acc
        
        

# open the two files 
Q2Seq = open('sequences.txt')
Q2Acc = open('AccessionNumbers.txt')

# read one line from each file since each line corresponds to a new sequence 
seq = Q2Seq.readline()
acc = Q2Acc.readline()

# loop to read one line from each file during each iteration and format those lines each time with function calls
# while loop terminates when there are no more lines to read in each file 
while(seq and acc):
    seq = formatSeq(seq)
    acc = formatAcc(acc)

    outfile = open(acc + '.txt', 'w') # match the file name with the accession number 
    outfile.write('>' + acc + '\n')
    outfile.write(seq)
    outfile.close()

    seq = Q2Seq.readline()
    acc = Q2Acc.readline()

Q2Seq.close()
Q2Acc.close()

print("Your output files have been created.")
print("\n\n")

#######################################################################################################################    
#######################################################################################################################

# Question 3: Predict size of organism population over a set time period

print("Question 3:")

# function for starting size input validation; do not accept a number less than 2  
def startNum():
    while True:
        startNum = int(input("What is your starting organism sample number? "))
        if startNum < 2:
            print("Your starting number must be greater than or equal to 2.")
        else:
            print("Your starting sample number is: " + str(startNum))
            break
    return startNum


# function for average daily % pop increase; do not accept a negative number 
def popIncrease():
    while True:
        popIncrease = float(input("What is the percent average daily population increase? "))
        if popIncrease < 0:
            print("You must enter a positive percentage value.")
        else:
            print("Your average percent daily population increase is: " + str(popIncrease) + "%")
            break
    return popIncrease

# function for number of days; do not accept a number less than 1
def daysMult():
    while True:
        daysMult = int(input("How many days will your population multiply? "))
        if daysMult < 1:
            print("You must enter at least one day.")
        else:
            print("Your entered number of days is: " + str(daysMult))
            break
    return daysMult



# function calls and catching return values in variables 
pop = startNum()
increase = popIncrease()
days = daysMult()
day = 1 # counter 

print("Day               Organisms")
print("---------------------------")
print(str(day) + "                   " + str(float(pop)))
for x in range(1, days, 1):                                 # start at 1 day, step by 1 up to total number of days entered by user 
    
    day = day + 1
    pop = (pop * (increase / 100)) + pop

    print(str(day) + "                   " + str(float(pop)))

print("\n\n")
###############################################################################################################################        
###############################################################################################################################        
    
# Question 4: Return the position of a sequence where an enzyme cuts ----------------------------------------------

print("Question 4:")

import re
import os.path 


# get the file and the sequence inside it 
while True:
    filename = input("Please enter your chosen FASTA file: ")
    fileExists = os.path.exists(filename)

    if fileExists == True:
        print ("Your file exists.")
        break
    else:
        print("Cannot find the file you specified. Please enter a new filename.")

# store the sequence (without the header) as a string 
Q4sequence = ""


# now open the file and check if it is in FASTA format 
while True:
    with open(filename) as myfile:
        if myfile.readline()[0] != '>':
            print("The file you entered is not in FASTA format.")
        else:
            for line in myfile:
                Q4sequence = Q4sequence + line.strip() # removes whitespace(s) from the beginning and end of string; getting rid of header 
            break

    while True:
        filename = input("Please enter your chosen FASTA file: ")
        fileExists = os.path.exists(filename)

        if fileExists == True:
            print ("Your file exists.")
            break
        else:
            print("Cannot find the file you specified. Please enter a new filename.")

print("Here is the sequence from your file: " + "\n" + Q4sequence + "\n\n")

# get the name of the restriction enzyme
REname = input("Enter the name of your restriction enzyme: ")

# open enzyme file and store content as a dictionary where keys are enzyme names and values are cut sites 
myfile = open("ResEnzymes.txt", 'r')
dictionary = {}
for line in myfile:
    acc, value = line.strip().split(',')
    val1 = value.strip().replace("'","").replace('N',"[A|T|G|C]") # replace IUB codes with their associated base pair possibilities (a list of possibilities)  
    val2 = val1.replace('W',"[A|T]")
    val3 = val2.replace('V',"[A|C|G]")
    val4 = val3.replace('H',"[A|C|T]")
    val5 = val4.replace('D',"[A|T|G]")
    val6 = val5.replace('B',"[C|G|T]")
    val7 = val6.replace('S',"[G|C]")
    val8 = val7.replace('M',"[A|C]")
    val9 = val8.replace('K',"[G|T]")
    val10 = val9.replace('Y',"[C|T]")
    val11 = val10.replace('R',"[A|G]")
    
    dictionary[acc.strip()] = val11 
    
myfile.close()
#print(dictionary)


# function to get all the cut sites associated with the start of the restriction fragment match  
def reSites(Q4sequence, REname):
    resites = []
    for site in re.finditer(REname, Q4sequence):
        resites.append(site.start() + 1) # +1 to make sure count is not starting at 0 

    return resites

# function call 
RE = reSites(Q4sequence, dictionary[REname])

print("Cut site(s) found at position(s): " + str(list(RE)))

print("\n\n")
################################################################################################################
################################################################################################################

# Question 5: Read in a whole genome and display the background codon frequencies for all codons accumulated from each reading frame ---------------------------------------------

import os.path 
print("Question 5: ")

# get the file and the sequence inside it 
while True:
    filename = input("Please enter your chosen genome file (in FASTA format): ")
    fileExists = os.path.exists(filename)

    if fileExists == True:
        print ("Your file exists.\n\n")
        break
    else:
        print("Cannot find the file you specified. Please enter a new filename.")

# store the sequence (without the header) as a string 
Q5sequence = ""


# now open the file and check if it is in FASTA format 
while True:
    with open(filename) as myfile:
        if myfile.readline()[0] != '>':
            print("The file you entered is not in FASTA format.")
        else:
            for line in myfile:
                Q5sequence = Q5sequence + line.strip() # removes whitespace(s) from the beginning and end of string; getting rid of header 
            break

    while True:
        filename = input("Please enter your chosen FASTA file: ")
        fileExists = os.path.exists(filename)

        if fileExists == True:
            print ("Your file exists.\n\n")
            break
        else:
            print("Cannot find the file you specified. Please enter a new filename.")

#print(str(Q5sequence))
print("The length of your genome is: " + str(len(Q5sequence)) + " base pairs\n\n")


# create a dictionary where the codons are the keys and the values are initialized counts starting at zero
# as we read through our sequence, counts will be accumulated for each key encountered 
dictionary = {'AAA': 0, 'AAT': 0, 'AAC': 0, 'AAG': 0,

             'AGA': 0, 'AGG': 0, 'AGC': 0, 'AGT': 0,

             'ATA': 0, 'ATG': 0, 'ATC': 0, 'ATT': 0,

             'ACA': 0, 'ACT': 0, 'ACG': 0, 'ACC': 0,


             'CCA': 0, 'CCT': 0, 'CCG': 0, 'CCC': 0,

             'CTA': 0, 'CTT': 0, 'CTG': 0, 'CTC': 0,

             'CAA': 0, 'CAT': 0, 'CAG': 0, 'CAC': 0,

             'CGA': 0, 'CGT': 0, 'CGG': 0, 'CGC': 0,

             
             'GGG': 0, 'GGA': 0, 'GGT': 0, 'GGC': 0,

             'GAA': 0, 'GAT': 0, 'GAG': 0, 'GAC': 0,

             'GCA': 0, 'GCT': 0, 'GCG': 0, 'GCC': 0,

             'GTG': 0, 'GTA': 0, 'GTC': 0, 'GTT': 0,

             
             'TTT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0,

             'TAA': 0, 'TAG': 0, 'TAC': 0, 'TAT': 0,

             'TGA': 0, 'TGT': 0, 'TGC': 0, 'TGG': 0,

             'TCC': 0, 'TCA': 0, 'TCT': 0, 'TCG': 0}


# Now get the reverse complement
flip = str.maketrans('ATGC','TACG')
revcomplement = Q5sequence.translate(flip)[::-1] # list slicing to get reverse of string
#print("The reverse complement of your DNA sequence is: " + "\n" + revcomplement)


count  = 0
while (count < 3):

    # use a 'for loop' to accumulate each codon match in the sequence
    for x in range(count, len(Q5sequence) - 2, 3):
        
        codon = Q5sequence[x: x + 3]             # iterate each key instance into codon 
    
        if codon in dictionary:
            dictionary[codon] = dictionary[codon] + 1 # accumulate the value into each key; the value is number of times each key appears in the sequence

    for x in range(count, len(revcomplement) - 2, 3):
    
        codon = revcomplement[x: x + 3]             # iterate each key instance into codon 
    
        if codon in dictionary:
            dictionary[codon] = dictionary[codon] + 1 # accumulate the value into each key; the value is number of times each key appears in the sequence

    count = count + 1      

#print(dictionary)


# calculate the total number of codons in each frame for both plus and minus strands 
total_codons = dictionary.values()
total = sum(total_codons)

print("The total number of codons in your genome, for all 6 reading frames, is: " + str(total))
print('\n\n')
print("List of background frequencies for each codon:\n")

# get the background frequencies 
for x in dictionary:
    codons = dictionary[x]
    frequ = (100 * codons)/total

    print(x + ": " + "{:.4f}".format(frequ))

print("\n\n")
# get the average for each codon 
av = 0
av2 = 0
for val in dictionary.values():
    av = av + val
    av2 = av2 + val

#print(av)
#print(len(dictionary))
av = av / len(dictionary)
   
print("Average number of times each codon appears in the genome: " + str(av))

# get the average for each reading frame
av2 = av2 / 6
print("Average number of codons appearing in each frame: " + str(av2))











       






 

    

   


 

    

 

    














