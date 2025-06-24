import os
from tokenize import String

CODON_SIZE = 3

#-------------------------------------------------------------------------------------------------
# Prompt User for filename with DNA sequences in FASTA format
def getFileName():
    while True:
        filename = input("Please enter your chosen filename: ")
        fileExists = os.path.exists(filename)

        if fileExists:
            print("Your file exists.")
            break
        else:
            print("The file you entered does not exist. Please enter a new filename.")
    return filename

#--------------------------------------------------------------------------------------------------
# Get the minimum ORF from the User as a number (default 50)

def getMinORF():
    minORF = ""

    while type(minORF) != int or minORF < 50:
        try:
            minORF = int(input("Enter a number for the minimum ORF you want to search for. The default is 50 bases "
                               "minimum: "))
        except ValueError:
            print("Please enter a valid number.")

    print("Your minimum ORF is " + str(minORF) + " bases.")
    return minORF

#---------------------------------------------------------------------------------------------------
# Formats the sequence by eliminating whitespaces and lowercase letters.
# Returns a formatted sequence

def formatSequence(sequence_format):
    full_string = "".join(sequence_format)
    words = full_string.split()
    return_sequence = ""
    for i in words:
        return_sequence += i
    return return_sequence.upper()

#----------------------------------------------------------------------------------------------------
# Formats the sequence by eliminating whitespaces and lowercase letters.
# Void Function

def addSequenceToDictionary(dictionary: dict, key_dict, sequence_dict: list):
    dictionary[key_dict] = formatSequence(sequence_dict)

#----------------------------------------------------------------------------------------------------
# Read all the sequences from the file and output a dictionary relating each header to its sequence

def readSequences(seq_file_name):
    myfile = open(seq_file_name, 'r')
    dic = {}
    key_ = ''
    seq_read = []
    for line in myfile:
        if line.startswith(">") and key_ == '':
            key_ = line.split(' ')[0]

        elif line.startswith(">") and key_ != '':
            addSequenceToDictionary(dic, key_, seq_read)
            key_ = line.split(' ')[0]
            seq_read = []

        else:
            seq_read.append(line.rstrip())

    addSequenceToDictionary(dic, key_, seq_read)
    myfile.close()
    return dic

#--------------------------------------------------------------------------------------------------
# Find the start sites given a certain codon offset.
# Returns a list of start sites in order
def findStartSites(sequence_input):
    index = 0
    start_codon = "ATG"
    return_list = []
    while index + CODON_SIZE - 1 < len(sequence_input):
        if sequence_input[index: index + CODON_SIZE] == start_codon:
            return_list.append(index)
        index += 1
    return return_list

#----------------------------------------------------------------------------------------------------
# Find the stop sites.
# Returns a list of start sites in order
def findStopSites(sequence_for_stop_sites):
    index = 0
    stop_codons = ["TAA", "TAG", "TGA"]
    return_list = []
    while index + CODON_SIZE - 1 < len(sequence_for_stop_sites):
        if sequence_for_stop_sites[index: index + CODON_SIZE] in stop_codons:
            return_list.append(index)
        index += 1
    return return_list

#---------------------------------------------------------------------------------------------------
# Prints out the ORFs in a formatted manner.
def printORF(tag: String, frame_number: int, position: int, length: int, sequence_print: String):
    if frame_number > 3:
        position *= -1
    print(tag.strip() + " | " "FRAME = " + str(frame_number) + " " + "POS = " + str(position) + " " + "LEN = " + str(
        length))
    current_index = 0
    #print_count = 0
    sequence_length = len(sequence_print)
    while current_index + 3 <= sequence_length:  # Has a full codon.
        print(sequence_print[current_index:current_index + 3], end=" ")
        current_index += 3
        if (current_index / 3) % 15 == 0:
            print()
    if current_index < sequence_length:
        print(sequence_print[current_index:sequence_length], end=" ")
    print()

#------------------------------------------------------------------------------------------------------
# Get all the forward ORFs from each forward strand frame and store the sequences in a list 
# Return the list
def ProcessORF(all_starts, all_stops, user_minORF_entered):
    ORFs = []
    for key_, starts_ in all_starts.items():
        minStart = 0
        for start_indexes in starts_:
            if start_indexes > minStart:
                stops_ = all_stops[key_]

                for stop_index_ in stops_:
                    iLen = (stop_index_ + 3) - start_indexes
                    if (stop_index_ > start_indexes) and (iLen % 3 == 0):
                        if iLen >= user_minORF_entered:
                            tmpORF = seq[start_indexes:stop_index_ + 3]
                            # print('ORFforward |' + str(start_index) + '|' + str(stop_index) + '|' + seq[start_index:stop_index+3])
                            foundORFs = tmpORF
                            # print(key +  '|' + str(start_index) + '|' + str(stop_index) + '|' + tmpORF)
                            # print(iLen)
                            iFrame = (start_indexes % 3) + 1
                            Pos = start_indexes + 1
                            # print(iLen)
                            # print(Pos)
                            # print(iFrame)
                            # print(key.strip() + " | " + str(iFrame) + " | " + str(Pos) + " | " + str(iLen) + "\n" + tmpORF)
                            printORF(key_, iFrame, Pos, iLen, tmpORF)
                            ORFs.append(foundORFs)
                            minStart = stop_index_ + 3
                            break
                        else:
                            break  # Avoids intermediate stop codons that are too short of an ORF

    return ORFs

#-----------------------------------------------------------------------------------------------------------
# Get the reverse complement of "sequences"
# Return the reverse complement 

def rev_Complement(sequences_to_be_rev):
    allrevsequences = {}
    revcomplement = []

    for key_rev, revsequence in sequences_to_be_rev.items():
        if revsequence:
            flip = str.maketrans('ATGC', 'TACG')
            revcomplement = revsequence.translate(flip)[::-1]  # list slicing to get reverse of string
            allrevsequences[key_rev] = revcomplement

    return allrevsequences

#----------------------------------------------------------------------------------------------------------
# Find the start sites in the reverse complement.
# Returns a list of start sites in order
def findRevStartSites(sequence_start):
    return findStartSites(sequence_start)

#-----------------------------------------------------------------------------------------------------------
# Find the stop sites in the reverse complement.
# Returns a list of start sites in order
def findRevStopSites(sequence_stop):
    return findStopSites(sequence_stop)

#----------------------------------------------------------------------------------------------------------
# Get all the ORFs from each reverse complement frame and store the sequences in a list
# Return the list
def ReverseORF(allrevstarts, allrevstops, user_minORF_entered_rev):
    reverseORFs = []
    for key_rev_orf, starts in allrevstarts.items():
        minStart = 0
        for start_index_rev in starts:
            if start_index_rev > minStart:
                stops_rev = allrevstops[key_rev_orf]

                for stop_index_rev in stops_rev:
                    iLen = (stop_index_rev + 3) - start_index_rev
                    if (stop_index_rev > start_index_rev) and (iLen % 3 == 0):
                        if iLen >= user_minORF_entered_rev:
                            tmpORF = seq[start_index_rev:stop_index_rev + 3]
                            foundORFs = tmpORF
                            #print(key + '|' + str(start_index) + '|' + str(stop_index) + '|' + seq[start_index:stop_index+3])
                            iFrame = (start_index_rev % 3) + 4

                            Pos = (start_index_rev + 1)

                            #print(key.strip() + " | " + str(iFrame) + " | " + "-" + str(Pos) + " | " + str(iLen) + "\n" + tmpORF)
                            printORF(key_rev_orf, iFrame, Pos, iLen, tmpORF)
                            reverseORFs.append(foundORFs)
                            minStart = stop_index_rev + 3
                            break


                            #ORFs = (seq[start_index:stop_index+3])
                            #reverseORFs.append(ORFs)
                        else:
                            break # Avoids intermediate stop codons that are too short of an ORF
    return reverseORFs

#-----------------------------------------------------------------------------------------------
# Homology max score function
def homology_search(seq_1, seq_2):

    #algo based on the following lecture: https://youtu.be/6Udqou3vmng.
    #this will get a homology score based on match and mismatches but will not account for gaps

    cum_score = 0  # cumulative score of the whole

    global_min_score = 0  # min score to compare too

    left_index = 0  # left index to access in query

    right_index = 0  # right index to access in query

    max_dif = 0   # tracks the greatest diff between two bounds

    current_diff = 0    # tracks cumulative score minus global min

    max_score = 0

    # Determines which is the subject and which is the query by length
    if len(seq_2) <= len(seq_1):
        subject = seq_1
        query = seq_2
    else:
        subject = seq_2
        query = seq_1

    # Want to overlay the query over the subject such that we reduce the amount of comparisons making it O(mn)
    # m being the length of the query n being the length of the subject
    for index_j in range(len(subject) - len(query)):
        cum_score = 0
        for index_i in range(len(query)):
            if subject[index_j + index_i] == query[index_i]:
                cum_score += 1
            else:
                cum_score -= 1
            if index_i == 0 and index_j == 0:
                global_min_score = cum_score
            else:
                if cum_score < global_min_score:
                    global_min_score = cum_score
                    left_index = index_j + index_i + 1
                current_diff = cum_score - global_min_score

                if current_diff >= max_dif:
                    right_index = index_j + index_i + 1
                    max_dif = current_diff
                    max_score = cum_score - global_min_score

    print("right index is " + str(right_index))
    print("left index is " + str(left_index))
    length = right_index - left_index

    print("max score is: " + str(max_score) + " with length of: " + str(length))


    print("homologous sequence is: " + subject[left_index:right_index])

#########################################################################################
#########################################################################################
# Function calls 

# call 1
sequenceFile = getFileName()
print(sequenceFile)

# call 2
user_minORF = getMinORF()
print(user_minORF)

# call 3
sequences = readSequences(sequenceFile)

# call 4
all_starts_find = {}
for key, sequence in sequences.items():
    start_sites = findStartSites(sequence)
    # print(key)
    # print(start_sites)

    all_starts_find[key] = start_sites

# call 5
all_stops_found = {}
for key, sequence in sequences.items():
    stop_sites = findStopSites(sequence)
    # print(key)
    # print(stop_sites)

    all_stops_found[key] = stop_sites

# call 6
# print("Indexes for start sites:")
for key, starts in all_starts_find.items():
    # print(key)

    for start_index in starts:
        seq = sequences[key]
        # print(seq[start_index:start_index+3:1])
        # print(start_index)


# call 7 
# print("Indexes for stop sites:")
for key, stops in all_stops_found.items():
    # print(key)

    for stop_index in stops:
        seq = sequences[key]
        # print(seq[stop_index:stop_index+3:1])
        # print(stop_index)

# call 8 
forwardORF = ProcessORF(all_starts_find, all_stops_found, user_minORF)

# call 9
revComplement = rev_Complement(sequences)

# call 10
allrevstarts_found = {}
for key, sequence in revComplement.items():
    start_sites_rev = findRevStartSites(sequence)
    # print(key)
    # print(start_sites)

    allrevstarts_found[key] = start_sites_rev

# call 11
allrevstops_found = {}
for key, sequence in revComplement.items():
    stop_sites_rev = findRevStopSites(sequence)
    # print(key)
    # print(stop_sites)
    allrevstops_found[key] = stop_sites_rev

# call 12
#print("Indexes for reverse complement start sites:")
for key, starts in allrevstarts_found.items():
    #print(key)

    for start_index in starts:
        seq = revComplement[key]
        #print(seq[start_index:start_index+3:1])
        #print(start_index)


# call 13 
#print("Indexes for reverse complement stop sites:")
for key, stops in allrevstops_found.items():
    # print(key)

    for stop_index in stops:
        seq = revComplement[key]
        #print(seq[stop_index:stop_index+3:1])
        #print(stop_index)

# call 14
reverseORF = ReverseORF(allrevstarts_found, allrevstops_found, user_minORF)

# call 15 
case = True

while case:
    choice = input("Do you wish to quit? enter: 'yes'. Otherwise enter 'no' and I will show you homology between your frames: ")
    if choice == "yes" or choice == "y":
        print("The program has ended.")
        break
    else:
        print("Reverse ORF is:\n" + reverseORF[0])
        print('\n\n')
        print("Forward ORF is:\n" + forwardORF[0])
        print('\n\n')
        homology_search(reverseORF[0], forwardORF[0])
        case = False

#############################################################################
#############################################################################
# END OF FILE 
