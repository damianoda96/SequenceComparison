# Bioinformatics Project 2
import sys;
import math

sys.setrecursionlimit(0x100000)


def read_in_sequences(message):
    
    sequence = " "
    not_read = True
    
    while (not_read):
        
        file_name = input(message)
        
        try:
            with open(file_name, "r") as file:
                for line in file:
                    
                    # Skip data other then the sequence itself in the Genebank Coding Sequences FASTA Nucleotide 
                    # file format so that it may be used without prior modification. 
                    if line[0] != ">":
                        sequence += line.strip()
        except:
        
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");
            
        else:
            #print(sequence)
            print("\nFile read successfully.\n")
            not_read = False
            file.close()
    
    return sequence.strip();


def read_in_protien_codon_table(pc_table):
    
    try:

        with open("protien_codon_table.txt", "r") as file:
            
            for line in file:
                
                entry = line.split()
                
                pc_table[entry[0]] = entry[1]
        
        #print(pc_table)
        
    except:
        print("\nCould not read in protien codon table file.\n")
    
    else:
        file.close()
    
    return 0;


def print_table_to_file(table, s, t):
    
    notSaved = True;
    
    while (notSaved):
    
        user_input = input("\n\nEnter a name for the alignment table file:\n")
        file_name = user_input + ".txt"
        value_output = ""
        
        try:
        
            file = open(file_name, "w")
            
            label_start_s = "  "
            file.write(label_start_s)
            label_s = "e" + s
            label_t = "e" + t
            
            for j in range(len(label_s)):
                file.write(label_s[j])
                file.write(" ")
            
            file.write(" \n")
            
            for i in range(len(table)):	
                
                file.write(label_t[i])
                file.write(" ")
                
                for j in range(len(table[i])):
                    file.write(str(table[i][j]) + " ")
                    
                file.write("\n")
        
        except:
        
            if ".txt" in input:
                print("\nError writing to file.\n")
            else:
                print("\nPlease do not include file extension/type in the entered name.\n")
        else:
            
            print("\nTable save successful.\n");
            notSaved = False;
    
    file.close();
    
    return 0;


def print_table(table, s, t):
    
    print(" ", " ", end = '')
    
    letter_label_s = "e" + s
    letter_label_t = "e" + t
    
    for j in range(len(letter_label_s)):
        print(letter_label_s[j], " ", end = '')
    
    for i in range(len(table)):	
        
        print()
        
        print(letter_label_t[i], " ", end = '')
        
        for j in range(len(table[i])):
            print(table[i][j], " ", end = '')
    
    print("\n")


def make_table(s, t):
    
    # local alignment Smith and Waterman

    vals = []

    for i in range(len(t) + 1): # 1 extra for empty row / column

        vals_y = []

        for j in range(len(s) + 1): # 1 extra for empty row / column
            
            vals_y.append(0)
    
        vals.append(vals_y)
        
    del vals_y

    print("\n")

    return vals
    

def calc_table_vals(table, s, t):
    
    # always choose best case, should never be negative
    
    match = 1
    mismatch = -1
    gap = -2
    
    for i in range(len(table)):
        for j in range(len(table[i])):
                
            # some primitive logic
                
            if i == 0 or j == 0: # empty row/column always 0
                table[i][j] = 0
            else: 
                
                val_a = 0 # for diagonal
                val_b = 0 # for up
                val_c = 0 # for left
                
                # here we pick our max value
                
                # for 'diagonal' option
                
                if(t[i - 1] == s[j - 1]):
                    val_a = 1 + table[i-1][j-1] # sets string matches to 1
                elif(t[i - 1 != s[j - 1]]):
                    val_a = -1 + table[i-1][j-1]
                
                # for 'up' option:
                
                val_b = -2 + table[i-1][j]
                
                # for 'left' option:
                
                val_c = -2 + table[i][j-1]
                
                
                # assign the max of the three options
                
                max_val = max(val_a, val_b, val_c)
                
                if(max_val >= 0): # we ignore negative values and just assign 0
                
                    table[i][j] = max_val
                
                else:
                    table[i][j] = 0

    # print out to verify
    #print_table(table, s, t)
    
    return table


def max(a, b, c): # takes three options, will return the max 
    
    max_val = a
    
    if(b > max_val):
        max_val = b
    
    if(c > max_val):
        max_val = c
    
    return max_val


def get_alignment(table, s, t): # traverse through values for best alignment
   
    starting_x = 0
    starting_y = 0
    max_val = table[0][0] # starting point
    best = ["", ""]
    
    # get max value in all of table as starting point
    
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] >= max_val:
                max_val = table[i][j]

                starting_x = i
                starting_y = j
    
    # print("\nMax: " + str(max_val))
    # print("Index: " + str(j) + "," + str(i))

    # Add starting base pair, the build the best aligned sequence from there
    
    best + traverse_table(table, starting_x, starting_y, s, t, best)
    
    best[0] = best[0][::-1].strip()
    best[1] = best[1][::-1].strip()
    
    return best

def traverse_table(table, x, y, s, t, best): # recursive table traversal function

    # print(table[x][y])

    if(table[x][y] > 0):

        #print(t[x-1])

        # print(x)

        val_a = table[x][y]
        val_b = table[x][y]
        val_c = table[x][y]

        next_x = 0
        next_y = 0

        # diagonal

        if(t[x - 1] == s[y - 1]):

            val_a += 1 + table[x-1][y-1]

        else:

            val_a += -1 + table[x-1][y-1]

        # up

        val_b += -2 + table[x-1][y]

        # left

        val_c += -2 + table[x][y-1]

        max_val = max(val_a, val_b, val_c)

        # Set where to move to next, and decide if the two current letters
        # should be matched, or only one choosen and a gap inserted. 

        if(max_val == val_a):

            next_x = x - 1
            next_y = y - 1
            
            if(table[x][y] != 0):
                best[0] += t[next_x]
                best[1] += s[next_y]

        elif(max_val == val_b):

            next_x = x - 1
            next_y = y
            
            if(table[x][y] != 0):
                best[0] += t[next_x]
                best[1] += "-"

        else:

            next_x = x
            next_y = y - 1
            
            if(table[x][y] != 0):
                best[0] += "-"
                best[1] += s[next_y]
            
        traverse_table(table, next_x, next_y, s, t, best)
        
    return best # we want reversed string as best alignment


def analyze_alignment_mutations(best):
    
    pc_table = {}
    
    read_in_protien_codon_table(pc_table)
    
    i = 0;
    length = len(best[0])
    
    synonymous_muts = 0
    non_synonymous_muts = 0
    t_indels = 0
    s_indels = 0
    
    while (i <= length - 3):
        
        indel = False
        
        #print(best[0][i : i + 3], " ", best[1][i : i + 3], '\n')
        
        codon1 = best[0][i : i + 3]
        codon2 = best[1][i : i + 3]
        gap_index1 = codon1.find("-")
        gap_index2 = codon2.find("-")
        
        # if a gap exists, and is new, not an extension of an already detected gap, increase indel count 
        if gap_index1 != -1 and best[0][i + gap_index1 - 1] != "-":
            
            t_indels += 1
            indel = True
            
        if gap_index2 != -1 and best[1][i + gap_index2 - 1] != "-":
            
            s_indels += 1
            indel = True
        
        if indel == False:
        
            # if the codons in this frame both code the same amino acid, synonymous mutation found
            if pc_table[codon1] == pc_table[codon2]:
                        
                synonymous_muts += 1
            
            else:
                non_synonymous_muts += 1
        
        # Move in increments of 3 since checking codons. 
        i += 3
    
    print("\nSynonymous mutations: ", synonymous_muts)
    print("\nNon-Synonymous mutations: ", non_synonymous_muts)
    print("\nT Strand Indels: ", t_indels)
    print("\nS Strand Indels: ", s_indels)
    
    return 0


def print_alignment(best):
    
    print("\nBest Alignment: \n\n")
    
    # setting line sizes to match blast for easy comparison
    chunk_size = 60 # aka chars per line
    alignment_str = ""
    
    chunk_frac = len(best[0]) / chunk_size
    chunk_count = math.ceil(chunk_frac) 
    i = 0
    
    while (i < chunk_count):
        
        current_pos = i * chunk_size
        str_current_pos = str(current_pos) + " "
        index = str_current_pos + " "
        top_strand = best[0][current_pos : current_pos + chunk_size]
        bottom_strand = best[1][current_pos : current_pos + chunk_size]
        
        if i == chunk_count - 1:
            chunk_size = len(best[0]) - current_pos
            
        bond_str = (chunk_size * "|")
        
        alignment_str += index + top_strand + "\n " + (len(str_current_pos) * " ") + bond_str + "\n" + index + bottom_strand + "\n\n"
        
        i += 1
    
    print(alignment_str)
    
    return 0


def main():

    # Swapping sequences to ones n which gaps must be considered to test that/.
    #s = "CATCACCT"
    s = "TGTTACGGG"
    #t = "GATACCC"
    t = "GGTTGACTA"
    userInput = ""
    
    # simpler to test with above assinments for now
    s = read_in_sequences("\nEnter the name of the first sequence file to load:\n")
    t = read_in_sequences("\nEnter the name of the seconde sequence file to load:\n")

    # make table of correct size, fill with zeros

    table = make_table(s, t)

    # sift through strings and calc the table's values

    table = calc_table_vals(table, s, t)

    # save table (optional for user)
    
    # user_input = input("Save generated table to file? (y/n)")
    
    # if user_input == "y" or user_input == "Y":
    #print_table_to_file(table, s, t)

    # traverse our table for the best alignment

    best = get_alignment(table, s, t)

    print_alignment(best)
    print("\nLength: ", len(best[0]))

    print("S/Shanghai length: ", str(len(s)))
    print("T/Ohio length: ", str(len(t)))

    # let's remove best from the s and see what's left...

    #s_without_t = s.replace(best[0], '')
    # t_without_s = t.replace(best, '')

    #print("Remainder: " + s_without_t)
    #print("Remainder Length: " + str(len(s_without_t)))
    
    # find mutations and output the count of each type
    # need s and t for indel count, will pass in params
    analyze_alignment_mutations(best)


    #print("\n")

main()



