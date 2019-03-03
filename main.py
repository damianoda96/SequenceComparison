# Bioinformatics Project 2
import sys;

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
            
    return sequence;


def print_table_to_file(table, s, t):
    
    notSaved = True;
    
    while (notSaved):
    
        user_input = input("\nEnter a name for the alignment table file:\n")
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
    
    # alwasys choose best case, should never be negative
    
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
    # print_table(table, s, t)
    
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
    best = ""
    
    # get max value in all of table as starting point
    
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] > max_val:
                max_val = table[i][j]

                starting_x = i
                starting_y = j
    
    # print("\nMax: " + str(max_val))
    # print("Index: " + str(j) + "," + str(i))

    # recursive function to get best aligned sequence

    best += traverse_table(table, starting_x, starting_y, s, t, best)
    
    return best

def traverse_table(table, x, y, s, t, best): # recursive table traversal function

    # print(table[x][y])

    if(table[x][y] != 0):

        #print(t[x-1])

        best += t[x-1]

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

        if(max_val == val_a):

            next_x = x - 1
            next_y = y - 1

        elif(max_val == val_b):

            next_x = x - 1
            next_y = y

        else:

            next_x = x
            next_y = y - 1


        best = traverse_table(table, next_x, next_y, s, t, best)
        best = best[::-1]

    return best # we want reversed string as best alignment

def main():

    #s = "CATCACCT"
    #t = "GATACCC"
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
    print_table_to_file(table, s, t)

    # traverse our table for the best alignment

    best = get_alignment(table, s, t)

    print("Best Alignment: " + best)
    print(len(alignment))

    print("\n")

main()



