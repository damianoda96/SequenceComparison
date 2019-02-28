# Bioinformatics Project 2
import sys;


def read_in_sequences():
    
    sequence1 = " "
    sequence2 = " "
    not_read = True
    
    while (not_read):
        
        file_name = input("\nEnter the name of the file to load:\n")
        
        try:
            print(os.path.abspath(file_name))
            file = open(file_name, "r")
            
            sequence1 = file.readline().strip() 
            sequence2 = file.readline().strip() 
            
        except:
        
            print("Trouble reading input file. Make sure name is correct and the file is in the right format.");
            
        else:
            
            print("\nFile read successfully.\n")
            not_read = False
            file.close()
        
    retVal = (sequence1, sequence2)
    return retVal;


def print_table_to_file(table):
    
    notSaved = True;
    
    while (notSaved):
    
        file_name = input("\nEnter a name for the alignment table file:\n")
        file_name = file_name + ".txt"
    
        try:
        
            file = open(file_name, "w")
            
            for i in range(len(table)):	
                for j in range(len(table[i])):
                    file.write(table[i][j], " ", end = '')
                file.write()
        
        except:
            
            if ".txt" in file_name:
                print("\nError writing to file.\n")
            else:
                print("\nPlease do not include file extension/type in the entered name.\n")
        
        else:
            
            print("\nTable save successful.\n");
            notSaved = False;
    
    file.close();
    
    return 0;
    
    

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

    for i in range(len(table)):	
        for j in range(len(table[i])):
            print(table[i][j], " ", end = '')
        print()
    
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
    
    print("\nMax: " + str(max_val))
    print("Index: " + str(j) + "," + str(i))

    ########### TODO: starting at max index, get best alignment ##################

    # use rules here

    # something recursive could be cool

    best += traverse_table(table, starting_x, starting_y, s, t, best)
    
    return best

def traverse_table(table, x, y, s, t, best): # recursive table traversal function

	print(table[x][y])

	if(table[x][y] != 0):

		print(t[x-1])

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

    s = "CATCACCT"
    t = "GATACCC"
    userInput = ""

    #### TODO: Allow reading in flu_0.txt to strings in correct format.. ####
    
    #file_data = read_in_sequences()
    
    # simpler to test with above assinments for now
    #s = file_data[0]
    #t = file_data[1]

    # make table of correct size, fill with zeros

    table = make_table(s, t)

    # sift through strings and calc the table's values

    table = calc_table_vals(table, s, t)

    # save table (optional for user)
    
    # user_input = input("Save generated table to file? (y/n)")
    
    # if user_input in ["y", "Y"]:
        # print_table_to_file(table)

    # traverse our table for the best alignment

    best = get_alignment(table, s, t)

    print("Best Alignment: " + best)

    print("\n")



main()



