# Bioinformatics Project 2

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

		print(table[i])

	return table

def max(a, b, c): # takes three options, will return the bext one 

	max_val = a

	if(b > max_val):
		max_val = b

	if(c > max_val):
		max_val = c

	return max_val


def get_alignment(table, s, t): # traverse through values for best alignment

	best = ""
	starting_x = 0
	starting_y = 0
	max_val = table[0][0] # starting point

	# get max value in all of table as starting point

	for i in range(len(table)):
		for j in range(len(table[i])):
			if table[i][j] > max_val:
				max_val = table[i][j]

	print("\nMax: " + str(max_val))
	print("Index: " + str(j) + "," + str(i))


def main():

	s = "CATCACCT"
	t = "GATACCC"
	best = ""

	table = make_table(s, t)

	table = calc_table_vals(table, s, t)

	best = get_alignment(table, s, t)

	print("\n")

main()
