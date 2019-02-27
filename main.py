def make_table(s, t):
	# local alignment Smith and Waterman

	vals = []

	match = 1
	mismatch = -1
	gap = -2

	for i in range(len(t) + 1): # 1 extra for empty row / column

		vals_y = []

		for j in range(len(s) + 1): # 1 extra for empty row / column
			
			vals_y.append(0)

		vals.append(vals_y)
		del vals_y

	# print(s)
	for i in range(len(t) + 1):
		'''if i > 0:
			print(t[i - 1] + " "  + str(vals_x[i]))
		else:
			print("  "  + str(vals_x[i]))'''

		print(vals[i])



def align(table, s, t):
	pass


def main():

	s = "GATCACCT"
	t = "GATACCC"

	table = make_table(s, t)

	align(table, s, t)

main()
