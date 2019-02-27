def align(s, t):
	# local alignment Smith and Waterman

	print(s)
	print(t)

	vals_x = []

	match = 1
	mismatch = -1
	gap = -2

	for i in range(len(t)):

		vals_y = []

		for j in range(len(s)):
			
			vals_y.append(0)

		vals_x.append(vals_y)
		del vals_y

			# do something here...

	print(s)
	for i in range(len(t)):
		print(t[i] + " "  + str(vals_x[i]))

def main():

	s = "GATCACCT"
	t = "GATACCC"

	align(s, t)

main()
