# lines = []
# with open('nltcs.txt', 'r') as f:
# 	for line in f:
# 		temp = line.split('\t')[0]
# 		lines.append(temp)

# with open('nltcs_new.txt', 'w') as f:
# 	for l in lines:
# 		f.write(l + '\n')

arrays = []
with open('nltcs_new.txt', 'r') as f:
	for line in f:
		bit_array = []
		for bit in line.rstrip():
			bit_array.append(int(bit))
		arrays.append(bit_array)

with open('nltcs_new_csv.txt', 'w') as f:
	for array in arrays:
		for bit in array:
			f.write(str(bit) + ', ')
		f.write('\n')
