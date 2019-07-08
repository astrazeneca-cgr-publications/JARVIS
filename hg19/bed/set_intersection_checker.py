from sys import argv

file1 = argv[1]
file2 = argv[2]

set1 = []

with open(file1) as f:
	content = f.readlines()
set1 = [x.strip() for x in content]
set1 = set(set1)

with open(file2) as f:
	content = f.readlines()
set2 = [x.strip() for x in content]
set2 = set(set2)

intersect = set.intersection(set1, set2)

print(intersect)
print('Number of intersection elements:', len(intersect))
