#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys   
import difflib
from itertools import permutations
import filecmp    

if(len(sys.argv) != 1):
	print("USAGE: compare.py")
	print("Compares all files in the current folder")
	exit()
		
files = os.listdir(".")

checked = ["compare.py"]
for filename_one in files:
	checked.append(filename_one)
	for filename_two in files:
		if(filename_one == "compare.py" or filename_two in checked):
			continue
		else:
			lines_one = []
			f_one = open(filename_one, "r")
			for line in f_one:
				if(line[0] != '#' and line[0] != '\n'):
					temp = line.split()
					lines_one.append(temp[0])
			lines_two = []
			f_two = open(filename_two, "r")
			for line in f_two:
				if(line[0] != '#' and line[0] != '\n'):
					temp = line.split()
					lines_two.append(temp[0])		
			diff = difflib.ndiff(lines_one, lines_two)
			changes = [l for l in diff if l.startswith('+ ') or l.startswith('- ')]
			remove = []
			check = []
			for change_one in changes:
				check.append(change_one)
				for change_two in changes:
					if(change_two in check):
						continue
					if(change_one.split()[1] == change_two.split()[1] and change_one != change_two):
						if(change_one in changes and change_two in changes):
							remove.append(change_one)
							remove.append(change_two)
					else:
						d = difflib.get_close_matches(change_one, changes)
						for found in d:
							if(found == change_one or len(found.split()[1]) != len(change_one.split()[1])):
								continue
							else:
								perms = [''.join(p) for p in permutations(found.split()[1])]
								if(change_one.split()[1] in perms):
									if(change_one in changes and found in changes and change_one != found):
										remove.append(change_one)
										remove.append(found)
			remove_clean = list(set(remove))
			for delete in remove_clean:
				changes.remove(delete)
			if(len(changes) == 0):
				sys.stdout.writelines(filename_one + "  and " + filename_two + " are the same.\n")

				
	
