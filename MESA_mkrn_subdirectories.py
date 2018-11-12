# I created this script in order to tell Python to enter all the subdirectories in the same folder, compile and run MESA

# coding: utf-8

import os
from os import walk
import subprocess

pathw=os.getcwd()


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


subdirect = []
for i in range(len([x[1] for x in os.walk(pathw)][0])):
	subdirect.append(i)
	subdirect[i] = [x[1] for x in os.walk(pathw)][0][i] 	
	# List of all subdirectories. Change the parameter [0] depending on the needs (i.e., on the depth of the nested grid of subdirectories)
	# It should always be x[1], which corresponds to the directories present in the current work path. x[2] comprises the files present in the current work path
	
subdirect = sorted(subdirect, key=str.lower) 	# Otherwise, capitalized strings go first and the vector becomes a mess

direct = []
for i in range(len(subdirect)):
	direct.append(i)
	direct[i] = pathw + '/' + subdirect[i]


# call_args = [] 	# OK, so this actually works
# call_args2 = []
# for i in range(len(direct)):
# 	call_args.append(i)
# 	call_args2.append(i)
# 	with cd(direct[i]):
# 		#call_args[i] = 'touch ' + subdirect[i] + '.txt'
# 		call_args[i] = './mk'
# 		#call_args[i] = call_args[i].split() # because call takes a list of strings 
# 		subprocess.call(call_args[i])
# 		#call_args2[i] = 'touch ' + subdirect[i] + '.dat'
# 		call_args2[i] = './rn'
# 		#call_args2[i] = call_args2[i].split() # because call takes a list of strings 
# 		subprocess.call(call_args2[i])



process_clean = [] 
process_mk = [] 	
process_rn = []
output_mk = []
#output_rn = []
for i in range(len(direct)):
	process_clean.append(i)
	process_mk.append(i)
	process_rn.append(i)
	output_mk.append(i)
	#output_rn.append(i)
	with cd(direct[i]):
		print direct[i]
		process_clean[i] = subprocess.Popen('. ' + direct[i] + '/clean', shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE) 
		process_clean[i].wait()
		process_mk[i] = subprocess.Popen('. ' + direct[i] + '/mk', shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE) 	# stderr avoids getting a broken pipe
		# IMPORTANT: the absolute path to the mk and the rn files are required. The shell commands ./mk and ./rn are not enough because ./mk will yield "Fortran runtime error: end of the file reached"
		output_mk[i] =  subprocess.check_output('. ' + direct[i] + '/mk', shell=True)
		print output_mk[i]
		#print process[i].returncode # This command always retrieves a bunch of zeros until it is completely done; use check_output instead
		process_mk[i].wait()
		process_rn[i] = subprocess.Popen('. ' + direct[i] + '/rn | tee run_results.log', shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		#process_rn[i] = subprocess.Popen('. ' + direct[i] + '/rn', shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		#output_rn[i] =  subprocess.check_output('. ' + direct[i] + '/rn', shell=True)
		#print output_rn[i]
		process_rn[i].wait()