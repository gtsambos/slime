import os

def simulate_recent_history(file, outFile = "recent-history.trees", logFile = "recent-history.log"):
	"""
	Simulates genomic history from the start of admixture until
	the present day with SLiM.
	"""
	command_line = "slim" + " " + file + " " + "&> " + logFile
	slim_run = os.system(command_line)
	if slim_run != 0:
	    log_file = open(logFile, "r")
	    for line in log_file:
	        print(line)