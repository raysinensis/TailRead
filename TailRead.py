import glob
import csv
import os
import re

##read all sequence files in folder
path = "c:/bio/data/"
os.listdir(path)
files = glob.glob(path + '*.seq')

##set standard sequence
standfile = "c:/bio/standard.txt"
standdata = open(standfile, 'r')
standseq = standdata.read()
standtxtseq = standseq.replace("\n","")

##prepare output list of results
finalresults = [['name','sequence', '3 end', 'tail sequence', 'tail length']]

##for each single sequence
for singlefile in files:
	singlefile = singlefile.replace('\\','/')
	singledata = open(singlefile, 'r')
	seq=singledata.read()

##reverse complement if needed

##trim to designated sequence, around the 3' end
	anchorseq = standtxtseq[:19]
	txtseq = seq.replace("\n","")
	trimmedseq = txtseq[(txtseq.find(anchorseq)):]

##finding the templated 3'end position
	x = 0
	match = 0
	mismatch = 0
	while x < len(trimmedseq) and x < len(standseq) and mismatch == 0:
		if trimmedseq[x] == standseq[x]:
			match +=1
		else:
			mismatch +=1
		x +=1
	if mismatch == 1:
		x -=1

##assess tail details
	tossout = 0
	errmsg = ""
	matchafter = 0
	mismatch =1
	inmatch = 0
	delmatch = 0
	taillen = 0
	y = 1
	z = 1
	w = 1
	
	while y < (len(trimmedseq)-x) and y < (len(standseq)-x) and tossout == 0:
		if trimmedseq[x+y] == standseq[x+y]:
			matchafter +=1
		y +=1
		if y > 3 and float((matchafter)/(y-1)) > 0.75:
			tossout = 1
			errmsg = 'mismatch'
		if y == 3 and matchafter == 2:
			tossout = 1
			errmsg = 'mismatch'
	while z < (len(trimmedseq)-x-1) and z < (len(standseq)-x) and tossout ==0:
		if trimmedseq[x+z+1] == standseq[x+z]:
			inmatch +=1
		z +=1
		if z > 3 and (float(inmatch)/(z-1)) > 0.75:
			tossout = 1
			errmsg = 'insert'
	while w < (len(trimmedseq)-x) and w < (len(standseq)-x-1) and tossout ==0:
		if trimmedseq[x+w] == standseq[x+w+1]:
			delmatch +=1
		w +=1
		if w > 3 and (float(delmatch)/(w-1)) > 0.75:
			tossout = 1
			errmsg = 'deletion'
	if x == len(standseq):
		y = 0
	if tossout == 0 and y > x:
		tossout = 1
		errmsg = 'match no confidence'
	if tossout == 0:
		if len(trimmedseq) <= len(standseq):
			taillength = y
		else:
			taillength = len(trimmedseq) - x
		tailseq = trimmedseq[x:]
	else:
		taillength = 'error'
		tailseq = errmsg

	finalresults.append([singlefile,trimmedseq,x,tailseq,taillength])

outputfile = path + 'output.csv'
outputwriter = open(outputfile, 'w')
csv.writer(outputwriter).writerows(finalresults)
outputwriter.close()
	
	