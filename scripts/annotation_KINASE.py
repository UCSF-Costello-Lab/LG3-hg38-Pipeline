#!/usr/bin/python

####
#This script annotates a filtered mutect output with kinase information
####
import sys


def usage():
	print "Usage:"
	print "annotation_KINASE.py example.mutations kinaseList"
	sys.exit()


def parseMuts(mutsFilename):
	"""Parses the annotated.muts.txt file and returns a header and mutation list"""
	header = []
	snps = []
	muts = open(mutsFilename)
	for line in muts:
		line = line.rstrip("\n")
		line = line.split("\t")
		if (line[0] == "#gene"):
			header = line
			continue
		snps.append(line)
	muts.close()
	return (header, snps)


def annotation(mutsFilename, kinaseList):
	"""This script annotates a mutect file"""
	(header, snps) = parseMuts(mutsFilename)

	#find column numbers
	gene = header.index("#gene")
	accession = header.index("accession")

	#parse the kinase list
	kinases = []
	x = open(kinaseList)
	for line in x:
		line = line.rstrip("\n")
		line = line.split("\t")
		kinases.append(line)
	x.close()
	
	#annotate the mutect calls with kinase info
	header.append("KINASE?")
	for x in snps:
		answer = "NO"
		ids = x[accession].split(",")
		for id in ids:
			for y in kinases:
				if (id == y[1]):
					answer = "YES"
				id = id.split(".")[0]
				if (id == y[1]):
					answer = "YES"
		ids = x[gene].split(",")
		for id in ids:
			for y in kinases:
				if (id == y[0]):
					answer = "YES"
		x.append(answer)

	#print the (more) annotated file
	print "\t".join(header)
	for line in snps:
		print "\t".join(line)
	return()


if (len(sys.argv) != 3):
	usage()
else:
	annotation(sys.argv[1], sys.argv[2])