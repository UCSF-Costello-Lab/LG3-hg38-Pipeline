#!/usr/bin/python

####
#This script annotates a filtered mutect output with the Sanger Cancer Gene Census information
####
import sys


def usage():
	print "Usage:"
	print "annotation_CANCER.py example.mutations SangerCancerGeneCensus.txt RefseqToEntrez.txt"
	sys.exit()


def parseMuts(mutsFilename):
	"""Parses the example.mutations file and returns a header and mutation list"""
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

def parseCancer(cancerList):
	"""Parses the SangerCancerGeneCensus file and returns a geneID list"""
	genes = []
	cancer = open(cancerList)
	for line in cancer:
		line = line.rstrip("\n")
		line = line.split("\t")
		if (line[0] == "#Symbol"):
			continue
		genes.append(line)
	cancer.close()
	return (genes)

def parseRefSeq(RefseqToEntrez):
	"""Parses the RefseqToEntrez file"""
	convert = []
	x = open(RefseqToEntrez)
	for line in x:
		line = line.rstrip("\n")
		line = line.split("\t")
		if (line[0] == "#RefSeq_ID"):
			continue
		convert.append(line)
	x.close()
	return (convert)

def annotation(mutsFilename, cancerList, RefseqToEntrez):
	"""This script annotates a filtered mutect file"""
	(header, snps) = parseMuts(mutsFilename)
	convert = parseRefSeq(RefseqToEntrez)
	genes = parseCancer(cancerList)

	#find column numbers
	gene = header.index("#gene")
	accession = header.index("accession")
	
	#annotate the mutect calls with SangerCancerGeneCensus info
	header.append("SangerCancerGeneCensus?")
	for x in snps:
		answer = "NO"
		ids = x[accession].split(",")
		for id in ids:
			for y in convert:
				if (y[0] == id):
					for z in genes:
						if (y[1] == z[1]):
							answer = "YES"
		ids = x[gene].split(",")
		for id in ids:
			for z in genes:
				if (id == z[0]):
					answer = "YES"
		x.append(answer)

	#print the (more) annotated file
	print "\t".join(header)
	for line in snps:
		print "\t".join(line)
	return()


if (len(sys.argv) != 4):
	usage()
else:
	annotation(sys.argv[1], sys.argv[2], sys.argv[3])