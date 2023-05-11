#!/usr/bin/python

####
#This script parses a vcf file, compares a tumor-normal pair, and spits out stats
#on the relatedness of the two samples
####
import sys


def usage():
	print "Usage:"
	print "vcf_germline.py vcfFile normalName tumorName"
	sys.exit()


def parsing(fileName, normalName, tumorName):
	"""This script parses a vcf file"""
	header = []
	snps = []
	dbcount= 0
	tumorcount = 0
	filteredSNPs = []

	#parse the vcf file
	vcf = open(fileName)
	for line in vcf:
		line = line.rstrip("\n")
		line = line.split("\t")
		#save the column header list
		if (line[0][0] == "#"):
			if (line[0] == "#CHROM"):
				header = line
			continue
		snps.append(line)
	vcf.close()
	
	#find column numbers
	normal = header.index(normalName)
	tumor = header.index(tumorName)
	CHROM = header.index("#CHROM")
	POS = header.index("POS")
	ID = header.index("ID")
	INFO = header.index("INFO")
	REF = header.index("REF")
	ALT = header.index("ALT")
	FILTER = header.index("FILTER")
	FORMAT = header.index("FORMAT")
	
	for line in snps:
		#parse the format column
		snpformat = line[FORMAT].split(":")
		GT = snpformat.index("GT")
		AD = snpformat.index("AD")
		DP = snpformat.index("DP")
		GQ = snpformat.index("GQ")
		PL = snpformat.index("PL")
		#parse the data columns
		normalData = line[normal].split(":")
		tumorData = line[tumor].split(":")
		#***BEGIN FILTERS***
		if (line[normal] == "./."):
			continue
		if (line[tumor] == "./."):
			continue
		if (normalData[GT] != "0/1"):
			if (normalData[GT] != "1/1"):
				continue
		if (line[FILTER] != "PASS"):
			continue
		#***END FILTERS***
		filteredSNPs.append(line)
		if (line[ID][0:2] == "rs"):
			dbcount += 1
		if (tumorData[GT] == "0/1") or (tumorData[GT] == "1/1"):
			tumorcount += 1
		
	#calculate some stats
	dbcount = 100 * float(dbcount) / len(filteredSNPs)
	retained = 100 * float(tumorcount) / len(filteredSNPs)
	print " %s (Normal): %i non-reference variants (%.2f%% in dbSNP)" % (normalName,
	len(filteredSNPs), dbcount)
	print " %s (Tumor): %i retained (%.2f%%)" % (tumorName, tumorcount, retained)

	return()


if (len(sys.argv) != 4):
	usage()
else:
	parsing(sys.argv[1], sys.argv[2], sys.argv[3])

