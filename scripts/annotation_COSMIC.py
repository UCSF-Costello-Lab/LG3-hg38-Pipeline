#!/usr/bin/python

####
#This script annotates a filtered mutect output with COSMIC data
####
import sys


def usage():
	print "Usage:"
	print "annotation_COSMIC.py example.mutations cosmicFilename"
	sys.exit()


def parseMuts(mutsFilename):
	"""Parses the mutations file and returns a header and mutation list"""
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


def parseCOSMIC(cosmicFilename):
	"""Parses the COSMIC file and returns a header and mutation list"""
	header = []
	cosmic = []
	cos = open(cosmicFilename)
	for line in cos:
		line = line.rstrip("\n")
		line = line.split("\t")
		if (line[0] == "Gene name"):
			header = line
			continue
		cosmic.append(line)
	cos.close()
	return (header, cosmic)


def annotation(mutsFilename, cosmicFilename):
	"""This script annotates a mutations file"""
	(mutectHeader, mutectSNPs) = parseMuts(mutsFilename)
	(cosmicHeader, cosmic) = parseCOSMIC(cosmicFilename)

	#find column numbers in mutect
	gene = mutectHeader.index("#gene")
	contig = mutectHeader.index("contig")
	position = mutectHeader.index("position")
	accession = mutectHeader.index("accession")

	#find column numbers in cosmic
	geneNAME = cosmicHeader.index("Gene name")
	accessionNUMBER = cosmicHeader.index("Accession Number")
	hg38POS = cosmicHeader.index("Mutation genome position")
	primarySITE = cosmicHeader.index("Primary site")
	subtypeSITE = cosmicHeader.index("Site subtype")
	primaryHIST = cosmicHeader.index("Primary histology")
	subtypeHIST = cosmicHeader.index("Histology subtype")
	PMID = cosmicHeader.index("Pubmed_PMID")

	#annotate the mutect calls with COSMIC info
	mutectHeader.append("COSMIC_mutation_frequency")
	mutectHeader.append("COSMIC_mutation_within_3bp_frequency")
	mutectHeader.append("COSMIC_gene_frequency")
	mutectHeader.append("COSMIC_tumor_types")
	mutectHeader.append("COSMIC_PMIDs")
	for x in mutectSNPs:
		a, b, c, d, e, f = 0, 0, 0, [], [], 0
		for y in cosmic:
			if (len(y[hg38POS]) >= 5):
				snpCHR = "chr" + y[hg38POS].split(":")[0]
				snpPOS = y[hg38POS].split("-")[1]
			else:
				snpCHR = "NONE"
			matching = 0
			ids = x[gene].split(",")
			for id in ids:
				if (id == y[geneNAME]):
					matching = 1
			ids = x[accession].split(",")
			for id in ids:
				if (id == y[accessionNUMBER]):
					matching = 1
			if (x[contig] == snpCHR):
				if (abs(int(x[position]) - int(snpPOS)) <= 3):
					b += 1
					matching = 1
				if (x[position] == snpPOS):
					a += 1
					matching = 1
			if (matching == 1):
				c += 1
				types = (y[primarySITE] + ":" + y[subtypeSITE] + ":" +
				y[primaryHIST] + ":" + y[subtypeHIST])
				if types not in d:
					d.append(types)
				e.append(y[PMID])
		x.append(str(a))
		x.append(str(b))
		x.append(str(c))
		if (len(d) > 10):
			theNUMBER = str(len(d)) + " different types"
			d = [theNUMBER]
		x.append(";".join(d))
		if (len(e) > 10):
			theNUMBER = str(len(e)) + " PMIDs"
			e = [theNUMBER]
		x.append(",".join(e))

	#print the (more) annotated file
	print "\t".join(mutectHeader)
	for line in mutectSNPs:
		print "\t".join(line)
	return()


if (len(sys.argv) != 3):
	usage()
else:
	annotation(sys.argv[1], sys.argv[2])
