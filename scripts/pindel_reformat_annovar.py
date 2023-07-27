#!/usr/bin/python
import sys

def pindel_reformat_annovar(fileheader, headerfile):
  ## open 2 files - for data & header
  data = open(fileheader + '.exonic_variant_function').readlines()
  splicing = open(fileheader + '.variant_function').readlines()
  header = open(headerfile).readline()

  ## adjust header line to account for annovar added columns
  oldh = header.strip().split('\t')
  if oldh[0] != '#CHROM':
    print 'error: no header'
    sys.exit(1)
  h = ['line', 'type', 'coding_details'] + oldh

  ## identify which columns contain the data that we want in the final filtered output file
  indel_type = h.index('type')
  coding = h.index('coding_details')
  chr = h.index('#CHROM'); h[chr] = "contig"
  start = h.index('POS'); h[start] = "position"
  end = h.index('END')
  ref = h.index('REF'); h[ref] = "ref_allele"
  alt = h.index('ALT'); h[alt] = "alt_allele"
  homlen = h.index('HOMLEN')
  homseq = h.index('HOMSEQ')
  svtype = h.index('SVTYPE')
  svlen = h.index('SVLEN')
  ntlen = h.index('NTLEN')
  #norm = h.index('normal_name')
  #tum = h.index('tumor_name')
  norm_alt = h.index('normal_AD')
  tum_alt = h.index('tumor_AD')
  pat = h.index('patient_ID')
  samp = h.index('sample_type')
  h.append('1000g2010_11'); g2010 = h.index('1000g2010_11')
  h.append('1000g2011_05'); g2011 = h.index('1000g2011_05')
  h.append('dbSNP'); dbSNP = h.index('dbSNP')
  h.append('#gene'); gene = h.index('#gene')
  h.append('accession'); transcript = h.index('accession')
  h.append('exon'); exon = h.index('exon')
  h.append('nucleotide'); nucleotide = h.index('nucleotide')
  h.append('protein'); protein = h.index('protein')
  h.append('context'); context = h.index('context')
  h.append('algorithm'); alg = h.index('algorithm')
  h.append('known_variant_status'); known = h.index('known_variant_status')
 
  ## add in splicing data from other file
  ## work through splicing data; make list of "exonic;splicing" genes
  es_genes = []; sp_list = []
  for line in splicing:
    l = line.strip('\n').split('\t')
    if l[0] == 'exonic;splicing':
      es_genes.append(l[1].split(';')[0])
    if l[0] == 'splicing':
      l = ['NA', 'NA', 'NA'] + l[2:] + [l[1], 'NA', 'NA', 'NA', 'NA', l[0], 'Pindel', ','.join(kvs(l, g2010-1, g2011-1, dbSNP-1))]
      sp_list.append(l)

  ## read data, parse coding_details into several new columns
  parsed = []
  #data = data[:20]
  for line in data:
    l = line.strip('\n').split('\t')

    g=[]; t=[]; e=[]; n=[]; p=[]
    if l[coding] == "UNKNOWN":
      g.append("UNKNOWN")
      t.append("NA")
      e.append("NA")
      n.append("NA")
      p.append("NA")
    else:
      splitme = l[coding].strip().split(',')[:-1]
      for s in splitme:
        x = s.split(':')
        g.append(x[0])
        t.append(x[1])
        e.append(x[2])
	if len(x) >= 4: n.append(x[3])
	else: n.append('NA')
        if len(x) == 5: p.append(x[4])
        else: p.append('NA')

    l = l + [g[0], ','.join(t), ','.join(e), ','.join(n), ','.join(p)]
    if g[0] in es_genes: l.append('exonic;splicing')
    else: l.append('exonic')
    l.append('Pindel')
    l.append(','.join(kvs(l, g2010, g2011, dbSNP)))
    parsed.append(l)

  ## reformat the data for final mutations file - combine exonic & splicing lists
  new_col_order = [gene, chr, start, ref, alt, nucleotide, protein, context, indel_type, pat, samp, norm_alt, tum_alt, svlen, ntlen, homlen, homseq, alg, known, transcript]
  ordered = []
  for l in parsed:
    ordered.append([l[x] for x in new_col_order])
  for l in sp_list:
    ordered.append([l[x] for x in new_col_order])
 
  ## output to file
  outfile = open(fileheader+'.muts', 'w')
  outfile.write('\t'.join([h[x] for x in new_col_order]) + '\n')
  for line in ordered:
    outfile.write('\t'.join(line) + '\n')
  outfile.close()

def kvs(line, g2010, g2011, dbSNP):
  toreturn = []
  if line[g2010] != '': toreturn.append(line[g2010])
  if line[g2011] != '': toreturn.append(line[g2011])
  if line[dbSNP] != '': toreturn.append(line[dbSNP])
  if toreturn == []: toreturn = ['novel']
  return toreturn

if __name__=="__main__":
  if len(sys.argv) != 3:
    print 'usage: %s annovar.file.header file.with.header.row' %(sys.argv[0])
    sys.exit(1)

  pindel_reformat_annovar(sys.argv[1], sys.argv[2].strip())

