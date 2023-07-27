#!/usr/bin/python
import sys, subprocess, os.path

def pindel_annotate_normal_coverage(filename, projectname):
  ## read data
  infile = open(filename)
  data = infile.readlines()
  infile.close()
  header = data[0]
  data = data[1:]

  ## parse header
  h = header.strip().split('\t')
  chr = h.index('contig')
  pos = h.index('position')
  ref = ['A', 'C', 'T', 'G']
  pat = h.index('patient_ID')
  h.append('normal_raw_coverage'); cov = h.index('normal_raw_coverage')

  ## prepare outfile
  annofile = open(filename + '.norm.txt', 'w')
  annofile.write('\t'.join(h) + '\n')

  ## make list of all patients
  pats = []
  for line in data:
    l = line.split('\t')
    if l[pat] not in pats: pats.append(l[pat])
  print pats

  ## for each patient, pull out all mutations, generate a pileup, and annotate
  for p in pats:
    print p
    muts = filter(lambda x:x.split('\t')[pat] == p, data)
    cfg_file = open('../' + p + '.pindel.cfg')
    norm_file = cfg_file.readline().split('\t')[0]
### Ivan
    cfg_file.close()
    ## confirm norm_file still exists in this location - if not, try /data
    no_bam = False
    if not os.path.isfile(norm_file):
      norm_file = norm_file.replace("home","data",1)
      if not os.path.isfile(norm_file):
        print "ERROR: normal bam file not found"
        no_bam = True

    tmp_file_header = filename + '_' + p 

    if not no_bam:
      ## generate bedfile for this patient
      print '  making bedfile'
      bedfile = tmp_file_header+ '.bed'
      outfile = open(bedfile, 'w')
      for m in muts:
        split = m.split('\t')
        outfile.write(split[chr] + '\t' + split[pos] + '\t' + str(int(split[pos])+1) + '\n')
      outfile.close()

      ## generate mpileup for this patient
      print '  making old pileup'
      mpilefile = tmp_file_header + '.pileup'
      #command = [os.environ["LG3_HOME"] + '/tools/samtools-0.1.12a/samtools', 'pileup', '-l', bedfile, norm_file]
      command = [os.environ["SAMTOOLS_OLD"], 'pileup', '-l', bedfile, norm_file]
      task=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      (stdout,stderr)=task.communicate()
      outfile = open(mpilefile, 'w')
      outfile.write(stdout)
      outfile.close()

    print '  annotating indels'
    for line in muts:
      l = line.strip('\n').split('\t')
      if no_bam: count="NA"
      else: count = getNormalCoverage(mpilefile, l[chr], l[pos], ref)
      l.append(str(count))
      annofile.write('\t'.join(l) + '\n')

    ## clean up
    print '  cleaning up'
    for f in [bedfile, mpilefile]:
      command = ['rm' , '-f', f]
      task=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      (stdout,stderr)=task.communicate()

  annofile.close()


def getNormalCoverage(rnaFilename, contig, position, ref_allele):
  """Parses a mpileup file and returns the normal coverage stats"""
  RNA_ref = 0
  rnaString = ""
  skip = 0
  rna = open(rnaFilename)
  for line in rna:
    line = line.rstrip("\n")
    line = line.split("\t")
    if (line[0] == contig):
      if (line[1] == position):
        rnaString = line[4].upper()
        break
  rna.close()
  #process the reads from the pileup
  for base in rnaString:
    if (skip > 0):
      skip -= 1
      continue
    if (base == "^"):
      skip = 1
      continue
    if (base.isdigit() == True):
      skip = int(base)
      continue
    if (base in ref_allele):
      RNA_ref += 1
  return (RNA_ref)



if __name__=="__main__":
  if len(sys.argv) != 3:
    print 'usage: %s inputfile projectname' %(sys.argv[0])
    sys.exit(1)

  pindel_annotate_normal_coverage(sys.argv[1], sys.argv[2].strip())

