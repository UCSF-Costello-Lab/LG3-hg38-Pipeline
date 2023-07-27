#!/usr/bin/python
import sys, optparse

parser=optparse.OptionParser(usage="%prog all.pindel.vcf.muts")
parser.add_option('-n', action='store', dest='normal', default=0, help='max number alt reads in normal [default: %default]')
parser.add_option('-t', action='store', dest='tumor', default=6, help='min number alt reads in tumor [default: %default]')
parser.add_option('-s', action='store', dest='size', default=50, help='max size of indel (in bp) [default: %default]')

def pindel_filter(filename, opts):
  ## open file
  data = open(filename).readlines()
  header = data[0]
  h = header.strip().split('\t')
  data = data[1:]

  ## make tumor dictionary for renaming samples
  #tumdict={"pri":"Primary", "priv2":"Primaryv2", "priv3":"Primaryv3","priv4":"Primaryv4", "priv5":"Primaryv5", "priv6":"Primaryv6", "rec1":"Recurrence1", "rec2":"Recurrence2", "rec3":"Recurrence3", "rec1.12":"Recurrence1v2", "rec1.13":"Recurrence1v3", "rec1.23":"Recurrence1v4", "rec1v2":"Recurrence1v2", "GBM":"Tumor", "GNS":"GNS", "ML":"ML"}

  # ['#CHROM', 'POS', 'REF', 'ALT', 'END', 'HOMLEN', 'HOMSEQ', 'SVLEN', 'SVTYPE', 'NTLEN', 'normal_name', 'tumor_name', 'normal_GT', 'normal_AD', 'tumor_GT', 'tumor_AD']

  ## identify which columns contain the data that we want in the final filtered output file
  chr = h.index('#CHROM')
  start = h.index('POS')
  end = h.index('END')
  ref = h.index('REF')
  alt = h.index('ALT')
  homlen = h.index('HOMLEN')
  homseq = h.index('HOMSEQ')
  svtype = h.index('SVTYPE')
  svlen = h.index('SVLEN')
  ntlen = h.index('NTLEN')
  norm = h.index('normal_name')
  tum = h.index('tumor_name')
  norm_alt = h.index('normal_AD')
  tum_alt = h.index('tumor_AD')
  h.append('patient_ID'); pat = h.index('patient_ID')
  h.append('sample_type'); samp = h.index('sample_type')

  ## read data, check each filter, build new table of maintained indels
  filterdata = []; pat_list = []
  for line in data:
    l = line.strip().split('\t')

    if l[svtype] not in ['INS', 'DEL', 'RPL']: continue

    ### BEGIN  edit by Ivan to handle newer pindel VCF with AD n_ref,n_alt
    #if int(l[norm_alt]) > opts.normal: continue
    #if int(l[tum_alt]) < opts.tumor: continue
    if ',' in l[norm_alt]:
        norm_altA = l[norm_alt].split(',')[1]
    else:
        norm_altA = l[norm_alt]
    if int(norm_altA) > opts.normal: continue

    if ',' in l[tum_alt]:
        tum_altA = l[tum_alt].split(',')[1]
    else:
        tum_altA = l[tum_alt]
    if int(tum_altA) < opts.tumor: continue
    #### END of Ivan's edit

    if abs(int(l[svlen])) > opts.size: continue
    if l[svtype]=='RPL' and  abs(int(l[svlen])) == abs(int(l[ntlen])):
      continue
    
    ## clean up the data that we keep
    ## fix ref/alt alleles and start positions
    if l[svtype]=='INS':
      l[ref] = '-'
      l[alt] = l[alt][1:]
    elif l[svtype]=='DEL':
      l[ref] = l[ref][1:]
      l[alt] = '-'
      l[start] = str(int(l[start])+1)
    elif l[svtype]=='RPL':
       del_len = abs(int(l[svlen]))
       nt_len = int(l[ntlen])
       if (len(l[ref]) == del_len + 1):
         if l[ref][:1] != l[alt][:1]: print 'WeIrD'
         l[ref] = l[ref][1:]
	 l[alt] = l[alt][1:]
         l[start] = str(int(l[start])+1)
       else:
         print 'this should never happen'
	 print l
    else: print 'invalid svtype: %s' %(l[svtype])
    ## add two columns: 'patient_ID' and 'sample_type'
    pp = l[tum].split("_")

    ## TROUBLESHOOTING: This where/why underscores in 'patient_ID' fails.
    ## The code tries to split up the 'tumor_name' value (e.g.
    ## 'Patient157_t10_underscore_Normal' and 'Patient157_t10_underscore_Primary')
    ## into 'patient_ID' and 'sample_type' based on the assumption that they
    ## are separated by an underscore ('_') and that underscore is the first
    ## one in the 'tumor_name' string.  The currently implementation causes
    ## 'patient_ID' to become 'Patient157' when we need 'Patient157_t10_underscore'.
    ## /HB 2018-10-16
    
    if (pp[0] not in pat_list): pat_list.append(pp[0])
    l.append(pp[0])
    l.append("_".join(pp[1:len(pp)]))
    filterdata.append('\t'.join(l))

    '''
    if l[tum][0] == "P":
      l.append('Patient' + l[tum][1:3]);
      if ('Patient' + l[tum][1:3] not in pat_list): pat_list.append('Patient' + l[tum][1:3])
      l.append(tumdict[l[tum][3:]])
      filterdata.append('\t'.join(l))

    elif l[tum][0]=="G" and l[tum][1]=="l":
      l.append("sam")
      if ("sam" not in pat_list): pat_list.append("sam")
      l.append(l[tum])
      filterdata.append('\t'.join(l))

    elif l[tum][0]=="G":
      l.append('GBM' + l[tum][3:5]);
      if ('GBM' + l[tum][3:5] not in pat_list): pat_list.append('GBM' + l[tum][3:5])
      l.append(tumdict[l[tum][6:]])
      filterdata.append('\t'.join(l))

    else:
      print "unknown sample format"
      sys.exit(1)
   '''
    
  ## check for multiple indels at the same position, select the one with higher coverage across all samples for that patient
  cutdata = []
  for p in pat_list:
    pat_data = filter(lambda x:x.split('\t')[pat] == p, filterdata)

    ## make a dictionary of all mutations - 
    mut_dict = {}  ## chr_pos : [ref_alt, #supporting_reads]
    for line in pat_data:
      l = line.split('\t')
      key = l[chr] + '_' + l[start]
      val1 = l[ref] + '_' + l[alt]

      ## find all mutations that match key & val1 - ie, all samples with the same mutation called
      mut_data = filter(lambda x:x.split('\t')[chr] == l[chr] and x.split('\t')[start] == l[start] and x.split('\t')[ref] == l[start] and x.split('\t')[alt] == l[alt], pat_data)
      ## get alt_coverage for this site in all samples
      total_coverage = 0
      for m in mut_data: total_coverage += int(m.split('\t')[tum_alt])

      ## add data to dictionary
      ## check if key is already present in dict
      if key in mut_dict:
        ## compare val1:
	if val1 == mut_dict[key][0]:
	  ## this mutation has already been addresssed - so just sanity check
	  if mut_dict[key][1] != total_coverage: print 'ERROR: coverage mismatch'
	else:
	  ## compare coverages - this is two mutations at the same position but with different ref/alt
	  ## if this new mutation has higher coverage, replace the ref/alt
	  if total_coverage > mut_dict[key][1]: mut_dict[key] = [val1, total_coverage]
      ## key is not in mut_dict; add it.
      else:
        mut_dict[key] = [val1, total_coverage]

    ## now we have mut_dict - a reference of which mutations to include. so now let's go through the patient mut_data and keep only the mutations in the dictionary
    for line in pat_data:
      l = line.split('\t')
      key = l[chr] + '_' + l[start]
      val1 = l[ref] + '_' + l[alt]
      if key in mut_dict and mut_dict[key][0] == val1: cutdata.append(line)
      
  '''
  cutdata = []
  for line in filterdata:
    out = True
    l = line.split('\t')
    for i in xrange(len(cutdata)):
      s = cutdata[i].split('\t')
      if l[start] == s[start] and l[chr] == s[chr] and l[norm] == s[norm] and l[tum] == s[tum]:
        out = False
        if int(l[tum_alt]) >= int(s[tum_alt]):
          cutdata[i] = line
    if out == True: cutdata.append(line)
  '''

  ## reformat the data - we only want some columns and we need them in a format that AnnoVar will accept
  # ['#CHROM', 'POS', 'REF', 'ALT', 'END', 'HOMLEN', 'HOMSEQ', 'SVLEN', 'SVTYPE', 'NTLEN', 'normal_name', 'tumor_name', 'normal_GT', 'normal_AD', 'tumor_GT', 'tumor_AD']
  new_col_order = [chr, start, end, ref, alt, pat, samp, norm_alt, tum_alt, svtype, svlen, ntlen, homlen, homseq]
  ordered = []
  for line in cutdata:
    l = line.split('\t')
    ordered.append([l[x] for x in new_col_order])

  ## output to file
  outfile = open(filename+'.filter', 'w')
  outfile.write('\t'.join([h[x] for x in new_col_order]) + '\tdummy\n')
  for line in ordered:
    outfile.write('\t'.join(line) + '\tdummy\n')
  outfile.close()


if __name__=="__main__":
  (opts,args) = parser.parse_args()
  if len(args) != 1:
    parser.print_help()
    sys.exit(1)

  pindel_filter(args[0], opts)

