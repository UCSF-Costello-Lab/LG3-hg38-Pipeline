#!/usr/bin/python
import sys, subprocess, re

def vcf2tdt(filename):
  print filename
  data = open(filename).readlines()

  ## pull out header line to determine number of columns/samples in the vcf
  info_names = []; format_names = []
  i = 0
  for line in data:
    if line.startswith('##'):
      ## check if this line defines an INFO field; add to list of INFO names
      if line.startswith('##INFO='):
        info_names.append(re.sub('##INFO=<ID=(.*?),.*', "\\1", line).strip())
      elif line.startswith('##FORMAT='):
        format_names.append(re.sub('##FORMAT=<ID=(.*?),.*', "\\1", line).strip())
      i += 1
    elif line[0] == '#':
      header = line
      i += 1
      break
    else: print 'error'
  data = data[i:]

  ## parse header line
  h = header.strip().split('\t')
  chr = h.index('#CHROM')    # 0
  pos = h.index('POS')       # 1
  ref = h.index('REF')       # 3
  alt = h.index('ALT')       # 4
  info = h.index('INFO')  # 7
  format = h.index('FORMAT') # 8

  samples = {}
  for i in xrange(format+1, len(h)):
    if "Normal" in  h[i]: norm = i
    else: samples[h[i]] = i
  print norm
  print samples

  ## prep the output file
  outfile = open(filename+'.tdt', 'w')
  outheader = [h[x] for x in [chr, pos, ref, alt]] + info_names + ['normal_name', 'tumor_name'] + ["normal_"+ x for x in format_names] + ["tumor_"+x for x in format_names]
  #print outheader
  outfile.write('\t'.join(outheader) + '\n')

  ## now parse through data - one line of vcf will become len(samples) lines in the output file
  for line in data:
    l = line.strip().split('\t')
    if len(l) != len(h): print 'invalid number of columns'

    info_dict = dict(zip(info_names, ['' for i in info_names]))
    info_parsed = []
    for n in info_names:
      if n in l[info]:
        #info_dict[n] = re.sub('.*'+n+'=(.*?)', '\\1', l[info]).split(';')[0]
        info_parsed.append(re.sub('.*'+n+'=(.*?)', '\\1', l[info]).split(';')[0])
      else:
        #info_dict[n] = ''
        info_parsed.append('')
    #print l[info]
    #print info_dict
    #print info_parsed
   
    format_dict = {} 
    for n in format_names: ## eg  GT:0, AD:1
      format_dict[n] = l[format].split(':').index(n)
    #print format_dict

    sharedline = [l[x] for x in [chr, pos, ref, alt]]
    sharedline += info_parsed
    #sharedline += [info_dict[x] for x in info_names]

    for (s,i) in samples.iteritems():
      sampleline = sharedline + [h[norm], s]
      sampleline += [l[norm].split(':')[format_dict[x]] for x in format_names]
      sampleline += [l[i].split(':')[format_dict[x]] for x in format_names]
      outfile.write('\t'.join(sampleline) + '\n')
      
  outfile.close()


if __name__=="__main__":
  if len(sys.argv)<2:
    print 'usage: %s file'%(sys.argv[0])
    sys.exit(1)

  files = [sys.argv[1]]
  for f in files:
    vcf2tdt(f)
