import sys, subprocess, os.path, os

def pindel_setup(patient_ID, projectname, patIDs):
  ## read ID conversion file
  data = open(patIDs).readlines()
  header = data[0]
  data = data[1:]
  ## parse header
  h = header.strip().split('\t')
  col_pat = h.index("patient_ID")
  col_lib = h.index("lib_ID")
  col_st = h.index("sample_type")
  #col_file = h.index("file_header")

  ## pull out patient specific info
  data_pat = filter(lambda x:x.strip().split('\t')[col_pat] == patient_ID, data)

  ### After Recal use this:
  fileheader1 = ".bwa.realigned.rmDups.recal"
  bamheader1 = ".bwa.realigned.rmDups.recal"
  fileheader2 = "_multi"
  bamheader2 = ".bwa.mrkDups.sort.recal"
  #fileheader = "-trim.bwa.realigned.rmDups.recal"

  testID = data_pat[0].strip().split('\t')[col_lib]
  if "norm" in patient_ID:
    patient_ID_folder = patient_ID.split("norm")[0]
  else:
    patient_ID_folder = patient_ID

  fullpath = os.environ["LG3_INPUT_ROOT"] + "/" + projectname + "/bwa-mem/" + patient_ID_folder + "/"

  if os.path.isfile(fullpath + testID + fileheader1 + ".insert_size_metrics"):
    fileheader = fileheader1
    bamheader = bamheader1
  elif os.path.isfile(fullpath + testID + fileheader2 + ".insert_size_metrics"):
    fileheader = fileheader2
    bamheader = bamheader2
  else:
    print "ERROR: insert_size_metrics file can not be found. Tried:"
    print fullpath + testID + fileheader1 + ".insert_size_metrics"
    print fullpath + testID + fileheader2 + ".insert_size_metrics"
    sys.exit(1)

  print "Final fileheader " + fileheader 

  ## prepare outfile
  cfg = open(patient_ID + '.pindel.cfg', 'w')

  ## for each line in setup, open *insertsizemetrics, save median_insert_size, write to outfile
  for line in data_pat:
    l = line.strip().split('\t')
    sample_type = l[col_st]
    lib_ID = l[col_lib]

    ism = open(fullpath + lib_ID + fileheader + ".insert_size_metrics").readlines()
    foundit = False
    for i in ism:
      if foundit == True: med = i.split('\t')[0]; break
      if i.split('\t')[0] == "MEDIAN_INSERT_SIZE": foundit = True
    print sample_type, lib_ID, med

    if fileheader == ".merged":
      cfg.write(fullpath + lib_ID + fileheader + ".sorted.bam" + "\t" + med + "\t" + patient_ID + "_" + sample_type + "\n")
    else:
      cfg.write(fullpath + lib_ID + bamheader + ".bam" + "\t" + med + "\t" + patient_ID + "_" + sample_type + "\n")
      #cfg.write(fullpath + lib_ID + fileheader + ".bam" + "\t" + med + "\t" + patient_ID + "_" + sample_type + "\n")

  cfg.close()
    


if __name__=="__main__":
  if len(sys.argv) != 4:
    print 'usage: %s patient_ID projectname patient_ID_conversions' %(sys.argv[0])
    sys.exit(1)

  pindel_setup(sys.argv[1], sys.argv[2], sys.argv[3].strip())

