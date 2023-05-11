
##########################################################################################
# UCSF
# Costello Lab
# Formats the input for C++ snp-pileup code for all patients we have exome for (designed to
#   work for our unique file structure)
# Author: srhilz
# Version: v1 (2018.03.28)

# Input:
#   1. path to patient_ID_conversions.txt - must be located in runs_demo directory
#
# Output:
#   1. snp_pileup_input.txt
#
#
##########################################################################################

import sys, os

def create_input_snp_pileup(conversion_path, bampath):

    print("Reading in conversion file")

    ## read mutation data
    data = open(conversion_path).readlines()
    header = data[0]
    data = data[1:]

    ## parse header
    h = header.strip().split('\t')
    pid = h.index('patient_ID')
    lib = h.index('lib_ID')
    typ = h.index('sample_type')

    ## empty dic to fill with sample info, indexed by patient, then "Tumor" and "Normal"
    dic = {}

    ## empty dic to link libID with type, indexed by libid, used to provide an extra user-friendly identifier for tumor:normal pairs, which uses the tumor sample name
    dic_type = {}
    
    ## populate dic with patient and sample info from conversion file
    for line in data:
        line = line.rstrip().split('\t')
        if line[pid] not in dic:
            dic[line[pid]] = {}
        if line[typ] == 'Normal':
            if "Normal" not in dic[line[pid]]:
                dic[line[pid]]["Normal"] = [line[lib]]
            else:
                dic[line[pid]]["Normal"].append(line[lib])
        else:
            dic_type[line[lib]] = line[typ]
            if "Tumor" not in dic[line[pid]]:
                dic[line[pid]]["Tumor"] = [line[lib]]
            else:
                dic[line[pid]]["Tumor"].append(line[lib])

    ## empty dic to fill with sample info, indexed by patient, now formatted as
    ##   tumor-normal pairs
    dic_paired = {}

    ## populate dic_paired with tumor-normal pairs from dic, excluding pairs with
    ##  no normal, no tumor, or multiple normals (will have to figure these out
    ##  and run manually)
    print("Curating data for patients...")
    for patientID in dic:
        print(patientID)
        if "Normal" not in dic[patientID]:
            print("Warning(1): Patient does not have a normal sample...skipping.")
        else:
            if len(dic[patientID]["Normal"]) > 1:
                print("Warning(2): Patient has more than one normal sample...skipping.")
            else:
                normal_sample = dic[patientID]["Normal"][0]
                if "Tumor" not in dic[patientID]:
                    print("Warning(3): Patient does not have a tumor sample...skipping.")
                else:
                    dic_paired[patientID] = []
                    for tumor_sample in dic[patientID]["Tumor"]:
                        dic_paired[patientID].append([normal_sample, tumor_sample])

    ## for each patient, run snp-pileup
    print("Creating final SNP pileup input file...")

    recal_bam_ext = os.environ["RECAL_BAM_EXT"]

    outfile = open('snp_pileup_input.txt', 'w')# header = patient, nlibID, tlibID, nbampath, tbampath
    for patientID in dic_paired:
        print(patientID)
        for pair in dic_paired[patientID]:
            normal_libID = pair[0]
            tumor_libID = pair[1]
            # tries to get a normal bam
            print(bampath+'/'+patientID+'/'+normal_libID+'.'+recal_bam_ext+'.bam')
            if os.path.isfile(bampath+'/'+patientID+'/'+normal_libID+'.'+recal_bam_ext+'.bam'):
                normal_bam = bampath+'/'+patientID+'/'+normal_libID+'.'+recal_bam_ext+'.bam'
            elif os.path.isfile(bampath+'/'+patientID+'/'+normal_libID+'-trim.'+recal_bam_ext+'.bam'):
                normal_bam = bampath+'/'+patientID+'/'+normal_libID+'-trim.'+recal_bam_ext+'.bam'
            else:
                print("Warning(4): Normal bam path for patient cannot be determined...skipping.")
                continue
            # tries to get a tumor bam
            if os.path.isfile(bampath+'/'+patientID+'/'+tumor_libID+'.'+recal_bam_ext+'.bam'):
                tumor_bam = bampath+'/'+patientID+'/'+tumor_libID+'.'+recal_bam_ext+'.bam'
            elif os.path.isfile(bampath+'/'+patientID+'/'+tumor_libID+'-trim.'+recal_bam_ext+'.bam'):
                tumor_bam = bampath+'/'+patientID+'/'+tumor_libID+'-trim.'+recal_bam_ext+'.bam'
            else:
                print("Warning(5): Tumor bam path for patient cannot be determined...skipping.")
                continue
            sample_type = dic_type[tumor_libID]
            outfile.write('\t'.join([patientID, normal_libID, tumor_libID, normal_bam, tumor_bam, sample_type]) + '\n')


if __name__=="__main__":
    if len(sys.argv) != 3:
        print 'usage: %s patient_ID_conversion.txt bampath' %(sys.argv[0])
        sys.exit(1)

    create_input_snp_pileup(sys.argv[1], sys.argv[2].strip())
