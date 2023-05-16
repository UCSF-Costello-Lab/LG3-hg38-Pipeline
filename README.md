## LG3-hg38-Pipeline

## HG38 version of the original LG3 Pipeline.

## Installation

Clone repositary in designated directory

Required software and resources are listed in `lg3.conf`.
Edit as needed. Most software is available on C4 as loadable modules.
Currenly using `gatk/4.3.0.0`

## Preparations in workspace directory, e.g. ./runs_demo

1. Create sample table with 4 colums and the following headers:
"[lib_ID]	[SF]	[patient_ID]	[sample_type]"
e.g. `patient_ID_conversions.tsv`

NOTE: every patient must have a single sample_type == "Normal" sample!

2. Create run-time config file to set some env variables,
e.g. `run.conf`

```sh
	export LG3_HOME=/PATH/TO/ROOT/PIPELINE/DIRECTORY
	export PATH="${LG3_HOME}/bin:${PATH}"
	export CONV=patient_ID_conversions.tsv
	export PATIENT=PatientXXX
	export PROJECT=LG3_hg38
	export SAMPLES=
	export EMAIL=first.last@ucsf.edu
```

3. Create `./rawdata` directory (or symbolic link)
	This will be the path to the raw fastq files.
	Naming convention for FASTQ files:
		Lib_ID_R1.fastq.gz
		Lib_ID_R2.fastq.gz

4. Create `./output` directory (or symbolic link) 
	This will be the root of the output directory tree.

5. Create `./logs` directory (or symbolic link)
	This is where logs will be saved.


## Running pipeline from workspace directory step-by-step

1. Source the run-time config file
	e.g. 
	`source run.conf`

2. start `_run_Trim`

3. start `_run_Align_mem`

4. start `_run_Germline`

5. start `_run_Mutect2`

6. start `_run_Facets`

Note: 4, 5 and 6 could be run in parallel


## Pipeline Modules

1. `trim_galore.sh`

2. `Align_mem.sh` (all samples in parallel)
	- bwa mem
	- RevertSam (gatk)
	- AddOrReplaceReadGroups (gatk)
	- MergeBamAlignment (gatk)
	- MarkDuplicates (gatk)
	- flagstat (samtools)
	- SortSam | SetNmMdAndUqTags (gatk)
	- BaseRecalibrator (gatk)
	- ApplyBQSR (gatk)
	- ValidateSamFile (gatk)
	- CollectMultipleMetrics (gatk)
	- CollectHsMetrics (gatk)
	- flagstat (samtools)

3. Germline.sh (all samples together)
	- HaplotypeCaller (gatk)
	- CombineGVCFs (gatk)
	- GenotypeGVCFs (gatk)
	- VariantAnnotator (gatk)
	- VariantFiltration (gatk)
	- SelectVariants - SNPs (gatk)
	- check relatedness (custom)

4. Mutect2_TvsN.sh
	- Mutect2 (gatk)
	- LearnReadOrientationModel (gatk)
	- GetPileupSummaries: Normal (gatk)
	- GetPileupSummaries: Tumor (gatk)
	- CalculateContamination (gatk)
	- FilterMutectCalls (gatk)
	- Funcotator (gatk)
	- VariantsToTable (gatk)
	
5. FACETS
	- create_input_snp_pileup.py (custom)
	- snp_pileup (htstools)
	- runFACETS.R (facets R package)

