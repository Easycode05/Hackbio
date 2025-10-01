# Clinical Case Study: Whole-Genome Sequencing

**Clinical Background**

A 6-year-old boy presents with chronic cough, recurrent lung infections, and poor weight gain despite adequate nutrition. 
His physician suspects cystic fibrosis (CF). A sweat chloride test has been ordered but yields borderline results (45 mmol/L). 
To confirm diagnosis and identify the causative genetic variant(s), the clinical team orders whole-genome sequencing (WGS).

The boy’s father has also been sequenced as part of the family study. The mother’s genomic data are not available.

**Task**

You are part of the bioinformatics team tasked with analyzing the WGS data to determine whether the child carries pathogenic variants in the CFTR gene.

**Objectives**
1) Perform quality control and alignment of raw WGS reads to the human reference genome.
2) Call variants in both father and child.
3) Focus analysis on the CFTR locus (chromosome 7q31.2).
4)Compare the child’s variants against the father’s to determine inheritance patterns.
5) Annotate and interpret variants to identify known pathogenic CFTR mutations (e.g., ΔF508, G542X, W1282X, etc.) using variant databases such as ClinVar.
6) Conclude whether the child is affected by CF based on genotype, and comment on whether the observed mutations are inherited or potentially de novo.

## Quality Control and Preprocessing

A directory "projects" was created using - mkdir projects
A sub-directory "raw_reads" was created and the reads (forward and reverse of the father and child is contained here)

**1) QC Report:** I used fastqc to run a QC report for the reads. The output file is to be stored in qc/
```
mkdir qc
Fastqc raw_reads/*.gz -o qc/
```

**2) Html reports(forward of father and child)**

child_1_fastqc.html
- 10/11 flagged green
- Warning for Per base GC content
- %GC is 44 (suggestive of an human genome)
- Total Sequence: 51,146,287 |  Sequence length: 100bp
- Average Quality Per Read(Phred score): 38
- % of Sequence remaning if deduplicated: 80.89
- No Overrepresented sequence

father_1_fastqc.html
- 10/11 flagged green
- Warning for Per base GC content
- %GC is 41 (suggestive of an human genome)
- Total Sequence: 65,914,442 |  Sequence length: 100bp
- Average Quality Per Read(Phred score): 38
- % of Sequence remaning if deduplicated: 94.01
- No Overrepresented sequence

**3) Trimming with Fastp**
   
A bash script "trim" was created to handle the trimming of the reads.
The bash script starts with creating a directory trim/ to hold the trimmed reads.

Tool(s) used
- fastp
```
#!/bin/bash

#making directory
mkdir -p trim

#Samples
Samples=(child father)

#access the raw reads files and run fastp, send the output to a trim directory
for s in "${Samples[@]}"; do
   echo ">>> Running fastp for $s..."
   fastp \
      -i "/raw_reads/${s}_1.fastq.gz" \
      -I "/raw_reads/${s}_2.fastq.gz" \
      -o "trim/${s}_1_trim.fastq.gz" \
      -O "trim/${s}_2_trim.fastq.gz" \
      -h "trim/${s}.fastp.html" \
      -j "trim/${s}.fastp.json"\
      --verbose
done
```
**4) Genome mapping**
This involves mapping the reference genome (hg38) to the child and father genome.
A bash script "align.sh" was used to automate the genome mapping process. 

Tools used
- Repair.sh
- bwa mem
- samtools 

```
#!/bin/bash

#Reference genome
REF="/home/a_adegite/a_adegite/projects/refdata/hg38/hg38.fasta"

#make directories
mkdir -p repaired
mkdir -p alignment

#Samples
Samples=(child father)

for s in "${Samples[@]}"; do
   echo ">>>Repairing paired reads for $s..."
   
   #Repair forward and reverse reads, output repaired pairs + singletons
   repair.sh in1="trim/${s}_1_trim.fastq.gz" in2="trim/${s}_2_trim.fastq.gz" \
             out1="repaired/${s}_1.repaired.fastq.gz" out2="repaired/${s}_2.repaired.fastq.gz" \
             outsingle="repaired/${s}_singletons.fastq.gz"
   
   echo ">>> Mapping $s to reference genome..."
   
   #Map repaired reads to reference, output BAM directly
   bwa mem -R "@RG\tID:${s}\tSM:${s}\tPL:ILLUMINA" $REF "repaired/${s}_1.repaired.fastq.gz" "repaired/${s}_1.repaired.fastq.gz" \
   | samtools view -b -o "alignment/${s}_sample.bam"
   
   echo ">>> Finished $s"

done
```

**5) Variant Calling**

*i) Sorting & Marking duplicates*
```
#!/bin/bash

# Directories
mkdir -p sorted marked

# Samples
Samples=(child father)

# Loop through samples
for s in "${Samples[@]}"; do
(
    echo ">>> Processing $s"

    # 1. Sort BAM
    gatk SortSam \
        -I "alignment/${s}_sample.bam" \
        -O "sorted/${s}_sorted.bam" \
        -SORT_ORDER coordinate

    # 2. Mark duplicates
    gatk MarkDuplicates \
        -I "sorted/${s}_sorted.bam" \
        -O "marked/${s}_marked.bam" \
        -M "marked/${s}_metrics.txt"

    # 3. Build BAM index
    gatk BuildBamIndex -I "marked/${s}_marked.bam"

    echo ">>> Finished $s"
) &
done

# Wait for all background jobs to finish
wait
```

*ii) Calibration using Base Score Quality Recalibration*
```
#!/bin/bash

# Directories
mkdir -p BQSR

# Samples
Samples=(child father)

#Reference genome
REF="/home/a_adegite/a_adegite/projects/refdata/hg38/hg38.fasta"

# Known sites (both dbSNP + known indels)
DBSNP="/home/a_adegite/a_adegite/projects/refdata/hg38/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INDELS="/home/a_adegite/a_adegite/projects/refdata/hg38/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"

# Loop through samples
for s in "${Samples[@]}"; do
    echo ">>> Performing BQSR for $s..."

    #Base recalibration (generate recalibration table)
    gatk BaseRecalibrator \
        -I "marked/${s}_marked.bam" \
        -R $REF \
        --known-sites $DBSNP \
        --known-sites $INDELS \
        -O "BQSR/${s}_recal_data.table"

    #Apply recalibration
    gatk ApplyBQSR \
        -I "marked/${s}_marked.bam" \
        -R $REF \
        --bqsr-recal-file "BQSR/${s}_recal_data.table" \
        -O "BQSR/${s}_recal.bam"
done
```
*iii) Haplotype caller*
```
#!/bin/bash

#Directories
mkdir -p gvcfs

#Samples
Samples=(child father)

#Reference genome
REF="/home/a_adegite/a_adegite/projects/refdata/hg38/hg38.fasta"

#Loop through samples
for s in "${Samples[@]}"; do
    echo ">>> Running HaplotypeCaller for $s"

    gatk HaplotypeCaller \
        -I "BQSR/${s}_recal.bam" \
        -R $REF \
        -O "gvcfs/${s}.g.vcf.gz" \
        -ERC GVCF
done
```
*iv) Running combinegvcfs*
```
#!/bin/bash

#Directory
mkdir -p joint_genotyping

#Reference genome
REF="/home/a_adegite/a_adegite/projects/refdata/hg38/hg38.fasta"

gatk CombineGVCFs \
    -R $REF \
    -V "gvcfs/child.g.vcf.gz" \
    -V "gvcfs/father.g.vcf.gz" \
    -O "joint_genotyping/combined.g.vcf.gz"
```
*v) Genotype the combinedgvcfs*

This will generate a joint raw multi-sample VCF that will be hard-filtered
```
#!/bin/bash

#Directory
mkdir -p VCF

#Reference genome
REF="/home/a_adegite/a_adegite/projects/refdata/hg38/hg38.fasta"

gatk GenotypeGVCFs \
    -R $REF \
    -V "joint_genotyping/combined.g.vcf.gz" \
    -O "VCF/combined_raw.vcf.gz"
```

*vi) Filtering of variants*

This helps remove false positives and ensure that high quality variants are selected.
The sample size is small, hence, the hard filtering method will be used.
```
#!/bin/bash
set -euo pipefail

#Reference genome
REF="/home/a_adegite/a_adegite/projects/refdata/hg38/hg38.fasta"

#the combined vcf files from the joint genotyping process
RAW_VCF="VCF/combined_raw.vcf.gz"

#Directory
mkdir -p filtered

echo ">>> Extracting SNPs"
gatk SelectVariants \
  -R "$REF" \
  -V "$RAW_VCF" \
  --select-type-to-include SNP \
  -O "filtered/combined_raw_snps.vcf.gz"

echo ">>> Applying SNP filters (GATK-style: each filter separately)"
gatk VariantFiltration \
  -R "$REF" \
  -V "filtered/combined_raw_snps.vcf.gz" \
  -O "filtered/combined_filtered_snps.vcf.gz" \
  --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
  --filter-name "QUAL_lt_30" --filter-expression "QUAL < 30.0" \
  --filter-name "SOR_gt_3" --filter-expression "SOR > 3.0" \
  --filter-name "FS_gt_60" --filter-expression "FS > 60.0" \
  --filter-name "MQ_lt_40" --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum_lt_-12.5" --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum_lt_-8" --filter-expression "ReadPosRankSum < -8.0"

echo ">>> Extracting INDELs"
gatk SelectVariants \
  -R "$REF" \
  -V "$RAW_VCF" \
  --select-type-to-include INDEL \
  -O "filtered/combined_raw_indels.vcf.gz"

echo ">>> Applying INDEL filters (GATK-style)"
gatk VariantFiltration \
  -R "$REF" \
  -V "filtered/combined_raw_indels.vcf.gz" \
  -O "filtered/combined_filtered_indels.vcf.gz" \
  --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
  --filter-name "FS_gt_200" --filter-expression "FS > 200.0" \
  --filter-name "ReadPosRankSum_lt_-20" --filter-expression "ReadPosRankSum < -20.0"

echo ">>> Merging filtered SNPs and INDELs into final VCF"
gatk MergeVcfs \
  -I "filtered/combined_filtered_snps.vcf.gz" \
  -I "filtered/combined_filtered_indels.vcf.gz" \
  -O "filtered/combined_filtered.vcf.gz"

echo ">>> Indexing final VCF"
tabix -p vcf "filtered/combined_filtered.vcf.gz"

echo ">>> Done. Final filtered VCF: filtered/combined_filtered.vcf.gz"
```

*vi) Extracting chromosome 7 for analysis*
```
bcftools view -r chr7 "filtered/combined_filtered.vcf.gz -Oz -o chr7_filtered3.vcf.gz
```
*indexing*
```
bcftools index chr7_filtered3.vcf.gz
```

*vi) Variant Annotation*

All variants were annotated using the Ensembl Variant Effect Predictor (VEP).
The chr7_filtered3.vcf.gz" file was uploaded to ensembl for analysis
A VEP table was gotten as a result
The following filters were applied to VEP table 
- Symbol=CFTR
- Consequence=missense_variant
The resulting table was exported in as *vep_input.txt*

The individual genotype (father and child) was extracted from *chr7_filtered3.vcf.gz* and exported as *genotypes.txt*.

0/0 - homozygous reference
0/1 or 1/0 - heterozygous
1/1 - homozygous alternate

```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' chr7_filtered3.vcf.gz \
| sed '1s/.*/CHROM\tPOS\tREF\tALT\tFather_GT\tChild_GT/' > genotypes.txt
```

*vi) Data Merging*

The annotated VEP output "vep_input" was merged with the extracted father and child genotypes to create a unified table. 
This allowed direct comparison of genotypes for each variant while retaining the clinical annotations.

A python script was used to achieve this
```
import pandas as pd
import re

# Load VEP output

vep = pd.read_csv("vep_input.txt", sep="\t")

# Keep only CFTR variants
vep_cftr = vep[vep["SYMBOL"] == "CFTR"]

# Extract useful columns
vep_cftr = vep_cftr[[
    "Location",
    "Allele",
    "SYMBOL",
    "Existing_variation",
    "Consequence",
    "IMPACT",
    "CLIN_SIG",
    "HGVSp"  # Needed to calculate protein length
]]

# Extract start position from Location
# Format: chr:start-end, e.g., 7:117509089-117509089
vep_cftr['POS'] = vep_cftr['Location'].str.split(':').str[1].str.split('-').str[0].astype(int)

# Calculate protein length from HGVSp if Protein_length not present
# Example HGVSp: "p.Arg74Trp" or "p.Gly1234_Val1236del"
def calc_protein_length(hgvsp):
    if pd.isna(hgvsp):
        return 0
    # Extract numbers from the string
    nums = re.findall(r'\d+', hgvsp)
    if len(nums) == 0:
        return 0
    # Use the largest number as protein position
    return max(int(n) for n in nums)

vep_cftr['Protein_length'] = vep_cftr['HGVSp'].apply(calc_protein_length)

# Pick longest transcript per position

vep_sorted = vep_cftr.sort_values(['POS', 'Protein_length'], ascending=[True, False])
vep_best = vep_sorted.groupby('POS', as_index=False).first()

# Load genotype calls
geno = pd.read_csv("genotypes.txt", sep="\t")
geno["POS"] = geno["POS"].astype(int)

# Merge on position
merged = pd.merge(vep_best, geno, on="POS", how="inner")

# Add inheritance labeling
def label_inheritance(row):
    if row["Child_GT"] == "0/0":
        return "Not carried by Child"
    elif row["Father_GT"] == "0/0" and row["Child_GT"] != "0/0":
        return "De novo (not in Father)"
    elif row["Father_GT"] != "0/0" and row["Child_GT"] != "0/0":
        return "Inherited from Father"
    else:
        return "Unknown"

merged["Inheritance"] = merged.apply(label_inheritance, axis=1)

# Save final inheritance table

merged.to_csv("cftr_inheritance_table_best_transcripts.csv", index=False)
```

Variants affecting multiple transcripts were resolved by slecting the transcript with the largest protein position derived from HGVSp. 
This ensured one representative consequence per genomic position for downstream analysis.

*The inheritance was labelled*
- Child 0/0: Not carried by Child
- Father 0/0 + Child non-zero: De novo (not in father)
- Father non-zero + Child non-zero: Inherited from Father

**The Unified Table**
<img width="1176" height="207" alt="chr7_unified" src="https://github.com/user-attachments/assets/4d7b0f19-9d15-40a8-a1f2-785134af5da4" />


**Inheritance Analysis**

For each variant, child genotype were compared fo the father genotypes to determine the inheritance pattern:
- Inherited from father: Child carries at least one alternate allele present in the father
- Not inherited from father (Not carried by Child): Child does not carry any alternate alleles present in father.
- De novo (not in Father): Child carries an alternate allele absent in father.

**Variant Interpretation**

*rs115545701*
- This variant has been associated with Cystic fibrosis but they have conflicting interpretation of its pathogenicity, hence it is likely benign. 
- The variant is present in the father but absent in the child - No inheritance seen in child
- ClinVar Classification, Review status - 1 star
- T=0.0005542 (776/1400248, GnomAD_exomes)
   
<img width="747" height="359" alt="rs115545701" src="https://github.com/user-attachments/assets/0ab70c2c-6f28-4167-8b45-958ecf3dc5c7" />

*rs11971167*
- Pathogenic for Congenital bilateral aplasia of vas deferens from CFTR mutation
- Conflicitng interpretations of pathogenicity for cystic fibrosis
- A=0.004169 (622/149212, GnomAD_genomes)
  
<img width="1179" height="343" alt="rs11971167" src="https://github.com/user-attachments/assets/7ac51330-1c2c-43fb-9bd8-77696c842610" />
<img width="1184" height="322" alt="rs119_2" src="https://github.com/user-attachments/assets/d2ad9cde-3faa-426b-9faf-314587069294" />



