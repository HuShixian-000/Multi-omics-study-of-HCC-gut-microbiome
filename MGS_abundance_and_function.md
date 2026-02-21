## Health:kneaddata and metaphlan
```bash
#!/bin/bash
#SBATCH --ntasks=96
#SBATCH --job-name=kneaddata+metaphlan
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1
exec 3<>/tmp/fd1
rm -rf /tmp/fd1
for ((i=1;i<=16;i++))
do 
echo >&3
done 
for ((n=1;n<=54;n++))
do 
read -u3
id=$(sed -n "${n}p" /groups/g5840105/home/liuyuyao/Health.txt| awk '{ print $1 }')
file=$(grep "${id}$" /groups/g5840105/home/share/BIGMetaG/sampleName_clientId.txt | awk '{ print $1 }')
file_1=/groups/g5840105/home/share/BIGMetaG/cleandata/${file}_good_1.fq.gz
file_2=$(echo ${file_1} | sed 's/1.fq/2.fq/g')
patient_name=$(basename ${file_1} | cut -d "_" -f2)
echo -e "${patient_name}\t${id}\tHealth">>/groups/g5840105/home/share/BIGMetaG_result/kneaddata/health/health_sample.txt
{
kneaddata \
 --input1 $file_1 \
 --input2 $file_2 \
 --output-prefix $patient_name \
 -db /groups/g5840087/home/share/refGenome/metagenome_database/databases/humanDatabase/  \
 --serial \
 --cat-final-output \
 --remove-intermediate-output \
 -v -t 6 \
 --output /groups/g5840105/home/share/BIGMetaG_result/kneaddata/health \
 --trimmomatic /groups/g5840123/home/share/trim/Trimmomatic-0.39\
 --run-fastqc-start \
 --run-fastqc-end \
 --fastqc /home/share/bin/FastQC && \
metaphlan \
/groups/g5840105/home/share/BIGMetaG_result/kneaddata/health/${patient_name}.fastq \
--bowtie2out /groups/g5840105/home/share/BIGMetaG_result/metaphlan/health/${patient_name}.bowtie2.bz2 \
--nproc 6 \
--input_type fastq \
--bowtie2db /groups/g5840087/home/share/metaphlan_databases/ \
-t rel_ab_w_read_stats \
-s /groups/g5840105/home/share/BIGMetaG_result/metaphlan/health/${patient_name}.sam \
-o /groups/g5840105/home/share/BIGMetaG_result/metaphlan/health/${patient_name}.profiled_metagenome.txt 
echo >&3
}&
done 
wait
exec 3<&-
````

## merge result
```bash
#!/bin/bash
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
kneaddata_read_count_table --input /groups/g5840105/home/share/BIGMetaG_result/kneaddata/health/ --output /groups/g5840105/home/share/BIGMetaG_result/kneaddata/health/health_kneaddata_table.tsv
merge_metaphlan_tables.py /groups/g5840105/home/share/BIGMetaG_result/metaphlan/health/*.profiled_metagenome.txt >>/groups/g5840105/home/share/BIGMetaG_result/metaphlan/health/Health_metagenome.txt
````


## Health:humann
```bash
#!/bin/bash
#SBATCH --ntasks=96
#SBATCH --job-name=humann
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1
exec 3<>/tmp/fd1
rm -rf /tmp/fd1
for ((i=1;i<=16;i++))
do 
echo >&3
done 
for ((n=1;n<=54;n++))
do 
read -u3
file=$(sed -n "${n}p" /groups/g5840105/home/share/BIGMetaG_result/kneaddata/health/health_sample.txt | awk '{ print $1 }')
file_1=/groups/g5840105/home/share/BIGMetaG_result/kneaddata/health/${file}.fastq
{
humann \
--input ${file_1} \
--output /groups/g5840105/home/share/BIGMetaG_result/humann/health \
--taxonomic-profile /groups/g5840105/home/share/BIGMetaG_result/metaphlan/health/${file}.profiled_metagenome.txt \
--threads 6 \
--o-log /groups/g5840105/home/share/BIGMetaG_result/humann/health/${file}.log \
--protein-database /users/liuyuyao/ref/uniref/ \
--nucleotide-database /groups/g5840087/home/share/humann_databases/chocophlan_v31/
echo >&3
}&
done 
wait
exec 3<&-
````

## merge result
```bash
#!/bin/bash
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
cp /groups/g5840105/home/share/BIGMetaG_result/metaphlan/health/*.profiled_metagenome.txt /groups/g5840105/home/liuyuyao/profiled_metagenome
merge_metaphlan_tables.py /groups/g5840105/home/liuyuyao/profiled_metagenome/*.profiled_metagenome.txt >>/groups/g5840105/home/share/BIGMetaG_result/metaphlan/health_metagenome.txt
````

## HCC: kneaddata and metaphlan
```bash
#!/bin/bash
#SBATCH --ntasks=192
#SBATCH --job-name=kneaddata+metaphlan
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1
exec 3<>/tmp/fd1
rm -rf /tmp/fd1
for ((i=1;i<=24;i++))
do 
echo >&3
done 
for ((n=1;n<=388;n++))
do 
read -u3
id=$(sed -n "${n}p" /groups/g5840105/home/liuyuyao/HCC.txt| awk '{ print $1 }')
file=$(grep "${id}$" /groups/g5840105/home/share/BIGMetaG/sampleName_clientId.txt | awk '{ print $1 }')
file_1=/groups/g5840105/home/share/BIGMetaG/cleandata/${file}_good_1.fq.gz
file_2=$(echo ${file_1} | sed 's/1.fq/2.fq/g')
patient_name=$(basename ${file_1} | cut -d "_" -f2)
echo -e "${patient_name}\t${id}\tHCC">>/groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC/HCC_sample.txt
{
kneaddata \
 --input1 $file_1 \
 --input2 $file_2 \
 --output-prefix $patient_name \
 -db /groups/g5840087/home/share/refGenome/metagenome_database/databases/humanDatabase/  \
 --serial \
 --cat-final-output \
 --remove-intermediate-output \
 -v -t 8 \
 --output /groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC \
 --trimmomatic /groups/g5840123/home/share/trim/Trimmomatic-0.39\
 --run-fastqc-start \
 --run-fastqc-end \
 --fastqc /home/share/bin/FastQC && \
metaphlan \
/groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC/${patient_name}.fastq \
--bowtie2out /groups/g5840105/home/share/BIGMetaG_result/metaphlan/HCC/${patient_name}.bowtie2.bz2 \
--nproc 8 \
--input_type fastq \
--bowtie2db /groups/g5840087/home/share/metaphlan_databases/ \
-t rel_ab_w_read_stats \
-s /groups/g5840105/home/share/BIGMetaG_result/metaphlan/HCC/${patient_name}.sam \
-o /groups/g5840105/home/share/BIGMetaG_result/metaphlan/HCC/${patient_name}.profiled_metagenome.txt 
echo >&3
}&
done 
wait
exec 3<&-
````

## merge result
```bash
#!/bin/bash
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
kneaddata_read_count_table --input  /groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC/ --output  /groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC/HCC_kneaddata_table.tsv
merge_metaphlan_tables.py /groups/g5840105/home/share/BIGMetaG_result/metaphlan/HCC/*.profiled_metagenome.txt >>/groups/g5840105/home/share/BIGMetaG_result/metaphlan/HCC/HCC_profiled_metagenome.txt
````

## HCC:humann
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --job-name=HCC_humann
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1
exec 3<>/tmp/fd1
rm -rf /tmp/fd1
for ((i=1;i<=16;i++))
do 
echo >&3
done 
for ((n=1;n<=388;n++))
do 
read -u3
file=$(sed -n "${n}p" /groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC/HCC_sample.txt | awk '{ print $1 }')
file_1=/groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC/${file}.fastq
{
humann \
--input ${file_1} \
--output /groups/g5840105/home/share/BIGMetaG_result/humann/HCC \
--taxonomic-profile /groups/g5840105/home/share/BIGMetaG_result/metaphlan/HCC/${file}.profiled_metagenome.txt \
--threads 6 \
--o-log /groups/g5840105/home/share/BIGMetaG_result/humann/HCC/${file}.log \
--protein-database /users/liuyuyao/ref/uniref/ \
--nucleotide-database /groups/g5840087/home/share/humann_databases/chocophlan_v31/
echo >&3
}&
done 
wait
exec 3<&-
````

## merge result
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --job-name=merge
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
cp /groups/g5840105/home/share/BIGMetaG_result/humann/HCC/*_pathabundance.tsv  /groups/g5840105/home/share/BIGMetaG_result/humann/pathabundance
cp /groups/g5840105/home/share/BIGMetaG_result/humann/health/*_pathabundance.tsv  /groups/g5840105/home/share/BIGMetaG_result/humann/pathabundance
humann_join_tables --input /groups/g5840105/home/share/BIGMetaG_result/humann/pathabundance --output pathabundance_all.tsv
````

## build CARD database
```bash
#!/bin/bash
source /opt/app/anaconda3/bin/activate
conda activate qiime2-2023.2
python AdjustFastaHeadersForShortBRED.py  <protein_fasta_protein_homolog_model.fasta> CARD_2023.fasta

shortbred_identify.py \
--goi /input/CARD_2023.fasta \
--ref /ref/uniref90.fasta.gz \
--markers /input/CARD_2023.faa \
--tmp /input/tmp
```

## quantify ARGs
### Health(health_ARGs.sh)
```bash
#!/bin/bash
#SBATCH --job-name=ARGs_health
for log in /input/*.log
do 
id=$(basename ${log} .log)
file=/input/${id}.fastq
shortbred_quantify.py \
--markers /ref/CARD_2023.faa \
--wgs ${file} \
--results /output/${id}.txt \
--tmp /output/tmp \
--threads 20
done 
```
```bash
#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --job-name=ARGs_health
singularity exec --bind /groups/g5840105/home/share/BIGMetaG_result/kneaddata/health:/input --bind /groups/g5840105/home/share/VFARDB:/ref --bind /groups/g5840105/home/share/BIGMetaG_result/ARGs/health:/output /opt/app/sif/shortbred.sif /bin/bash /output/health_ARGs.sh
```

### HCC(HCC_ARGs.sh)
```bash
#!/bin/bash
#SBATCH --job-name=ARGs_HCC
for log in /input/*.log
do 
id=$(basename ${log} .log)
file=/input/${id}.fastq
shortbred_quantify.py \
--markers /ref/CARD_2023.faa \
--wgs ${file} \
--results /output/${id}.txt \
--tmp /output/tmp \
--threads 20
done 
```
```bash
#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --job-name=ARGs_HCC
#SBATCH --nodelist=node[1]
singularity exec --bind /groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC:/input --bind /groups/g5840105/home/share/VFARDB:/ref --bind /groups/g5840105/home/share/BIGMetaG_result/ARGs/HCC:/output /opt/app/sif/shortbred.sif /bin/bash /output/HCC_ARGs.sh
```

## build VF database(build_faa.sh)
```bash
#!/bin/bash
shortbred_identify.py \
--goi /input/VFDB_2022.fasta \
--ref /ref/uniref90.fasta \
--markers /input/VFDB_2022.faa \
--tmp /input/tmp
```
```bash
#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --job-name=VFDB
#SBATCH --nodelist=cpu1
singularity exec --bind /groups/g5840105/home/share/VFARDB/VFDB:/input --bind /users/liuyuyao/ref/uniref:/ref /opt/app/sif/shortbred.sif /bin/bash /input/build_faa.sh
```

## quantify VF
### HCC (VF_HCC.sh)
```bash
#!/bin/bash
#SBATCH --job-name=VF_HCC
for log in /input/*.log
do 
id=$(basename ${log} .log)
file=/input/${id}.fastq
shortbred_quantify.py \
--markers /ref/VFDB_2022.faa \
--wgs ${file} \
--results /output/${id}.txt \
--tmp /output/tmp \
--threads 30
done 
```
```bash
#!/bin/bash
#SBATCH --ntasks=30
#SBATCH --job-name=VF_HCC
singularity exec --bind /groups/g5840105/home/share/BIGMetaG_result/kneaddata/HCC:/input --bind /groups/g5840105/home/share/VFARDB/:/ref --bind /groups/g5840105/home/share/BIGMetaG_result/VF/HCC:/output /opt/app/sif/shortbred.sif /bin/bash /output/VF_HCC.sh
```

### Health (health_VF.sh)
```bash
#!/bin/bash
#SBATCH --job-name=VF_Health
for log in /input/*.log
do 
id=$(basename ${log} .log)
file=/input/${id}.fastq
shortbred_quantify.py \
--markers /ref/VFDB_2022.faa \
--wgs ${file} \
--results /output/${id}.txt \
--tmp /output/tmp \
--threads 20
done 
```
```bash
#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --job-name=VF_Health
singularity exec --bind /groups/g5840105/home/share/BIGMetaG_result/kneaddata/health:/input --bind /groups/g5840105/home/share/VFARDB/:/ref --bind /groups/g5840105/home/share/BIGMetaG_result/VF/health:/output /opt/app/sif/shortbred.sif /bin/bash /output/health_VF.sh
```