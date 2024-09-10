### 1.liftover chain file
sourceFA=MS_chr5.fa.gz
targetFA=chr5.fa
nextflow run /disk222/miaoj/nf-LO --source ${sourceFA} --target ${targetFA} --distance near --aligner GSAlign --outdir GSAlign --max_memory 100.GB

### 2.mask the unliftable region
bedtools maskfasta -fi MS_chr5.fa -bed unliftable.bed -fo MS_chr5.mask.fa

### 3. simulate reads
/disk191/miaoj/software/mason2/bin/mason_simulator -ir MS_chr5.mask.fa -o MS5_REP1_R1.fq.gz -or MS5_REP1_R2.fq.gz -oa MS5_REP1.bam \
-n 3800000 --num-threads 10 --illumina-read-length 150 --seed 111

### 4. simulate reads containing known variants
# reference
samtools faidx Sscrofa11.fa chr5 > chr5.fa
# call variants
REF=chr5.fa
QUERY=MS_chr5.mask.fa
NAME=MS5
minimap2 -cx asm5 -t 32 --cs $REF $QUERY > ${NAME}.paf
sort -k6,6 -k8,8n ${NAME}.paf > ${NAME}.sort.paf
paftools.js call -l 10000 -L 50000 -q 5 -f $REF -s $QUERY ${NAME}.sort.paf > ${NAME}.vcf 2>${NAME}.stat &
bgzip MS5.vcf
tabix MS5.vcf.gz
# simulate
/disk191/miaoj/software/mason2/bin/mason_simulator -ir chr5.fa -iv MS5.vcf.gz -o REP1_R1.fq.gz -or REP1_R2.fq.gz -oa REP1.bam \
-n 3800000 --num-threads 10 --illumina-read-length 150 --seed 111 &