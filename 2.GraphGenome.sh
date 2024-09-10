### 1. construct the graph genome
# reference graph
vg autoindex --workflow giraffe -r chr5.fa -p REFg

# graph genome with known variants
vg autoindex --workflow giraffe -r chr5.fa -v variants.vcf.gz -p Var

# Pangenome constructed using MC pipeline
for i in BBH BM CH DNXE NCMD NX SWT TC LR MS
do
echo -e "${i}\t${i}.fa.gz" >> samples.txt
done

cactus-minigraph ./jobstore samples.txt Sus.gfa.gz --reference Sscrofa11 --mapCores 20 --workDir ./temp &> Minigraph.map.log
cactus-graphmap ./jobstore samples.txt Sus.gfa.gz Sus.paf --outputFasta Sus.gfa.fa --reference Sscrofa11 --mapCores 30 --workDir ./temp &> Minigraph.remap.log
cactus-align ./jobstore samples.txt Sus.paf Sus.hal \
--reference Sscrofa11 --pangenome --outVG --maxLen 10000 --workDir ./temp &> Cactus.align.log
cactus-graphmap-join ./jobstore --vg Sus.vg --hal Sus.hal \
--outDir MC/ --outName Sus-pan --reference Sscrofa11 --giraffe clip filter --filter 2 --vcf --gbz --gfa &> Cactus.index1.log


### 2. mapping performance evaluation
# index
gatk CreateSequenceDictionary -R chr5.fa -O chr5.dict
bwa index chr5.fa

# linear reference genome
SAMPLE_ID=simulation
bwa mem -M -t 20 -R \"@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tLB:${SAMPLE_ID}\tPL:illumina\" ${REF} ${fqDir}/MS5_${rep}_R1.fq.gz ${fqDir}/MS5_${rep}_R2.fq.gz | samtools view -bS - > ${SAMPLE_ID}.bam
samtools sort -@ 20 -o ${SAMPLE_ID}.sorted.bam ${SAMPLE_ID}.bam
samtools index ${SAMPLE_ID}.sorted.bam
gatk HaplotypeCaller --native-pair-hmm-threads 50 -R ${REF} -I ${SAMPLE_ID}.sorted.bam -O ${SAMPLE_ID}.vcf.gz
rm ${SAMPLE_ID}.bam

# graph genome
vg giraffe -Z ${indexDir}${prefix}.giraffe.gbz -m ${indexDir}${prefix}.min -d ${indexDir}${prefix}.dist -f ${fqDir}/MS5_${rep}_R1.fq.gz -f ${fqDir}/MS5_${rep}_R2.fq.gz -p -o GAM -t 20 > ${SAMPLE_ID}.gam
samtools sort -@ 20 -o ${SAMPLE_ID}.sorted.bam ${SAMPLE_ID}.bam
samtools addreplacerg -r \"@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tLB:${SAMPLE_ID}\tPL:illumina\" ${SAMPLE_ID}.sorted.bam -o ${SAMPLE_ID}.sorted.rg.bam
samtools index ${SAMPLE_ID}.sorted.rg.bam
rm ${SAMPLE_ID}.bam ${SAMPLE_ID}.sorted.bam
gatk HaplotypeCaller --native-pair-hmm-threads 50 -R ${REF} -I ${SAMPLE_ID}.sorted.rg.bam -O ${SAMPLE_ID}.vcf.gz
