
## 2病毒组原始表格获得

1.coverM：计算中位TPM、MSP_TPM、abundance、count
https://github.com/wwood/CoverM
1.1 数据准备
#根据cd-hit得到的聚类结果，提取vOTU序列，
cd ~/virus2
mkdir -p temp/coverm/cdhit
seqkit stat temp/cdhit/HQMQLQ_dereplicated.fna #8829
#拆分8829个contig为单独文件到指定目录
awk '/^>/{filename="temp/coverm/cdhit/"substr($0,2) ".fasta"}; {print >> filename; close(filename)}' temp/cdhit/HQMQLQ_dereplicated.fna
ls -l temp/coverm/cdhit/ | grep ^- | wc -l  #确定文件数 8829

#根据MGV得到的聚类结果，提取vOTU序列，
cd ~/virus2
mkdir -p temp/coverm/dereplicated_genomes result/coverm
seqkit stat result/vOTUs.fna # 7153
#拆分7153个contig为单独文件到指定目录
awk '/^>/{filename="temp/coverm/dereplicated_genomes/"substr($0,2) ".fasta"}; {print >> filename; close(filename)}' result/vOTUs.fna
ls -l temp/coverm/dereplicated_genomes/ | grep ^- | wc -l  #确定文件数 7153

#拆分genecat，40262个 3.11Gb为单独文件到指定目录
#cd ~/virus2
#seqkit stat ~/virus2/temp/cdhit/HQMQLQ_dereplicated.fna # 40262
#awk '/^>/{filename="temp/coverm/genecat_dereplicated_genomes/"substr($0,2) ".fasta"}; {print >> filename; close(filename)}' temp/cdhit/HQMQLQ_dereplicated.fna
#ls -l temp/coverm//genecat_dereplicated_genomes/ | grep ^- | wc -l  #确定文件数 40262
1.2 计算TPM 
######  1.计算中位TPM  ######
#并行计算
conda activate coverm
cd ~/virus2
mkdir -p temp/coverm/TPM
#需要用到原始数据
tail -n+2 result/metadata.txt|cut -f1|rush -j 6 \
"time coverm genome --coupled /data7/public/lyx/age/meta/temp/hr/{}_1.fastq /data7/public/lyx/age/meta/temp/hr/{}_2.fastq \
--genome-fasta-directory ~/virus2/temp/coverm/dereplicated_genomes/ -x fasta \
-m tpm \
--min-covered-fraction 0 \
-t 10 \
-o temp/coverm/TPM/{}.txt"
# 结果合并
#    conda activate humann3
    sed -i 's/_1.fastq TPM//' temp/coverm/TPM/*.txt
#    humann_join_tables --input temp/coverm/TPM \
#      --file_name txt \
#      --output result/coverm/TPM.txt   会报错
mkdir -p result/coverm/TPM     
for file in temp/coverm/TPM/*.txt; do  
    outfile="result/coverm/TPM/$(basename "$file")"  
    awk '{$1=""; sub(/^[ \t]+/, ""); print}' "$file" > "$outfile"
done

#合并
find result/coverm/TPM/ -type f -name "*.txt" -print0 | xargs -0 paste -d'\t' > temp/coverm/TPM.txt
#添加表头
awk 'NR==FNR {a[NR]=$1; next}   
    {print a[FNR] "\t" $0}  
' temp/coverm/TPM/A01.txt temp/coverm/TPM.txt > result/coverm/TPM.txt
rm -rf result/coverm/TPM
rm -rf temp/coverm/TPM.txt 

1.3 计算genecat_TPM 
######  1.3计算genecat_TPM  ######
#并行计算
conda activate coverm
cd ~/virus2
mkdir -p temp/coverm/genecat_TPM
#需要用到原始数据
tail -n+2 result/metadata.txt|cut -f1|rush -j 6 \
"time coverm genome --coupled /data7/public/lyx/age/meta/temp/hr/{}_1.fastq /data7/public/lyx/age/meta/temp/hr/{}_2.fastq \
--genome-fasta-directory ~/virus2/temp/coverm/genecat_dereplicated_genomes/ -x fasta \
-m tpm \
--min-covered-fraction 0 \
-t 10 \
-o temp/coverm/genecat_TPM/{}.txt"
# 结果合并
#    conda activate humann3
    sed -i 's/_1.fastq TPM//' temp/coverm/genecat_TPM/*.txt
#    humann_join_tables --input temp/coverm/TPM \
#      --file_name txt \
#      --output result/coverm/TPM.txt   会报错
mkdir -p result/coverm/genecat_TPM     
for file in temp/coverm/genecat_TPM/*.txt; do  
    outfile="result/coverm/genecat_TPM/$(basename "$file")"  
    awk '{$1=""; sub(/^[ \t]+/, ""); print}' "$file" > "$outfile"
done

#合并
find result/coverm/genecat_TPM/ -type f -name "*.txt" -print0 | xargs -0 paste -d'\t' > temp/coverm/genecat_TPM.txt
#添加表头
awk 'NR==FNR {a[NR]=$1; next}   
    {print a[FNR] "\t" $0}  
' temp/coverm/genecat_TPM/A01.txt temp/coverm/genecat_TPM.txt > result/coverm/genecat_TPM.txt
rm -rf result/coverm/genecat_TPM
rm -rf temp/coverm/genecat_TPM.txt 

1.4计算 MSP_TPM和Pangenomes
###如果宏基因组流程中用salmon基因定量计算出了gene.TPM，则可直接跳到1.3.4步骤运行MSPminer
1.3.1 metaProdigal基因预测Gene prediction
# 输入文件：去冗余得到的基因目录 temp/cdhit/HQMQLQ_dereplicated.fna
# 输出文件：prodigal预测的基因序列 temp/prodigal/gene.fa
# 基因大，可参考附录prodigal拆分基因文件，并行计算
conda activate megahit
cd ~/virus2
mkdir -p temp/prodigal
# prodigal的meta模式预测基因，>和2>&1记录分析过程至gene.log。1.8G1.5h
time prodigal -i temp/cdhit/HQMQLQ_dereplicated.fna \
    -d temp/prodigal/gene.fa \
    -o temp/prodigal/gene.gff \
    -p meta -f gff > temp/prodigal/gene.log 2>&1 
# 查看日志是否运行完成，有无错误
tail temp/prodigal/gene.log
# 统计基因数量,6G18s3M
seqkit stat temp/prodigal/gene.fa 
# 统计完整基因数量，数据量大可只用完整基因部分
grep -c 'partial=00' temp/prodigal/gene.fa 
# 提取完整基因(完整片段获得的基因全为完整，如成环的细菌基因组)
grep 'partial=00' temp/prodigal/gene.fa | cut -f1 -d ' '| sed 's/>//' > temp/prodigal/full_length.id
seqkit grep -f temp/prodigal/full_length.id temp/prodigal/gene.fa > temp/prodigal/full_length.fa
seqkit stat temp/prodigal/full_length.fa

1.3.2 cd-hit基因聚类/去冗余cluster & redundancy
# 输入文件：prodigal预测的基因序列 temp/prodigal/gene.fa
# 输出文件：去冗余后的基因和蛋白序列：result/NR/nucleotide.fa, result/NR/protein.fa

mkdir -p result/NR
# aS覆盖度，c相似度，G局部比对，g最优解，T多线程，M内存0不限制，n是word size，蛋白用5，核酸用10
# 2万基因2m，3M384p15m，2千万需要2000h，多线程可加速
cd-hit-est -i temp/prodigal/gene.fa \
    -o result/NR/nucleotide.fa \
    -aS 0.9 -c 0.95 -G 0 -g 0 -T 8 -M 0
# 统计非冗余基因数量，单次拼接结果数量下降不大，如3M-2M，多批拼接冗余度高
grep -c '>' result/NR/nucleotide.fa #5772506
# 翻译核酸为对应蛋白序列, --trim去除结尾的*
seqkit translate --trim result/NR/nucleotide.fa \
    > result/NR/protein.fa 
# 两批数据去冗余使用cd-hit-est-2d加速，见附录

1.3.3 salmon基因定量quantitfy
# 输入文件：去冗余后的基因序列：result/NR/nucleotide.fa
# 输出文件：Salmon定量：result/salmon/gene.count, gene.TPM

mkdir -p temp/salmon
salmon -v # 1.8.0

# 建索引, -t序列, -i 索引，10s,#num 5772506
salmon index -t result/NR/nucleotide.fa \
  -p 3 -i temp/salmon/index 

# 定量，l文库类型自动选择，p线程，--meta宏基因组
# 多个任务并行, 18s30m
time tail -n+2 result/metadata.txt | cut -f1 | rush -j 3 \
  "salmon quant -i temp/salmon/index -l A -p 6 --meta \
    -1 /data7/public/lyx/age/meta/temp/hr/{1}_1.fastq -2 /data7/public/lyx/age/meta/temp/hr/{1}_2.fastq \
    -o temp/salmon/{1}.quant"

# 合并
mkdir -p result/salmon
salmon quantmerge --quants temp/salmon/*.quant \
    -o result/salmon/gene.TPM
salmon quantmerge --quants temp/salmon/*.quant \
    --column NumReads -o result/salmon/gene.count
sed -i '1 s/.quant//g' result/salmon/gene.*

# 预览结果表格
head -n3 result/salmon/gene.*
1.3.4 MSPminer计算pangenomes
#需更换settings.ini文件内：计数矩阵文件路径、输出文件路径
count_matrix_file = /data4/machuang/virus2/result/salmon/gene.TPM
output_dir = ./output
#计算MSP中的pangenomes（core and accessory genes）
cd ~/virus2/temp/MSPminer 
mkdir -p output
#自行更换 settings.ini文件内的 计数矩阵文件路径
#count_matrix_file= /data4/machuang/virus2/result/salmon/gene.TPM
#output_dir=./output
./mspminer settings.ini   

##最终汇总结果在 all_msps.tsv 文件内
cp output/all_msps.tsv ~/virus2/result/salmon/

#core 164010 :指在所有被分析的基因组中都存在的基因集合。核心基因组代表了所有个体（或菌株）共有的功能。
#accessory 21243 :这部分基因只存在于部分个体（或菌株）中，而不是所有个体中。附属基因组代表了不同个体间的可变性和特异性，常常与环境适应性和菌株特异性功能有关。
#shared_core 8483 :在某个特定群体或子集的基因组中共同存在的基因集合。共享核心基因组表示某一特定群体内所有个体共有的功能。
#shared_accessory 428 :在某个特定群体或子集的基因组中部分存在的基因集合。共享附属基因组表示特定群体内的一些个体所特有的功能。
1.5 计算abundance
######2.abundance计算######
#并行计算相对丰度
conda activate coverm
cd ~/virus2
mkdir -p temp/coverm/abundance
tail -n+2 result/metadata.txt|cut -f1|rush -j 6 \
"time coverm genome --coupled /data7/public/lyx/age/meta/temp/hr/{}_1.fastq /data7/public/lyx/age/meta/temp/hr/{}_2.fastq \
--genome-fasta-directory temp/coverm/dereplicated_genomes/ -x fasta \
-o temp/coverm/abundance/{}.txt"
# 结果合并
#    conda activate humann3
    sed -i 's/_1.fastq Relative Abundance (%)//' temp/coverm/abundance/*.txt
#    humann_join_tables --input temp/coverm/abundance \
#     --file_name txt \
#     --output result/coverm/abundance.txt
      
mkdir -p result/coverm/abundance
for file in temp/coverm/abundance/*.txt; do  
    outfile="result/coverm/abundance/$(basename "$file")"  
    awk '{$1=""; sub(/^[ \t]+/, ""); print}' "$file" > "$outfile"
done
#合并
find result/coverm/abundance/ -type f -name "*.txt" -print0 | xargs -0 paste -d'\t' > temp/coverm/abundance.txt
#添加第一列
awk 'NR==FNR {a[NR]=$1; next}   
    {print a[FNR] "\t" $0}  
' temp/coverm/abundance/A01.txt temp/coverm/abundance.txt > result/coverm/abundance.txt
rm -rf result/coverm/abundance
rm -rf temp/coverm/abundance.txt    

1.6 计算count
######3.计算count######

#并行计算
conda activate coverm
cd ~/virus2
mkdir -p temp/coverm/count
tail -n+2 result/metadata.txt|cut -f1|rush -j 6 \
"time coverm genome --coupled /data7/public/lyx/age/meta/temp/hr/{}_1.fastq /data7/public/lyx/age/meta/temp/hr/{}_2.fastq \
--genome-fasta-directory ~/virus2/temp/coverm/dereplicated_genomes/ -x fasta \
-m count \
--min-covered-fraction 0 \
-t 10 \
-o temp/coverm/count/{}.txt"
# 结果合并
#    conda activate humann3
    sed -i 's/_1.fastq Read Count//' temp/coverm/count/*.txt
#    humann_join_tables --input temp/coverm/count \
#      --file_name txt \
#      --output result/coverm/count.txt
mkdir -p result/coverm/count
for file in temp/coverm/count/*.txt; do  
    outfile="result/coverm/count/$(basename "$file")"  
    awk '{$1=""; sub(/^[ \t]+/, ""); print}' "$file" > "$outfile"
done
#合并
find result/coverm/count/ -type f -name "*.txt" -print0 | xargs -0 paste -d'\t' > temp/coverm/count.txt
#添加第一列
awk 'NR==FNR {a[NR]=$1; next}   
    {print a[FNR] "\t" $0}    
' temp/coverm/count/A01.txt temp/coverm/count.txt > result/coverm/count.txt
rm -rf result/coverm/count
rm -rf temp/coverm/count.txt

#提取出iqtree用MGV做的的1600个votu
awk 'NR==FNR {ids[$1]; next} FNR==1 || $1 in ids' result/iqtree/vOTUs.concat.txt result/coverm/count.txt > result/coverm/count1600.txt

#提取出iqtree的1245个votu
awk 'NR==FNR {ids[$1]; next} FNR==1 || $1 in ids' result/iqtree2/vOTUs.concat.txt result/coverm/count.txt > result/coverm/count1245.txt
1.7 计算rpkm
######3.计算rpkm######

#并行计算
conda activate coverm
cd ~/virus2
mkdir -p temp/coverm/rpkm
tail -n+2 result/metadata.txt|cut -f1|rush -j 6 \
"time coverm genome --coupled /data7/public/lyx/age/meta/temp/hr/{}_1.fastq /data7/public/lyx/age/meta/temp/hr/{}_2.fastq \
--genome-fasta-directory ~/virus2/temp/coverm/dereplicated_genomes/ -x fasta \
-m rpkm \
--min-covered-fraction 0 \
-t 10 \
-o temp/coverm/rpkm/{}.txt"
# 结果合并
sed -i 's/_1.fastq RPKM//' temp/coverm/rpkm/*.txt
mkdir -p result/coverm/rpkm
for file in temp/coverm/rpkm/*.txt; do  
    outfile="result/coverm/rpkm/$(basename "$file")"  
    awk '{$1=""; sub(/^[ \t]+/, ""); print}' "$file" > "$outfile"
done
#合并
find result/coverm/rpkm/ -type f -name "*.txt" -print0 | xargs -0 paste -d'\t' > temp/coverm/rpkm.txt
#添加第一列
awk 'NR==FNR {a[NR]=$1; next}   
    {print a[FNR] "\t" $0}    
' temp/coverm/rpkm/A01.txt temp/coverm/rpkm.txt > result/coverm/rpkm.txt
rm -rf result/coverm/rpkm
rm -rf temp/coverm/rpkm.txt
2.与数据库blastn比对，找Novel (9h) 

2.1 对数据库进行CheckV质控（只保留HQ、MQ）
     https://note.youdao.com/s/MvKLetyA
2.2 合并MQ及以上质量的7个数据库，找Novel
        取交集很严格，因此找到的Novel很少
------------------------------------------
#合并CHVD、DEVoC、GPD、GVD、NCBI、CHGV、MGV七个数据库CheckV质控后的病毒序列
ref_name   原始     >MQ        Novel
#MGV      54118  → 52537     463/6495   
#DEVoC    12986  → 1850      4296/6495   
#CHVD     45033  → 27356     529/6495    
#GPD      142809 → 82436     367/6495   
#GVD      33242  → 8895      1026/6495  
#NCBI     18723  → 14077     4850/6495  
#CHGV     21646  → 11428     1728/6495   

#split_five_merge            2056/6495

#seven             198579    175/6495    

#合并
cat ~/db/checkv/CHVD/CHVD_HQMQ.fna ~/db/checkv/DEVoC/DEVoC_HQMQ.fna ~/db/checkv/GPD/GPD_HQMQ.fna ~/db/checkv/GVD/GVD_HQMQ.fna ~/db/checkv/NCBI/NCBI_HQMQ.fna ~/db/checkv/CHGV/CHGV_HQMQ.fna ~/db/checkv/MGV/MGV_HQMQ.fna > ~/db/merge_ref.fna
seqkit stat ~/db/merge_ref.fna  #198,579
#过滤重复contig，filter.py脚本有改动，请更新
~/virus1/temp/scripts/filter.py ~/db/merge_ref.fna ~/db/seven_ref.fna
seqkit stat ~/db/seven_ref.fna #198,579
rm ~/db/merge_ref.fna
-------------------------------------------#上述步骤只用运行一次。

cd ~/virus2/temp/MGV/mgv/ani_cluster
conda activate blast

#用/data4/machuang/db/MGV/MGV_v1.0_2021_07_08/mgv_contigs.fna数据建库，用~/virus1/result/checkv/HQMQLQ.fna进行all vs all比对
#合并后使用all vs all比对

makeblastdb -in ~/db/seven_ref.fna  -out blastdb -dbtype nucl
#输入6338
#~/db/software/ncbi-blast-2.16.0+/bin/blastn -query ~/virus2/result/vOTUs.fna  -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
blastn -query ~/virus2/result/vOTUs.fna  -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90

#Compute ANI from BLAST results
#python cluster.py -i blast.tsv -o ani1.tsv
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥10的序列，做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_merge_seven1.tsv
wc -l ref_merge_seven1.tsv #带表头5480
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_cluster/ref_merge_seven1.tsv result/vOTUs.list > temp/MGV/mgv/ani_cluster/ref_merge1.tsv
# 5479/6338 ref_merge
# 859/6338 novel_merge  13.45%
#文献Novel 1746/4422   39%
                               
2.旧方法找novel
##结果在ani.tsv文件,第一列我们的病毒，第二列MGV，第四列ani值
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $4="ani";print $1, $4; next} $4>=70 {print $1,$4}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2, "ref"}' OFS='\t' > ref_merge_seven.tsv
#wc -l ref_merge_seven.tsv # 带表头6138

cd ~/virus2
#得到表格，标明6338个病毒哪些是MGVref哪些是Novel
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_cluster/ref_merge_seven.tsv result/vOTUs.list > temp/MGV/mgv/ani_cluster/ref_merge.tsv

# 6137/6338  ref_merge
# 201/6338 novel_merge
#文献Novel 1746/4422   39%
2.3 单独与MGV比较
mkdir -p ~/virus2/temp/MGV/mgv/ani_mgv ~/virus2/temp/MGV/mgv/new_old
cd ~/virus2/temp/MGV/mgv/ani_mgv
cp ../ani_cluster/blastani.py ./
cp ../ani_cluster/cluster.py ./
conda activate blast

#用质控后的MGV数据库建库，~/db/checkv/MGV/MGV_HQMQ.fna #52537
makeblastdb -in ~/db/checkv/MGV/MGV_HQMQ.fna  -out blastdb -dbtype nucl
#输入6338
#~/db/software/ncbi-blast-2.16.0+/bin/blastn -query ~/virus2/result/vOTUs.fna  -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
blastn -query ~/virus2/result/vOTUs.fna -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥10的序列（因为我们的病毒length都是≥5000的，所以在该覆盖度下，我们命中的的序列长度是≥500，此标准参考陈卫华老师2024advance-science），做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_mgv1.tsv
wc -l ref_mgv1.tsv #带表头4617
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_mgv/ref_mgv1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_mgv1.tsv

2，NM文献方法，比对上的是已知，剩余是未知
cd ~/virus2/temp/MGV/mgv/ani_mgv
awk 'BEGIN {FS=OFS="\t"} NR == 1 {$1="contig_id";$2="ref";print $1,$2; next} !seen[$1]++ {print $1, "ref"}' ani.tsv > ref_mgv.tsv
#awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $4="ani"; $6="tcov";print $1, $4, $6; next} $4>=70 {print $1,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,"ref"}' OFS='\t' > ref_mgv3.tsv
wc -l ref_mgv.tsv #带表头5830

cd ~/virus2
#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_mgv/ref_mgv.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_mgv.tsv



# 5829/6338 ref_mgv
# 509/6338  novel_mgv
2.4 (可选)单独与DEVoC比较
mkdir -p ~/virus2/temp/MGV/mgv/ani_devoc
cd ~/virus2/temp/MGV/mgv/ani_devoc
cp ../ani_cluster/blastani.py ./
cp ../ani_cluster/cluster.py ./
conda activate blast

#用质控后的DEVoC数据库建库，~/db/checkv/DEVoC/DEVoC_HQMQ.fna #1850
makeblastdb -in ~/db/checkv/DEVoC/DEVoC_HQMQ.fna  -out blastdb -dbtype nucl
#输入6338
blastn -query ~/virus2/result/vOTUs.fna   -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000  -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥10的序列，做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_devoc1.tsv
wc -l ref_devoc1.tsv #带表头800
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_devoc/ref_devoc1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_devoc1.tsv

2，NM文献方法，比对上的是已知，剩余是未知
cd ~/virus2/temp/MGV/mgv/ani_devoc
awk 'BEGIN {FS=OFS="\t"} NR == 1 {$1="contig_id";$2="ref";print $1,$2; next} !seen[$1]++ {print $1, "ref"}' ani.tsv > ref_devoc.tsv
wc -l ref_devoc.tsv #带表头2142

cd ~/virus2
#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_devoc/ref_devoc.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_devoc.tsv

# 2141/6338 ref_devoc
# 4197/6338 novel_devoc
2.5 单独与CHVD比较
mkdir -p ~/virus2/temp/MGV/mgv/ani_chvd
cd ~/virus2/temp/MGV/mgv/ani_chvd
cp ../ani_cluster/blastani.py ./
cp ../ani_cluster/cluster.py ./
conda activate blast

#用质控后的CHVD数据库建库，~/db/checkv/CHVD/CHVD_HQMQ.fna #27356
makeblastdb -in ~/db/checkv/CHVD/CHVD_HQMQ.fna  -out blastdb -dbtype nucl
#输入6338
blastn -query ~/virus2/result/vOTUs.fna   -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000  -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥70的序列，做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_chvd1.tsv
wc -l ref_chvd1.tsv #带表头4509
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_chvd/ref_chvd1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_chvd1.tsv

2，NM文献方法，比对上的是已知，剩余是未知
cd ~/virus2/temp/MGV/mgv/ani_chvd
awk 'BEGIN {FS=OFS="\t"} NR == 1 {$1="contig_id";$2="ref";print $1,$2; next} !seen[$1]++ {print $1, "ref"}' ani.tsv > ref_chvd.tsv
wc -l ref_chvd.tsv #带表头5775

cd ~/virus2
#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_chvd/ref_chvd.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_chvd.tsv

# 5774/6338 ref_chvd
# 564/18898 novel_chvd
2.6 单独与GPD比较
mkdir -p ~/virus2/temp/MGV/mgv/ani_gpd
cd ~/virus2/temp/MGV/mgv/ani_gpd
cp ../ani_cluster/blastani.py ./
cp ../ani_cluster/cluster.py ./
conda activate blast

#用质控后的GPD数据库建库，~/db/checkv/GPD/GPD_HQMQ.fna #82436
makeblastdb -in ~/db/checkv/GPD/GPD_HQMQ.fna  -out blastdb -dbtype nucl
#输入6338
blastn -query ~/virus2/result/vOTUs.fna   -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000  -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥70的序列，做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_gpd1.tsv
wc -l ref_gpd1.tsv #带表头4896
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_gpd/ref_gpd1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_gpd1.tsv

2，NM文献方法，比对上的是已知，剩余是未知
cd ~/virus2/temp/MGV/mgv/ani_gpd
awk 'BEGIN {FS=OFS="\t"} NR == 1 {$1="contig_id";$2="ref";print $1,$2; next} !seen[$1]++ {print $1, "ref"}' ani.tsv > ref_gpd.tsv
wc -l ref_gpd.tsv #带表头5932

cd ~/virus2
#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_gpd/ref_gpd.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_gpd.tsv

# 5931/6338 ref_gpd
# 407/6495 novel_gpd
2.7 单独与GVD比较
mkdir -p ~/virus2/temp/MGV/mgv/ani_gvd
cd ~/virus2/temp/MGV/mgv/ani_gvd
cp ../ani_cluster/blastani.py ./
cp ../ani_cluster/cluster.py ./
conda activate blast

#用质控后的GVD数据库建库，~/db/checkv/GVD/GVD_HQMQ.fna #8895
makeblastdb -in ~/db/checkv/GVD/GVD_HQMQ.fna  -out blastdb -dbtype nucl
#输入6338
blastn -query ~/virus2/result/vOTUs.fna   -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000  -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥70的序列，做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_gvd1.tsv
wc -l ref_gvd1.tsv #带表头3665
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_gvd/ref_gvd1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_gvd1.tsv

2，NM文献方法，比对上的是已知，剩余是未知
cd ~/virus2/temp/MGV/mgv/ani_gvd
awk 'BEGIN {FS=OFS="\t"} NR == 1 {$1="contig_id";$2="ref";print $1,$2; next} !seen[$1]++ {print $1, "ref"}' ani.tsv > ref_gvd.tsv
wc -l ref_gvd.tsv #带表头5287

cd ~/virus2
#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_gvd/ref_gvd.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_gvd.tsv

# 5286/6338 ref_gvd
# 1052/18898 novel_gvd
2.8 （可选）单独与NCBI比较
mkdir -p ~/virus2/temp/MGV/mgv/ani_ncbi
cd ~/virus2/temp/MGV/mgv/ani_ncbi
cp ../ani_cluster/blastani.py ./
cp ../ani_cluster/cluster.py ./
conda activate blast

#用质控后的NCBI数据库建库，~/db/checkv/GVD/GVD_HQMQ.fna #14077
makeblastdb -in ~/db/checkv/NCBI/NCBI_HQMQ.fna  -out blastdb -dbtype nucl
#输入6338
blastn -query ~/virus2/result/vOTUs.fna   -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000  -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥70的序列，做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_ncbi1.tsv
wc -l ref_ncbi1.tsv #带表头487
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_ncbi/ref_ncbi1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_ncbi1.tsv

2，NM文献方法，比对上的是已知，剩余是未知
cd ~/virus2/temp/MGV/mgv/ani_ncbi
awk 'BEGIN {FS=OFS="\t"} NR == 1 {$1="contig_id";$2="ref";print $1,$2; next} !seen[$1]++ {print $1, "ref"}' ani.tsv > ref_ncbi.tsv
wc -l ref_ncbi.tsv #带表头1617

cd ~/virus2
#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_ncbi/ref_ncbi.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_ncbi.tsv

# 1616/6338 ref_ncbi
# 4722/6338 novel_ncbi
2.9 单独与CHGV比较
mkdir -p ~/virus2/temp/MGV/mgv/ani_chgv
cd ~/virus2/temp/MGV/mgv/ani_chgv
cp ../ani_cluster/blastani.py ./
cp ../ani_cluster/cluster.py ./
conda activate blast

#用质控后的CHGV数据库建库，~/db/checkv/CHGV/CHGV_HQMQ.fna  #11428
makeblastdb -in ~/db/checkv/CHGV/CHGV_HQMQ.fna  -out blastdb -dbtype nucl
#输入18898
blastn -query ~/virus2/result/vOTUs.fna   -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000  -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#筛选ani值≥70，覆盖度(tcov)≥10的序列，做为已知病毒
awk 'BEGIN {FS=OFS="\t"} NR==1 {$1="contig_id"; $2="tname"; $4="ani"; $6="tcov"; print $1, $2, $4, $6; next} $4>=70 && $6>=10 {print $1,$2,$4,$6}' ani.tsv | sort -k1,1 | awk '!seen[$1]++ {print $1,$2,$3,$4, "ref"}' OFS='\t' > ref_chgv1.tsv
wc -l ref_chgv1.tsv #带表头5308 
cd ~/virus2
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_chgv/ref_chgv1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_chgv1.tsv

2，NM文献方法，比对上的是已知，剩余是未知
awk 'BEGIN {FS=OFS="\t"} NR == 1 {$1="contig_id";$2="ref";print $1,$2; next} !seen[$1]++ {print $1, "ref"}' ani.tsv > ref_chgv.tsv
wc -l ref_chgv.tsv #带表头11615

cd ~/virus2
#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "ref" : "Novel");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/ani_chgv/ref_chgv.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_chgv.tsv

# 11614/18898 ref_chgv
# 7284/18898 novel_chgv
2.10 合并单独比较结果
MGV、CHVD、GPD、GVD 、CHGV  取并集 （DEVoC和NCB的Novel太多，可能是因为数据库不全，因此不合并）
#合并Novel的，并去重
cd  ~/virus2
#取交集     
awk '
    FNR == 1 && NR == FNR {next}
    FNR == 1 && NR != FNR {header = $0; next}
    FNR > 1 && $2 == "Novel" {novel[$1]++}
    END {
        for (id in novel) {
            if (novel[id] == ARGC - 1) {
                print id, "Novel"
            }
        }
    }
' temp/MGV/mgv/new_old/ref_mgv1.tsv temp/MGV/mgv/new_old/ref_chvd1.tsv temp/MGV/mgv/new_old/ref_gpd1.tsv temp/MGV/mgv/new_old/ref_gvd1.tsv temp/MGV/mgv/new_old/ref_chgv1.tsv > temp/MGV/mgv/new_old/ref_split_five1.tsv

wc -l temp/MGV/mgv/new_old/ref_split_five1.tsv # Novel带表头5228

#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "Novel" : "ref");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/new_old/ref_split_five1.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_five1.tsv

2.旧方法
#合并Novel的，并去重
cd  ~/virus2
#取并集
awk 'FNR==1 && NR!=1 {next} 
     FNR==1 {print; next} 
     $2 != "ref" && !seen[$1]++ {print}' temp/MGV/mgv/new_old/ref_mgv.tsv temp/MGV/mgv/new_old/ref_chvd.tsv temp/MGV/mgv/new_old/ref_gpd.tsv temp/MGV/mgv/new_old/ref_gvd.tsv temp/MGV/mgv/new_old/ref_chgv.tsv > temp/MGV/mgv/new_old/ref_split_five.tsv
wc -l temp/MGV/mgv/new_old/ref_split_five.tsv # Novel带表头8553

#标注新旧病毒
awk 'NR==FNR {
    unique[$1];
    next
}    
FNR==1 && !header_printed {
    print "contig_id\ttype";
    header_printed=1;
}    
{
    printf "%s\t%s", $1, ($1 in unique ? "Novel" : "ref");
    for (i=2; i<=NF; i++) printf "\t%s", $i;
    print "";
}' temp/MGV/mgv/new_old/ref_split_five.tsv result/vOTUs.list > temp/MGV/mgv/new_old/ref_five.tsv

cp temp/MGV/mgv/new_old/ref_five.tsv result/MGV/
#ref  10345/6495 
#Novel  8553/18898   45.26%
#文献Novel 1746/4422   39%

3.宿主预测 iPHoP
https://bitbucket.org/srouxjgi/iphop/src/main/
24/6/25 17:11  24/6/30 4:56
#运行iphop，记得删掉文件夹中旧的运行文件
conda activate iphop
cd ~/virus2
mkdir -p  temp/iphop result/iphop

iphop predict --fa_file ~/virus2/result/vOTUs.fna --db_dir /data4/machuang/db/viwrap/iPHoP_db/iPHoP_db/ --out_dir ~/virus2/temp/iphop/ -t 10
#保留最高Confidence score的宿主预测
awk -F, '!seen[$1]++' temp/iphop/Host_prediction_to_genus_m90.csv > result/iphop/Host_prediction_to_genus_m90.csv
 #该文件不包含种水平
cp temp/iphop/Host_prediction_to_genome_m90.csv result/iphop/ #该文件包含种水平，但要自己根据得分筛选
wc -l result/iphop/Host_prediction_to_genus_m90.csv # 13782/18898 具有宿主预测
#用不包含种水平的Host_prediction_to_genus_m90.csv文件提取出phylum水平
awk -F, '{print $1","$3}' temp/iphop/Host_prediction_to_genus_m90.csv | sed '1s/Virus/contig_id/; 1s/Host genus/phylum/' | awk -F, 'BEGIN {OFS=","} {sub(/\*p__[^,]*|;c[^,]*$/, "", $2); print}' | awk -F, 'BEGIN {OFS=","} {sub(/^.*p__/, "", $2); print}' | awk -F, '!seen[$1]++' > result/iphop/Host_prediction_phylum.csv

#awk -F, '{print $1","$3}' temp/iphop/Host_prediction_to_genus_m90.csv > temp/iphop/Host_prediction_columns1_3.csv
#sed -i '1s/Virus/contig_id/; 1s/Host genus/phylum/' temp/iphop/Host_prediction_columns1_3.csv
#awk -F, 'BEGIN {OFS=","} {sub(/\*p__[^,]*|;c[^,]*$/, "", $2); print}' temp/iphop/Host_prediction_columns1_3.csv > temp/iphop/Host_prediction_modified.csv
#awk -F, 'BEGIN {OFS=","} {sub(/^.*p__/, "", $2); print}' temp/iphop/Host_prediction_modified.csv > temp/iphop/Host_prediction_phylum.csv
#awk -F, '!seen[$1]++' temp/iphop/Host_prediction_phylum.csv > result/iphop/Host_prediction_phylum.csv
#rm temp/iphop/Host_prediction_columns1_3.csv
#rm temp/iphop/Host_prediction_modified.csv

#删除中间文件
cd ~/virus2
rm -rf temp/iphop/Wdir
4.物种注释
4.1geNomed 
https://bitbucket.org/srouxjgi/iphop/src/main/
大部分注释到纲水平
########## geNomed ########
#genomed的使用,annotate进行注释
conda activate genomad
cd ~/virus2
mkdir -p result/taxonomy/genomad temp/taxonomy/genomad
genomad annotate result/vOTUs.fna temp/taxonomy/genomad ~/db/genomad_db 

cp temp/taxonomy/genomad/vOTUs_annotate/vOTUs_taxonomy.tsv result/taxonomy/genomad/
#个性化运行
#genomad end-to-end metagenome.fna genomad_output ~/db/genomad_db
#genomad annotate metagenome.fna genomad_output ~/db/genomad_db
#genomad find-proviruses metagenome.fna genomad_output ~/db/genomad_db
#genomad marker-classification metagenome.fna genomad_output ~/db/genomad_db
#genomad nn-classification metagenome.fna genomad_output
#genomad aggregated-classification metagenome.fna genomad_output
#genomad score-calibration metagenome.fna genomad_output
#genomad summary metagenome.fna genomad_output
4.2  VirusTaxo (Genus level :7079/18898)
https://github.com/omics-lab/VirusTaxo
#属水平
cd ~/virus2/temp/virustaxo
conda activate virustaxo
source ./environment/bin/activate
python3 predict.py -h
mkdir -p ~/virus2/result/taxonomy/virustaxo ~/virus2/temp/taxonomy/virustaxo
#运行
 python3 predict.py \
   --model_path ~/db/virustaxo/vt_db_jan21_2024/DNA_RNA_18451_k20.pkl \
   --seq ~/virus2/result/vOTUs.fna  > ~/virus2/temp/taxonomy/virustaxo/virustaxo_genus.txt
#查看Genus多少未注释
awk '$3 == "Unknown" {count++} END {print count}' ~/virus2/temp/taxonomy/virustaxo/virustaxo_genus.txt
# 11819/18898  62.54% Unknown genus
cp ~/virus2/temp/taxonomy/virustaxo/virustaxo_genus.txt ~/virus2/result/taxonomy/virustaxo/
#退出./environment/bin/activate环境
deactivate
4.3  VirusTaxo_Hierarchical (order, family and genus level ) （测试）
https://github.com/omics-lab/VirusTaxo_Hierarchical
#order, family and genus
conda activate VH
source ./environment/bin/activate
cd ~/virus2/temp/VirusTaxo_Hierarchical

# Train DNA model
python3 train.py \
  --data ./dataset/DNA/seq_data \
  --data_metainfo ./dataset/DNA/train.csv \
  --model_dir ./model/DNA   
# Train RNA model
python3 train.py \
  --data ./dataset/RNA/seq_data \
  --data_metainfo ./dataset/RNA/train.csv \
  --model_dir ./model/RNA
  
## Predict taxonomy of DNA viruses. Please change the k-mer lengh to 21 in the `config.py` file.
python3 predict.py --input /data4/machuang/virus2/result/vOTUs.fna \
--model_dir ./model/DNA/ > ./output.txt
##test
python3 test.py \
  --type single \
  --input /data4/machuang/virus2/result/vOTUs.fna \
  --model_dir ./model/DNA \
  --output ./output.txt
#报错，提issue
#FileNotFoundError: [Errno 2] No such file or directory: './dataset/DNA/root'
#因为没有运行模型训练脚本    
#退出./environment/bin/activate环境
deactivate
4.4  DemoVir（Family level）
https://github.com/omics-lab/VirusTaxo_Hierarchical
# （作者提示该软件自2018年来就没有维护，数据库过旧，不能再用于注释了）
# 近几年发表的文章中，还在用此软件
#于是对数据库及注释文件进行更新
#安装依赖
cd ~/virus2/temp
git clone https://github.com/feargalr/Demovir.git
mv Demovir demovir
chmod -R 755 demovir

#运行,环境中需要有R，hmmer中装的有
conda activate hmmer
cd ~/virus2/temp/demovir
#加载prodigal
export PATH="$HOME/virus1/temp/Prodigal/:$PATH"
#加载usearch
export PATH="$HOME/db/software/:$PATH"
#数据库添加权限
chmod 755 ~/virus2/temp/demovir/uniprot_trembl.viral.udb
#demovir.sh：第11行改为
#usearch -ublast AA.fasta -db $DIR/uniprot_trembl.viral.udb -evalue 1e-5 -trunclabels -blast6out trembl_ublast.viral.txt -threads $2 &> /dev/null
#demovir.R：第85、86行改为
# taxassO[sums_o>1] = "Unassigned"
# taxassF[sums_f>1] = "Unassigned"
#demovir.R：第88行
# towrite = data.frame(levels(factor(results3$contig)),taxassO,margin_o,taxassF,margin_f)
#demovir.R：第10行
# taxa_file = paste(script.basename,"uniprot_trembl_viral.RDS",sep="/")
#确保自己temp/demovir里面有数据库uniprot_trembl.viral.udb和注释文件uniprot_trembl_viral.RDS ，没有的话去/data4/machuang/virus2/temp/demovir里面复制，或者参考 ‘0病毒组软件安装’ 下载并处理
./demovir.sh ~/virus2/result/vOTUs.fna 8
#结果文件DemoVir_assignments.txt
#order : 123
#family : 820
4.5  vConTACT2
https://bitbucket.org/MAVERICLab/vcontact2/src/master/
# (kingdom、phylum、class、order、family、genus都有)
mkdir -p ~/virus2/temp/vcontact2 
cd ~/virus2/temp/vcontact2
conda activate vcontact2

#使用自带函数构建[gene-to-genome mapping file]
vcontact2_gene2genome -p ~/virus2/result/vOTUs.faa -o ../taxonomy/vcontact2/vOTUs_g2g.tsv -s Prodigal-FAA
#出现警告，不影响结果Warning: 'U' mode is deprecated
 
#运行
vcontact2 --raw-proteins ~/virus2/result/vOTUs.faa  --proteins-fp ../taxonomy/vcontact2/vOTUs_g2g.tsv --db 'ProkaryoticViralRefSeq211-Merged' --threads 20 --output-dir ../taxonomy/vcontact2/ --c1-bin cluster_one-1.0.jar
#运行过程中会有WARNING，但不是错误，请忽略
cp ../taxonomy/vcontact2/genome_by_genome_overview.csv ../taxonomy/vcontact2/c1.ntw ~/virus2/result/taxonomy/vcontact2/
4.6  PhaGCN和PhaGCN2.0
PhaGCN 24/9/24 11:26  22:04 10.5h
https://github.com/KennthShang/PhaGCN

# 科水平 13096/18898
conda activate phagcn
mkdir -p  ~/virus2/result/taxonomy/phagcn
cd ~/virus2/temp/phagcn
#运行
python run_Speed_up.py --contigs ~/virus2/result/vOTUs.fna --len 5000

cp final_prediction.csv out/gene_to_genome.csv out/contig_gene_to_genome.csv out/network.ntw ~/virus2/result/taxonomy/phagcn/
PhaGCN2.0   9545/18898   
v 2.3.0 ：界门纲目科水平
https://github.com/KennthShang/PhaGCN2.0
cd ~/virus2/temp/phagcn2.0
mkdir -p  ~/virus2/result/taxonomy/phagcn2.0
conda activate phagcn2
export MKL_SERVICE_FORCE_INTEL=1
#使用方法
python run_Speed_up.py --contigs ~/virus2/result/vOTUs.fna --len 5000
#参数设置--contigs #需要输入的contigs文件，通过megahit组装--len #对contigs长度进行阈值设置，默认值为8000 bp
#注意事项
1.确认输入的contigs为病毒contigs,可通过Virsorter或DeepVirFinder获得
2.congtigs序列中的id不能含有空格等中文字符
#保存重要文件
cp result/*.{ntw,csv} ~/virus2/result/taxonomy/phagcn2.0/
#下载某文件到本地
scp machuang@192.168.60.214:/data4/machuang/virus2/temp/deephage/output/vOTUs_lifestyle.csv /Users/14346/downloads/
5.生活方式预测（只保留完整病毒的生活预测结果）
5.1 Bacphlip
https://github.com/adamhockenberry/bacphlip
#输入的噬菌体必须是完整的
conda activate bacphlip
cd ~/virus2/
mkdir -p temp/bacphlip result/bacphlip
cd ~/virus2/temp/bacphlip

bacphlip -i ~/virus2/result/vOTUs.fna --multi_fasta
#结果在~/virus2/result/vOTUs.fna.bacphlip文件内
cd ~/virus2
#提取结果，标明病毒哪些是Temperate，哪些是Virulent
awk 'NR==1 {next}  
BEGIN {  
    print "contig_id\tlifestyle";  
}  
{  
    lifestyle = ($2 > $3) ? "Virulent" : ($2 < $3) ? "Temperate" : "Unknown";
    print $1 "\t" lifestyle;
}' result/vOTUs.fna.bacphlip > temp/bacphlip/bacphlip.tsv

#7908 Temperate;10974 Virulent; 16 Unknown
cp temp/bacphlip/bacphlip.tsv result/bacphlip/
5.2 Deephage
https://github.com/shufangwu/DeePhage
#根据DeePhage评分，
温和型（score ≤ 0.3）
不确定温和型（0.3 < score ≤ 0.5）
不确定裂解型（0.5 < score ≤ 0.7）
裂解型（score > 0.7）
conda activate deephage
mkdir -p ~/virus2/result/deephage
cd ~/virus2/temp/deephage
#加载环境
export LD_LIBRARY_PATH="$HOME/db/software/MCR/v94/runtime/glnxa64:$HOME/db/software/MCR/v94/bin/glnxa64:$HOME/db/software/MCR/v94/sys/os/glnxa64:$HOME/db/software/MCR/v94/extern/bin/glnxa64:$LD_LIBRARY_PATH"
#运行命令
#输入文件不能太大,当输入为900M时会出现下列报错，大文件可分开跑
#拆分
cd ~/virus2
mkdir -p temp/deephage/input temp/deephage/output result/deephage
awk -v n=$(grep -c "^>" ~/virus2/result/vOTUs.fna) 'BEGIN {s=int((n+2)/6)} /^>/ {if (++c > s) {c=1; f++}} {print > "temp/deephage/input/input" f+1 ".fna"}' ~/virus2/result/vOTUs.fna
cd temp/deephage/input/
seqkit stat input1.fna # 3150 353Mb
seqkit stat input2.fna # 3150 346Mb
seqkit stat input3.fna # 3150 383Mb
seqkit stat input4.fna # 3150 699Mb
seqkit stat input5.fna # 3150 562Mb
seqkit stat input6.fna # 3148 383Mb

#中间文件夹会重复，因此需单独依次运行
cd ~/virus2/temp/deephage
./DeePhage input/input1.fna output/vOTUs_lifestyle1.csv
./DeePhage input/input2.fna output/vOTUs_lifestyle2.csv
./DeePhage input/input3.fna output/vOTUs_lifestyle3.csv
./DeePhage input/input4.fna output/vOTUs_lifestyle4.csv
./DeePhage input/input5.fna output/vOTUs_lifestyle5.csv
./DeePhage input/input6.fna output/vOTUs_lifestyle6.csv

#合并
cp output/vOTUs_lifestyle1.csv output/vOTUs_lifestyle.csv

for file in output/vOTUs_lifestyle2.csv output/vOTUs_lifestyle3.csv output/vOTUs_lifestyle4.csv output/vOTUs_lifestyle5.csv output/vOTUs_lifestyle6.csv; do
tail -n +2 "$file" >> output/vOTUs_lifestyle.csv
done
wc -l output/vOTUs_lifestyle.csv #带表头18899
#根据得分划分四个等级
cd ~/virus2/temp/deephage/
awk -F, 'NR==1 {print $1","$2","$3",lifestyle"; next} {
    if ($3 <= 0.3) {
        $4 = "Temperate";
    } else if ($3 > 0.3 && $3 <= 0.5) {
        $4 = "Uncertain_Temperate";
    } else if ($3 > 0.5 && $3 <= 0.7) {
        $4 = "Uncertain_Virulent";
    } else if ($3 > 0.7) {
        $4 = "Virulent";
    }
    print $1","$2","$3","$4
}' output/vOTUs_lifestyle.csv > ~/virus2/result/deephage/vOTUs_lifestyle.csv



# 10143/18898 temperate
# 5528/18898 Uncertain_Temperate
# 893/18898 virulent
# 2234/18898 Uncertain_Virulent

#输入文件不能太大,当输入为900M时会出现下列报错，大文件可分开跑

#报错
#Warning:Variable 'test matrix'was not savedFor variables larger than 2GB use MAT-file vrsion 7.3 or later.
#> In DeePhage(line 51)
#Using TensorFlow backend
#data5/baidefeng/miniconda3/envs/deephage/libpython3.6/importlib/ bootstrap.py:219:Runtim:Warning: compiletime version 3.5 of modul'tensorflow.python.framework.fast tensor utidoes not match runtime version 3.6
#return f(*args, **kwds)
#Traceback(most recent call last):File "software model prediction all train.pyline 12, in <module>test matrix= file['test matrix'l:l.rave.
#KeyError:'test matrix
#Index exceeds array bounds.
#Error in DeePhage:(line 78)
#MATLAB:badsubscript

#conda activate deephage
#python -c "import tensorflow as tf; print(tf.__version__)" #1.4.0
#save('filename.mat', '-v7.3')
#pip install --upgrade tensorflow

5.3 PhaTYP
https://github.com/KennthShang/PhaTYP
cd ~/virus2/temp/PhaTYP
conda activate phatyp
#输出结果默认在当前文件夹下，结果文件为prediction.csv
python preprocessing.py --contigs ~/virus2/result/vOTUs.fna
#python PhaTYP.py --out prediction.csv

# 12458/18681 temperate
# 6222/18681 virulent
6.流行率计算prevalence
#准备abundance.txt文件，去除原始文件中第二行unmapped
cd virus2
mkdir -p result/iqtree/plan
mkdir -p temp/iqtree/plan
sed '2d' result/coverm/abundance.txt > temp/iqtree/plan/abundance.txt
#把非0的数都变成1
awk -v OFS="\t" 'NR==1 {print; next} { $1 = $1; for(i=2; i<=NF; i++) if ($i != 0) $i = 1; print }' temp/iqtree/plan/abundance.txt > abundance.tmp && mv abundance.tmp temp/iqtree/plan/abundance.txt

#再计算所有votu的流行率，代码中 avg = sum / 279需根据自己样本数量更改
awk '  
    NR==1 {          
        gsub("Genome", "contig_id", $0);  
        print $0 "\tprevalence";  
        next  
    }  
    {      
        sum = 0;  
        for(i=2; i<=NF; i++) {  
            sum += $i;  
        }  
        avg = sum / 279;  
        print $0 "\t" avg;  
    }  
' temp/iqtree/plan/abundance.txt > temp/iqtree/plan/prevalence.txt

awk '{print $1 "\t" $NF}' temp/iqtree/plan/prevalence.txt > result/iqtree/plan/prevalence.txt
#制作core_individual箱线图原始文件core_individual.csv
sed '2d' ~/virus2/result/coverm/abundance.txt > ~/virus2/result/coverm/abundance_tmp.txt
cut -f2 ~/virus2/result/iqtree/plan/prevalence.txt > ~/virus2/result/iqtree/plan/prevalence_col2.txt
paste ~/virus2/result/coverm/abundance_tmp.txt ~/virus2/result/iqtree/plan/prevalence_col2.txt | awk 'BEGIN {OFS="\t"} {$2=$NF FS $2; NF--; print}' > ~/virus2/result/iqtree/plan/core_individual.csv

rm ~/virus2/result/coverm/abundance_tmp.txt ~/virus2/result/iqtree/plan/prevalence_col2.txt

#根据每组样本筛选，分为四个组CELY,根据自己分组情况
#######计算C组流行率#######
#######先得到C组原始文件
awk -v FS="\t" -v OFS="\t" '
NR==FNR {
    cols[$1];
    next
}
FNR==1 {
    # 初始化header数组
    for (i=1; i<=NF; i++) {
        if (i == 1 || $i in cols) {
            header[i] = 1
        }
    }
    # 打印标题行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
    next
}
{
    # 打印数据行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
}
' result/metadataC.txt temp/iqtree/plan/abundance.txt > temp/iqtree/plan/abundanceC.txt  #得出的结果如果第二行跑到表头上了，需手动调整temp/iqtree/plan/abundanceC.txt，按个回车，并确保表头最后一行只有一个制表符，此处正在测试，先手动调整。

#再计算C组流行率，代码中 avg = sum / 100需根据自己样本数量更改
awk '  
    NR==1 {          
        gsub("Genome", "contig_id", $0);  
        print $0 "\tprevalence";  
        next  
    }  
    {      
        sum = 0;  
        for(i=2; i<=NF; i++) {  
            sum += $i;  
        }  
        avg = sum / 100;  
        print $0 "\t" avg;  
    }  
' temp/iqtree/plan/abundanceC.txt > temp/iqtree/plan/prevalenceC.txt

awk '{print $1 "\t" $NF}' temp/iqtree/plan/prevalenceC.txt > result/iqtree/plan/prevalenceC.txt

#######计算E组流行率，#######
#####先得到E组原始文件
awk -v FS="\t" -v OFS="\t" '
NR==FNR {
    cols[$1];
    next
}
FNR==1 {
    # 初始化header数组
    for (i=1; i<=NF; i++) {
        if (i == 1 || $i in cols) {
            header[i] = 1
        }
    }
    # 打印标题行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
    next
}
{
    # 打印数据行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
}
' result/metadataE.txt temp/iqtree/plan/abundance.txt > temp/iqtree/plan/abundanceE.txt
#得出的结果如果第二行跑到表头上了，需手动调整temp/iqtree/plan/abundanceE.txt，按个回车，并确保表头最后一行只有一个制表符，此处正在测试，先手动调整。

#再计算E组流行率，代码中 avg = sum / 54需根据自己样本数量更改
awk '  
    NR==1 {          
        gsub("Genome", "contig_id", $0);  
        print $0 "\tprevalence";  
        next  
    }  
    {      
        sum = 0;  
        for(i=2; i<=NF; i++) {  
            sum += $i;  
        }  
        avg = sum / 54;  
        print $0 "\t" avg;  
    }  
' temp/iqtree/plan/abundanceE.txt > temp/iqtree/plan/prevalenceE.txt

awk '{print $1 "\t" $NF}' temp/iqtree/plan/prevalenceE.txt > result/iqtree/plan/prevalenceE.txt
#########计算L组流行率##########
#####先得到L组原始文件
awk -v FS="\t" -v OFS="\t" '
NR==FNR {
    cols[$1];
    next
}
FNR==1 {
    # 初始化header数组
    for (i=1; i<=NF; i++) {
        if (i == 1 || $i in cols) {
            header[i] = 1
        }
    }
    # 打印标题行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
    next
}
{
    # 打印数据行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
}
' result/metadataL.txt temp/iqtree/plan/abundance.txt > temp/iqtree/plan/abundanceL.txt   #得出的结果如果第二行跑到表头上了，需手动调整temp/iqtree/plan/abundanceL.txt，按个回车，并确保表头最后一行只有一个制表符，此处正在测试，先手动调整。

#再计算L组流行率，代码中 avg = sum / 60需根据自己样本数量更改
awk '  
    NR==1 {          
        gsub("Genome", "contig_id", $0);  
        print $0 "\tprevalence";  
        next  
    }  
    {      
        sum = 0;  
        for(i=2; i<=NF; i++) {  
            sum += $i;  
        }  
        avg = sum / 60;  
        print $0 "\t" avg;  
    }  
' temp/iqtree/plan/abundanceL.txt > temp/iqtree/plan/prevalenceL.txt

awk '{print $1 "\t" $NF}' temp/iqtree/plan/prevalenceL.txt > result/iqtree/plan/prevalenceL.txt

##########计算Y组流行率，代码中 avg = sum / 65需根据自己样本数量更改#########
#####先得到Y组原始文件
awk -v FS="\t" -v OFS="\t" '
NR==FNR {
    cols[$1];
    next
}
FNR==1 {
    # 初始化header数组
    for (i=1; i<=NF; i++) {
        if (i == 1 || $i in cols) {
            header[i] = 1
        }
    }
    # 打印标题行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
    next
}
{
    # 打印数据行
    first=1
    for (i=1; i<=NF; i++) {
        if (header[i]) {
            if (first) {
                printf "%s", $i
                first=0
            } else {
                printf OFS "%s", $i
            }
        }
    }
    print ""  # 打印换行符
}
' result/metadataY.txt temp/iqtree/plan/abundance.txt > temp/iqtree/plan/abundanceY.txt   #得出的结果如果第二行跑到表头上了，需手动调整temp/iqtree/plan/abundanceY.txt，按个回车，并确保表头最后一行只有一个制表符，此处正在测试，先手动调整。

######再计算Y组流行率，代码中 avg = sum / 65需根据自己样本数量更改
awk '  
    NR==1 {          
        gsub("Genome", "contig_id", $0);  
        print $0 "\tprevalence";  
        next  
    }  
    {      
        sum = 0;  
        for(i=2; i<=NF; i++) {  
            sum += $i;  
        }  
        avg = sum / 65;  
        print $0 "\t" avg;  
    }  
' temp/iqtree/plan/abundanceY.txt > temp/iqtree/plan/prevalenceY.txt

awk '{print $1 "\t" $NF}' temp/iqtree/plan/prevalenceY.txt > result/iqtree/plan/prevalenceY.txt

######计算总样本各病毒流行率，代码中 avg = sum / 279需根据自己样本数量更改#####
#awk '  
#    NR==1 {          
#        gsub("Genome", "contig_id", $0);  
#        print $0 "\tprevalence";  
#        next  
#    }  
#    {      
#        sum = 0;  
#        for(i=2; i<=NF; i++) {  
#            sum += $i;  
#        }  
#        avg = sum / 279;  
#        print $0 "\t" avg;  
#    }  
#' temp/iqtree/plan/abundance.txt > temp/iqtree/plan/prevalence.txt

#awk '{print $1 "\t" $NF}' temp/iqtree/plan/prevalence.txt > result/iqtree/plan/prevalence.txt

7.中位TPM的vOTU是否在组间富集

8.建树
8.1 IQtree
#使用vOTUs.fna 18898
方法1 生成iqtree的输入vOTUs.concat.faa文件（只用MGV数据库蛋白和我们的蛋白）

conda activate hmmer
cd ~/virus2
mkdir -p temp/iqtree/input result/iqtree
---------------------------------------------- 虚线内的只用运行1次
#预测MGV中代表性votu的蛋白
# export PATH="$HOME/virus1/temp/Prodigal/:$PATH"
prodigal -i ~/db/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.fna -o ~/db/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.genes -a ~/db/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.faa -p meta
seqkit stat ~/db/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.faa
-------------------------------------------- 虚线内的只用运行1次
#合并自己数据和MGV数据
cat result/vOTUs.faa ~/db/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.faa > temp/iqtree/input/vOTUs.mgv.faa
#加载trimal环境
export PATH="$HOME/db/software/trimal/source:$PATH"

#要确保markerout是空文件夹，如果以前的结果在里面，需要删除，否则会报错KeyError:
cp ~/virus1/VOG77.hmm temp/iqtree ./

python ~/virus1/temp/vcentenarian/code/06_evaluation/ViralMarkerTree77.py --in_faa temp/iqtree/input/vOTUs.mgv.faa --out_dir temp/iqtree/markerout --threads 10 #生成concat.faa、hmmsearch.txt文件
#提取MGV的concat.faa
awk '/^>MGV/{f=1} f; /^>/ && f{f=0}' temp/iqtree/markerout/concat.faa > temp/iqtree/markerout/mgv.concat.txt
 wc -l temp/iqtree/markerout/mgv.concat.txt #7792
 
awk '/^>MGV/{f=1; print; next} f && !/^>/{print} /^>/ && !/^>MGV/{f=0}' temp/iqtree/markerout/concat.faa > temp/iqtree/markerout/mgv.concat.faa
seqkit stat temp/iqtree/markerout/mgv.concat.faa #7792
#提取出自己数据的concat.faa
 awk '/^>k141/{f=1} f; /^>/ && f{f=0}' temp/iqtree/markerout/concat.faa > temp/iqtree/markerout/vOTUs.concat.txt
 wc -l temp/iqtree/markerout/vOTUs.concat.txt #1600
cp temp/iqtree/markerout/vOTUs.concat.txt result/iqtree
sed -i 's/>//g; s/ percent_gaps=.*//g' result/iqtree/vOTUs.concat.txt 
 
#awk '/^>/ && /::/ {print}' temp/iqtree/markerout/concat.faa > temp/iqtree/markerout/vOTUs.concat.txt
 
awk '/^>k141/{f=1; print; next} f && !/^>/{print} /^>/ && !/^>k141/{f=0}' temp/iqtree/markerout/concat.faa > temp/iqtree/markerout/vOTUs.concat.faa
seqkit stat temp/iqtree/markerout/vOTUs.concat.faa #1600/18898

#awk '/^>/ && /::/{f=1; print; next} f && !/^>/{print} /^>/ && (!/^>/ || $0 !~ /::/){f=0}' temp/iqtree/markerout/concat.faa > temp/iqtree/markerout/vOTUs.concat.faa
#seqkit stat temp/iqtree/markerout/vOTUs.concat.faa  #597

#运行iqtree
conda activate iqtree
mkdir -p ~/virus2/temp/iqtree/output
#-m MFP 寻找最佳模型，并建树   
iqtree -s ~/virus2/temp/iqtree/markerout/vOTUs.concat.faa -m MFP -n 10 -nt 8 -pre ~/virus2/temp/iqtree/output/vOTUs
#会自己寻找模型并建树，如果想看用的什么模型，结果在~/virus1/temp/iqtree/markerout1/vbvfvs2.concat.faa.log中，显示用的是Q.pfam+I+R4
# -n 10 迭代次数
# -nt 线程数


#保存树文件到result
mkdir -p ~/virus2/result/iqtree
cp ~/virus2/temp/iqtree/output/vOTUs.treefile ~/virus2/result/iqtree/

方法2 生成iqtree的输入vOTUs.concat.faa文件（用5个数据库蛋白和我们的蛋白，MGV、CHVD、GPD、GVD 、CHGV  ）
conda activate hmmer
cd ~/virus2
mkdir -p temp/iqtree2/input  result/iqtree2
#合并五个数据库的蛋白和我们的蛋白
cat result/vOTUs.faa ~/db/checkv/MGV/MGV_HQMQ.faa ~/db/checkv/CHVD/CHVD_HQMQ.faa ~/db/checkv/GPD/GPD_HQMQ.faa ~/db/checkv/GVD/GVD_HQMQ.faa ~/db/checkv/CHGV/CHGV_HQMQ.faa > temp/iqtree2/input/vOTUs.five.faa 
seqkit stat temp/iqtree2/input/vOTUs.five.faa # 14050915 4.67Gb

#加载trimal环境
export PATH="$HOME/db/software/trimal/source:$PATH"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH

#要确保markerout是空文件夹，如果以前的结果在里面，需要删除，否则会报错KeyError:
cp ~/virus1/VOG77.hmm temp/iqtree2

python ~/virus1/temp/vcentenarian/code/06_evaluation/ViralMarkerTree77.py --in_faa temp/iqtree2/input/vOTUs.five.faa --out_dir temp/iqtree2/markerout --threads 10 #生成concat.faa、hmmsearch.txt文件

#提取出自己数据的concat.faa
awk '/^>k141/{f=1} f; /^>/ && f{f=0}' temp/iqtree2/markerout/concat.faa > temp/iqtree2/markerout/vOTUs.concat.txt

wc -l temp/iqtree2/markerout/vOTUs.concat.txt #1245
cp temp/iqtree2/markerout/vOTUs.concat.txt result/iqtree2
sed -i 's/>//g; s/ percent_gaps=.*//g' result/iqtree2/vOTUs.concat.txt 
  
awk '/^>k141/{f=1; print; next} f && !/^>/{print} /^>/ && !/^>k141/{f=0}' temp/iqtree2/markerout/concat.faa > temp/iqtree2/markerout/vOTUs.concat.faa
seqkit stat temp/iqtree2/markerout/vOTUs.concat.faa #1245/18898

#运行iqtree
conda activate iqtree
mkdir -p temp/iqtree2/output
#-m MFP 寻找最佳模型，并建树   
iqtree -s ~/virus2/temp/iqtree2/markerout/vOTUs.concat.faa -m MFP -n 2 -nt 100 -pre ~/virus2/temp/iqtree2/output/vOTUs
cp temp/iqtree2/output/vOTUs/vOTUs.treefile result/iqtree2/
#会自己寻找模型并建树，如果想看用的什么模型，结果在~/virus1/temp/iqtree/markerout1/vbvfvs2.concat.faa.log中，显示用的是Q.pfam+I+R4
# -n 10 迭代次数
# -nt 线程数


方法3 直接mafft对齐6495个蛋白，作为iqtree输入（直接被kill）
mkdir -p ~/virus2/temp/iqtree1  ~/virus2/temp/iqtree1/output
cd ~/virus2
conda activate iqtree
export PATH="$HOME/virus1/temp/Prodigal/:$PATH"
#6495个vOTU,生成蛋白质文件final-viral-combined.faa
prodigal -i result/vOTUs.fna -o result//vOTUs.genes -a result//vOTUs.faa -p meta
#对齐蛋白文件
mafft --auto --thread 8 result//vOTUs.faa > temp/iqtree1/vOTUs.faa
iqtree -s temp/iqtree1/vOTUs.faa -m MFP -n 4 -nt 16 -pre ~/virus2/temp/iqtree1/output/vOTUs


8.2 nwk格式树文件(测试中)
#蛋白预测，使用Prodigal和fasta36
cd ~/virus3
mkdir -p temp/tree_nwk result/tree_nwk
export PATH="$HOME/virus1/temp/Prodigal/:$PATH"
#输入是最终得到的votu文件HQMQLQ.fna
cat ~/virus3/result/HQMQLQ.fna | prodigal -a ~/virus3/temp/tree_nwk/hqmqlq.faa -p meta > ~/virus3/temp/tree_nwk/HQMQLQ.gbk
#.faa文件就是预测的蛋白文件，.gbk是预测的genus文件
#删除faa文件head中多余信息，只保留id
#sed 's/^[^>]*>\([^#]*\).*$/>\1/' ~/virus3/temp/tree_nwk/hqmqlq.faa > ~/virus3/temp/tree_nwk/HQMQLQ.faa
sed -e 's/>\([^ ]*\).*/\>\1/' -e 's/length_.*_//g' ~/virus3/temp/tree_nwk/hqmqlq.faa > ~/virus3/temp/tree_nwk/HQMQLQ.faa
 

 
export PATH="$HOME/virus1/temp/fasta36/bin/:$PATH"
fasta36 -h #version: 36.3.8i May, 2023
fasta36 ~/virus3/temp/tree_nwk/HQMQLQ.faa ~/virus3/temp/tree_nwk/HQMQLQ.faa -m 8 > ~/virus3/temp/tree_nwk/HQMQLQ.fasta36
#出现*** Warning: ktup = 0 out of range [1..3], reset to 3  

conda activate hmmer #我们加载hmmer环境，里面安装过mcl下面会用到
export PATH="$HOME/virus1/temp/VFCs/:$PATH"  #加载命令所在环境
cat ~/virus3/temp/tree_nwk/HQMQLQ.faa | f2s | seqlengths > ~/virus3/temp/tree_nwk/HQMQLQ.faa.lengths

#cat ~/virus3/temp/tree_nwk/HQMQLQ.fasta36 | joincol ~/virus3/temp/tree_nwk/HQMQLQ.faa.lengths | joincol ~/virus3/temp/tree_nwk/HQMQLQ.faa.lengths 2 | awk '{print $1 "\t" $2 "\t" $11 "\t" $13/$14 "\t" ($8-$7)/(2*$13)+($10-$9)/(2*$14) "\t" ($7+$8-$9-$10)/($13+$14)}' | awk '{if ($3 <= 0.05) print}' | awk '{if ($5 >= 0.4) print}' | awk '{if (sqrt(($4-1)^2) - (sqrt(sqrt($5))-.8) + sqrt($6^2) <= 0.1) print $1 "\t" $2}' | mcl - -o - --abc | awk '{j++; for (i = 1; i <= NF; i++) {print $i "\t" j}}' > ~/virus3/temp/tree_nwk/HQMQLQ.VOGs.tsv

cat ~/virus3/temp/tree_nwk/HQMQLQ.fasta36 | joincol ~/virus3/temp/tree_nwk/HQMQLQ.faa.lengths | joincol ~/virus3/temp/tree_nwk/HQMQLQ.faa.lengths 2 | awk '{if ($13 < $14) {s = ($3*$4)/($13*100)-.75; if (s < 0) {s = 0}} else {s = 0} print $1 "\t" $2 "\t" $11 "\t" $13/$14 "\t" ($8-$7)/(2*$13)+($10-$9)/(2*$14) "\t" ($7+$8-$13)/$13-($9+$10-$14)/$14 "\t" s}' | awk '{if ($3 <= 0.05) print}' | awk '{if ($5 >= 0.4) print}' | awk '{if (sqrt(log($4)^2) - (-0.0181/($5-0.32)+0.23) + sqrt($6^2) <= 0.15 + $7) print $1 "\t" $2}' | mcl - -o - --abc  | awk '{j++; for (i = 1; i <= NF; i++) {print $i "\t" j}}' > ~/virus3/temp/tree_nwk/HQMQLQ.VOGs.tsv
#cat test_2removed.fastab | ./joincol test_2removed.lengths | ./joincol test_2removed.lengths 2 | awk '{print $1 "\t" $2 "\t" $11 "\t" $13/$14 "\t" ($8-$7)/(2*$13)+($10-$9)/(2*$14) "\t" ($7+$8-$9-$10)/($13+$14)}' | awk '{if ($3 <= 0.05) print}' | awk '{if ($5 >= 0.4) print}' | awk '{if (sqrt(($4-1)^2) - (sqrt(sqrt($5))-.8) + sqrt($6^2) <= 0.1) print $1 "\t" $2}' > test_2removed.tsv

awk '{if ($11 <= 0.05) print $1 "\t" $2 "\t" $12}' ~/virus3/temp/tree_nwk/HQMQLQ.fasta36 | rev | sed 's/\t[[:digit:]]\+_/\t/' | rev | sed 's/_[[:digit:]]\+\t/\t/' | sort | hashsums | tree_bray > ~/virus3/temp/tree_nwk/HQMQLQ_votus.mat 
#生成.mat文件报错，Argument "" isn't numeric in addition (+) at /data4/machuang/virus1/temp/VFCs/hashsums line 20, <STDIN> line 3217755. Illegal division by zero at /data4/machuang/virus1/temp/VFCs/tree_bray line 22, <STDIN> line 516127.  提issue

export PATH="$HOME/virus1/temp/rapidnj/bin/:$PATH"
rapidnj -i pd ~/virus3/temp/tree_nwk/HQMQLQ_votus.mat >~/virus3/result/tree_nwk/HQMQLQ_votus.nwk
9.IQtree外圈美化
cd ~/virus2
mkdir -p result/iqtree/plan
#制作宿主预测外圈
awk -F, 'NR==FNR{a[$1]; next} $1 in a' result/iqtree/vOTUs.concat.txt result/iphop/Host_prediction_phylum.csv > result/iqtree/plan/iTOL_HostFamily.csv
#根据宿主，添加颜色
awk -F, -v OFS=, '  
    BEGIN {  
        color["Firmicutes_A"] = "#8ECFC9"  
        color["Firmicutes"] = "#FA7F6F"  
        color["Firmicutes_C"] = "#BEB8DC"  
        color["Proteobacteria"] = "#82B0D2"  
        color["Actinobacteriota"] = "#999999"  
        color["Bacteroidota"] = "#FFBE7A"  
        color["Fusobacteriota"] = "#E7DAD2"  
    }  
    NR=1 {
        if ($2 in color) {  
            $0 = $1 OFS color[$2] OFS substr($0, index($0, $2));  
        }  
        print  
    }  
' result/iqtree/plan/iTOL_HostFamily.csv > result/iqtree/plan/iTOL_HostFamily_colored.csv
#先手动将iTOL_HostFamily_colored.csv内容复制到树美化文件中，后续再优化

#制作Novel预测外圈
cd ~/virus2
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a {print $1 "," $2}' result/iqtree/vOTUs.concat.txt temp/MGV/mgv/new_old/ref_five1.tsv > result/iqtree/plan/iTOL_Novel.csv
#查看第二列有多少Novel
awk -F, '{print $2}' result/iqtree/plan/iTOL_Novel.csv | grep -c "Novel" #169/1600

#添加颜色
awk -F, -v OFS=, '  
    BEGIN {  
        color["Novel"] = "#FF0000"  
        color["ref"] = "#000000"            
    }  
    NR=1 {
        if ($2 in color) {  
            $0 = $1 OFS color[$2] OFS substr($0, index($0, $2));  
        }  
        print  
    }  
' result/iqtree/plan/iTOL_Novel.csv > result/iqtree/plan/iTOL_Novel_colored.csv
#先手动将iTOL_Novel_colored.csv内容复制到树美化文件中，后续再优化

#制作病毒Family预测外圈
cd ~/virus2
awk -F',' 'NR==FNR {a[$1]; next} $1 in a {print $1 "," $3}' result/iqtree/vOTUs.concat.txt result/taxonomy/phagcn/final_prediction.csv > result/iqtree/plan/iTOL_Family.csv


#添加颜色
awk -F, -v OFS=, '  
    BEGIN {  
        color["Siphoviridae"] = "#1f77b4"  
        color["Podoviridae"] = "#ff7f0e" 
        color["Myoviridae"] = "#2ca02c" 
        color["Herelleviridae"] = "#d62728"
        color["Demerecviridae"] = "#9467bd"
        color["Autographiviridae"] = "#8c564b"          
    }  
    NR=1 {
        if ($2 in color) {  
            $0 = $1 OFS color[$2] OFS substr($0, index($0, $2));  
        }  
        print  
    }  
' result/iqtree/plan/iTOL_Family.csv > result/iqtree/plan/iTOL_Family_colored.csv
#先手动将iTOL_Genus_colored.csv内容复制到树美化文件中，后续再优化

#添加颜色
awk -F, -v OFS=, '  
    BEGIN {  
        color["Siphoviridae"] = "#1f77b4"  
        color["Podoviridae"] = "#ff7f0e" 
        color["Myoviridae"] = "#2ca02c" 
        #color["crAss-phage"] = "#d62728"
        color["Microviridae"] = "#9467bd"
        color["Inoviridae"] = "#8c564b"          
    }  
    NR=1 {
        if ($2 in color) {  
            $0 = $1 OFS color[$2] OFS substr($0, index($0, $2));  
        }  
        print  
    }  
' result/iqtree/plan/paper_Family.csv > result/iqtree/plan/paper_Family_colored.csv




#制作病毒Genus预测外圈
cd ~/virus2
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a {print $1 "," $3}' result/iqtree/vOTUs.concat.txt result/taxonomy/virustaxo/virustaxo_genus.txt > result/iqtree/plan/iTOL_Genus.csv


#添加颜色
awk -F, -v OFS=, '  
    BEGIN {  
        color["Toutatisvirus"] = "#1f77b4"  
        color["Taranisvirus"] = "#ff7f0e" 
        color["Moineauvirus"] = "#2ca02c" 
        color["Brigitvirus"] = "#d62728"
        color["Glaedevirus"] = "#9467bd"
        color["Lagaffevirus"] = "#8c564b"          
    }  
    NR=1 {
        if ($2 in color) {  
            $0 = $1 OFS color[$2] OFS substr($0, index($0, $2));  
        }  
        print  
    }  
' result/iqtree/plan/iTOL_Genus.csv > result/iqtree/plan/iTOL_Genus_colored.csv
#先手动将iTOL_Genus_colored.csv内容复制到树美化文件中，后续再优化
10.功能注释
10.1 CARD耐药基因 
10.1.1 rgi （Resistance Gene Identifier v.5.1.0）
识别 CARD 耐药基因
#使用
conda activate rgi
cd ~/virus2
#加载数据库
rgi load --card_json ~/db/card/localDB/card.json --local

# 简化蛋白ID
cut -f 1 -d ' ' ~/virus2/result/vOTUs.faa > temp/rgi/protein.fa
# 这个错误忽略即可，不是报错，没有任何影响  grep: 写错误: 断开的管道
#grep '>' result/NR/protein.fa | head -n 3
grep '>' temp/rgi/protein.fa | head -n 3

## 蛋白层面注释ARG
time rgi main -i temp/rgi/protein.fa -t protein \
  -n 9 -a DIAMOND --include_loose --clean --local --low_quality\
  -o temp/rgi/protein
head -n3 temp/rgi/protein.txt

# contig层面注释ARG
time rgi main -i ~/virus2/result/vOTUs.fna -t contig \
  -n 9 -a DIAMOND --include_loose --clean --local --low_quality\
  -o temp/rgi/contig
head -n3 temp/rgi/contig.txt
10.1.2 Resfams database
识别 CARD 耐药基因
conda activate hmmer
cd ~/virus2
mkdir -p temp/resfams
hmmsearch --tblout temp/resfams/output.tbl --domtblout temp/resfams/output.domtbl ~/db/Resfams/Resfams.hmm result/vOTUs.faa
#output.tbl 文件：简洁的表格格式，包含基本的匹配信息，适合快速查看和分析。
#output.domtbl 文件：详细的域比对信息，适合对比对结果进行深入分析。

#过滤bit-score >=50的，并且根据得分去重，选择更高得分的
conda activate vr  #环境中要有pandas
python ~/db/scripts/filter_hmmsearch.py --domtbl temp/resfams/output.domtbl --output temp/resfams/filtered_output.domtbl --bit-score-threshold 50
# $1 contig名称 $4 耐药基因名称  $8 score

10.2 AMG辅助代谢基因
10.2.1 AMRFinder v.3.12.8
识别 AMG 辅助代谢基因
https://github.com/michaelwoodworth/AMRFinder_scripts
cd ~/virus2
conda activate amrfinder
#输入为核苷酸

amrfinder -n result/vOTUs.fna --plus -o temp/amrfinder/vOTUs_fna_amrfinder.tsv
#输入为蛋白
amrfinder -p result/vOTUs.faa --plus -o temp/amrfinder/vOTUs_faa_amrfinder.tsv

#作图 
# 00_amrfinder_filter.py
cd temp/amrfinder
00_amrfinder_filter.py -i vOTUs_faa_amrfinder.tsv -o  filter.tsv -m complete
# 01_amrfinder_binary_matrix.py
01_amrfinder_binary_matrix.py -i filter_amrfinder.tsv -o matrix.tsv
# 02_amrfinder_validate_and_summarize_RPKM.py 
python 02_amrfinder_validate_and_summarize_RPKM.py -a filter.tsv -m ${coverage_magic_path} -o amrfinder_RPKM.tsv -v -V
10.2.2 VIBRANT
识别 AMG 辅助代谢基因
https://github.com/AnantharamanLab/VIBRANT
cd ~/virus2
conda activate VIBRANT
python3 temp/VIBRANT/VIBRANT_run.py -i ~/virus2/result/vOTUs.fna -l 5000 -t 8 -folder ~/virus2/temp/VIBRANT/AMG -d ~/virus2/temp/VIBRANT/databases
10.2.3 DRAM
识别 AMG 辅助代谢基因
https://github.com/WrightonLabCSU/DRAM
conda activate DRAM
cd ~/virus2
mkdir -p result/DRAM
#使用
#annotate some MAGs
temp/DRAM/scripts/DRAM-v.py annotate -i result/vOTUs.fna -o annotation --bit_score_threshold 50 --threads 8
#summarize annotations
~/virus2/temp/DRAM/scripts/DRAM-v.py distill -i annotation/annotations.tsv -o genome_summaries --trna_path annotation/trnas.tsv --rrna_path annotation/rrnas.tsv
#结果文件
DRAM-v liquor是对已在带注释的病毒重叠群中检测到的潜在 AMG (pAMG) 的总结，HTML 文件的形式。
annotations.tsv包含所有预测的开放阅读框架的所有注释
vMAG_stats.tsv提供有关每个病毒重叠群的详细信息
10.3 EGGnog功能注释
识别 COG、KO、CAZy
在线分析http://eggnog-mapper.embl.de/
10.3.1 prodigal预测基因序列
cd ~/virus2
conda activate megahit

time prodigal -i result/vOTUs.fna \
    -d result/gene.fa \
    -o result/gene.gff \
    -p meta -f gff > result/gene.log 2>&1 
10.3.2 cdhit
# 输入文件：prodigal预测的基因序列 result/gene.fa
# 输出文件：去冗余后的基因和蛋白序列：result/eggnog/nucleotide.fa, result/eggnog/protein.fa
cd ~/virus2
conda activate megahit

mkdir -p result/eggnog
# aS覆盖度，c相似度，G局部比对，g最优解，T多线程，M内存0不限制
# 2万基因2m，3M384p15m，2千万需要2000h，多线程可加速
cd-hit-est -i result/gene.fa \
    -o result/eggnog/nucleotide.fa \
    -aS 0.9 -c 0.95 -G 0 -g 0 -T 8 -M 0
# 统计非冗余基因数量，单次拼接结果数量下降不大，如3M-2M，多批拼接冗余度高
grep -c '>' result/eggnog/nucleotide.fa
# 翻译核酸为对应蛋白序列, --trim去除结尾的*
seqkit translate --trim result/eggnog/nucleotide.fa \
    > result/eggnog/protein.fa 
10.3.3 salmon定量
# 输入文件：去冗余后的基因序列：result/eggnog/nucleotide.fa
# 输出文件：Salmon定量：result/vsalmon/gene.count, gene.TPM

mkdir -p temp/vsalmon
salmon -v # 1.8.0

# 建索引, -t序列, -i 索引，10s
salmon index -t result/eggnog/nucleotide.fa \
  -p 3 -i temp/vsalmon/index 

# 定量，l文库类型自动选择，p线程，--meta宏基因组
# 多个任务并行, 18s30m
time tail -n+2 result/metadata.txt | cut -f1 | rush -j 3 \
  "salmon quant -i temp/vsalmon/index -l A -p 6 --meta \
    -1 /data7/public/lyx/age/meta/temp/hr/{1}_1.fastq -2 /data7/public/lyx/age/meta/temp/hr/{1}_2.fastq \
    -o temp/vsalmon/{1}.quant"

# 合并
mkdir -p result/vsalmon
salmon quantmerge --quants temp/vsalmon/*.quant \
    -o result/vsalmon/gene.TPM
salmon quantmerge --quants temp/vsalmon/*.quant \
    --column NumReads -o result/vsalmon/gene.count
sed -i '1 s/.quant//g' result/vsalmon/gene.*
#计算FPKM
#python ~/db/scripts/calculate_RPKM.py --quants "temp/vsalmon/*.quant/quant.sf" --tpm "result/vsalmon/gene.TPM" --counts "result/vsalmon/gene.count" --output "result/vsalmon/gene.FPKM"

# 预览结果表格
head -n3 result/vsalmon/gene.*
10.3.4 COG/KO/CAZy丰度汇总表
#运行
conda activate eggnog
db=/data/meta/db
cd ~/virus2
#数据库宏基因组流程里有，在/db/eggnog
time emapper.py --data_dir /db/eggnog \
  -i result/eggnog/protein.fa --cpu 6 -m diamond --override \
  -o temp/eggnog/output

# 格式化结果并显示表头
grep -v '^##' temp/eggnog/output.emapper.annotations | sed '1 s/^#//' \
  > temp/eggnog/output
csvtk -t headers -v temp/eggnog/output

# 生成COG/KO/CAZy丰度汇总表
mkdir -p result/eggnog
# 显示帮助
summarizeAbundance.py -h
# 汇总，7列COG_category按字母分隔，12列KEGG_ko和19列CAZy按逗号分隔，原始值累加
summarizeAbundance.py \
  -i result/vsalmon/gene.TPM \
  -m temp/eggnog/output --dropkeycolumn \
  -c '7,12,19' -s '*+,+,' -n raw \
  -o result/eggnog/eggnog
sed -i 's#^ko:##' result/eggnog/eggnog.KEGG_ko.raw.txt
sed -i '/^-/d' result/eggnog/eggnog*
head -n3 result/eggnog/eggnog*
# eggnog.CAZy.raw.txt  eggnog.COG_category.raw.txt  eggnog.KEGG_ko.raw.txt

# 添加注释生成STAMP的spf格式
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
  ${db}/EasyMicrobiome/kegg/KO_description.txt \
  result/eggnog/eggnog.KEGG_ko.raw.txt | \
  sed 's/^\t/Unannotated\t/' \
  > result/eggnog/eggnog.KEGG_ko.TPM.spf
head -n 5 result/eggnog/eggnog.KEGG_ko.TPM.spf
# KO to level 1/2/3
summarizeAbundance.py \
  -i result/eggnog/eggnog.KEGG_ko.raw.txt \
  -m ${db}/EasyMicrobiome/kegg/KO1-4.txt \
  -c 2,3,4 -s ',+,+,' -n raw --dropkeycolumn \
  -o result/eggnog/KEGG
head -n3 result/eggnog/KEGG*

# CAZy
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
   ${db}/EasyMicrobiome/dbcan2/CAZy_description.txt result/eggnog/eggnog.CAZy.raw.txt | \
  sed 's/^\t/Unannotated\t/' > result/eggnog/eggnog.CAZy.TPM.spf
head -n 3 result/eggnog/eggnog.CAZy.TPM.spf
# CAZy to Description
summarizeAbundance.py \
  -i result/eggnog/eggnog.CAZy.raw.txt \
  -m ${db}/EasyMicrobiome/dbcan2/CAZy_description.txt \
  -c 1,2 -s ',+' -n raw --dropkeycolumn \
  -o result/eggnog/CAZy
head -n3 result/eggnog/CAZy*

# COG
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
  ${db}/EasyMicrobiome/eggnog/COG.anno result/eggnog/eggnog.COG_category.raw.txt > \
  result/eggnog/eggnog.COG_category.TPM.spf
head -n 3 result/eggnog/eggnog.COG_category.TPM.spf
10.4 prokka基因注释
#生成gff文件  --outdir指定输出目录，--prefix 指定输出文件前缀 ，--centre X --compliant自动生成符合条件的contig ID 长度
cd ~/virus2
conda activate prokka

prokka --force --outdir temp/prokka/gff --prefix vOTUs --centre X --compliant --kingdom Viruses --gffver 3 result/vOTUs.fna
#--force            Force overwriting existing output folder (default OFF)
#--gffver 3 输出的gff格式为9列

#生成的文件
gff：This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV.
gbk：This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence.
fna：Nucleotide FASTA file of the input contig sequences.
faa： Protein FASTA file of the translated CDS sequences.
ffn：Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
sqn：An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc.
fsa：Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines.
tbl：Feature Table file, used by "tbl2asn" to create the .sqn file.
err：Unacceptable annotations - the NCBI discrepancy report.
tsv：Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product
txt：Statistics relating to the annotated features found.
log：Contains all the output that Prokka produced during its run. This is a record of what settings you used, even if the --quiet option was enabled.

11.汇总表
#add checkv
echo -e "contig_id\tcheckv_quality" > result/master_table.tsv && awk 'BEGIN { FS=OFS="\t" } NR==FNR {vOTUs[$1]=1; next} $1 in vOTUs {print $1, $8}' result/vOTUs.list temp/checkv/HQMQLQ.txt >> result/master_table.tsv

#add Genus
echo -e "contig_id\tcheckv_quality\tGenus" > result/master_table_new.tsv && awk 'BEGIN { FS=OFS="\t" } FNR==NR {master_table[$1]=$2; next} FNR==1 {next} {print $1, master_table[$1], $3}' result/master_table.tsv result/taxonomy/virustaxo/output.txt >> result/master_table_new.tsv && mv result/master_table_new.tsv result/master_table.tsv

#add Kingdom、Phylum、Class、Order、Family

# 创建新文件并添加表头
(echo -e "contig_id\tcheckv_quality\tGenus\tKingdom\tPhylum\tClass\tOrder\tFamily" && awk 'BEGIN { FS=OFS="\t" }
NR==FNR {ictv[$1]=$0; next}
FNR==1 {next}
{
    if ($3 in ictv) {
        split(ictv[$3], a, "\t");
        print $0, a[2], a[3], a[4], a[5], a[6]
    } else {
        print $0, "NA", "NA", "NA", "NA", "NA"
    }
}
' result/taxonomy/ICTV.txt result/master_table.tsv) > result/master_table_new.tsv && mv result/master_table_new.tsv result/master_table.tsv
#补充空白为NA
awk 'BEGIN {FS=OFS="\t"} {for(i=1; i<=NF; i++) if($i == "") $i="NA"; print}'  result/master_table.tsv > result/master_table_new.tsv && mv result/master_table_new.tsv result/master_table.tsv

