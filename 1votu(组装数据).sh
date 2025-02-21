
## 1votu(组装数据)

1.准备数据
#新建文件夹
mkdir -p /data4/machuang/virus2
cd /data4/machuang/virus2
mkdir -p seq temp result
# result  储存重要结果，及metadata文件
# temp    存放临时文件
# seq     存放宏基因组分析得到的组装结果
### C:90  Y:34 L:62 E:28 NA: 3  共计217个样本

#如果是单个组装，需要将单个组装得到的final_assembly.fasta合并。在temp/megahit文件夹下
find /data1/liuyongxin/age/meta/temp/megahit/ -type f -name "final_assembly.fasta" -exec cat {} \; > ~/virus2/seq/final_assembly.fasta
seqkit stat ~/virus2/seq/final_assembly.fasta #num_seq : 8646893  33Gb
#num_seqs    sum_len     min_len  avg_len    max_len
#8646893   34786222477     202    4023      1526733

#py脚本去重（不去重genomad会报错）
~/db/scripts/filter.py ~/virus2/seq/final_assembly.fasta ~/virus2/seq/no_duplicates1.fna
seqkit stat ~/virus2/seq/no_duplicates1.fna #num_seq : 8644224  32.7Gb
#num_seqs   sum_len     min_len  avg_len    max_len
#8644224   34786222477    202    4024.2     1526733

#seqkit stat统计的是所有唯一序列的总长度，所以sum_len不变

rm ~/virus2/seq/final_assembly.fasta
2.使用六种病毒鉴定软件
2.1 geNomad  (57.5h)
https://github.com/apcamargo/genomad
conda activate genomad
mkdir -p ~/virus2/temp/genomad_high ~/virus2/temp/genomad_moderate ~/virus2/result/genomad 
cd ~/virus2
#genomad的使用
#生成High confidence病毒,命令为--conservative
genomad end-to-end --cleanup -t 8 seq/no_duplicates.fna temp/genomad_high/ ~/db/genomad_db --conservative
wc -l temp/genomad_high/no_duplicates_summary/no_duplicates_virus_summary.tsv #108635
#去除多余字符，且length ≥ 2500
awk '$2 >= 2500 {
    print $1, $2
}' OFS="\t" temp/genomad_high/no_duplicates_summary/no_duplicates_virus_summary.tsv > result/genomad/gn_high.tsv



wc -l result/genomad/gn_high.tsv #60082 
seqkit stat temp/genomad_high/no_duplicates_summary/no_duplicates_virus.fna  #108634 772Mb
# num_seqs      sum_len     min_len  avg_len  max_len
# 108,634     792,686,122      444   7,296.9  412,909
#提取High confidence病毒
awk '{print $1}' result/genomad/gn_high.tsv | seqkit grep -f - temp/genomad_high/no_duplicates_summary/no_duplicates_virus.fna > temp/genomad_high/no_duplicates_summary/gn_high.fna
seqkit stat temp/genomad_high/no_duplicates_summary/gn_high.fna #60081 695Mb
#num_seqs      sum_len    min_len   avg_len   max_len
#60,081      715284535    2,500     11,905.3  412,909

#生成Moderate confidence病毒，命令为--relaxed
genomad end-to-end --cleanup -t 8 seq/no_duplicates.fna temp/genomad_moderate/ ~/db/genomad_db --relaxed
#去除多余字符，且length ≥ 2500
awk '$2 >= 2500 {

    print $1, $2
}' OFS="\t" temp/genomad_moderate/no_duplicates_summary/no_duplicates_virus_summary.tsv > result/genomad/gn_moderate.tsv
wc -l result/genomad/gn_moderate.tsv #235386 
seqkit stat temp/genomad_moderate/no_duplicates_summary/no_duplicates_virus.fna  #1064933 2.7Gb

#提取Moderate confidence病毒
awk '{print $1}' result/genomad/gn_moderate.tsv | seqkit grep -f - temp/genomad_moderate/no_duplicates_summary/no_duplicates_virus.fna > temp/genomad_moderate/no_duplicates_summary/gn_moderate.fna
seqkit stat temp/genomad_moderate/no_duplicates_summary/gn_moderate.fna #235385 607Mb
#num_seqs      sum_len    min_len   avg_len   max_len
#235,385     1646081597    2500     6993.1   412909
2.2 VIBRANT (110.5h)
https://github.com/AnantharamanLab/VIBRANT
conda activate VIBRANT
mkdir -p ~/virus2/temp/VIBRANT
#更换工作环境需复制这些，如果第一次安装，请参考 0病毒组软件安装
cp -r ~/virus1/temp/VIBRANT/databases ~/virus2/temp/VIBRANT/
cp -r ~/virus1/temp/VIBRANT/scripts ~/virus2/temp/VIBRANT/
cp -r ~/virus1/temp/VIBRANT/files ~/virus2/temp/VIBRANT/
cp -r ~/virus1/temp/VIBRANT/.git ~/virus2/temp/VIBRANT/
cp ~/virus1/temp/VIBRANT/VIBRANT_run.py ~/virus2/temp/VIBRANT/
cd temp/VIBRANT/databases
#下载完成的运行
./VIBRANT_setup.py -test
#有乱码是因为浮点问题，脚本代码未更新。不影响结果
#安装完成提示VIBRANT v1.2.1 is good to go!
#See example_data/ for quick test files.rm

#输入8644224个votu
cd ~/virus2
python3 temp/VIBRANT/VIBRANT_run.py -i ~/virus2/seq/no_duplicates.fna -l 2500 -t 8 -folder ~/virus2/temp/VIBRANT/result -d ~/virus2/temp/VIBRANT/databases
wc -l temp/VIBRANT/result/VIBRANT_no_duplicates/VIBRANT_phages_no_duplicates/no_duplicates.phages_combined.txt #45497带表头
seqkit stat temp/VIBRANT/result/VIBRANT_no_duplicates/VIBRANT_phages_no_duplicates/no_duplicates.phages_combined.fna #45496 669Mb
#num_seqs      sum_len   min_len   avg_len   max_len
#45,496    699,882,315    3,010    15,383.4  305,274
#输出45496个病毒

#### High confidence 
echo -e "scaffold" > result/VIBRANT/vb_high.tsv && cat temp/VIBRANT/result/VIBRANT_no_duplicates/VIBRANT_phages_no_duplicates/no_duplicates.phages_combined.txt >> result/VIBRANT/vb_high.tsv

wc -l result/VIBRANT/vb_high.tsv #45497 带表头

#删除中间大文件
cd ~/virus2
rm -rf temp/VIBRANT/result/VIBRANT_no_duplicates/VIBRANT_HMM_tables_parsed_no_duplicates
rm -rf temp/VIBRANT/result/VIBRANT_no_duplicates/VIBRANT_HMM_tables_unformatted_no_duplicates
rm -rf temp/VIBRANT/result/VIBRANT_no_duplicates/*.{ffn,faa,gff}

2.3 VirFinder（ 120h） 
https://github.com/jessieren/VirFinder
conda activate hmmer
cd ~/virus2
mkdir -p result/vf

#输入文件太大，拆分运行,~/virus2/seq/no_duplicates.fna #num_seq : 8644224 
awk -v n=$(grep -c "^>" ~/virus2/seq/no_duplicates.fna) 'BEGIN {s=int((n+2)/3)} /^>/ {if (++c > s) {c=1; f++}} {print > "temp/MGV/mgv/viral_detection_pipeline/input/virfinder_input" f+1 ".fna"}' ~/virus2/seq/no_duplicates.fna
cd temp/MGV/mgv/viral_detection_pipeline/input/
seqkit stat virfinder_input1.fna # 2881408 10.9Gb
seqkit stat virfinder_input2.fna # 2881408 11.1Gb
seqkit stat virfinder_input3.fna # 2881408 10.5Gb
#拆分Run VirFinder 24/8/28/14:57 - -24/9/5 13:01 120h
cd ~/virus2
Rscript temp/MGV/mgv/viral_detection_pipeline/virfinder.R temp/MGV/mgv/viral_detection_pipeline/input/virfinder_input1.fna temp/MGV/mgv/viral_detection_pipeline/output/virfinder1.tsv
wc -l temp/MGV/mgv/viral_detection_pipeline/output/virfinder1.tsv #算上表头2881409
Rscript temp/MGV/mgv/viral_detection_pipeline/virfinder.R temp/MGV/mgv/viral_detection_pipeline/input/virfinder_input2.fna temp/MGV/mgv/viral_detection_pipeline/output/virfinder2.tsv
wc -l temp/MGV/mgv/viral_detection_pipeline/output/virfinder2.tsv #算上表头2881409
Rscript temp/MGV/mgv/viral_detection_pipeline/virfinder.R temp/MGV/mgv/viral_detection_pipeline/input/virfinder_input3.fna temp/MGV/mgv/viral_detection_pipeline/output/virfinder3.tsv
wc -l temp/MGV/mgv/viral_detection_pipeline/output/virfinder3.tsv #算上表头2881409
#合并
cd temp/MGV/mgv/viral_detection_pipeline/output/
awk 'FNR==1 && NR!=1 {next} {print}' virfinder1.tsv virfinder2.tsv virfinder3.tsv > vf.tsv
wc -l vf.tsv #8644225

#结果的表头、内容不规范
cd ~/virus2
awk 'BEGIN {print "name\tlength\tscore\tpvalue"} NR>1 {gsub(/\"/, "", $2); print $2"\t"$3"\t"$4"\t"$5}' OFS="\t" <(sed '1s/.*/name\tlength\tscore\tpvalue/' temp/MGV/mgv/viral_detection_pipeline/output/vf.tsv) > tmp && mv tmp temp/MGV/mgv/viral_detection_pipeline/output/vf_modified.tsv

#High confidence :length >=2500 score ≥ 0.9 p-value < 0.05
awk 'NR==1 || ($2 >= 2500 && $3 >= 0.9 && $4 < 0.05 ) NR>1 {print $1"\t"$2"\t"$3"\t"$4}' OFS="\t" temp/MGV/mgv/viral_detection_pipeline/output/vf_modified.tsv > result/vf/vf_high.tsv
wc -l result/vf/vf_high.tsv  #带表头 46486    

#Moderate confidence :length >=5000   0.9> score ≥ 0.7   p-value < 0.05
awk 'NR==1 || ($2 >= 2500 && $3 < 0.9 && $3 >= 0.7 && $4 < 0.05 ) NR>1 {print $1}' OFS="\t"  temp/MGV/mgv/viral_detection_pipeline/output/vf_modified.tsv > result/vf/vf_moderate.tsv
wc -l result/vf/vf_moderate.tsv  #带表头  112484   
#virfinder：/  %

#提取序列
cut -f1 result/vf/vf_high.tsv | seqkit grep -f - ~/virus2/seq/no_duplicates.fna > temp/MGV/mgv/viral_detection_pipeline/output/vf_high.fna
  
seqkit stat temp/MGV/mgv/viral_detection_pipeline/output/vf_high.fna # 13851  281Mb
#num_seqs      sum_len    min_len   avg_len  max_len
#46,485     288,701,038    2,500    6,210.6  412,909
cut -f1 result/vf/vf_moderate.tsv | seqkit grep -f - ~/virus2/seq/no_duplicates.fna > temp/MGV/mgv/viral_detection_pipeline/output/vf_moderate.fna
 
seqkit stat temp/MGV/mgv/viral_detection_pipeline/output/vf_moderate.fna # 112483 704Mb
#num_seqs      sum_len    min_len   avg_len   max_len
#112,483    722,806,829    2,500    6,425.9  437,137
2.4 DeepVirFinder (100h)
https://github.com/jessieren/DeepVirFinder
conda activate dvf
cd ~/virus2
mkdir -p result/dvf
temp/dvf/dvf.py -i ~/virus2/seq/no_duplicates.fna -m temp/dvf/models -o temp/dvf/result/ -l 2500 -c 1
wc -l temp/dvf/result/no_duplicates.fasta_gt5000bp_dvfpred.txt #1264617

#High confidence :length >=5000 score ≥ 0.9 p-value < 0.05
awk 'NR==1 || ($3 >= 0.9 && $4 < 0.05 )' temp/dvf/result/no_duplicates.fasta_gt5000bp_dvfpred.txt > result/dvf/dvf_high.tsv
wc -l result/dvf/dvf_high.tsv  #带表头71544，无重复

#Moderate confidence :length >=5000   0.9> score ≥ 0.7   p-value < 0.05
awk 'NR==1 || ($3 < 0.9 && $3 >= 0.7 && $4 < 0.05 )' temp/dvf/result/no_duplicates.fasta_gt5000bp_dvfpred.txt > result/dvf/dvf_moderate.tsv
wc -l result/dvf/dvf_moderate.tsv  #带表头99633，无重复

#提取序列
cut -f1 result/dvf/dvf_high.tsv | seqkit grep -f - ~/virus2/seq/no_duplicates.fna > temp/dvf/result/dvf_high.fna
seqkit stat temp/dvf/result/dvf_high.fna  #71543 1.13Gb
#num_seqs        sum_len    min_len    avg_len    max_len
#71,543      1,191,600,250    5,000   16,655.7   1,526,733
cut -f1 result/dvf/dvf_moderate.tsv | seqkit grep -f - ~/virus2/seq/no_duplicates.fna > temp/dvf/result/dvf_moderate.fna  
seqkit stat temp/dvf/result/dvf_moderate.fna  #99632 2.12Gb
#num_seqs        sum_len    min_len   avg_len   max_len
# 99,632     2,241,981,748    5,000   22,502.6  972,405
2.5 Virsorter（h）（可选）
https://github.com/simroux/VirSorter
24/7/29 11:33
#运行
conda activate virsorter
cd ~/virus2
mkdir -p result/vs
temp/virsorter/wrapper_phage_contigs_sorter_iPlant.pl -f ~/virus2/seq/no_duplicates.fna --db 1 --wdir temp/virsorter --ncpu 8 --data-dir ~/db/virsorter/virsorter-data
#--db 1 : RefSeq 数据库更加标准化和经过仔细的注释.
#--db 2 : Virome 数据库更注重广泛性和多样性，适合发现新的或未充分注释的病毒。
#处理结果文件
#python ~/db/scripts/process_virsorter.py /path/to/your/virsorterfile.csv

#High confidence :(categories 1, 2, 4, 5)

#Moderate confidence :  (categorized under cat3, cat6, circular) 

2.6 Virsorter2 (443.25h)
https://github.com/jiarong/VirSorter2
conda activate vs2
cd ~/virus2
mkdir -p result/vs2
#输入文件太大会报错，控制输入文件在30G左右
virsorter run -w temp/vs2/result/no_duplicates -i ~/virus2/seq/no_duplicates.fna --include-groups "dsDNAphage,ssDNA"  --min-score 0.5 --min-length 5000 --keep-original-seq  --seqname-suffix-off -j 16 all
wc -l temp/vs2/result/no_duplicates/final-viral-score.tsv #70819带表头
seqkit stat temp/vs2/result/no_duplicates/final-viral-combined.fa #70818 1.05Gb
#--keep-original-seq保留原始序列
#--seqname-suffix-off保留原始名称
seqkit stat temp/vs2/result/no_duplicates/final-viral-combined.fa #70818 1.05Gb
#num_seqs        sum_len    min_len   avg_len   max_len
#70,818     1,132,179,453      711    15,987.2  508,340

#High confidence :  score >= 0.9  hallmark gene >=1  
awk 'NR > 1 && ($4 >= 0.9 && $7 >= 1 && $6 >= 2500)' temp/vs2/result/no_duplicates/final-viral-score.tsv > result/vs2/vs2_high.tsv
wc -l result/vs2/vs2_high.tsv  #带表头30770 ，里面有重复，我们放到合并步骤进行了去重
#提取High confidence病毒
awk '{print $1}' result/vs2/vs2_high.tsv | seqkit grep -f - temp/vs2/result/no_duplicates/final-viral-combined.fa > temp/vs2/result/no_duplicates/vs2_high.fna 
seqkit stat temp/vs2/result/no_duplicates/vs2_high.fna #30799
#num_seqs      sum_len    min_len   avg_len   max_len
#30,799      667,381,462    1,394  21,668.9  483,234

#Moderate confidence : 0.9 > score >= 0.5  hallmark gene >=1  
awk 'NR==1 || ($4 < 0.9 && $4 >= 0.5 && $7 >=1 && $6 >= 2500)' temp/vs2/result/no_duplicates/final-viral-score.tsv > result/vs2/vs2_moderate.tsv
wc -l result/vs2/vs2_moderate.tsv  #带表头6670 ，里面有重复，我们放到合并步骤进行了去重
#提取Moderate confidence病毒
awk '{print $1}' result/vs2/vs2_moderate.tsv | seqkit grep -f - temp/vs2/result/no_duplicates/final-viral-combined.fa > temp/vs2/result/no_duplicates/vs2_moderate.fna
seqkit stat temp/vs2/result/no_duplicates/vs2_moderate.fna#6690
#num_seqs      sum_len    min_len   avg_len   max_len
#6,690     179,472,186    1,394    26,826.9  508,340

2.7 ViralVerify （277h）
https://github.com/ablab/viralVerify
#使用
cd ~/virus2
conda activate viralverify
mkdir -p result/vv
temp/viralverify/bin/viralverify -f ~/virus2/seq/no_duplicates.fna -o temp/viralverify/output -t 8 --hmm ~/db/Pfam_A/Pfam-A.hmm 

seqkit stat temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus.fasta #35249 387Mb
seqkit stat temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus_uncertain.fasta #609749 2.31Gb
#两个数据库都行，pfam更大、更新一些，先用这个
#~/db/viralverify/nbc_hmms.hmm 
#~/db/Pfam_A/Pfam-A.hmm
#提取结果
cut -d',' -f1-3 temp/viralverify/output/no_duplicates_result_table.csv > temp/viralverify/output/vv.csv

#High confidence : classified as virus
awk -F',' 'BEGIN {print "Contig name,Prediction,Length"} $2 == "Virus" && $3 >= 2500' temp/viralverify/output/vv.csv > result/vv/vv_high.tsv
#更换分隔符为\t
awk '{gsub(/,/, "\t"); print}' result/vv/vv_high.tsv > temp/vv_high_tmp.tsv && mv temp/vv_high_tmp.tsv result/vv/vv_high.tsv




wc -l result/vv/vv_high.tsv  #带表头29791
seqkit stat temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus.fasta #35249 387Mb
#num_seqs      sum_len  min_len   avg_len  max_len
#35,249  405,087,943    1,000  11,492.2  256,155 

#Moderate confidence : classified as uncertain virus
awk -F',' 'BEGIN {print "Contig name,Prediction,Length"} $2 == "Uncertain - viral or bacterial" && $3 >= 2500' temp/viralverify/output/vv.csv > result/vv/vv_moderate.tsv
#更换分隔符为\t
awk '{gsub(/,/, "\t"); print}' result/vv/vv_moderate.tsv > temp/vv_moderate.tsv && mv temp/vv_moderate.tsv result/vv/vv_moderate.tsv


wc -l result/vv/vv_moderate.tsv  #带表头489758
seqkit stat temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus_uncertain.fasta#609749 2.31Gb
#num_seqs      sum_len    min_len   avg_len  max_len
#609,749    2,466,063,472    1,000  4,044.4  210,659 
#删除大文件
rm temp/viralverify/output/*.fa temp/viralverify/output/*.fasta temp/viralverify/output/no_duplicates_domtblout temp/viralverify/output/no_duplicates_out_pfam

seqkit stat temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus.fasta
seqkit stat temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus_uncertain.fasta
3.合并六种软件结果（高质量、中等质量）
#合并所有鉴定软件的序列fna文件
#合并所有软件的鉴定结果fna文件
cd ~/virus2
cat temp/genomad_high/no_duplicates_summary/gn_high.fna temp/genomad_moderate/no_duplicates_summary/gn_moderate.fna temp/VIBRANT/result/VIBRANT_no_duplicates/VIBRANT_phages_no_duplicates/no_duplicates.phages_combined.fna temp/MGV/mgv/viral_detection_pipeline/output/vf_high.fna temp/MGV/mgv/viral_detection_pipeline/output/vf_moderate.fna temp/dvf/result/dvf_high.fna temp/dvf/result/dvf_moderate.fna temp/vs2/result/no_duplicates/vs2_high.fna temp/vs2/result/no_duplicates/vs2_moderate.fna temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus.fasta temp/viralverify/output/Prediction_results_fasta/no_duplicates_virus_uncertain.fasta > seq/six.fna
seqkit stat seq/six.fna #1076029  9.69Gb
# num_seqs  sum_len       min_len   avg_len    max_len
#1,076,029  10,255,717,944    1,000  9,531.1  1,526,733
#去除length小于5000的，方法1
seqkit seq -m 5000 -g seq/six.fna > seq/six_5000.fna
seqkit stat seq/six_5000.fna #568,854 8.21Gb
#num_seqs        sum_len     min_len   avg_len    max_len
# 568,854     8,649,910,590    5,000  15,205.9  1,526,733

#方法2
#awk '/^>/ {if (seq && length(seq) >= 5000) print header"\n"seq; header=$0; seq=""; next} {seq = seq $0} END {if (length(seq) >= 5000) print header"\n"seq}' seq/six.fna > seq/six_5000.fna
#seqkit stat seq/six_5000.fna # 8.07Gb
#num_seqs        sum_len     min_len   avg_len    max_len
# 568,854     8,649,910,590    5,000  15,205.9  1,526,733

#去重           
~/db/scripts/clean.py seq/six_5000.fna seq/gnvbvfdvfvs2vv.fna
seqkit stat seq/gnvbvfdvfvs2vv.fna #309231 , 4.73Gb
#num_seqs        sum_len    min_len   avg_len    max_len
#309,231     5,068,022,810    5,000  16,389.1  1,526,733

rm -rf seq/six.fna
rm -rf seq/six_5000.fna
3.1 # High confidence viruses（两个以上软件鉴定出来的高质量）（先三个取交集，再取并集）
#使用7种鉴定方法（gn\vb\vf\dvf\vs\vs2\vv），三三取交集，共35种交集组合
#使用6种鉴定方法（gn\vb\vf\dvf\vs2\vv），三三取交集，共20种交集组合
#py脚本在/data4/machuang/db/scripts/Find.py，先复制到自己的scripts文件夹内
#脚本运行,需要环境中有pandas
conda activate vr
cd ~/virus2
python ~/db/scripts/Find.py -i result/genomad/gn_high.tsv result/VIBRANT/vb_high.tsv result/vf/vf_high.tsv result/dvf/dvf_high.tsv result/vs2/vs2_high.tsv result/vv/vv_high.tsv -o temp/high.tsv
wc -l temp/high.tsv # 23494

#提取高质量序列，从病毒鉴定软件的输出结果中提取
seqkit grep -f temp/high.tsv ~/virus2/seq/gnvbvfdvfvs2vv.fna > temp/high.fna  
seqkit stat temp/high.fna #23494 390Mb
# num_seqs        sum_len     min_len    avg_len   max_len
# 23,494       402,003,172    5,000  17,110.9  412,909

#脚本正确与否，验证过程，使用4种鉴定方法（gn\vb\vf\dvf），三三取交集，共4种交集组合，得到的最终结果一致，所以脚本正确
#方法1：脚本得到 #10634
#cd ~/virus2
#conda activate vr
#python ~/db/scripts/Find.py -i result/genomad/gn_high.tsv result/VIBRANT/vb_high.tsv result/vf/vf_high.tsv result/dvf/dvf_high.tsv -o temp/high.txt
#wc -l temp/high.txt #10634

#方法2 ：Linux命令得到 #10634
#cd ~/virus2
#取gnvb交集

 
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/genomad/gn_high.tsv result/VIBRANT/vb_high.tsv > temp/gnvb_high.txt

#wc -l temp/gnvb_high.txt  #23199
#将gnvb与vf再取交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} $1 in a {print $1}' result/vf/vf_high.tsv temp/gnvb_high.txt > temp/gnvbvf_high.txt
#wc -l temp/gnvbvf_high.txt #910

#取gnvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/genomad/gn_high.tsv result/vf/vf_high.tsv > temp/gnvf_high.txt


#wc -l temp/gnvf_high.txt  # 1029
#将gnvf与dvf再取交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} $1 in a {print $1}' result/dvf/dvf_high.tsv temp/gnvf_high.txt > temp/gnvfdvf_high.txt
#wc -l temp/gnvfdvf_high.txt   #893

#取vbdvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/VIBRANT/vb_high.tsv result/dvf/dvf_high.tsv > temp/vbdvf_high.txt

#wc -l temp/vbdvf_high.txt  #14879
#将vbdvf与gn再取交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} $1 in a {print $1}' result/genomad/gn_high.tsv temp/vbdvf_high.txt > temp/gnvbdvf_high.txt
#wc -l temp/gnvbdvf_high.txt   #10102

#取vbvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/VIBRANT/vb_high.tsv result/vf/vf_high.tsv > temp/vbvf_high.txt

#wc -l temp/vbvf_high.txt  #1302
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} $1 in a {print $1}' result/dvf/dvf_high.tsv temp/vbvf_high.txt > temp/vbvfdvf_high.txt
#wc -l temp/vbvfdvf_high.txt   #1129

#取并集，去除重复项
#cat temp/gnvbvf_high.txt temp/gnvfdvf_high.txt temp/gnvbdvf_high.txt temp/vbvfdvf_high.txt  | sort | uniq > temp/high.txt 
#wc -l temp/high.txt  #10634


#脚本正确与否，验证过程，使用4种鉴定方法（gn\vb\vf\dvf），两两取交集，共6种交集组合，得到的最终结果一致，所以脚本正确
#方法1：脚本得到 #20950
#cd ~/virus2
#conda activate vr
#python ~/db/scripts/find.py result/genomad/gn_high.tsv result/VIBRANT/vb_high.tsv result/vf/vf_high.tsv result/dvf/dvf_high.tsv --output temp/high.tsv
#wc -l temp/high.tsv #20950

#方法2 ：Linux命令得到 #20950
#cd ~/virus2

#取gnvb交集

 
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/genomad/gn_high.tsv result/VIBRANT/vb_high.tsv > temp/gnvb_high.txt

#wc -l temp/gnvb_high.txt  #2067

#取gnvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/genomad/gn_high.tsv result/vf/vf_high.tsv > temp/gnvf_high.txt


#wc -l temp/gnvf_high.txt  # 1038  

#取gndvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/genomad/gn_high.tsv result/dvf/dvf_high.tsv > temp/gndvf_high.txt


#wc -l temp/gndvf_high.txt  #18926

#取vbdvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/VIBRANT/vb_high.tsv result/dvf/dvf_high.tsv > temp/vbdvf_high.txt

#wc -l temp/vbdvf_high.txt  #652

#取vbvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/VIBRANT/vb_high.tsv result/vf/vf_high.tsv > temp/vbvf_high.txt

#wc -l temp/vbvf_high.txt  #69

#取vfdvf交集
#awk -F'\t' 'NR==FNR {if (NR>1) a[$1]; next} FNR>1 && $1 in a {print $1}' result/dvf/dvf_high.tsv result/vf/vf_high.tsv > temp/vfdvf_high.txt


#wc -l temp/vfdvf_high.txt  #2471

#取并集，去除重复项
#cat temp/gnvb_high.txt temp/gnvf_high.txt temp/gndvf_high.txt temp/vbdvf_high.txt temp/vbvf_high.txt temp/vfdvf_high.txt | sort | uniq > temp/High.txt 
#wc -l temp/High.txt #20950
3.2 # Complete viruses（需同时满足高质量和CheckV鉴定出的完整基因）
https://bitbucket.org/berkeleylab/checkv/src/master/
#运行checkv，记录版本号
conda activate checkv
checkv -h #v1.0.3
cd ~/virus2
#设置数据库，运行checkv
export CHECKVDB=/data4/machuang/db/viwrap/CheckV_db
#重复运行不同数据，不会覆盖原结果，若不删除旧文件则会出现key_error
#输入为高质量基因组，23494
checkv end_to_end temp/high.fna temp/checkv_complete -t 10

#挑出完整基因组
awk -F'\t' 'NR==1 {print; next} $10 ~ /^(100.0)$/ {print}' temp/checkv_complete/quality_summary.tsv > temp/checkv_complete/complete.tsv
wc -l temp/checkv_complete/complete.tsv #带表头1252
seqkit stat temp/checkv_complete/viruses.fna #22462 358Mb


rm -rf temp/checkv_complete/tmp
3.3 #  Moderate confidence viruses（两个及以上软件鉴定出的中等质量  and  只有一个软件鉴定出的高质量）
#脚本运行，两个软件以上鉴定得到的中等质量
cd ~/virus2
conda activate vr
python ~/db/scripts/find.py result/genomad/gn_moderate.tsv result/vf/vf_moderate.tsv result/dvf/dvf_moderate.tsv result/vs2/vs2_moderate.tsv result/vv/vv_moderate.tsv --output temp/moderate1.tsv
wc -l temp/moderate1.tsv #71074

#只有一个软件鉴定得到的高质量，也视为中等质量
#取并集去重
awk -F'\t' 'FNR>1 {print $1}' result/genomad/gn_high.tsv result/VIBRANT/vb_high.tsv result/vf/vf_high.tsv result/dvf/dvf_high.tsv result/vs2/vs2_high.tsv result/vv/vv_high.tsv | sort | uniq > temp/six_u.tsv
wc -l temp/six_u.tsv # 122817 
#找单个软件鉴定出是高质量，我们认为是中等质量
comm -23 <(sort temp/six_u.tsv) <(sort temp/high.tsv) > temp/moderate2.tsv
wc -l temp/moderate2.tsv # 99323
#应满足 moderate2.tsv + high.tsv (23494)= six_u.tsv

#合并两种情况的moderate
cat temp/moderate1.tsv temp/moderate2.tsv | sort | uniq > temp/moderate.tsv
wc -l temp/moderate.tsv # 147186 

#提取
seqkit grep -f temp/moderate.tsv  ~/virus2/seq/gnvbvfdvfvs2vv.fna > temp/moderate.fna
seqkit stat temp/moderate.fna #146,958 2.22Gb
#num_seqs        sum_len     min_len   avg_len    max_len
#146,958     2,340,310,906    5,000     15,925   1,526,733

3.4 # Low confidence viruses （只有一个软件鉴定出的中等质量）
#只有一个软件鉴定得到的中等质量，被视为低质量
#取并集去重
awk -F'\t' 'FNR>1 {print $1}' result/genomad/gn_moderate.tsv result/vf/vf_moderate.tsv result/dvf/dvf_moderate.tsv result/vs2/vs2_moderate.tsv result/vv/vv_moderate.tsv | sort | uniq > temp/low_u.tsv
wc -l temp/low_u.tsv # 265503
#找单个软件鉴定出是中等质量，我们认为是低质量
comm -23 <(sort temp/low_u.tsv) <(sort temp/moderate.tsv) > temp/low.tsv
wc -l temp/low.tsv # 162256
3.5 # 合并High 和 Moderate
cd ~/virus2
#取并集，high.tsv+moderate.tsv=170452,去除重复项
cat temp/high.tsv temp/moderate.tsv | sort | uniq > temp/high_moderate.tsv
wc -l temp/high_moderate.tsv # 156852
#提取序列
seqkit grep -f temp/high_moderate.tsv  ~/virus2/seq/gnvbvfdvfvs2vv.fna > temp/high_moderate.fna
seqkit stat temp/high_moderate.fna # 156624 2.4Gb
#因为有些length小于5000，虽然命令运行过程中设置了长度，但是结果还是有一些length小于5000，导致156624略低于156852
# num_seqs        sum_len      min_len      avg_len    max_len
#156624        2,535,697,220    5,000      16,189.7  1,526,733

#rm temp/high_moderate.fna
4.CheckV (6h)
https://bitbucket.org/berkeleylab/checkv/src/master/
#按照下列要求筛选
#(1)CheckV用AAI模型确定的至少MQ的votu
#(2)LQ病毒中,基因数≥10、基因含量≥40%的病毒基因和<10%的宿主基因
#(3)排除HQMQ中HMM预测的

#运行checkv，记录版本号
conda activate checkv
checkv -h #v1.0.3
cd ~/virus2
#设置数据库
export CHECKVDB=/data4/machuang/db/viwrap/CheckV_db
#重复运行不同数据，不会覆盖原结果，若不删除旧文件则会出现key_error
checkv end_to_end ~/virus2/temp/high_moderate.fna temp/checkv -t 10 
seqkit stat temp/checkv/viruses.fna  #144445 1.94Gb
# num_seqs        sum_len      min_len      avg_len    max_len
#144,445       2,078,276,335    5,000        14,388  1,526,733
seqkit stat temp/checkv/proviruses.fna  #12226  224Mb
# num_seqs        sum_len      min_len      avg_len    max_len
#12,226          235,184,300      159       19,236.4  260,516
wc -l temp/checkv/quality_summary.tsv #带表头156625

#对原噬菌体名称进行修改
#cp temp/checkv/proviruses.fna temp/checkv/prov.fna #先复制
#awk '/^>k141/{f=1} f; /^>/ && f{f=0}' temp/checkv/prov.fna > temp/checkv/prov.txt
#sed 's/_1 [0-9]*-[0-9]*\/[0-9]*//' temp/checkv/prov.fna > temp/checkv/prov_cleaned.fna
#awk '/^>k141/{f=1} f; /^>/ && f{f=0}' temp/checkv/prov_cleaned.fna > temp/checkv/prov1.txt

#对结果进行筛选，除去HQ、MQ中HMM预测的, LQ需满足 (geme_count≥10; virus_gene≥40%; host_gene＜10%)
# HQMQ.txt ：除去HMM预测的  
awk -F'\t' 'NR==1 {print; next} $8 ~ /^(Complete|High-quality|Medium-quality)$/ && $11 != "HMM-based (lower-bound)" && $2 >= 5000 {print}' temp/checkv/quality_summary.tsv > temp/checkv/HQMQ.txt

wc -l temp/checkv/HQMQ.txt #12301

# LQ.txt ：geme_count≥10; virus_gene≥40%; host_gene＜10% 
#v_genes=$6; h_genes=$7; g_count=$5
awk -F'\t' 'NR==1 {print > "temp/checkv/LQ.txt"; next} \
            $8 ~ /^(Low-quality)$/ && $5 >= 10 && $2 >=5000 && $6/$5 >= 0.4 && $7/$5 < 0.1 {print > "temp/checkv/LQ.txt"}' temp/checkv/quality_summary.tsv

wc -l temp/checkv/LQ.txt  #6891

# HQMQLQ.txt ： 合并 
{ cat temp/checkv/HQMQ.txt; tail -n +2 temp/checkv/LQ.txt; } > temp/checkv/HQMQLQ.txt

wc -l temp/checkv/HQMQLQ.txt #检查行数19191

#根据结果筛选序列文件,checkv会对proviruses序列进行切割输出，并添加额外信息，对viruses病毒则不会，所以我们提取病毒序列可以直接从temp/high_moderate.fna进行提取，如果直接从checkv的输出结果temp/checkv/viruses.fna提取，则会少掉proviruses的序列
#先提取出contig_id文件
awk '{print $1}' temp/checkv/HQMQLQ.txt > temp/checkv/contig_id.txt
wc -l temp/checkv/contig_id.txt #检查行数带表头，19191
seqkit grep -f temp/checkv/contig_id.txt temp/high_moderate.fna > temp/checkv/HQMQLQ.fna 
seqkit stat temp/checkv/HQMQLQ.fna #19190,625Mb
# num_seqs      sum_len      min_len   avg_len   max_len
# 19,190       643,858,952    5,000   33,551.8   634,451

#判断从哪个文件中提取
#先提取出contig_id文件
#awk '{print $1}' temp/checkv/HQMQLQ.txt > temp/checkv/contig_id.txt
#wc -l temp/checkv/contig_id.txt #检查行数带表头，13497
#seqkit grep -f temp/checkv/contig_id.txt temp/checkv/viruses.fna > temp/checkv/HQMQLQ.fna 
#seqkit stat temp/checkv/HQMQLQ.fna #13198 257MB

#seqkit grep -f temp/checkv/contig_id.txt seq/no_duplicates.fna > temp/checkv/HQMQLQ_old.fna
#seqkit stat temp/checkv/HQMQLQ_old.fna  #13496 433MB
#保存最终结果到result/
cp temp/checkv/HQMQLQ.txt result/
cp temp/checkv/HQMQLQ.fna result/


#删除中间文件
cd ~/virus2
rm -rf temp/checkv/tmp
5.去冗余 
1. 用cd-hit 去冗余
####cd-hit去冗余#### 10/25/14:00-
cd ~/virus2
mkdir -p temp/cdhit
export PATH="$HOME/db/software/cdhit/:$PATH"
#输入19190
#95% identity and 90% coverage
cd-hit-est -i result/HQMQLQ.fna \
        -o temp/cdhit/HQMQLQ_dereplicated.fna \
        -aS 0.90 -c 0.95 -G 0 -g 1 -T 16 -M 0 
        #-d 0 -n 10 默认的
seqkit stat temp/cdhit/HQMQLQ_dereplicated.fna #8829 328Mb
# num_seqs      sum_len     min_len   avg_len   max_len
#   8,829     338,274,892    5,013  38,314.1  634,451
2. 用 blast聚类,去冗余 
https://github.com/snayfach/MGV/tree/master/ani_cluster
mkdir -p ~/virus2/temp/MGV/mgv/ani_cluster/REP
cd ~/virus2/temp/MGV/mgv/ani_cluster/REP
cp ../blastani.py ./
cp ../cluster.py ./
conda activate blast

#用自身~/virus2/result/ref_votu.fna建库，并all vs all比对
#合并后使用all vs all比对
makeblastdb -in ~/virus2/result/HQMQLQ.fna  -out blastdb -dbtype nucl
blastn -query ~/virus2/result/HQMQLQ.fna  -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90

#Compute ANI from BLAST results
python blastani.py -i blast.tsv -o ani.tsv

#Perform centroid-based clustering    聚类去重复（去冗余）
python cluster.py --fna ~/virus2/result/HQMQLQ.fna --ani ani.tsv --out clusters.tsv --min_ani 95 --min_qcov 0 --min_tcov 85
  
#其中符合要求的contig名字就在clusters.tsv文件中第一列
cut -f1 clusters.tsv > ~/virus2/result/vOTUs.list
wc -l ~/virus2/result/vOTUs.list # 7153
#从 FASTA 中提取有代表性的基因组
seqkit grep -f ~/virus2/result/vOTUs.list ~/virus2/result/HQMQLQ.fna > ~/virus2/result/vOTUs.fna  
seqkit stat ~/virus2/result/vOTUs.fna  #7153 263Mb
# num_seqs      sum_len     min_len   avg_len   max_len
#  7153       271,934,623    5,016    38,016.9  634,451


#蛋白预测
cd ~/virus2/result
export PATH="$HOME/db/software/Prodigal/:$PATH"
#7153个vOTU,生成蛋白质文件final-viral-combined.faa
prodigal -i vOTUs.fna -o ./vOTUs.genes -a ./vOTUs.faa -p meta
seqkit stat vOTUs.faa

#根据non-redundant gene catalogue 获取  HQMQLQ.genecat.m6
#04_viral_postprocessing.wdl中MGX.genecat.fna就是 = temp/cdhit/HQMQLQ_dereplicated.fna
mkdir -p ~/virus2/temp/MGV/mgv/ani_cluster/genecat
cd ~/virus2
conda activate blast
#用得到的18898个vOTUs，建库
cat ~/virus2/result/vOTUs.fna > temp/viruses.fna
makeblastdb -in temp/viruses.fna -dbtype nucl -out temp/MGV/mgv/ani_cluster/genecat/HQMQLQdbreps
#将非冗余基因目录temp/cdhit/HQMQLQ_dereplicated.fna与建的库比对，生成HQMQLQ.genecat.m6，用于后续median_TPM分析
blastn -task megablast -db temp/MGV/mgv/ani_cluster/genecat/HQMQLQdbreps -query temp/cdhit/HQMQLQ_dereplicated.fna -out temp/MGV/mgv/ani_cluster/genecat/HQMQLQ.genecat.m6 -outfmt 6 -num_threads 16 -evalue 0.01 -qcov_hsp_perc 80 -perc_identity 90 -word_size 20 -reward 2 -penalty -3 -max_target_seqs 10 -max_hsps 1

#用18898和MGV合并的vOTUs建库
# mkdir -p ~/virus2/temp/MGV/mgv/ani_cluster/genecat_merge
# cat ~/virus2/result/vOTUs.fna ~/db/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.fna > temp/viruses_merge.fna
# makeblastdb -in temp/viruses_merge.fna -dbtype nucl -out temp/MGV/mgv/ani_cluster/genecat_merge/HQMQLQdbreps
#将非冗余基因目录temp/cdhit/HQMQLQ_dereplicated.fna与建的库比对，生成HQMQLQ.genecat.m6，用于后续median_TPM分析
# blastn -task megablast -db temp/MGV/mgv/ani_cluster/genecat_merge/HQMQLQdbreps -query temp/cdhit/HQMQLQ_dereplicated.fna -out temp/MGV/mgv/ani_cluster/genecat_merge/HQMQLQ.genecat_merge.m6 -outfmt 6 -num_threads 16 -evalue 0.01 -qcov_hsp_perc 80 -perc_identity 90 -word_size 20 -reward 2 -penalty -3 -max_target_seqs 10 -max_hsps 1

6.用VirRep验证得到的病毒序列
https://github.com/Dongyq815/VirRep
#用最终votu数据运行VirRep
conda activate vr
cd ~/virus2
mkdir -p temp/VirRep/output
#-n 是将fna分为多少片段
python temp/VirRep/src/VirRep.py -i ~/virus2/result/vOTUs.fna -o temp/VirRep/output -l 5000 -n 100 -m temp/VirRep/src/model/VirRep.pth --conservative  --use-amp --cpu --cpu-threads 8
