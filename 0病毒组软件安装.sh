
## 0病毒组软件安装

#新建文件夹
mkdir -p /data4/machuang/virus1
cd /data4/machuang/virus1
mkdir -p seq temp result
# result  储存重要结果，及metadata文件
# seq     存放宏基因组分析得到的组装结果
# temp    存放临时文件
1.安装checkv
https://bitbucket.org/berkeleylab/checkv/src/master/
#conda安装
conda create -n checkv -c conda-forge -c bioconda checkv
#数据库下载
wget -c https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
#解压数据库到指定文件夹
tar -zxvf checkv-db-v1.5.tar.gz -C ~/db
#最新版数据库删除了DIAMOND以减少空间，缺少checkv_reps.dmnd文件，因此需要下载（https://portal.nersc.gov/CheckV/）老版本数据库（v1.2及之后版本都没有该文件）后本地构建。
#我们下载v1.1版本 https://portal.nersc.gov/CheckV/checkv-db-v1.1.tar.gz下载解压后，将checkv_reps.dmnd文件放入~/db/checkv-db-v1.5/genome_db文件夹中
wget -c https://portal.nersc.gov/CheckV/checkv-db-v1.1.tar.gz
tar -zxvf checkv-db-v1.1.tar.gz -C ~/db
cp ~/db/checkv-db-v1.1/genome_db/checkv_reps.dmnd ~/db/checkv-db-v1.5/genome_db/

#进入数据库目录
   cd ~/db/viwrap/CheckV_db/genome_db
 #手动建立数据库，不然会报错DIAMOND task failed. Program should be rerun.查看diamond.log文件The log says “Error: Database was built with a different version of diamond as is incompatible.”
 diamond makedb --in checkv_reps.faa --db checkv_reps
#最后我们设置一下数据库位置就可以使用了
export CHECKVDB=~/db/viwrap/CheckV_db

 #输入测试文件
 checkv end_to_end /data4/machuang/virus/vs2_result2/final-viral-combined.fa /data4/machuang/virus/vs2_result2/checkv -t 4
2.安装MGV、CD-HIT、meghit
2.1MGV
https://github.com/snayfach/MGV
######MGV安装
conda create -n blast -c bioconda blast pip wget unzip
conda activate blast
pip install biopython
mkdir -p /data4/machuang/virus1/temp/MGV
cd /data4/machuang/virus1/temp/MGV
#下载MGV聚类所需脚本
wget -c https://codeload.github.com/snayfach/MGV/zip/refs/heads/master
unzip master
mv MGV-master MGV
chmod -R 755 MGV


2.2    CD-HIT
https://github.com/weizhongli/cdhit
######CD-HIT安装
cd ~/virus1/temp
wget -c https://codeload.github.com/weizhongli/cdhit/zip/refs/heads/master
unzip master
mv cdhit-master cdhit
chmod -R 755 cdhit
cd cdhit
make
export PATH="$HOME/virus1/temp/cdhit/:$PATH"

#使用
#cd ~/virus2
#mkdir -p temp/cdhit result/cdhit
#export PATH="$HOME/virus1/temp/cdhit/:$PATH"
#cd-hit-est -i result/HQMQLQ.fna -o temp/cdhit/nucleotide.fna -c 0.99 -g 1 -M 0 -d 0 -n 10 
#8963
2.3 MEGAHIT
https://github.com/voutcn/megahit
conda creat -n megahit -c bioconda megahit
conda activate megahit
3.安装Prodigal、hmmsearch
3.1 Prodigal
https://github.com/hyattpd/Prodigal
（推荐）方法1：#安装Prodigal
mkdir -p ~/virus1/temp/Prodigal
cd ~/virus1/temp/Prodigal
wget -c https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux
mv prodigal.linux prodigal
chmod +x prodigal
export PATH="$HOME/virus1/temp/Prodigal/:$PATH"
prodigal -h

Input/Output Parameters

  -i, --input_file:     Specify input file (default stdin).
  -o, --output_file:    Specify output file (default stdout).
  -a, --protein_file:   Specify protein translations file.
  -d, --mrna_file:      Specify nucleotide sequences file.
  -s, --start_file:     Specify complete starts file.
  -w, --summ_file:      Specify summary statistics file.
  -f, --output_format:  Specify output format.
                          gbk:  Genbank-like format (Default)
                          gff:  GFF format
                          sqn:  Sequin feature table format
                          sco:  Simple coordinate output
  -q, --quiet:          Run quietly (suppress logging output).

#方法2 ##病毒蛋白预测安装prodigal
#cd ~/virus1/temp
#先下载https://codeload.github.com/hyattpd/Prodigal/zip/refs/heads/GoogleImport压缩包，解压并修改权限
#wget -c https://codeload.github.com/hyattpd/Prodigal/zip/refs/heads/GoogleImport
#unzip GoogleImport
#mv Prodigal-GoogleImport prodigal
#rm -rf GoogleImport
#chmod -R 755 prodigal
#cd prodigal
#make install INSTALLDIR=/data4/machuang/virus1/temp/Prodigal/
#需要改成自己的路径
#export PATH="$HOME/virus1/temp/Prodigal/:$PATH"
#prodigal -h  #确定一下是否安装成功
#cd ..
#rm -rf prodigal
3.2 hmmsearch
https://github.com/EddyRivasLab/hmmer
#安装hmmsearch
conda create -n hmmer -c biocore -c bioconda -c conda-forge r-base=4.3.1 hmmer diamond bioconda::mcl pandas
conda activate hmmer
hmmsearch -h#确认是否安装成功


##Download & decompress reference databases of HMMs
#wget -O input/imgvr.hmm.gz https://img.jgi.doe.gov//docs/final_list.hmms.gz
#wget -O input/pfam.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz 
#gunzip input/imgvr.hmm.gz
#gunzip input/pfam.hmm.gz


4.病毒鉴定

Cenote-Taker3 https://github.com/mtisza1/Cenote-Taker3
4.1 VIBRANT
https://github.com/EddyRivasLab/hmmer
######安装VIBRANT
conda create -n VIBRANT -c conda-forge -c bioconda vibrant==1.2.1
conda activate VIBRANT
#下载脚本，设置数据库
cd ~/virus1/temp
git clone https://github.com/AnantharamanLab/VIBRANT
chmod -R 755 VIBRANT
cd VIBRANT/databases
#代码下载数据库运行 ，数据库会下载到~/virus1/temp/VIBRANT/databases
 ./VIBRANT_setup.py
#下载完成的运行
./VIBRANT_setup.py -test
#安装完成提示VIBRANT v1.2.1 is good to go!
#See example_data/ for quick test files.
4.2 VirSorter、VirSorter2
4.2.1 VirSorter
https://github.com/simroux/VirSorter
#下载数据库
mkdir -p ~/db/virsorter
cd ~/db/virsorter
wget -c https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
tar -xvzf virsorter-data-v2.tar.gz
#安装
conda create --name virsorter -y -c bioconda mcl=14.137 muscle blast pandas perl-bioperl perl-file-which hmmer=3.1b2 perl-parallel-forkmanager perl-list-moreutils diamond=0.9.14 metagene_annotator
conda activate virsorter
cd ~/virus2/temp
git clone https://github.com/simroux/VirSorter.git
mv VirSorter virsorter
chmod -R 755 virsorter
cd virsorter/Scripts
make clean
make

#安装所需脚本
cd ~/miniconda3/envs/virsorter/bin
wget -c http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
tar -xvzf mga_x86_64.tar.gz

#运行，(categories 1, 2, 4, 5)
conda activate virsorter
cd ~/virus2
temp/virsorter/wrapper_phage_contigs_sorter_iPlant.pl -f ~/virus2/seq/no_duplicates.fna --db 1 --wdir temp/virsorter --ncpu 16 --data-dir ~/db/virsorter/virsorter-data
4.2.2 VirSorter2
https://github.com/jiarong/VirSorter2
######安装VirSorter2
conda create -n vs2 -c conda-forge -c bioconda virsorter=2
conda activate vs2
#我选择手动下载数据库https://osf.io/v46sc/download ；
#然后解压
tar -xzf db.tgz -C /data/machuang/db/viwrap/
mv db VirSorter2_db
chmod -R 755 VirSorter2_db
#手动下载依赖项，并解压手动将依赖项中virsorter文件夹中文件放入db内
cd ~/virus1/temp
wget -c https://codeload.github.com/jiarong/VirSorter2/zip/refs/heads/master
unzip master
mv VirSorter2-master vs2
chmod -R 755 vs2
rm -rf master
cp -R ~/virus1/temp/vs2/virsorter/* ~/db/viwrap/VirSorter2_db/
chmod -R 755 ~/db/viwrap/VirSorter2_db
#设置环境
virsorter config --init-source --db-dir=/data4/machuang/db/viwrap/VirSorter2_db  #VirSorter 2.2.4
4.3 Virfinder、Deepvirfinder
4.3.1 Virfinder
conda activate hmmer
mkdir -p /data4/machuang/R/x86_64-pc-linux-gnu-library/4.3
#进入Linux的R版本
R
#可以先用命令行下载，若报错再手动下载，virfinder最好去github上下载，手动安装，注意所有依赖装完才能装virfinder
install.packages('shape')
install.packages('iterators')
install.packages('foreach')
install.packages('Rcpp')
install.packages("lattice")
#install.packages("matrix")
install.packages("matrix", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

install.packages("survival")
install.packages('RcppEigen')
install.packages('qvalue')
install.packages('glmnet')


install.packages("codetools", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("shape", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("survival", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
install.packages("devtools", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/RcppEigen_0.3.4.0.0.tar.gz", repos = NULL, type="source")
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/Rcpp_1.0.12.tar.gz", repos = NULL, type="source")
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/iterators_1.0.14.tar.gz", repos = NULL, type="source")
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/foreach_1.5.2.tar.gz", repos = NULL, type="source")
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/glmnet_4.1-8.tar.gz", repos = NULL, type="source")
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/qvalue.tar.gz", repos = NULL, type="source")
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/VirFinder_1.1.tar.gz", repos = NULL, type="source")
#依赖包Rcpp,需手动下载，然后放入~/virus1/temp/MGV/mgv/viral_detection_pipeline  下载地址https://cran.rstudio.com/bin/windows/contrib/4.3/Rcpp_1.0.12.zip
unzip Rcpp_1.0.12.zip -d Rcpp_1.0.12
mv Rcpp_1.0.12/Rcpp ./
rm -rf  Rcpp_1.0.12
tar -czvf Rcpp_1.0.12.tar.gz Rcpp/
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/Rcpp_1.0.12.tar.gz", repos = NULL, type="source")

#依赖包RcppEigen,需手动下载，然后放入~/virus1/temp/MGV/mgv/viral_detection_pipeline https://cran.rstudio.com/bin/windows/contrib/4.3/RcppEigen_0.3.4.0.0.zip
unzip RcppEigen_0.3.4.0.0.zip -d RcppEigen_0.3.4.0.0
mv RcppEigen_0.3.4.0.0/RcppEigen ./
rm -rf  RcppEigen_0.3.4.0.0
tar -czvf RcppEigen_0.3.4.0.0.tar.gz RcppEigen/
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/RcppEigen_0.3.4.0.0.tar.gz", repos = NULL, type="source")



#依赖包iterators,需手动下载，然后放入~/virus1/temp/MGV/mgv/viral_detection_pipeline  https://cran.rstudio.com/bin/windows/contrib/4.3/iterators_1.0.14.zip
unzip iterators_1.0.14.zip -d iterators_1.0.14
mv iterators_1.0.14/iterators ./
rm -rf  iterators_1.0.14
tar -czvf iterators_1.0.14.tar.gz iterators/
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/iterators_1.0.14.tar.gz", repos = NULL, type="source")

#依赖包foreach,需手动下载，然后放入~/virus1/temp/MGV/mgv/viral_detection_pipeline  下载地址https://cran.rstudio.com/bin/windows/contrib/4.3/foreach_1.5.2.zip
unzip foreach_1.5.2.zip -d foreach_1.5.2
mv foreach_1.5.2/foreach ./
rm -rf  foreach_1.5.2
tar -czvf foreach_1.5.2.tar.gz foreach/
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/foreach_1.5.2.tar.gz", repos = NULL, type="source")

##依赖包glmnet,需手动下载，然后放入~/virus1/temp/MGV/mgv/viral_detection_pipeline  下载地址https://cran.rstudio.com/bin/windows/contrib/4.3/glmnet_4.1-8.zip
unzip glmnet_4.1-8.zip -d glmnet_4.1-8
mv glmnet_4.1-8/glmnet ./
rm -rf glmnet_4.1-8
tar -czvf glmnet_4.1-8.tar.gz glmnet/
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/glmnet_4.1-8.tar.gz", repos = NULL, type="source")

#依赖包qvalue，R安装后在library文件夹内找到并压缩，然后再移入~/virus1/temp/MGV/mgv/viral_detection_pipeline
unzip qvalue.zip -d qvalue1
mv qvalue1/qvalue ./
rm -rf  qvalue1
tar -czvf qvalue.tar.gz qvalue/
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/qvalue.tar.gz", repos = NULL, type="source")


#安装virfinder,需手动下载，然后放入~/virus1/temp/MGV/mgv/viral_detection_pipeline  下载地址https://github.com/jessieren/VirFinder/blob/master/windows/VirFinder_1.1.zip
unzip VirFinder_1.1.zip -d VirFinder_1.1
mv VirFinder_1.1/VirFinder ./VirFinder1
rm -rf  VirFinder_1.1
tar -czvf VirFinder_1.1.tar.gz VirFinder/
install.packages("/data4/machuang/virus1/temp/MGV/mgv/viral_detection_pipeline/VirFinder_1.1.tar.gz", repos = NULL, type="source")
4.3.2 Deepvirfinder
https://github.com/jessieren/DeepVirFinder
######安装Deepvirfinder
conda create -n dvf python=3.6 numpy theano=1.0.3 keras=2.2.4 scikit-learn Biopython h5py=2.10.0
conda activate dvf
conda install mkl-service
cd ~/virus1/temp
git clone https://github.com/jessieren/DeepVirFinder
mv DeepVirFinder dvf

chmod -R 755 dvf
4.4  geNomad
 https://portal.nersc.gov/genomad/index.html
https://github.com/apcamargo/genomad
conda create -n genomad -c conda-forge -c bioconda genomad
conda activate genomad

cd ~/db
#安装数据库https://portal.nersc.gov/genomad/__data__/genomad_db_v1.7.tar.gz
#genomad download-database ~/db
#代码下载容易断开出错，直接下载
wget -c https://portal.nersc.gov/genomad/__data__/genomad_db_v1.7.tar.gz
tar -xzvf genomad_db_v1.7.tar.gz
chmod 755 -R genomad_db
#genomad的使用
genomad end-to-end -t 8 metagenome.fna genomad_output ~/db/genomad_db #进行完整的pipline运行

#个性化运行
genomad annotate metagenome.fna genomad_output ~/db/genomad_db
genomad find-proviruses metagenome.fna genomad_output ~/db/genomad_db
genomad marker-classification metagenome.fna genomad_output ~/db/genomad_db
genomad nn-classification metagenome.fna genomad_output
genomad aggregated-classification metagenome.fna genomad_output
genomad score-calibration metagenome.fna genomad_output
genomad summary metagenome.fna genomad_output
4.5  VirRep
https://github.com/Dongyq815/VirRep
conda create -n vr python=3.8
conda activate vr
conda install -c bioconda -c conda-forge biopython numpy pandas tqdm scipy scikit-learn packaging
#如果有NVIDIA GPU and supprots CUDA 11.8可以安装以下
conda install -c pytorch -c nvidia pytorch torchvision torchaudio pytorch-cuda=11.8 

#安装脚本
cd ~/virus2/temp
git clone https://github.com/Dongyq815/VirRep.git
chmod -R 755 VirRep

#运行,报错 AssertionError: Torch not compiled with CUDA enabled（添加命令--cpu）

2024/7/10 19:55 -2024/7/

#用原始组装数据运行VirRep
conda activate vr
cd ~/virus2
mkdir -p temp/VirRep/vr

python temp/VirRep/src/VirRep.py -i ~/virus2/seq/no_duplicates.fna -o temp/VirRep/vr -l 5000 -n 100 -m temp/VirRep/src/model/VirRep.pth --use-amp --conservative --cpu --cpu-threads 8
#--conservative：Apply conservative settings. This is equivalent to the following criteria: baseline 0.8; min_score 1000-8000:0.9,8001-10000:0.85,10001-Inf:0.8. (default: False)
#--use-amp：Use automatic mixed precision to accelerate computational effeciency. (default: False)
#--provirus-off：Skip the iterative segment extension procedure (default: False)

#用最终votu数据运行VirRep
conda activate vr
cd ~/virus2
mkdir -p temp/VirRep/VR

python temp/VirRep/src/VirRep.py -i ~/virus2/result/vOTUs.fna -o temp/VirRep/VR -l 5000 -n 100 -m temp/VirRep/src/model/VirRep.pth --conservative  --use-amp --cpu --cpu-threads 8

4.6  ViralVerify
https://github.com/ablab/viralVerify
#viralVerify is a Python script, thus, installation is not required. However, it has the following dependencies:
#Python 3.6+,
#Prodigal (https://github.com/hyattpd/Prodigal
#hmmsearch (from the hmmer package, http://hmmer.org/download.html
#provided decompressed database of virus/chromosome-specific HMMs (https://figshare.com/s/f897d463b31a35ad7bf0
#or
#recent release of the Pfam-A database (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/).

#下载数据库
mkdir -p ~/db/viralverify ~/db/Pfam_A
cd ~/db/viralverify
wget -c https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/17904323/nbc_hmms.hmm.gz
gunzip nbc_hmms.hmm.gz
cd ~/db/Pfam_A
wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
#conda 安装
conda create -n viralverify -y -c bioconda -c conda-forge python=3.6 hmmer prodigal
conda activate viralverify
#下载脚本
cd ~/virus2/temp
git clone https://github.com/ablab/viralVerify.git
mv viralVerify viralverify
chmod -R 755 viralverify
cd viralverify


#使用
cd ~/virus2/temp/viralverify
conda activate viralverify
export PATH="$HOME/virus1/temp/Prodigal/:$PATH"

./bin/viralverify -f ~/virus2/seq/no_duplicates.fna -o output -t 16 --hmm ~/db/viralverify  

        Optional arguments:
        -h, --help  Show the help message and exit
        --db DB     Run BLAST on input contigs against provided db
        -t          Number of threads
        -thr THR    Sensitivity threshold (minimal absolute score to classify sequence, default = 7)
        -p          Output predicted plasmidic contigs separately
5.病毒定量、树图
5.1coverm
https://github.com/wwood/CoverM
######安装converm
conda create -n coverm -c conda-forge -c bioconda coverm 
5.2 humann3
https://github.com/biobakery/humann
######安装humann3
conda create -n humann3 -c bioconda -c conda-forge humann=3.7
5.3 iqtree
https://github.com/iqtree/iqtree2
######安装iqtree及对齐软件mafft
conda create -n iqtree -c bioconda -c conda-forge iqtree mafft
conda activate iqtree
#下载VOG数据库，获得77个vog_marker
mkdir ~/db/VOG/2017_04_02hmm
cd ~/db/VOG/2017_04_02hmm
wget -c https://fileshare.lisc.univie.ac.at/vog/vog81/vog.hmm.tar.gz
tar -xzvf vog.hmm.tar.gz 
#下载77个vog_id,下载的文件为表格且有多个页面，已经整理到了/data4/machuang/db/VOG/VOGmarker_id.txt
#wget -c https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-019-0448-z/MediaObjects/41564_2019_448_MOESM3_ESM.xlsx
#筛选77个vog.hmm
mkdir -p ~/db/VOG/vog77  
#筛选VOG_hmm数据库里我们需要的77个模型
grep -oP '^\s*VOG\d+' ~/db/VOG/VOGmarker_id.txt | while read vog_id; do  
    find ~/db/VOG/2017_04_02hmm -type f -name "*$vog_id*.hmm" -exec cp -t ~/db/VOG/vog77 {} +  
done
#合并需要的77个marker模型
cat ~/db/VOG/vog77/*.hmm > ~/db/VOG/vog77/VOG77.hmm 
#复制到
cp ~/db/VOG/vog77/VOG77.hmm ~/virus1/
5.4 trimal
https://github.com/inab/trimal
#######安装trimal
#使用trimAI，对marker进行裁剪
cd ~/virus1/temp/iqtree
wget -c https://codeload.github.com/inab/trimal/zip/refs/heads/trimAl
unzip trimAl
rm -rf trimAl
mv trimal-trimAl trimal
cd trimal/source
make
#将该文件夹添加到环境中
export PATH="$HOME/virus1/temp/iqtree/trimal/source:$PATH"
5.5 fasta36
https://github.com/wrpearson/fasta36
#######安装fasta36，及依赖命令脚本
cd ~/virus1/temp
wget -c https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8i.tar.gz
tar -xzvf fasta-36.3.8i.tar.gz
mv fasta-36.3.8i fasta36
chmod -R 755 fasta36
cd ~/virus1/temp/fasta36/src
make -f ../make/Makefile.linux64_sse2 all #会出现warning不用管
#所有的blast命令脚本都在cd ~/virus1/temp/fasta36/bin文件夹内，加载环境后就可用了
export PATH="$HOME/virus1/temp/fasta36/bin/:$PATH"
fasta36 -h

#下载用到的命令脚本（f2s、joincol、hashsums、s2f、seqlengths、tree_bray）
cd ~/virus1/temp
wget -c https://codeload.github.com/shiraz-shah/VFCs/zip/refs/heads/main
unzip main
rm -rf main
mv VFCs-main VFCs
chmod -R 755 VFCs
#下载用到的命令脚本rapidnj
wget -c https://codeload.github.com/somme89/rapidNJ/zip/refs/heads/master
unzip master
rm -rf master
mv rapidNJ-master rapidnj
chmod -R 755 rapidnj
cd rapidnj
make
5.6 MSPminer
cd ~/virus2/temp
# MSPminer v 1.1.3
wget -c https://www.enterome.com/wp-content/uploads/2023/04/MSPminer.zip
unzip MSPminer.zip
mv MSPminer_1_1_3 MSPminer
chmod -R 755 MSPminer
cd MSPminer

mkdir -p output
#自行更换 settings.ini文件内的 计数矩阵
#count_matrix_file=
#output_dir=
./mspminer settings.ini
6.生活方式预测，安装bacphlip、Deephage
6.1bacphlip
https://github.com/adamhockenberry/bacphlip
#需按照此方法安装，否则会报错
#numpy要<1.24.0，不然会出现AttributeError: module 'numpy' has no attribute 'float'.
conda create -n bacphlip -c bioconda -c conda-forge  python biopython=1.77 pandas joblib scikit-learn hmmer numpy=1.23.5
conda activate bacphlip
pip install bacphlip
6.2Deephage
https://github.com/shufangwu/DeePhage
######conda 安装依赖项
conda create -n deephage -c bioconda -c conda-forge  python==3.6.7 numpy==1.16.4 h5py==2.9.0 tensorflow==1.4.0 keras==2.1.3
conda activate deephage

######下载脚本
cd ~/virus2/temp
wget -c https://codeload.github.com/shufangwu/DeePhage/zip/refs/heads/master
unzip master
rm master
mv DeePhage-master deephage
chmod -R 755 deephage

######下载依赖软件MCR
mkdir -p ~/db/software/mcr 
cd ~/db/software/mcr
#安装包1.6G 位置~/db/software/mcr/MCR_R2018a_glnxa64_installer.zip


wget -c https://ssd.mathworks.com/supportfiles/downloads/R2018a/deployment_files/R2018a/installers/glnxa64/MCR_R2018a_glnxa64_installer.zip
#解压


unzip MCR_R2018a_glnxa64_installer.zip
#安装


./install -mode silent -agreeToLicense yes -destinationFolder  /data4/machuang/db/software/MCR
#加载环境
export LD_LIBRARY_PATH="$HOME/db/software/MCR/v94/runtime/glnxa64:$HOME/db/software/MCR/v94/bin/glnxa64:$HOME/db/software/MCR/v94/sys/os/glnxa64:$HOME/db/software/MCR/v94/extern/bin/glnxa64:$LD_LIBRARY_PATH"




######使用测试
cd ~/virus2
temp/DeePhage/DeePhage temp/DeePhage/example.fna temp/DeePhage/example.csv

6.3PhaTYP
https://github.com/KennthShang/PhaTYP
cd ~/virus2/temp
git clone https://github.com/KennthShang/PhaTYP.git
cd PhaTYP/
conda env create -f phatyp.yaml -n phatyp #有时候会报出来HTTPs连接错误，多试几遍
conda activate phatyp

fileid="1tsUArctGf9Fd3xa-0sEcp6ykwxTy9uxG"
filename="model.zip"
html=`curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}"`
curl -Lb ./cookie "https://drive.google.com/uc?export=download&`echo ${html}|grep -Po '(confirm=[a-zA-Z0-9\-_]+)'`&id=${fileid}" -o ${filename}
#命令行国内不能下载，需要vpn，下载网址 https://drive.google.com/uc?export=download&id=1tsUArctGf9Fd3xa-0sEcp6ykwxTy9uxG
#该文件服务器地址/data4/machuang/db/PhaTYP/model.zip
 
unzip model.zip
pip install .  
#Successfully installed PhaTYP-0.3.0

7.下载百岁老人文件脚本
https://github.com/RasmussenLab/vCentenarian
cd ~/virus1/temp
wget -c https://codeload.github.com/RasmussenLab/vCentenarian/zip/refs/heads/master
unzip master
tar -xzvf vcentenarian/code/06_evaluation/mgvrepo.tar.gz -C vcentenarian/code/06_evaluation/
mv vCentenarian-master vcentenarian
chmod -R 755 vcentenarian
rm -rf master
8.宿主预测 iphop  
 https://bitbucket.org/srouxjgi/iphop/src/main/
conda create -n iphop python=3.8
conda activate iphop
conda install -c conda-forge -c bioconda iphop
iphop -h #v1.3.3
iphop predict -h

#查看数据库版本为Sept_2021_pub_rw     iphop需>= 1.3.0
#tar -tf /data/meta/db/viwrap/iPHoP.latest_rw.tar.gz
#安装数据库108G  /data/meta/db/viwrap/iPHoP.latest_rw.tar.gz
#直接解压，替换代码中path_to_iPHoP_db 为自己db位置
tar -xzf /data/meta/db/viwrap/iPHoP.latest_rw.tar.gz -C path_to_iPHoP_db
chmod -R 755 path_to_iPHoP_db
9.病毒物种注释
9.1 安装DemoVir
https://github.com/feargalr/Demovir（作者提示该软件自2018年来就没有维护，数据库过旧，不能再用于注释了）
#安装依赖
cd ~/virus2/temp
git clone https://github.com/feargalr/Demovir.git
mv Demovir demovir
chmod -R 755 demovir

#装所需软件
cd ~/db/software
#wget -c http://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
#gunzip usearch11.0.667_i86linux32.gz
#mv usearch11.0.667_i86linux32 usearch
#chmod 755 usearch

wget -c https://github.com/rcedgar/usearch_old_binaries/raw/main/bin/usearch11.0.667_i86linux64
mv usearch11.0.667_i86linux64 usearch
chmod +x usearch

#下载数据库（软件自带，比较旧）
cd ~/virus2/temp/demovir
wget -c https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/10442241/nr.95.fasta.bz2
bunzip2 nr.95.fasta.bz2
#建库
~/db/software/usearch -makeudb_ublast nr.95.fasta -output uniprot_trembl.viral.udb &> usearch_database.log
chmod 755 uniprot_trembl.viral.udb

#下载新版uniprot病毒数据库 2024-04
#https://www.uniprot.org/uniprotkb?query=(reviewed:false)%20AND%20(taxonomy_id:10239)
mkdir -p ~/db/demovir
cd ~/db/demovir
wget -O uniprotkb_10239_2024_07_29.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28reviewed%3Afalse%29+AND+%28taxonomy_id%3A10239%29%29"
gunzip uniprotkb_10239_2024_07_29.fasta.gz 
mv uniprotkb_10239_2024_07_29.fasta uniprot_trembl_viral.fasta
seqkit stat uniprot_trembl_viral.fasta #5,671,882
#改序列名称,只保留>和||中间的内容
awk '/^>/ {split($0, a, "\\|"); print ">"a[2]; next} {print}' uniprot_trembl_viral.fasta > output.fasta
#去冗余
export PATH="$HOME/virus1/temp/cdhit/:$PATH"
cd-hit-est -i output.fasta -o nr.95.fasta -c 0.95 -g 0 -M 0 -d 0 -n 10 
seqkit stat nr.95.fasta
 #5,652,936
#建库
~/db/software/usearch --makeudb_usearch nr.95.fasta --output uniprot_trembl.viral.udb &> usearch_database.log
chmod 755 uniprot_trembl.viral.udb
cp uniprot_trembl.viral.udb ~/virus2/temp/demovir/

#下载对应注释
wget -O uniprotkb_10239_2024_07_29.tsv.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clineage_ids%2Clineage&format=tsv&query=%28%28reviewed%3Afalse%29+AND+%28reviewed%3Afalse%29+AND+%28taxonomy_id%3A10239%29%29"
gunzip  uniprotkb_10239_2024_07_29.tsv.gz
mv uniprotkb_10239_2024_07_29.tsv  uniprot_trembl_viral.tsv
wc -l uniprot_trembl_viral.tsv #5,671,883
#处理对应注释
cut -f 1,5,7 uniprot_trembl_viral.tsv | \
awk -F '\t' 'BEGIN {OFS="\t"} 
NR==1 {print "id", "order", "family", "organism"; next} 
{
  order = ""; family = ""; 
  split($3, arr, ", ");
  for (i in arr) {
    if (arr[i] ~ /order/) {
      order = arr[i];
      gsub(/ \(order\)| \(suborder\)/, "", order);
    }
    if (arr[i] ~ /family/) {
      family = arr[i];
      gsub(/ \(family\)| \(subfamily\)/, "", family);
    }
  }
  print $1, order, family, $2
}' > uniprot_trembl_viral_1.tsv

wc -l uniprot_trembl_viral_1.tsv #5,671,883
cp uniprot_trembl_viral_1.tsv ~/virus2/temp/demovir/
 
#制作注释的RDS文件
conda activate hmmer
cd ~/virus2/temp/demovir/
R
# 安装和加载 readr 包
install.packages("readr")
library(readr)
# 读取 TSV 文件
tsv_data <- read_tsv("uniprot_trembl_viral_1.tsv")
# 保存为 RDS 文件
saveRDS(tsv_data, file ="./uniprot_trembl_viral.RDS")
# 安装和加载tibble包
install.packages("tibble")
library(tibble)
#tibble转化为普通数据
taxa <- readRDS("uniprot_trembl_viral.RDS")
taxa <- as.data.frame(taxa)

9.2 安装vConTACT2
https://bitbucket.org/MAVERICLab/vcontact2/src/master/
#conda 安装
conda create -n vcontact2 python=3.8 -y
conda activate vcontact2
conda install -c conda-forge networkx=2.2 numpy=1.20.1 scipy=1.6.0 pandas=1.1.5 scikit-learn=0.24.1 biopython=1.78 hdf5=1.10.4 pyparsing=2.4.6 psutil=5.8.0 -y
conda install -c conda-forge pytables=3.6.1 -y


conda install -c bioconda mcl=14.137 diamond=0.9.14 -y
conda install -c bioconda blast=2.10.1 -y
conda install -c conda-forge openjdk #安装java


#下载脚本
cd ~/virus2/temp
#直接网页下载，~/db/MAVERICLab-vcontact2-c0413a6c92e8.tar.gz
wget -c https://bitbucket.org/MAVERICLab/vcontact2/get/master.tar.gz #wget下载报错

tar -xzvf ~/db/MAVERICLab-vcontact2-c0413a6c92e8.tar.gz -C ~/virus2/temp/

mv MAVERICLab-vcontact2-c0413a6c92e8 vcontact2
chmod -R 755 vcontact2
cd ~/virus2/temp/vcontact2
wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -O cluster_one-1.0.jar
python3 -m pip install .

9.3 安装VirusTaxo
https://github.com/omics-lab/VirusTaxo
cd  ~/virus2/temp
git clone https://github.com/omics-lab/VirusTaxo
mv VirusTaxo virustaxo
chmod -R 755 virustaxo
conda create -n virustaxo python=3.10
conda activate virustaxo
cd virustaxo
python3 -m venv environment
source ./environment/bin/activate
pip install -r requirements.txt #有点慢
#biopython==1.83  numpy==2.0.0 pandas==2.2.2 python-dateutil==2.9.0 pytz==2024.1  six==1.16.0 tqdm==4.66.4 tzdata==2024.1





python3 predict.py -h
####  ~/db/virustaxo/vt_db_jan21_2024.tar.gz
mkdir -p ~/db/virustaxo
cd ~/db/virustaxo
#数据库下载，网页下载需要挂vpn，已下载至/data4/machuang/db/virustaxo/vt_db_jan21_2024.tar.gz
wget -c https://drive.google.com/file/d/1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs/view?pli=1 #wget下载会报错
tar -xzvf vt_db_jan21_2024.tar.gz
#DB Building cmd:
#python3 build.py \
#   --meta ./Dataset/RNA_meta.csv \
#   --seq ./Dataset/RNA_seq.fasta \
#   --k 17 \
#   --saving_path ./model/RNA.pkl

#共有三个，使用~/db/virustaxo/vt_db_jan21_2024/DNA_RNA_18451_k20.pkl
#~/db/virustaxo/vt_db_jan21_2024/DNA_9384_k21.pkl
#~/db/virustaxo/vt_db_jan21_2024/DNA_RNA_18451_k20.pkl
#~/db/virustaxo/vt_db_jan21_2024/RNA_9067_k17.pkl

#运行
cd ~/virus2/temp/virustaxo
conda activate virustaxo
source ./environment/bin/activate
python3 predict.py -h
mkdir -p ~/virus2/result/taxonomy/virustaxo ~/virus2/temp/taxonomy/virustaxo
python3 predict.py \
   --model_path ~/db/virustaxo/vt_db_jan21_2024/DNA_RNA_18451_k20.pkl \
   --seq ~/virus2/result/vOTUs.fna  > ~/virus2/temp/taxonomy/virustaxo/output.txt
#Entropy<0.5的被认为可信
deactivate
9.3 安装VirusTaxo_Hierarchical （测试）
https://github.com/omics-lab/VirusTaxo_Hierarchical
cd ~/virus2/temp/
git clone https://github.com/omics-lab/VirusTaxo_Hierarchical.git
chmod -R 755 VirusTaxo_Hierarchical
cd VirusTaxo_Hierarchical
conda create -n VH python=3.8 pandas
conda activate VH
python3 -m venv environment
source ./environment/bin/activate
pip install -r requirements.txt
#pip版本 pip install --upgrade pip
#将requirements.txt中的依赖提取出来单独下载 #可能会下载超时，多试几次
#pip cache purge 清理缓存，重新安装；
#设置timeout时间，避免pip下载报错 --default-timeout=100

#模型训练
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
9.4 安装PhaGCN
https://github.com/KennthShang/PhaGCN
cd ~/virus2/temp
#git clone  https://github.com/KennthShang/PhaGCN
wget -c https://codeload.github.com/KennthShang/PhaGCN/zip/refs/heads/main
unzip main
mv PhaGCN-main phagcn
chmod -R 755 phagcn
cd phagcn
conda env create -f environment.yaml -n phagcn
conda activate phagcn

#运行
python run_Speed_up.py --contigs ~/virus2/result/vOTUs.fna --len 5000
9.5 安装PhaGCN2
https://github.com/KennthShang/PhaGCN2.0
#PhaGCN2是一个基于GCN的模型，它可以通过深度学习分类器学习物种掩蔽特征，用于新的病毒分类学分类。要使用PhaGCN 2，只需将您的Contigs输入到程序中。
#推荐使用Conda安装，environment.yaml从作者github链接下载：https://github.com/KennthShang/PhaGCN2.0
cd ~/virus2/temp
#wget -c https://codeload.github.com/KennthShang/PhaGCN2.0/zip/refs/heads/main
#unzip main
git clone https://github.com/KennthShang/PhaGCN2.0.git  #不稳定多试
chmod -R 755 PhaGCN2.0-main
mv PhaGCN2.0-main phagcn2.0
cd phagcn2.0
conda env create -f environment.yaml -n phagcn2
#在使用之前准备数据库，从作者github链接下载，链接参上
cd database
tar -zxvf ALL_protein.tar.gz
cd ..
#加载环境
conda activate phagcn2
export MKL_SERVICE_FORCE_INTEL=1
#使用方法
python run_Speed_up.py --contigs ~/virus2/result/vOTUs.fna --len 5000
#参数设置--contigs #需要输入的contigs文件，通过megahit组装--len #对contigs长度进行阈值设置，默认值为8000 bp
#注意事项
1.确认输入的contigs为病毒contigs,可通过Virsorter或DeepVirFinder获得
2.congtigs序列中的id不能含有空格等中文字符
（可以不安装）8.6 安装PhaBOX
https://github.com/KennthShang/PhaBOX
cd ~/virus2/temp
wget -c https://codeload.github.com/KennthShang/PhaBOX/zip/refs/heads/main
unzip main
mv PhaBOX-main phabox
chmod -R 755 phabox
cd phabox
pip install aiohttp==3.8.1
conda env create -f webserver.yml -n phabox
conda activate phabox


#下载数据库
mkdir -p ~/db/phabox
cd ~/db/phabox
#提供的有百度云下载
Link: https://pan.baidu.com/s/18gx_p-Y4g22W5LcXvIyO_A pwd: uran
Link: https://pan.baidu.com/s/1QJQAIr89xbt4e3pJr_QhaQ pwd: 2gjb

#使用
cd ~/virus2
conda activate phabox
#run all pipline
python main.py --contigs ~/virus2/result/vOTUs.fna  --threads 8 --len 5000 --rootpth simple_test --out out/ --dbdir database/ --parampth parameters/ --scriptpth scripts/
# run PhaMer
python PhaMer_single.py --contigs test_contigs.fa --threads 8 --len 3000 --rootpth simple_test --out out/ --dbdir database/ --parampth parameters/ --scriptpth scripts/

# run PhaTYP
python PhaTYP_single.py --contigs test_contigs.fa --threads 8 --len 3000 --rootpth simple_test --out out/ --dbdir database/ --parampth parameters/ --scriptpth scripts/

# run PhaGCN
python PhaGCN_single.py --contigs test_contigs.fa --threads 8 --len 3000 --rootpth simple_test --out out/ --dbdir database/ --parampth parameters/ --scriptpth scripts/

# run CHERRY
python Cherry_single.py --contigs test_contigs.fa --threads 8 --len 3000 --rootpth simple_test --out out/ --dbdir database/ --parampth parameters/ --scriptpth scripts/
10.病毒功能注释
10.1 数据库下载(找主页中的FTP)：
KEGG  https://www.kegg.jp/kegg/            在线注释https://www.kegg.jp/blastkoala/                /data/meta/db/kegg/2101/genes/viruses
TIGRFAM https://www.jcvi.org/research/tigrfams  FTP： https://ftp.ncbi.nlm.nih.gov/hmm/current/
Pfam-A http://pfam.xfam.org/        FTP：  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/   
VOGDB https://vogdb.org/      FTP：https://fileshare.lisc.univie.ac.at/vog/vog224/
Earth’s Virome viral protein families database https://www.nature.com/articles/nature19094  FTP: https://portal.nersc.gov/dna/microbial/prokpubs/EarthVirome_DP/
碳水化合物活性酶CAZy  https://bcb.unl.edu/dbCAN2/blast.php
毒力基因注释数据库 VFDB http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
耐药基因CARD https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
mkdir -p ~/db/KEGG/virus ~/db/VOG ~/db/Pfam_A ~/db/TIGRFAM
#KEGG
cd ~/db/KEGG/virus
cp /data/meta/db/kegg/2101/genes/viruses/*.gz ./


#VOGDB
cd ~/db/VOG
wget -c  https://fileshare.lisc.univie.ac.at/vog/vog224/vog.virusonly.tsv.gz
gunzip vog.virusonly.tsv.gz
wget -c https://fileshare.lisc.univie.ac.at/vog/vog224/vog.faa.tar.gz
tar -xzvf vog.faa.tar.gz
mv faa 2024_06_09faa  #共60250个
cd 2024_06_09faa
cat *.faa > ../24_6_9merge.faa

#Pfam-A,蛋白
cd ~/db/Pfam_A
wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz
gunzip Pfam-A.fasta.gz 
wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
#Pfam-A.fasta： Protein 56937050


#TIGRFAM,隐马尔可夫模型（HMM）
cd ~/db/TIGRFAM
wget -c https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB

#Earth’s Virome viral protein families database
mkdir -p ~/db/Earth
cd ~/db/Earth
#隐马尔可夫模型（HMM）
wget -c https://portal.nersc.gov/dna/microbial/prokpubs/EarthVirome_DP/final_list.hmms
#病毒fna
wget -c https://portal.nersc.gov/dna/microbial/prokpubs/EarthVirome_DP/mVGs_sequences_v2.fna

#数据库注释，、CAZy、VFDB、CARD
# Swiss-Prot
mkdir -p ~/db/SwissProt
cd ~/db/SwissProt
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz 
#构建索引   --in接刚才下载的fasta文件，-d接数据库的名字
diamond makedb --in uniprot_sprot.fasta -d swissprot
#-d指定数据库, -q指定输入文件, -k指定hits数，也就是获取比对结果的多少条，这里我们做注释其实一条就可，-e就是e-value啦，最后是输出文件，根据实际情况命名
diamond blastx -d swissprot -q /mnt/d/WorkSpace/100Doing/102Work/xxx/Seq/xxx.fas -k 1 -e 0.00001 -o xxx.swiss_dia_matches.m8

#Pfam数据库注释  版本：2024-05-28
conda activate hmmer
conda install -c bioconda -c conda-forge pfam_scan
#建索引
cd ~/db/Pfam_A
hmmpress Pfam-A.hmm
#比对
pfam_scan.pl -fasta ../xxx.faa -dir . -outfile xxx.pfam_matches.txt
#克隆pfam_scan的脚本到本地并运行
cd ~/virus2/temp
git clone https://github.com/aziele/pfam_scancd pfam_scan
./pfam_scan.py --help
#运行命令
./pfam_scan.py test/test.fasta pfamdb/
#test.fasta为输入文件，pfamdb/为数据库位置，可通过-out添加输出文件位置

#hmmsearch
#conda activate hmmer
#hmmsearch -h#3.1b2
#hmmbuild /data4/machuang/virus/iqtree/hmmer/aligend_all.hmm /data4/machuang/virus/iqtree/hmmer/aligend_all.faa
#hmmsearch --noali -E 1e-5 merge.hmm /data4/machuang/virus/prodigal_result/1118.proteins.faa > /data4/machuang/virus/iqtree/hmmer/final_virus.txt
10.2 安装prokka
#生成gff文件，需先安装prokka
cd ~/virus2
conda create -y -n prokka -c conda-forge -c bioconda prokka blast=2.9
conda activate prokka

cd ~/virus2/temp
wget -c https://codeload.github.com/tseemann/prokka/zip/refs/heads/master
unzip master

mv prokka-master prokka
chmod -R 755 prokka
cd prokka

#运行
#先生成gff文件  --outdir指定输出目录，--prefix 指定输出文件前缀 ，--centre X --compliant自动生成符合条件的contig ID 长度
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
10.3 安装Resistance Gene Identifier v.6.0.3
识别CARD耐药基因
https://github.com/arpcard/rgi
card https://card.mcmaster.ca/download
cd ~/virus2/temp
git clone https://github.com/arpcard/rgi
chmod -R 755 rgi
cd rgi
conda env create -f conda_env.yml
conda activate rgi
pip install --use-pep517 -r requirements.txt
#pip install lxml==4.9.1

python setup.py build
python setup.py test #只是测试，出现error，找不到文件可忽略
python setup.py install #如果出现Download error，多试几次

rgi -h #6.0.3

#下载配置CARD数据库
mkdir -p ~/db/card/3.2.9
cd ~/db/card

# 清理之前的加载
rgi clean --local

# 下载并解压CARD数据
wget -c https://card.mcmaster.ca/download/0/broadstreet-v3.2.9.tar.bz2
tar -xvf broadstreet-v3.2.9.tar.bz2 -C ./3.2.9


# 创建注释文件
rgi card_annotation -i ./3.2.9/card.json > card_annotation.log 2>&1


# 加载数据到RGI
rgi load \
  --card_json ./3.2.9/card.json \
  --debug --local \
  --card_annotation card_database_v3.2.9.fasta \
  --card_annotation_all_models card_database_v3.2.9_all.fasta
#~/db/card/localDB


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
10.4 安装NCBI AMRFinder tool v.3.12.8
识别AMG辅助代谢基因
https://github.com/michaelwoodworth/AMRFinder_scripts
 #blast版本要≥2.2，amrfinder=3.12
conda create -y -c conda-forge -c bioconda -n amrfinder ncbi-amrfinderplus
conda activate amrfinder
#数据库只能下载到默认位置，不能指定，https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/
#Downloading AMRFinder database version 2024-05-02.2 into '/data4/machuang/miniconda3/envs/amrfinder/share/amrfinderplus/data/2024-05-02.2/'
amrfinder -u    # Database version: 2024-05-02.2

#下载脚本
cd ~/virus2/temp
wget -c https://codeload.github.com/michaelwoodworth/AMRFinder_scripts/zip/refs/heads/master
unzip master
mv AMRFinder_scripts-master amrfinder
chmod -R 755 amrfinder
 
#下载测试文件
mkdir -p ~/virus2/temp/amrfinder/test
cd ~/virus2/temp/amrfinder/test
BASE_URL=https://raw.githubusercontent.com/ncbi/amr/master
curl --silent -L \
     -O ${BASE_URL}/test_dna.fa \
     -O ${BASE_URL}/test_prot.fa \
     -O ${BASE_URL}/test_prot.gff \
     -O ${BASE_URL}/test_both.expected \
     -O ${BASE_URL}/test_dna.expected \
     -O ${BASE_URL}/test_dna_mut_all.expected \
     -O ${BASE_URL}/test_prot.expected \
     -O ${BASE_URL}/test_amrfinder.sh
#测试
bash test_amrfinder.sh path  #最后出现Success!即为成功

cd ~/virus2
conda activate amrfinder
#输入为核苷酸

amrfinder -n result/vOTUs.fna --plus -o temp/amrfinder/vOTUs_fna_amrfinder.tsv
#输入为蛋白
amrfinder -p result/vOTUs.faa --plus -o temp/amrfinder/vOTUs_faa_amrfinder.tsv
10.5 安装Resfams database
耐药基因，核心数据库来自CARD、LacED、 β-内酰胺酶蛋白的独特抗生素耐药蛋白序列
http://www.dantaslab.org/resfams
mkdir -p ~/db/Resfams
cd ~/db/Resfams
#Resfams HMM Database (Core) - v1.2, updated 2015-01-27
#Database version for annotation of microbial proteins in the absence of any functional confirmation for antibiotic resistance.
wget -c https://dantaslab.wustl.edu/resfams/Resfams.hmm.gz
gunzip Resfams.hmm.gz
wget -c https://static1.squarespace.com/static/5402b1a0e4b02a7794494453/t/5a8e1150085229db68b03444/1519259984513/180102_resfams_metadata_updated_v1.2.2.xlsx
mv 180102_resfams_metadata_updated_v1.2.2.xlsx Resfams_profile_HMM_Metadata_v1.2.2.xlsx
#Resfams HMM Database (Full) - v1.2, updated 2015-01-27
#Database version for annotation of microbial proteins when functional confirmation for antibiotic resistance is available (such as functional metagenomic selections).
wget -c https://dantaslab.wustl.edu/resfams/Resfams-full.hmm.gz
gunzip Resfams-full.hmm.gz

# Using hmmsearch
conda activate hmmer

10.6 安装eggNOG
https://github.com/eggnogdb/eggnog-mapper
在线分析http://eggnog-mapper.embl.de/
#安装eggnog-emapper
conda create -n eggnog -c bioconda python=3.7 biopython=1.76 psutil=5.7.0 xlsxwriter=1.4.3
conda activate eggnog
conda install -c bioconda -y eggnog-mapper
cd ~/virus2/temp
wget -c https://codeload.github.com/eggnogdb/eggnog-mapper/zip/refs/heads/master
unzip master
mv eggnog-mapper-master eggnog
chmod -R 755 eggnog


#mkdir ~/db/eggnog
#./download_eggnog_data.py -y --data_dir ~/db/eggnog
#镜像数据下载
#ftp://download.nmdc.cn/tools/eggnog/eggnog.db.gz
#ftp://download.nmdc.cn/tools/eggnog/eggnog_proteins.dmnd.gz
10.7 安装DRAM/DRAM-V 
病毒注释和识别AMG
https://github.com/WrightonLabCSU/DRAM
cd ~/virus2/temp
wget -c https://codeload.github.com/WrightonLabCSU/DRAM/zip/refs/heads/master
unzip master
mv DRAM-master DRAM
chmod -R 755 DRAM
cd DRAM
conda env create -f environment.yaml -n DRAM
conda activate DRAM

#已有数据库
DRAM-setup.py prepare_databases --output_dir DRAM_data --kegg_loc kegg.pep

#没有数据库，设置DRAM可能需要很长时间（最多5个小时），并且默认情况下使用大量内存（512 gb）。
#profiles.tar.gz、ko_list.gz、uniref90.fasta.gz、Pfam-A.full.gz、Pfam-A.hmm.dat.gz、dbCAN-HMMdb-V9.txt、CAZyDB.07302020.fam-activities.txt、vog.hmm.tar.gz、vog.annotations.tsv.gz、viral.x.protein.faa.gz、
#mag_annotator/database_setup.py中376行改为  merge_files(glob(path.join(hmm_dir, 'hmm/VOG*.hmm')), vog_hmms)
mkdir -p ~/db/DRAM
# --skip_uniref :reduce memory usage to ~64 Gb
#https://github.com/WrightonLabCSU/DRAM/issues/305
#https://github.com/WrightonLabCSU/DRAM/issues/182
./scripts/DRAM-setup.py prepare_databases --output_dir ~/db/DRAM --threads 6 --skip_uniref

#使用
#annotate some MAGs
~/virus2/temp/DRAM/scripts/DRAM-v.py annotate -i 'my_bins/*.fa' -o annotation
#summarize annotations
~/virus2/temp/DRAM/scripts/DRAM-v.py distill -i annotation/annotations.tsv -o genome_summaries --trna_path annotation/trnas.tsv --rrna_path annotation/rrnas.tsv
#结果文件
DRAM-v liquor是对已在带注释的病毒重叠群中检测到的潜在 AMG (pAMG) 的总结，HTML 文件的形式。
annotations.tsv包含所有预测的开放阅读框架的所有注释
vMAG_stats.tsv提供有关每个病毒重叠群的详细信息
11.蛋白聚类
11.1 安装HMMScan
https://github.com/josh-wilde/hmmscan


11.2 安装DGRscan
Diversity-Generating Retroelements
https://github.com/YuzhenYe/DGRscan

11.3安装MMseqs2 
sequence search and clustering
https://github.com/soedinglab/MMseqs2
