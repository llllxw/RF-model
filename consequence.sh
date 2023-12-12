##切换至工作目录
cd /c/amplicon/22pipeline

##建立临时文件
mkdir -p temp

##建立最终结果的文件
mkdir -p result

##解压原始测序数据
gunzip T1seq/*.gz

##less按页查看，空格翻页、q退出；head查看前10行
head seq/SRR7584004_1.fastq

# 查看序列1和2谁为前段，谁为reverse  
less T1seq/SRR21402622_1.fastq
grep -i "GTGCCAGC" T1seq/SRR21402622_1.fastq


###pipeline.sh 分析主流程及依赖注释数据库
#db/*。fa.gz 16S数据库：GreenGene(GG),RDP,SILVA and UNIT
#usearch数据库database下载页：http:www.drive5.com/sintax/
gunzip database/rdp_16s_v16_sp.fa.gz
gunzip database/SILVA_138_SSURef_NR99_tax_silva.fasta.gz
head -n2 database/rdp_16s_v16_sp.fa
#QIIME greengene 有参数据库：ftp://greengenes,microbio.me/

##2.合并双端序列与样品标签
# 已WT1单样品合并为例，
time vsearch -fastq_mergepairs seq/SRR7584005_1.fastq \
-reverse seq/SRR7584005_2.fastq \
-fastqout temp/SRR7584005.merged.fq \
-relabel SRR7584005

#依照实验设计批量处理并合并
for i in `tail -n+2 t1metadata.txt | cut -f 1`;do
vsearch --fastq_mergepairs T1seq/${i}_2.fastq \
--reverse T1seq/${i}_1.fastq \
--fastqout t1temp/${i}.merged.fq --relabel ${i}.
done &
 
#将单端序列改名  
# # 单个序列改名示例
i=WT1
gunzip -c seq/${i}_1.fq.gz > seq/${i}.fq
usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.

# 批量改名，需要有单端fastq文件，且解压(usearch不支持压缩格式)
gunzip seq/*.gz
time for i in `tail -n+2 t1metadata.txt|cut -f1`;do
usearch -fastx_relabel T1seq/${i}.fastq -fastqout t1temp/${i}.merged.fq -prefix ${i}.
done &
  
  #合并所有样品至同一文件
cat t1temp/*.merged.fq > t1temp/all.fq
ls -l temp/all3/all3.fq
head temp/all3/all3.fq 

#3.切除引物与质控，根据文献中的左右序列长短以及相应的barcode位置进行jianqie
time vsearch --fastx_filter t1temp/all.fq \
--fastq_stripleft 20 --fastq_stripright 21 \
--fastq_maxee_rate 0.01 \
--fastaout t1temp/filtered.fa

#可选，每个样本切除接头和引物，纯净序列用于下游分析
mkdir -p seq/clean
for i in `tail -n+2 sample-metadata.tsv | cut -f 1`;do
vsearch --fastx_filter temp/${i}.merged.fq \
--fastq_strupleft 19 --fastq_stripright 20 \
--fastqout seq/clean/${i}.fq
vsearch --fastq_mergepairs seq/${i}_1.fq --reverse seq/${i}_2.fq \
--fastqout temp/${i}.merged.fq --relabel ${i}.

#4.去冗余与生成OTUS
#4.1序列去冗余

#去冗余Find unique read sequences and abundances
#并添加miniuniqusize最小为8或1/1M，去除低丰度，增加计算速度，
#-sizeout输出丰度， --relabel必须，否则文件头不正常
#输出非冗余序列数据小于1万条为宜
time vsearch --derep_fulllength t1temp/filtered.fa \
  --output t1temp/uniques.fa --relabel Uni --minuniquesize 120 --sizeout

#Ss,告丰度非冗余序列非常小(<2Mb)
ls -l t1temp/uniques.fa

#4.2 生成OTU/ASV

#有两种方法：unoise3推荐（种/株水平精度），传统的97%聚类OTU（属水平精度）
#usearch自带de novo去嵌合体

#方法1. 97%聚类OTU，不推荐，除非reviewer要求或数据量太大
./usearch -cluster_otus temp/uniques.fa \
-otus temp/otus.fa \
-relabel OTU_

#方法2.ASV非聚类去噪法
usearch -unoise3 t1temp/uniques.fa \
-zotus t1temp/zotus.fa

#修改序列名：格式调整 format OTU prefix
sed 's/Zotu/OTU_/g' t1temp/zotus.fa > t1temp/otus.fa

#方法3.vsearch生成OTU，无自动de novo去冗余功能
#权限usearch免费版本受限时（通过提高minisize减少数据量）使用，不推荐
#重命名、相似97%，不屏蔽，输入和输出count
vsearch --cluster_size temp/uniques.fa \
--centroids temp/otus_raw.fa \
--relabel OTU_ --id 0.97 --qmask none --sizein --sizeout
#vsearch还需连用
vsearch --uchime_denovo temp/otus_raw.fa \
--nonchimeras result/otus.fa
#4.3 基于参考去嵌合

#方法1.usearch+rdp去嵌合（快但容易假阴性）
#官方模式sensitive容易假阴性高，推荐balanced，high_confidence
./usearch -uchime2_ref temp/otus.fa -db database/rdp_16s_v16_sp.fa \
-chimeras temp/otus_chimeras.fa -strand plus -mode balanced

#获得非嵌合体序列ID
cat temp/otus.fa temp/otus_chimeras.fa| grep '>' | sort | uniq -u \
| sed 's/>//' > temp/non_chimeras.id
#筛选非嵌合体
./usearch -fastx_getseqs temp/otus.fa -labels temp/non_chimeras.id \
-fastaout result/otus.fa

#方法2.vsearch+silva去嵌合（慢15min-3d但略好）
#注释数据库和统计
silva=database/SILVA_138_NR99.fasta
gunzip ${silva}.gz
ls -l ${silva}
grep -c '>' ${silva}
time vsearch --uchime_ref t1temp/otus.fa \
--db ${silva} \
--nonchimeras t1temp/otus_silva.fa
gzip ${silva} &
  ##确定使用，修改名称进入主流程
  mv t1temp/otus_silva.fa t1result/otus.fa 

#方法3.不去嵌合，执行如下命令
cp -f temp/otus.fa result/otus.fa

#5.1生成OTU表
time vsearch --usearch_global t1temp/filtered.fa --db t1result/otus.fa \
--otutabout t1result/otutab.txt --id 0.97 --threads 30

# SILVA物种注释
gunzip database/SILVA_138_NR99.fasta.gz
vsearch --usearch_global result/otus.fa --db database/SILVA_138_NR99.fasta --biomout out_tax.txt --id 0.97

time usearch --sintax result/otus.fa --db database/SILVA_138_NR99.fasta \
-strand both --tabbedout result/SILVAotus.sintax --sintax_cutoff 0.6

vsearch --usearch_global t1result/otus.fa -db database/SILVA_138_NR99.fasta --id 0.97 --query_cov 0.97 --strand both --biomout t1result/SILVAotus.txt --alnout t1result/SILVAotus.aln --blast6out t1result/SILVAotus_blast.xls --fastapairs t1result/SILVAotus_pairs.fa --notmatched t1result/SILVAotus_notmatched.fa --userfields query+target+id+qcov+tcov --userout t1result/SILVAotus_stat.xls


# RDP物种注释
gunzip database/rdp_16s_v16.fa.gz
vsearch --usearch_global result/otus.fa --db database/rdp_16s_v16.fa --biomout result/RDPout_tax.txt --id 0.97

time usearch --sintax result/otus.fa --db database/rdp_16s_v16_sp.fa \
-strand both --tabbedout result/RDPotus.sintax --sintax_cutoff 0.6

vsearch --usearch_global result/otus.fa -db database/rdp_16s_v16.fa --id 0.97 --query_cov 0.97 --strand both --biomout result/RDPotus.txt --alnout result/RDPotus.aln --blast6out result/RDPotus_blast.xls --fastapairs result/RDPotus_pairs.fa --notmatched result/RDPotus_notmatched.fa --userfields query+target+id+qcov+tcov --userout result/RDPotus_stat.xls

# GG物种注释
vsearch --usearch_global result/otus.fa --db database/gg_16s_13.5.fa --biomout gg/GGotus.txt --id 0.97
vsearch --usearch_global result/otus.fa -db database/gg_16s_13.5.fa --id 0.97 --query_cov 0.97 --strand both --biomout gg/GGotus.txt --alnout gg/GGotus.aln --blast6out gg/GGotus_blast.xls --fastapairs gg/GGotus_pairs.fa --notmatched gg/GGotus_notmatched.fa --userfields query+target+id+qcov+tcov --userout gg/GGotus_stat.xls

time usearch --sintax result/otus.fa --db database/rdp_16s_v16_sp.fa \
-strand both --tabbedout result/RDPotus.sintax --sintax_cutoff 0.6


#usearch-32位运行内存为4G，数量大时不能使用
# time usearch --sintax result/otus.fa --db database/gg_16s_13.5.fa \
# -strand both --tabbedout gg/GGotus.sintax --sintax_cutoff 0.6


#5.2OTU表简单统计
time usearch -otutab_stats result/otutab.txt \
-output result/otutab.stat
cat result/otutab.stat

# #5.3等量抽样标准化(可见文件呢OTU_standanrdization,R)
# usearch -otutab_norm result/otutab.txt \
#     -sample_size 80000 \
#     -output result/otutab_norm.txt
# usearch1 -otutab_stats result/otutab_norm.txt \
#     -otuput result/otutab_norm_report.txt
# cat result/outtab_norm_report.txt


#6.与已知物种比较更像谁()SILVA
usearch -sintax result/otus.fa -db database/SILVA_138_NR99.fasta \
-strand both -tabbedout result/out_tax.txt -sintax_cutoff 0.6
cut -f 1,4 out_tax.txt|sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/taxonomy-silva.txt

#6.与已知物种比较更像谁()RDP
usearch -sintax result/otus.fa -db database/rdp_16s_v16_sp.fa \
-strand both -tabbedout result/RDPotus_tax.txt -sintax_cutoff 0.6
cut -f 1,4 result/RDPotus_tax.txt|sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/taxonomy-RDP.txt

#6.与已知物种比较更像谁()GG
usearch -sintax gg/otus.fa -db database/gg_16s_13.5.fa \
-strand both -tabbedout gg/sintax.txt -sintax_cutoff 0.6
cut -f 1,4 gg/GGotus.txt|sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > gg/taxonomy-GG.txt


##比对greengene数据库生成的OTU表,基于它的OTUs表可进行功能预测（PICRUSt）和形态学预测（Bugbase）
# gg=database/gg_16s_13.5.fa
# gunzip ${gg}.gz
# ls -l ${gg}
# grep -c '>' ${gg}
# time vsearch --uchime_ref gg/otus.fa \
#   --db ${gg} \
#   --nonchimeras gg/otus_gg.fa
# gzip ${gg} &
# 
# mv gg/otus_gg.fa gg/otus1_PICRUST.fa
# 
# vsearch --usearch_global temp/filtered.fa --db gg/otus1_PICRUST.fa \
#   --otutabout gg/otutab.txt --id 0.97 --threads 14

###9.0有参比对——功能预测，如Greengenes,可用于picurst，bugbase分析
mkdir -p gg/
  gunzip database/gg_16s_13.5.fa.gz
#与gg所有97%OTUs比对，用于功能预测
#方法一
usearch -otutab temp/filtered.fa -otus gg/97_otus.fasta \
--otutabout gg/otutabs.txt -threads 14

#方法2，vsearch比对
time vsearch --usearch_global temp/filtered.fa --db gg/otus.fa \
--otutabout gg/otutab.txt --id 0.97 --threads 14

#统计
usearch -otutab_stats gg/otutab.txt -output gg/otutab.stat
cat gg/otutab,stat


###使用PICRUST2预测菌落功能
#使用 bioconda 安装 PICRUSt2 环境
conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.3.0_b
source activate picrust2

#将扩增子序列插入参考进化树中，结果为epa_result_parsed.newick
mkdir -p picrust2/

#序列不合格，将otu.fa反向设置
seqtk seq -r otus.fa >otus_rc.fa

place_seqs.py -s otus_rc.fa \
  -o picrust2/ -p 24 \
  --intermediate placement_rcworking

#预测16S的拷贝数，-n计算NSTI，预测基因组的E.C.编号，KO相度
hsp.py -i 16S -t picrust2/epa_result_parsed.newick -o picrust2/marker_nsti_predicted.tsv.gz -p 8 -n

#预测16S的酶学委员会EC编号
hsp.py -i EC -t picrust2/epa_result_parsed.newick -o picrust2/EC_predicted.tsv.gz -p 8

#预测16S的基因同源簇KO编号(时间长)
hsp.py -i KO -t picrust2/epa_result_parsed.newick -o picrust2/KO_predicted.tsv.gz -p 8

#预测E.C.和KO的丰度，—strat_out输出分层结果(极大的增加计算时间)
metagenome_pipeline.py -i otutab.txt \
   -m picrust2/marker_nsti_predicted.tsv.gz \
   -f picrust2/EC_predicted.tsv.gz \
   -o picrust2/EC_metagenome
 
metagenome_pipeline.py -i otutab.txt \
   -m picrust2/marker_nsti_predicted.tsv.gz \
   -f picrust2/KO_predicted.tsv.gz \
   -o picrust2/KO_metagenome
   
#MetaCyc通路丰度，基于EC结果汇总
  
pathway_pipeline.py -i picrust2/EC_metagenome/pred_metagenome_unstrat.tsv.gz \
    -o picrust2/pathways \
    --intermediate pathways_working \
    -p 8

#KEGG通路，基于KO
pathway_pipeline.py -i picrust2/KO_metagenome/pred_metagenome_unstrat.tsv.gz -o picrust2/KEGG_pathways --no_regroup --map ~/.conda/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv

#结果注释
cd picrust2

#添加EC的注释
add_descriptions.py -i EC_metagenome/pred_metagenome_unstrat.tsv.gz -m EC \
  -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

#KO添加注释
add_descriptions.py -i KO_metagenome/pred_metagenome_unstrat.tsv.gz -m KO \
  -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz


#pathway添加注释
add_descriptions.py -i pathways/path_abun_unstrat.tsv.gz -m METACYC \
  -o pathways/path_abun_unstrat_descrip.tsv.gz

#7.0计算样品内的丰富度、均匀度
#alpha_div命令基于标准化OUT表计算14种指数
mkdir -p alpha/
usearch -alpha_div t1result/otu_Flattening.txt \
-output alpha/alpha.txt

#结果有多个文件，需要目录

#稀释抽样：1%-100%抽样一百次（richness）
usearch -alpha_div_rare t1result/otu_Flattening.txt \
-output alpha/alpha_rare.txt -method without_replacement



##附件；Reads count整数值如何标准化相对丰度
# 求取各个OTU在对应样品的丰度频率
usearch -otutab_counts2freqs result/otutab_mean_Genotype.txt \
-output result/otutab_rare_Genotype.txt

##筛选各组高丰度菌用于比较

#按组求均值，需根据实验设计metadata.txt修改组列名
#输入文件为feautre表result/otutab.txt,实验设计metadata.txt
#输出为特征表按组的均值-一个实验可能有多种分组方式
Rscript script/otu_mean.R

# 筛选各组高丰度菌用于比较
awk 'BEGIN{OFS=FS="\t"}{if (FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
      else {for(i=2; i<=NF;i++) if($i>0.05) print $1, a[i]; }}' \
result/otutab_mean_Genotype.txt > alpha/otu_group_exist.txt
head -n9 alpha/otu_group_exist.txt 

#基于OTU构建进化树 Make OTU tree
usearch -cluster_agg t1result/otus.fa \
-treeout t1result/otus.tree
#结果有多个文件，需要目录
mkdir -p beta/
  # 生成5种距离矩阵：bray——curtis，euclidean，jaccard,manhatten,unifrac,又有非权重版本（_binary）
  usearch -beta_div t1result/otu_Flattening.txt \
-tree t1result/otus.tree \
-filename_prefix beta/
  
  
  ## 8. 物种注释分类汇总
  
  #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
  cut -f 1,4 result/RDPotus_tax.txt \
|sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
> result/taxonomy.txt
head -n3 result/taxonomy.txt

#OTU对应物种8列格式：注意注释是非整齐
#生成物种表格OTU/ASV中空白补齐为Unassigned
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
result/taxonomy-RDP.txt > temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> result/taxonomy.txt
head -n3 result/taxonomy.txt

#统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
mkdir -p result/tax
for i in p c o f g;do
usearch -sintax_summary result/otus.sintax \
-otutabin result/otutab_rare.txt -rank ${i} \
-output result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
# 列出所有文件
wc -l result/tax/sum_*.txt
head -n3 result/tax/sum_g.txt


###9.0物种注释结果分类汇总
#筛选对应序列
cut -f 1 result/otutab.txt | tail -n+2 > temp/otutab.id
usearch -fastx_getseqs result/otus.fa \
-labels temp/otutab.id \
-fastaout result/otus.fa


#筛选与OTU表顺序对应的注释
awk 'ARGIND==1{a[$1]=$0}NR>ARGIND==2{print a[$1]}' \
result/otus.sintax temp/otutab.id \
> result/final/sintax,txt


awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' \
result/otus.sintax temp/otutab.id \
> result/final/sintax,txt

#补齐末尾列
# sed -i 's/\t$/\td:Unassigned/' result/final/otus.sintax
head -n2 result/final/sintax.txt 


#报错处理删除末注释行
head -n2 result/sintax.txt

#统计门纲目科属，使用rank参数p c o f g
mkdir -p-tax
for i in p c o f g;do
usearch -sintax_summary result/final/sintax.txt \
-otutabin result/otu_Flattening.txt \
-rank ${i} \
-output tax/sum_${i}.txt
done 
sed -i 'S/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' tax/sum_*.txt

#OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释
cut -f 1,4 result/RDPotus_tax.txt \
|sed 's/\td/\tk/;s/:/_/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
> result/RDPtaxonomy2.txt
head -n1 result/RDPtaxonomy2.txt

# #23. 用Unassign补齐上面taxonomy中的空白
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned"; \
 split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
 print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}'
result/taxonomy-RDP.txt >temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> result/final/RDPtaxonomy.txt
head -n3 result/RDPtaxonomy.txt



#10.项目空间清理

###10.1删除中间大文件

#删除双端合并文件
rm temp/*.merged.fq
#删除合并文件
rm temp/all.fq

###10.2压缩原始数据和数据库
#注：原始数据及时压缩并上传数据中心备份
#压缩原始文件节省空间
gzip seq/*
  #数据库经常用，保存至指定位置
  #压缩数据库
  time gzip database/*.fa


#23,LEfSe,STAMP统计分析和可视化
mkdir -p result/lefse

Rscript ${ea}/script/format2lefse.R 
--input result/otutab.txt \
--taxonomy result/taxonomy.txt --design result/metadata.txt \
--group Group --threshold 0.1 \
--output result/lefse/LEfSe



#Bugbase细菌表型预测
git clone https://github.com/knights-lab/BugBase
cd BugBase
    #安装依赖包
    export BUGBASE_PATH=`pwd`
    export PATH=$PATH:`pwd`/bin
    #安装了所有依赖包
    run.bugbase.r -h
    #测试数据
    run.bugbase.r -i doc/data/HMP_s15.txt -m doc/data/HMP_map.txt -c HMPBODYSUBSITE -o output
    
 ### 2. 准备输入文件
cd ~/amplicon/result
#输入文件：基于greengene OTU表的biom格式(本地分析支持txt格式无需转换)和mapping file(metadata.txt首行添加#)
#上传实验设计+刚才生成的otutab_gg.txt
#生成在线分析使用的biom1.0格式
biom convert -i otutabs.txt -o otutab_gg.biom --table-type="OTU table" --to-json
sed '1 s/^/#/' metadata.txt > MappingFile.txt
#下载otutab_gg.biom 和 MappingFile.txt用于在线分析

### 3. 本地分析
export BUGBASE_PATH=`pwd`
export PATH=$PATH:`pwd`/bin
run.bugbase.r -i otutab_gg.txt -m MappingFile.txt -c Group -o phenotype/
