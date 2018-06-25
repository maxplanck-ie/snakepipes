#$ -cwd
#$ -j y
#$ -V
#$ -S /bin/bash
#$ -o /dev/null

INPREFIX=$1
OUTFILE=$2

#Assumes a that the prefix is a set of files trimmed with trim_galore ==> ending with val_1.fq.gz

#make use of that C and G are complementary; in the CT read, the number of Gs are an estimate of the number of Cs and vice versa in the GA read. 

for s in ${INPREFIX}'_R1.fastq.gz'
do
 gunzip -c $s |awk -vRS=`gunzip -c $s|head -1 |awk '{print substr($0,0,7)}'` -vFS='\n' ' {seq=$2} NR>1 {startG+=gsub(/^G/,"G",seq);G+=gsub(/G/,"G",seq);CG+=gsub(/CG/,"CG",seq);C+=gsub(/C/,"C",seq);endC+=gsub(/C$/,"C",seq)} END {print startG,G,CG,C,endC}'
done |awk '{startG+=$1;G+=$2;CG+=$3;C+=$4;endC+=$5} END{sumG=G-CG-startG;sumC=C-CG-endC; if(sumG>sumC){print "READ1\tCT",100*(1-sumC/sumG)}else{print "READ1\tGA",100*(1-sumG/sumC)}}' >> $OUTFILE

for s in ${INPREFIX}'_R2.fastq.gz'
do
 gunzip -c $s |awk -vRS=`gunzip -c $s|head -1 |awk '{print substr($0,0,7)}'` -vFS='\n' ' {seq=$2} NR>1 {startG+=gsub(/^G/,"G",seq);G+=gsub(/G/,"G",seq);CG+=gsub(/CG/,"CG",seq);C+=gsub(/C/,"C",seq);endC+=gsub(/C$/,"C",seq)} END {print startG,G,CG,C,endC}'
done |awk '{startG+=$1;G+=$2;CG+=$3;C+=$4;endC+=$5} END{sumG=G-CG-startG;sumC=C-CG-endC; if(sumG>sumC){print "READ2\tCT",100*(1-sumC/sumG)}else{print "READ2\tGA",100*(1-sumG/sumC)}}' >> $OUTFILE

