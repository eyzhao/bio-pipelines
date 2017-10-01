file=$1
regionFile=$2
output=$3

while read file
do
    regionSetName=$(basename $regionFile)
    while read regionFile
    do
        filename=$(basename $file)
        echo 'Running mpileup on file $file' > ./logs/mpileup.log
        echo "samtools mpileup -g -f /projects/ezhao_prj/data/references/hg19/GRCh37-lite.fa -l $regionFile $file | bcftools call -c -v > '$output/$filename.$regionSetName.vcf'"
    done < $regionFile
done < $file
