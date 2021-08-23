####### Parse EasySFS output
wd=/net/harris/vol1/home/beichman/vaquita/sfs/easySFS/projection_preview
categories="missense synonymous"

for category in $categories 
do
easyOut=${category}.easySFS.projPreview.txt
pops="GoC"  # just one pop
for pop in $pops
do
# put a "$" at the end of pop to signify that it should be at the end of the line
# get the line from the easy SFS output; get rid of of the parentheses and make a column
echo "projection,snps" > $wd/$pop.${easyOut%.txt}.R.format.txt
grep -A1 "$pop$" $wd/$easyOut | tail -1 | sed 's/(//g' | sed 's/)//g' | sed 's/, /,/g' |  tr '\t' '\n' >> $wd/$pop.${easyOut%.txt}.R.format.txt
done
done