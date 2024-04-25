
cd /fs/ess/PAS1444/orthology_Vl_Vv

for inputfile in *.blast; do



cut -f1,2 "$inputfile" |
sort -u -k1,1 > checkpointfile1.txt 
sort -u -k2,2 checkpointfile1.txt > checkpointfile2.txt
#Do the above to get a line count of 'inputfile' when columns 1 and 2 when all duplicates are removed


#Actual work starts here \/
awk -F '\t' -v OFS='\t' '{ gsub("-mRNA*", "", $1) ; print }' "$inputfile" |
awk -F '\t' -v OFS='\t' '{ gsub(".m*", "", $1) ; print }' |
awk -F '\t' -v OFS='\t' '{ gsub("-mRNA*", "", $2) ; print }' |
awk -F '\t' -v OFS='\t' '{ gsub(".m*", "", $2) ; print }' > intermediatefile1.txt
#replace any locations where a '_t' and then anything after it exists with 'nothing' in column 1
#replace any locations where a '-T' and then anything after it exists with 'nothing' in column 1
#replace any locations where a '_t' and then anything after it exists with 'nothing' in column 2
#replace any locations where a '-T' and then anything after it exists with 'nothing' in column 2

awk -F'\t' 'match($1, "VvCab") {print NR}' intermediatefile1.txt > col1wVitvinumbers.txt
awk -F'\t' 'match($2, "VvCab") {print NR}' intermediatefile1.txt > col2wVitvinumbers.txt
#find the rows which contain 'Vitvi' in the 1st column, then print those line numbers to a new file
#find the rows which contain 'Vitvi' in the 2nd column, then print those line numbers to a new file
#(I checked the number of lines in each of the two paired files (with and without 'Vitvi' for a column)
#			in excel by sorting and filtering and they match exactly what I was expecting)

cat col1wVitvinumbers.txt col2wVitvinumbers.txt > intermediatefile2.txt
sort intermediatefile2.txt > intermediatefile3.txt
uniq -u intermediatefile3.txt > intermediatefile4.txt
#concatentate the file which stated which lines contained 'Vitvi' in col 1 with the file which stated
#			which lines contained 'Vitvi' in col 2
#sort the newly concatenated list
#keep only those row values which are uniq 

awk 'NR==FNR{data[$1]; next}FNR in data' intermediatefile4.txt intermediatefile1.txt > intermediatefile5.txt
#pull out the rows of data which contained a GREM4 to PN collinear gene association
#			(FNR = folder number record; tells program to use those numbers as input for action
#			the NR tells it that the FNR numbers are actually NR (row number) numbers
#			so the above script is just telling it to pull those rows)

sort -k1,1 -k3,3nr intermediatefile5.txt > intermediatefile6.txt
#sort by the first column, then by the third (the percent identity column); in the third column, the data is
#		numeric and it should be sorted in reverse order (checks out, I checked)

sort -u -k1,1 intermediatefile6.txt > intermediatefile7.txt 
#remove duplicate gene names from the first column (they do exist, I checked).
#the sort+uniq command is down here because if we simply took the first hit up at the beginning it would
#		simply take the perfect match of the gene which would just say that gene GREM123 matched gene
#		GREM123 and then we would never see any orthologous relationships come up (GREM123 matches PN152)
#		By waiting until now to do the sort+uniq, we have already gotten rid of all of the instances
#		where column one and two contain genes from the same species.  (checked on excel and checked out)

sort -k2,2 -k3,3nr intermediatefile7.txt > intermediatefile8.txt
#sort by the second column, then by the third (the percent identity column); in the third column, the data is
#		numeric and it should be sorted in reverse order (checks out, I checked)

sort -u -k2,2 intermediatefile8.txt > intermediatefile9.txt
#remove duplicate gene names from the first column (they do exist, I checked).
#then did a sort+uniq on the 2nd column to get rid of duplicates there (checked in excel and it checked out)
#		It is important to remember that this file still contains two listings of most all genes though,
#		they have one listing in column 1 and one listing in column 2 (checks out, I checked)

awk -F '\t' 'match($1, "VvCab") {print $2 "\t" $1} !match($1, "VvCab") {print $1 "\t" $2}' intermediatefile9.txt > intermediatefile10.txt
#flip flop the 1st and 2nd columns if 'Vitvi' is found in the first column

sort -u -k1,1 intermediatefile10.txt > intermediatefile11.txt
#remove duplicates from first column (as the flip flop has introduced a second of copy of each name now)

sort -u -k2,2 intermediatefile11.txt > "$inputfile".orthology_list.txt
#remove duplicates from second column (as the flip flop has introduced a second of copy of each name now)

done