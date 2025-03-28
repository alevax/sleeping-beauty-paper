# This script is called by TAP2.pl

# This script calls list2tab.pl

# Version Fisher_clean_TKS_v1

# Sarver annotation
# this is a script to automate data workup following Fishers exact test association analyses 
# This script takes the R_result and the names file and then generates
# 4 4xls spreadsheets
# 3 images

@name;
@sub_name;
@result;
@sub_result;
$count = 0;
unless (-d "results/Assoc") {
    mkdir("results/Assoc");
}
open OUT1, "> results/Assoc/Associations.xls";
open OUT2, "> results/Assoc/Ann_cis_list.txt";
open OUT3, "> results/Assoc/Ann_cis_table.xls";
print OUT3 "Event1\tEvent2\tFET Pvalue\tBonferroni Significance\tBH 20%FDR Significancei\n";
open OUT4, "> results/Assoc/Cis_cis_list.txt";
open OUT5, "> results/Assoc/Cis_cis_table.xls";
print OUT5 "Event1\tEvent2\tFET Pvalue\tBonferroni Significance\tBH 20%FDR Significance\n";
open OUT6, "> results/Assoc/Ann_ann_list.txt";
open OUT7, "> results/Assoc/Ann_Ann_table.xls";
print OUT7 "Event1\tEvent2\tFET Pvalue\tBonferroni Significance\tBH 20%FDR Significance\n";
open SOURCE, "< results/Assoc/R_result.txt";
@tableAnnCis;
@tableAnnAnn;
@tableCisCis;
while (defined($line = <SOURCE>)) {
$count++;
chomp $line;
@result= split(/\t/, $line);
@sub_result= split(/\s=\s|\s<\s/, $result[1]);
print $sub_result[1];
$pval[$count]= $sub_result[1];
}
$count = 0;
open SOURCE, "< results/Assoc/Fisher_pre_named.txt";
while (defined($line = <SOURCE>)) {
$count++;
chomp $line;
@name= split(/\t/, $line);
@sub_name= split(/__/, $name[0]);
print OUT1 "$line\t$pval[$count]\n";
if (($sub_name[1] =~m/descriptor/) and ($sub_name[0] !~ m/descriptor/)) {
$viewpval= -(log($pval[$count])/log(10));
print OUT2 "$sub_name[1]\t$sub_name[0]\t$viewpval\n";
if ($pval[$count] < 0.05) {
$data4table = "$sub_name[1]\t$sub_name[0]\t$pval[$count]\n";
push(@tableAnnCis, $data4table);
}
}


if (($sub_name[1] !~m/descriptor/) and ($sub_name[0] !~ m/descriptor/)) {
$viewpval= -log($pval[$count]);
print OUT4 "$sub_name[1]\t$sub_name[0]\t$viewpval\n";
print OUT4 "$sub_name[0]\t$sub_name[1]\t$viewpval\n";
if ($pval[$count] < 0.05) {
$data4table = "$sub_name[1]\t$sub_name[0]\t$pval[$count]\n";
push(@tableCisCis, $data4table);
}
}

if (($sub_name[1] =~m/descriptor/) and ($sub_name[0] =~ m/descriptor/)) {
$viewpval= -log($pval[$count]);
print OUT6 "$sub_name[1]\t$sub_name[0]\t$viewpval\n";
print OUT6 "$sub_name[0]\t$sub_name[1]\t$viewpval\n";
if ($pval[$count] < 0.05) {
$data4table = "$sub_name[1]\t$sub_name[0]\t$pval[$count]\n";
push(@tableAnnAnn, $data4table);
}
}


}
@sortedAnnCis = sort { (split '\t', $a)[2] <=> (split '\t', $b)[2] } @tableAnnCis;
@sortedCisCis = sort { (split '\t', $a)[2] <=> (split '\t', $b)[2] } @tableCisCis;
@sortedAnnAnn = sort { (split '\t', $a)[2] <=> (split '\t', $b)[2] } @tableAnnAnn;

open SOURCE, "< FISH/multi_test.txt";
while (defined($line = <SOURCE>)) {
chomp $line;
@name=split(/\t/, $line);
$corrAnnCis = $name[2];
$corrCisCis = $name[3];
$corrAnnAnn = $name[4];
}

$pos = 1;
foreach (@sortedAnnCis) {
$bonf = 0.05/$corrAnnCis;
$fdr = ($pos/$corrAnnCis)*0.20;
$pos++;
chomp $_;
print OUT3 "$_\t$bonf\t$fdr\n";
}
$pos = 1;

foreach (@sortedCisCis) {
$bonf = 0.05/$corrCisCis;
$fdr = ($pos/$corrCisCis)*0.20;
$pos++;
chomp $_;
print OUT5 "$_\t$bonf\t$fdr\n";
}
$pos =1;
foreach (@sortedAnnAnn) {
$bonf = 0.05/$corrAnnAnn;
$fdr = ($pos/$corrAnnAnn)*0.20;
$pos++;
chomp $_;
print OUT7 "$_\t$bonf\t$fdr\n";
}


close OUT7;
close OUT6;
close OUT5;
close OUT4;
close OUT3;
close OUT2;
close OUT1;
system ("perl lib/list2tab.pl 1 2 3 0 results/Assoc/Ann_cis_list.txt > results/Assoc/Ann_cis_matrix.txt");
system ("perl lib/list2tab.pl 1 2 3 0 results/Assoc/Cis_cis_list.txt > results/Assoc/Cis_cis_matrix.txt");
system ("perl lib/list2tab.pl 1 2 3 0 results/Assoc/Ann_ann_list.txt > results/Assoc/Ann_ann_matrix.txt");
