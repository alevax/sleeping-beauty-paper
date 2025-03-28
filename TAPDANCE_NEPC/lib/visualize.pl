# This script is called by TAP2.pl

# This script calls list2tab.pl

require 'config.pl';
open OUT, "> FISH/library_cis_list.txt";
open OUT2, ">FISH/multi_test.txt";

$query1 = "select count(name) from cis_result_final_$proj where pvalueregion < $cocis_threshold and type = 'all_$library_percent' and gene_name not like '%BAD%' and region > 1 order by pvalueregion limit 100";
$sth = $dbh->prepare("$query1");

$sth->execute;

while ((@row) = $sth->fetchrow_array) {
$cis_number = $row[0];
}

$query2="select name,gene_name, library_name from cis_result_final_$proj where pvalueregion < $cocis_threshold and type = 'all_$library_percent' and gene_name not like '%BAD%' and region > 1 order by pvalueregion limit 100";

$sth = $dbh->prepare("$query2");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
#print $row[2];
@chunks = split(/::/m, $row[2]);
foreach $lib(@chunks) {
if ($lib ne '') {
print OUT "$row[0]::$row[1]\t$lib\t1\n";
}
} 
}
close OUT;

$sth = $dbh->prepare("select count(distinct descriptor) from metadata_$proj where descriptor != 'all';");

$sth->execute;

while ((@row) = $sth->fetchrow_array) {
$desc_number = $row[0];
}
$ann_cis_testing = $desc_number*$cis_number;
$cis_cis_testing = $cis_number*($cis_number-1)/2;
$ann_ann_testing = $desc_number*($desc_number-1)/2;
print OUT2 "$cis_number\t$desc_number\t$ann_cis_testing\t$cis_cis_testing\t$ann_ann_testing\n";
close OUT2;
$sth = $dbh->prepare("select 'descriptor', descriptor, library from metadata_$proj;");
$sth->execute;

open OUT, "> FISH/library_descriptor_list.txt";

while ((@row) = $sth->fetchrow_array) {
print OUT "$row[0]::$row[1]\t$row[2]\t1\n";
}

close OUT;

open (my $list, ">", "FISH/list.txt");
open (my $cis_list, "<", "FISH/library_cis_list.txt");
while(<$cis_list>) {
    print $list $_;
}
close($cis_list);
open (my $desc_list, "<", "FISH/library_descriptor_list.txt");
while (<$desc_list>) {
    print $list $_;
}
close($desc_list);
close($list);
#system ("cat FISH/library_cis_list.txt FISH/library_descriptor_list.txt > FISH/list.txt");
system ("perl lib/list2tab.pl 1 2 3 0 FISH/list.txt  > FISH/cis4cluster.txt");
system ("perl lib/list2tab.pl 1 2 3 0 FISH/library_cis_list.txt > results/Assoc/Cis.txt");
system ("perl lib/list2tab.pl 1 2 3 0 FISH/library_descriptor_list.txt > results/Assoc/Ann.txt");

print "generated the matrix file cis4cluster.txt";


