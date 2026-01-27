# This is the TAP4.pl script 

# The goal of this script is to combine several TAPDANCE projects
#	These projects will have overlap in their libraries and TAP4 will combine all the 
#	insertions from multiple projects that map to a single library.
#	TAP4 will then identify all the CISs present in a NR set

# Version: TKS v1

# Version history
# one file is required nr.txt sorted in chromosome order!
# Aaron Lyman Sarver Oct 21, 2010
# This script takes over following the bowtie mapping and generates CIS.
# Update NOV15,2010 to pass library information through 
#$Fully modified december 1 2010i

require 'config.pl';
require 'lib/feature_finder_methods.pl';
my $path = $0;
$path =~ s/\/\w*\.pl$//g;

# I don't think this following line is needed anymore
# require "$path/project_man.pl";

unless (-d "CIS") {
    mkdir("CIS");
}
unless (-d "results") {
    mkdir("results");
}

$proj_all = $proj."_all";

$sth = $dbh->prepare("drop table if exists lib_mappable_$proj_all");
$sth->execute;

$sth = $dbh->prepare("drop table if exists lib_mappable_$proj");
$sth->execute;

$sth = $dbh->prepare("drop table if exists illumina_hist_$proj_all");
$sth->execute;

$sth = $dbh->prepare("drop table if exists illumina_hist_$proj");
$sth->execute;

$sth = $dbh->prepare("drop table if exists bowtie_lib_$proj_all");
$sth->execute;

$sth = $dbh->prepare("drop table if exists bowtie_lib_$proj");
$sth->execute;

$sth = $dbh->prepare("drop table if exists metadata_$proj_all");
$sth->execute;


open DIN, "< data/meta.tab";
@met;
while (defined($line = <DIN>)) {
    chomp $line;
    @field= split(/\t/, $line);
    push (@met,$line);
}

$proj_all = $proj."_all";
$var=shift(@met);

$sth = $dbh->prepare("create table lib_mappable_$proj_all select * from lib_mappable_$var");
$sth->execute;
$sth = $dbh->prepare("create table illumina_hist_$proj_all select * from illumina_hist_$var");
$sth->execute;
$sth = $dbh->prepare("create table bowtie_lib_$proj_all select * from bowtie_lib_$var");
$sth->execute;
$sth = $dbh->prepare("create table metadata_$proj_all select * from metadata_$var");
$sth->execute;

foreach $var(@met) {
$sth = $dbh->prepare("insert into lib_mappable_$proj_all select * from lib_mappable_$var");
$sth->execute;
$sth = $dbh->prepare("insert into  illumina_hist_$proj_all select * from illumina_hist_$var");
$sth->execute;
$sth = $dbh->prepare("insert into  bowtie_lib_$proj_all select * from bowtie_lib_$var");
$sth->execute;
$sth = $dbh->prepare("insert into metadata_$proj_all select * from metadata_$var");
$sth->execute;
}

unless(-e "data/metadata.tab") {
    $sth = $dbh->prepare("select library, descriptor, type from metadata_$proj_all");
    $sth->execute;
    open(my $meta_out, ">", "data/metadata.tab") || die "Unable to open data/metadata.tab for output. $!\n";
    while ((@row) = $sth->fetchrow_array) {
	print $meta_out join("\t", @row) . "\n";
    }
    close($meta_out);
}


$sth = $dbh->prepare("create table lib_mappable_$proj select distinct library, sum(total) as total from lib_mappable_$proj_all group by library;");
$sth->execute;

$sth = $dbh->prepare("create table illumina_hist_$proj select distinct library, chromo, start, sum(count) as count,orient from illumina_hist_$proj_all group by library,chromo,start,orient;");
$sth->execute;

$sth = $dbh->prepare("create table bowtie_lib_$proj select distinct library,chromo, start,stop, pos, sum(count) as count,strand, orient from bowtie_lib_$proj_all group by library,chromo,start,stop,pos,strand order by chromo,start+0,orient ;");
$sth->execute;

$sth = $dbh->prepare("select chromo, start, stop, concat(library,'-',count),count,strand,start,stop from bowtie_lib_$proj");
$sth->execute;
#7.2
open OUT, "> results/raw_$proj.BED";
print OUT "track type ='BED' name='raw$proj' description='all$proj' visibility=2 itemRgb='On'";
while ((@row) = $sth->fetchrow_array) {
 print OUT "\n$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]";}
close OUT;
#7.3
$sth = $dbh->prepare("drop table if exists convert_2_hex_$proj");
$sth->execute;
$sth = $dbh->prepare("create table convert_2_hex_$proj (orient varchar(5), hex varchar(20))");
$sth->execute;
$sth = $dbh->prepare("insert into convert_2_hex_$proj VALUES(\"A\", \"255,0,0\")");
$sth->execute;
$sth = $dbh->prepare("insert into convert_2_hex_$proj VALUES(\"B\", \"0,0,255\")");
$sth->execute;
$sth = $dbh->prepare("insert into convert_2_hex_$proj VALUES(\"U\", \"0,0,0\")");
$sth->execute;

$sth = $dbh->prepare("select chromo, start, stop, concat(library,'-',count),count,strand,start,stop from bowtie_lib_$proj where count > 9");
### JJE
#$sth = $dbh->prepare("select chromo, start, start+1, concat(library,count),count,'+',start,start+1, hex from bowtie_lib_$proj A, convert_2_hex_$proj B where A.orient = B.orient and count > 9");
###
$sth->execute;
open OUT, "> results/raw10_$proj.BED";
print OUT "track type ='BED' name='raw10$proj' description='10$proj' visibility=2 itemRgb='On'";
while ((@row) = $sth->fetchrow_array) {
 print OUT "\n$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]";}
close OUT;

&set_project_status($proj, 1, 0);

