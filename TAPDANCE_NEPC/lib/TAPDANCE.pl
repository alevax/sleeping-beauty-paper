# Description:  Processing of transposon insertion sequences.
# Script must be executed from within the project folder with the command "perl lib/TAPDANCE.pl"
# For detailed instructions see the file "Tim's Guide to TAPDANCE".

# Version history:
# 	This script conducts the PAPS steps 1
# 	Originally written by Aaron Lyman Sarver July 12 2010
# 	Updated to version 3.1 December 1 2010
# 	Updated to version 3.2 using module load bowtie/0.12.8 Nuri Alpay Temiz May 9 2013
# 	Updated to version TKS1 November 21, 2013. This script works after loading bowtie/1.0.0
#		and the bowtie index files are described below

# MySQL DB: requires perl DBI read write accessible MySQL database. The MySQL DB name and password
#	must be entered in the config.pl file.

# Bowtie: 
#	As of 11/23/13 using bowtie/1.0.0. 
#	Currently running using MSI's bowtie/1.0.0.
#	Before running tapdance, I load the bowtie module from MSI system using the following command: 
#		module load bowtie/1.0.0
#	Bowtie index files (.ebwt) reside in the following directory:
#		/home/starrt2/shared/genomes/mm10/bowtieindex/
#	The index files are named like: genome.1.ebwt
#	Bowtie is executed five times in tapdance. If you are using a different genome index
#		you need to change the directory in five system call lines.  The lines look like this:
#			system ('bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f /home/starrt2/shared/genomes/mm10/bowtieindex/genome mapping/4map_1.txt mapping/out.txt');
#		these lines reside in section 5.1, 5.2, 5.3, 5.4, and 5.5
#	If you change the version of bowtie, make sure the options used in this script will work with
#		the new version.  If not, change them on the same lines as above.

# Project folder must contain the config.pl file and two folders, data/ and lib/
#	The data/ folder must contain three files: barcode.txt, chromo.tab, and seqs.tab
#	The lib/ folder must contain the following 17 perl scripts and the mm9.bed file: 
#		cis_finder.pl
#		convert_nr_2_annot.pl
#		feature_finder_methods.pl
#		feature_finder.pl
#		Fisher_clean_nocis.pl
#		Fisher_clean.pl
#		generate_cis_matrix.pl
#		generate_genome_location.pl
#		invert.pl
#		list2tab.pl
#		matrix2fisher_nocis.pl
#		matrix2fisher.pl
#		mm9.bed
#		shuffle.pl
#		TAP2.pl
#		TAP4.pl
#		TAPDANCE.pl
#		visualize.pl

# TAPDANCE.pl will work on fiLes of less than 40 million sequences. 

# TAPDANCE identifies and removes barcodes, IRDR sequences and linker sequences from sequencing data, 
#	then compiles a list of non redundant sequences. TAPDANCE creates and uses MySQL databases which
#	all contain the Project name, which is specified in the config.pl file.

# After running TAPDANCE, use the script TAP2.pl to calculate common insertion sites.

###################################################################################################
# Print the starting time of TAPDANCE
# Use the time function and localtime(time) function to create a time variable
# Set the abbreviations
my @days_week = qw(Sun Mon Tue Wed Thu Fri Sat);
my @months = qw (Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
# Set the initiating time using the time function
my $start_time = time;
# Convert the start time to the local time
my @time=localtime($start_time);
my $local_sec = $time[0];
my $local_min = $time[1];
my $local_hr = $time[2];
my $local_day_of_month = $time[3];
my $local_mon = $time[4]+1;
my $local_yr = $time[5]+1900;
my $local_weekday = $time[6];
print "TAPDANCE v.TKS1\n";
# To print out time as 7/15/2013 14:15:32
print "Started:\t$local_mon/$local_day_of_month/$local_yr\t$local_hr:$local_min:$local_sec\n";
###################################################################################################
# 1.0 Set up directories and MySQL tables
require 'config.pl';
system ('mkdir results');
print "Created results/ directory\n";

#Remove old data if it exists
$sth = $dbh->prepare("drop table if exists illumina_raw_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists barcode_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists illumina_decoded_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists illumina_without_IRDR_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists illumina_nr_pre_map_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists illumina4dec_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists illumina4dec2_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists illumina_blastout_$proj");
$sth->execute;
print "Removed pre-existing tables\n";

#################################################################################################
#1.1  Part 1 load database with seqs.tab data and barcode2lib.txt data

$sth = $dbh->prepare("create table illumina_raw_$proj ( id varchar(20), description varchar(20), sequence varchar(100) ) ");
$sth->execute;

print "Generated illumina_raw_$proj table\n";

$sth = $dbh->prepare("load DATA local INFILE 'data/seqs.tab' INTO TABLE illumina_raw_$proj ");
$sth->execute;

print "1.1 complete - Loaded illumina_raw_$proj data/seqs.tab file \n";

#1.2
$sth = $dbh->prepare("create table barcode_$proj ( seq varchar(20), library varchar(30), direction varcHar(20) )");
$sth->execute;

$sth = $dbh->prepare("load DATA local INFILE 'data/barcode2lib.txt' INTO TABLE barcode_$proj");
$sth->execute;

print "1.2 complete - Loaded barcode_$proj from data/barcode2lib.txt\n";

#1.3
$sth = $dbh->prepare("select count(*) from illumina_raw_$proj");
$sth->execute;

print "finished lines 129-130\n";

while ((@row) = $sth->fetchrow_array) {
       print $row[0]."\tNumber of lines in illumina_raw_$proj file\n";
}

print "finished lines 134-136\n";

$sth = $dbh->prepare("select count(distinct id) from illumina_raw_$proj");
$sth->execute;

print "finished lines 140-141\n";

while ((@row) = $sth->fetchrow_array) {
       print $row[0]."\tNumber of unique ids in illumina raw\n";
}

print "finished lines 145-146\n";

$sth = $dbh->prepare("select count(*) from barcode_$proj");
$sth->execute;

print "finished lines 151-152\n";

while ((@row) = $sth->fetchrow_array) {
       print $row[0]."\tNumber of lines in barcode_$proj\n";
}

print "1.3 complete - Calculated size of barcode and sequence files\n";

#################################################################################################
#2.1 map the barcodes to the sequences using subroutine "resolve barcodes" in config.pl.
resolve_barcodes();
print "2.1 complete - Resolved barcodes using config.pl subroutine and Mapped barcodes to sequences in table illumina_decoded_$proj\n";

#2.2 Print out the list of unique libraries and the number of sequences that match the barcode
$sth = $dbh->prepare("select distinct library, count(id) from illumina_decoded_$proj group by library;");
$sth->execute;

print "Library\t# of inserts\n";
while ((@row) = $sth->fetchrow_array) {
       print $row[0]."\t".$row[1]."\n";
}
print "2.2 complete - Printed list of libraries and number of insertions\n";

#2.3 Resolve IRDRs using the config.pl subroutine resolve_IRDR.
resolve_IRDR();
print "2.3 complete - Mapped and removed IRDR sequences, created table illumina_without_IRDR_$proj\n"; 

#################################################################################################
# 3.1 Find and Remove internal linker

$sth = $dbh->prepare("update illumina_without_IRDR_$proj set type = 'idR'  where insertion_sequence like '%GTCCCTTAAGCGGAGCC%'");
$sth->execute;
$sth = $dbh->prepare("update illumina_without_IRDR_$proj set insertion_sequence =  substring_index(insertion_sequence,'GTCCCTTAAGCGGAGCC',1) where type = 'idR';");

print "3.1 complete - Identified Right linker in remaining seqeunces and renamed to idR.\n";

$sth = $dbh->prepare("update illumina_without_IRDR_$proj set type = 'idL'  where insertion_sequence like '%TAGTCCCTTAAGCGGAG%'");
$sth->execute;
$sth = $dbh->prepare("update illumina_without_IRDR_$proj set insertion_sequence =  substring_index(insertion_sequence,'TAGTCCCTTAAGCGGAG',1) where type = 'idL';");
$sth->execute;
print "3.2 complete - Identified Left linker in remaining seqeunces and renamed to idL.\n";


# Part 3.3 Find and Remove GGATCC sequences
$sth = $dbh->prepare("update illumina_without_IRDR_$proj set type = 'bam'  where insertion_sequence like '_____GGATCC%'");
$sth->execute;
print "3.3 complete - Identified seqeunce with GGATCC and renamed.\n";

# Part 3.4 Find and remove any remaining seqeunces that do not start with TA seqeunce
$sth = $dbh->prepare("update illumina_without_IRDR_$proj set type = 'noTA'  where insertion_sequence not like 'TA%'");
$sth->execute;
print "3.4 complete - Identified seqeunce without TA and renamed to noTA\n";

# 3.5 Find and identify any sequences that are two short for meaningful mapping
$sth = $dbh->prepare("update illumina_without_IRDR_$proj set type = 'min'  where (type = 'good' or type = 'idR' or type = 'idL') and char_length(insertion_sequence)<24");

$sth->execute;
print "3.5 complete - Identified short seqeunces and label them min\n";

#################################################################################################
#4.1 
$sth = $dbh->prepare("alter table illumina_without_IRDR_$proj add index (type)");
$sth->execute;

#necessary variable type alteration
$sth = $dbh->prepare("alter table illumina_without_IRDR_$proj modify insertion_sequence varchar(290)"); 
$sth->execute;

$sth = $dbh->prepare("alter table illumina_without_IRDR_$proj add index (insertion_sequence)"); 
$sth->execute;

$sth = $dbh->prepare("alter table illumina_without_IRDR_$proj add index (library)");
$sth->execute;
print "4.1 complete - indexed table illumina_without_IRDR\n";

#4.2 Remove duplicate sequences maintaining count of all

$sth = $dbh->prepare("create table illumina_nr_pre_map_$proj select distinct count(type) as number, insertion_sequence, library from illumina_without_IRDR_$proj where (type = 'good' or type = 'idR' or type = 'idL')  group by library, insertion_sequence");
$sth->execute;

print "4.2 complete - Created illumina_nr_pre_map table\n";

#4.3
$sth = $dbh->prepare("create table illumina4dec_$proj (id INT unsigned not null auto_increment, primary key(ID), number varchar(20), library varchar(30),sequence varchar(100) )");
$sth->execute;

print "4.3 complete - Created illumina4dec table\n";

$sth = $dbh->prepare("insert into illumina4dec_$proj select NULL,number,library, insertion_sequence from illumina_nr_pre_map_$proj");
$sth->execute;

print "4.4 complete - populated nr_id table\n";

$sth = $dbh->prepare("create table illumina4dec2_$proj select * from illumina4dec_$proj");
$sth->execute;
print "4.5 complete - Copied illumina4dec table to illumina4dec2 for later use\n";

#################################################################################################
#5 This section of the script creates a file for use with bowtie from MySQL tables, executes
#	bowtie, and places the output files in /mapping. Bowtie will be executed five times.

system ('mkdir mapping');

# Note: the following is a test command I used to see if bowtie was working
# system ('bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -q /home/starrt2/shared/genomes/e_coli/bowtieindex/e_coli /home/starrt2/starrt2/tapdance/projects/test1/mapping/e_coli_1000.fq mapping/out.txt');

# 5.1 bowtie iteration #1
$sth = $dbh->prepare("select concat(id,':',library,':',number,'-',char_length(sequence)),sequence from illumina4dec_$proj where length(sequence) > 33");
$sth->execute;

open OUT, "> mapping/4map_1.txt";
while ((@row) = $sth->fetchrow_array) {
       print OUT ">$row[0]\n$row[1]\n";
}
close OUT;
print "5.1 Created first mapping file (4map_1.txt) and copied to /mapping directory\n";
print "5.1 Started bowtie mapping iteration #1\n";
#need to specify this relative to bowtie position
my $cmd = "bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_1.txt mapping/mapping1.txt";
system($cmd);
#system ('bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_1.txt mapping/mapping1.txt');
print "5.1 Finished bowtie iteration #1\n";

# Load bowtie iteration #1 into the MySQL database

$sth = $dbh->prepare("create table illumina_blastout_$proj ( id varchar(50),strand varchar(20),chromo varchar(20),start varchar(20),mismatch varchar(30))");
$sth->execute;
print "5.1 Generated illumina_blastout table\n";

$sth = $dbh->prepare("load DATA local INFILE 'mapping/mapping1.txt' INTO TABLE illumina_blastout_$proj ");
$sth->execute;

$sth = $dbh->prepare("drop table if exists temp");
$sth->execute;

$sth = $dbh->prepare("create table temp select substring_index(B.id,':',1) as id, substring_index(substring_index(B.id,'-',-2),'-',1) as number,substring_index(B.id,'-',-1) as size, strand, chromo, start, mismatch from illumina_blastout_$proj B; ");
$sth->execute;

$sth = $dbh->prepare("delete illumina4dec_$proj from illumina4dec_$proj,temp where illumina4dec_$proj.id = temp.id ");
$sth->execute;
print "5.1 - complete\n";

# 5.2 bowtie iteration #2

$sth = $dbh->prepare("select concat(id,':',library,':',number,'-','33'),left(sequence,33) from illumina4dec_$proj where length(sequence) > 32");
$sth->execute;

open OUT, "> mapping/4map_2.txt";
while ((@row) = $sth->fetchrow_array) {
       print OUT ">$row[0]\n$row[1]\n";
}
close OUT;

print "5.2 Created second mapping file (4map_2.txt) and copied to /mapping directory\n";
print "5.2 Started bowtie mapping iteration #2\n";
#need to specify this relative to bowtie position
my $cmd = "bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_2.txt mapping/mapping2.txt";
system($cmd);
#system ('bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_2.txt mapping/mapping2.txt');
print "5.2 Finished bowtie iteration #2\n";

$sth = $dbh->prepare("load DATA local INFILE 'mapping/mapping2.txt' INTO TABLE illumina_blastout_$proj ");
$sth->execute;

$sth = $dbh->prepare("drop table temp ");
$sth->execute;

$sth = $dbh->prepare("create table temp select substring_index(B.id,':',1) as id, substring_index(substring_index(B.id,'-',-2),'-',1) as number,substring_index(B.id,'-',-1) as size, strand, chromo, start, mismatch from illumina_blastout_$proj B; ");
$sth->execute;

$sth = $dbh->prepare("delete illumina4dec_$proj from illumina4dec_$proj,temp where illumina4dec_$proj.id = temp.id ");
$sth->execute;
print "5.2 - complete\n";

# 5.3 bowtie iteration #3 
$sth = $dbh->prepare("select concat(id,':',library,':',number,'-','30'),left(sequence,30) from illumina4dec_$proj where length(sequence) > 29");

$sth->execute;

open OUT, "> mapping/4map_3.txt";
while ((@row) = $sth->fetchrow_array) {
       print OUT ">$row[0]\n$row[1]\n";
}
close OUT;

print "5.3 Created third mapping file (4map_3.txt) and copied to /mapping directory\n";
print "5.3 Started bowtie mapping iteration #3\n";
#need to specify this relative to bowtie position
my $cmd = "bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_3.txt mapping/mapping3.txt";
system($cmd);
#system ('bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_3.txt mapping/mapping3.txt');
print "5.3 Finished bowtie iteration #3\n";

$sth = $dbh->prepare("load DATA local INFILE 'mapping/mapping3.txt' INTO TABLE illumina_blastout_$proj ");
$sth->execute;

$sth = $dbh->prepare("drop table temp ");
$sth->execute;

$sth = $dbh->prepare("create table temp select substring_index(B.id,':',1) as id, substring_index(substring_index(B.id,'-',-2),'-',1) as number,substring_index(B.id,'-',-1) as size, strand, chromo, start, mismatch from illumina_blastout_$proj B; ");
$sth->execute;

$sth = $dbh->prepare("delete illumina4dec_$proj from illumina4dec_$proj,temp where illumina4dec_$proj.id = temp.id ");
$sth->execute;

print "5.3 - complete\n";

# 5.4 bowtie iteration #4

$sth = $dbh->prepare("select concat(id,':',library,':',number,'-','28'),left(sequence,28) from illumina4dec_$proj where length(sequence) > 27");
$sth->execute;

open OUT, "> mapping/4map_4.txt";
while ((@row) = $sth->fetchrow_array) {
       print OUT ">$row[0]\n$row[1]\n";
}
close OUT;
print "5.4 Created fourth mapping file (4map_4.txt) and copied to /mapping directory\n";
print "5.4 Started bowtie mapping iteration #4\n";
#need to specify this relative to bowtie position
my $cmd = "bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_4.txt mapping/mapping4.txt";
system($cmd);
#system ('bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_4.txt mapping/mapping4.txt');
print "5.4 Finished bowtie iteration #4\n";

$sth = $dbh->prepare("load DATA local INFILE 'mapping/mapping4.txt' INTO TABLE illumina_blastout_$proj ");
$sth->execute;

$sth = $dbh->prepare("drop table temp ");
$sth->execute;

$sth = $dbh->prepare("create table temp select substring_index(B.id,':',1) as id, substring_index(substring_index(B.id,'-',-2),'-',1) as number,substring_index(B.id,'-',-1) as size, strand, chromo, start, mismatch from illumina_blastout_$proj B; ");
$sth->execute;

$sth = $dbh->prepare("delete illumina4dec_$proj from illumina4dec_$proj,temp where illumina4dec_$proj.id = temp.id ");
$sth->execute;

print "5.4 - complete\n";

# 5.5 bowtie iteration #5

$sth = $dbh->prepare("select concat(id,':',library,':',number,'-','24'),left(sequence,24) from illumina4dec_$proj where length(sequence) > 23");
$sth->execute;

open OUT, "> mapping/4map_5.txt";
while ((@row) = $sth->fetchrow_array) {
       print OUT ">$row[0]\n$row[1]\n";
}
close OUT;
print "5.5 Created fifth mapping file (4map_5.txt) and copied to /mapping directory\n";
print "5.5 Started bowtie mapping iteration #5\n";
#need to specify this relative to bowtie position
my $cmd = "bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_5.txt mapping/mapping5.txt";
system($cmd);
#system ('bowtie -t -a --best --strata -v 3 -m 1 --suppress 5,6,7 -f $bowtie_path mapping/4map_5.txt mapping/mapping5.txt');
print "5.5 Finished bowtie iteration #5\n";

$sth = $dbh->prepare("load DATA local INFILE 'mapping/mapping5.txt' INTO TABLE illumina_blastout_$proj ");
$sth->execute;

$sth = $dbh->prepare("drop table if exists bowtie_$proj ");
$sth->execute;
$sth = $dbh->prepare("drop table if exists bowtie_lib_$proj ");
$sth->execute;

print "5.5 - complete\n";

###################################################################################################
# 6 Reduce to sequences to Non-Redundant (NR) set and load library and location mapping.

$sth = $dbh->prepare("create table bowtie_$proj select substring_index(B.id,':',1) as id, substring_index(substring_index(B.id,':',-2),':',1) as library, substring_index(substring_index(B.id,':',-1),'-',1) as number,substring_index(B.id,'-',-1) as size, strand, chromo, start, mismatch from illumina_blastout_$proj B");
$sth->execute;

$sth = $dbh->prepare("alter table bowtie_$proj add pos varchar(20)");
$sth->execute;

$sth = $dbh->prepare("alter table bowtie_$proj add stop varchar(20)");
$sth->execute;

$sth = $dbh->prepare("update bowtie_$proj set pos = size + start where strand = '-'");
$sth->execute;

$sth = $dbh->prepare("update bowtie_$proj set pos = start + 1 where strand = '+'");
$sth->execute;

$sth = $dbh->prepare("update bowtie_$proj set stop = start + size");
$sth->execute;

print "6.1 6.2 6.3 6.4 - created bowtie table and adjusted \n";

#6.5 determine transposon orientation (is it in the A or B orientation?)

$sth = $dbh->prepare("alter table bowtie_$proj add orient varchar(5)");
$sth->execute;

$sth = $dbh->prepare("update bowtie_$proj set orient = '-' where library like '%-L' and strand = '+'");
$sth->execute;

$sth = $dbh->prepare("update bowtie_$proj set orient = '+' where library like '%-L' and strand = '-'");
$sth->execute;

$sth = $dbh->prepare("update bowtie_$proj set orient = '-' where library like '%-R' and strand = '-'");
$sth->execute;

$sth = $dbh->prepare("update bowtie_$proj set orient = '+' where library like '%-R' and strand = '+'");
$sth->execute;

print "6.5 complete - determined mapping orientation\n";

#6.6

$sth = $dbh->prepare("select count(*) from temp");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print "6.6 - ".$row[0]."   Number of unique seqeunces that mapped\n";
}

print "6.6 complete\n";

#6.7

$sth = $dbh->prepare("select count(distinct chromo, start) from temp order by chromo, start+0; ");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print "6.7 - ".$row[0]."  Number of distinct genomic regions that were mapped to\n";
$nr_count = $row[0];
}

print "6.7 complete\n";

###################################################################################################
# 7.1 Create and print .BED files
$sth = $dbh->prepare("create table bowtie_lib_$proj select distinct library,chromo, start,stop, pos, sum(number) as count,strand, orient from bowtie_$proj group by library,chromo,start,stop,pos,strand order by chromo,start+0,orient");
$sth->execute;

print "7.1 - created bowtie_lib table\n";

$sth = $dbh->prepare("select chromo, start, stop, concat(library,'-',count),count,strand,start,stop from bowtie_lib_$proj");
$sth->execute;

# 7.2

open OUT, "> results/raw_$proj.BED";
print OUT "track type ='BED' name='raw$proj' description='all$proj' visibility=2 itemRgb='On'";
while ((@row) = $sth->fetchrow_array) {
 print OUT "\n$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]";}
close OUT;

# 7.3
$sth = $dbh->prepare("select chromo, start, stop, concat(library,'-',count),count,strand,start,stop from bowtie_lib_$proj where count > 9");
$sth->execute;

# select chromo, start, start+1, concat(library,count),count,'+',start,start+1, hex from bowtie_lib_$proj A, convert_2_hex B where A.orient = B.orient and count > 9");
# $sth->execute;

open OUT, "> results/raw10_$proj.BED";
print OUT "track type ='BED' name='raw10$proj' description='10$proj' visibility=2 itemRgb='On'";
while ((@row) = $sth->fetchrow_array) {
 print OUT "\n$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]";}
close OUT;

print "7.1, 7.2 complete - placed raw_project.BED files in results directory\n";

###################################################################################################
# 8.1 Mapping

$sth = $dbh->prepare("drop table if exists illumina_hist_$proj");
$sth->execute;

$sth = $dbh->prepare("drop table if exists illumina_hist_raw_$proj");
$sth->execute;

$sth = $dbh->prepare("create table illumina_hist_raw_$proj select library, chromo,round(pos,-2) as start, count,orient from bowtie_lib_$proj;");
$sth->execute;

$sth = $dbh->prepare("update illumina_hist_raw_$proj set library = left(library,LENGTH(library)-2) where library like '%-L'");
$sth->execute;
$sth = $dbh->prepare("update illumina_hist_raw_$proj set library = left(library,LENGTH(library)-2) where library like '%-R'");
$sth->execute;

$sth = $dbh->prepare("create table illumina_hist_$proj select distinct library, chromo,start, sum(count) as count,orient from illumina_hist_raw_$proj group by library,chromo,start,orient;");
$sth->execute;

$sth = $dbh->prepare("drop table if exists illumina_hist_raw_$proj");
$sth->execute;
$sth = $dbh->prepare("drop table if exists lib_count_$proj");
$sth->execute;

#8.2

$sth = $dbh->prepare("create table lib_count_$proj select distinct library, sum(count) as total from illumina_hist_$proj group by library;");
$sth->execute;

#8.3 the following was added to improve process. Rather than using the total mapped as the denominator in the nr ratio, the total mappable should be used. This will minimize the effect of the crossmapping problem by increasing the requirements for inclusion into the library nr lists sfor librarys wiht low percentage mappiung such as woud exist form mapping to the wrong species

$sth = $dbh->prepare("drop table if exists lib_m_$proj");
$sth->execute;
$sth = $dbh->prepare("create table lib_m_$proj select distinct library, sum(number) as total from illumina4dec2_$proj group by library");
$sth->execute;
$sth = $dbh->prepare("update lib_m_$proj set library = left(library,LENGTH(library)-2) where library like '%-L'");
$sth->execute;
$sth = $dbh->prepare("update lib_m_$proj set library = left(library,LENGTH(library)-2) where library like '%-R'");
$sth->execute;
$sth = $dbh->prepare("drop table if exists lib_mappable_$proj");
$sth->execute;
$sth = $dbh->prepare("create table lib_mappable_$proj select library, sum(total) as total from lib_m_$proj group by library");
$sth->execute;

print "8 complete - mapping is finished!\n";

###################################################################################################
#9 Generate a report describing the mapping.

open OUT, "> results/summary_$proj.txt";
print OUT "Summary of project $proj\n\nA.Report of project based counts \n";

$sth = $dbh->prepare("select count(*) from illumina_raw_$proj");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of Seqeunces\n";
}

$sth = $dbh->prepare("select count(*) from illumina_decoded_$proj;");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of Sequences that map to a barcode\n";
}

$sth = $dbh->prepare("select count(*) from illumina_without_IRDR_$proj");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of Sequences with a barcode that also have an IRDR\n";
}

$sth = $dbh->prepare("select count(*) from illumina_without_IRDR_$proj where type = 'idR';");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of Ideal Right linker sequences\n";
}

$sth = $dbh->prepare("select count(*) from illumina_without_IRDR_$proj where type = 'idL';");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of Ideal left linkersequences\n";
}

$sth = $dbh->prepare("select count(*) from illumina_without_IRDR_$proj where type = 'good';");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of good sequences\n";
}

$sth = $dbh->prepare("select count(*) from illumina_without_IRDR_$proj where type = 'bam';");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of Bamh1 sequence Artifact (removed)\n";
}

$sth = $dbh->prepare("select count(*) from illumina_without_IRDR_$proj where type = 'noTA';");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of sequences removed because of lack of TA sequence\n";
}

$sth = $dbh->prepare("select count(*) from illumina_without_IRDR_$proj where type = 'min';");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal number of sequences removed because of sequence length < 24\n";
}     

$sth = $dbh->prepare("select count(*) from illumina4dec2_$proj;");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
      print OUT $row[0]."\tTotal number of unique seqeunces pryor to mapping\n";
}

$sth = $dbh->prepare("select sum(number) from bowtie_$proj;");
$sth->execute;
while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tTotal SUM  of seqeunces that mapped\n";
}

$sth = $dbh->prepare("select count(*) from bowtie_$proj;");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tNumber of unique sequences that mapped\n";
}

$sth = $dbh->prepare("select count(*) from bowtie_lib_$proj;");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tNumber of distinct locations mapped\n";
}

$sth = $dbh->prepare("select count(*) from illumina_hist_$proj;");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       print OUT $row[0]."\tNumber of distinct regions mapped\n";
}
close OUT;

# report part B

open OUT, ">> results/summary_$proj.txt";
print OUT"\n\nB. Library counts associated with project $proj\n";

%sequence='';
%barcode_count='';
%IRDR_count='';
%IRDR_good='';
%IRDR_unique='';
%map_count='';
%map_total='';
%nr_count='';
@library='';

$sth = $dbh->prepare("select * from barcode_$proj");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $sequence{$row[1]} = $row[0];
       push(@library,$row[1]);
}

$sth = $dbh->prepare("select library, count(id) from illumina_decoded_$proj group by library");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $barcode_count{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select library, count(id) from illumina_without_IRDR_$proj group by library");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $IRDR_count{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select library, sum(number) from illumina4dec2_$proj group by library;");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $IRDR_good{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select library, count(number) from illumina4dec2_$proj group by library;");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $IRDR_unique{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select library, sum(number) from bowtie_$proj group by library");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $map_count{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select library, count(id) from bowtie_$proj group by library");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $map_total{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select library, count(chromo) from bowtie_lib_$proj group by library");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $map_nr{$row[0]} = $row[1];
}

print OUT "Library\tBarcode Sequence\tBarcode\tIRDR\tMappable\tUnique mappable\tTotal map\tUnique map\tNR map";
foreach $item (@library) {
print OUT "$item\t$sequence{$item}\t$barcode_count{$item}\t$IRDR_count{$item}\t$IRDR_good{$item}\t$IRDR_unique{$item}\t$map_count{$item}\t$map_total{$item}\t$map_nr{$item}\n";
}
close OUT;

open OUT, ">> results/summary_$proj.txt";
print OUT "\n\nC. Regions associated with project $proj\n";

%mappable='';
%name='';
%count01='';
%count001='';
%count0001='';
%count0='';
@library='';

$sth = $dbh->prepare("select * from lib_count_$proj");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $name{$row[0]} = $row[1];
       push(@library,$row[0]);
}

$sth = $dbh->prepare("select * from lib_mappable_$proj");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $mappable{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select A.library, count(chromo) from illumina_hist_$proj A, lib_mappable_$proj B where A.library = B.library and count/total > .01  group by A.library order by chromo,start+0");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $count01{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select A.library, count(chromo) from illumina_hist_$proj A, lib_mappable_$proj B where A.library = B.library and count/total > .001  group by A.library order by chromo,start+0");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $count001{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select A.library, count(chromo) from illumina_hist_$proj A, lib_mappable_$proj B where A.library = B.library and count/total > .0001  group by A.library order by chromo,start+0");
$sth->execute;

while ((@row) = $sth->fetchrow_array) {
       $count0001{$row[0]} = $row[1];
}

$sth = $dbh->prepare("select A.library, count(chromo) from illumina_hist_$proj A, lib_mappable_$proj B where A.library = B.library and count/total > 0  group by A.library order by chromo,start+0");
$sth->execute;

print OUT "Library\tTotal mappable\ttotal map\tLocations >0.01\tLocations >0.001\tLocations >0.0001\t All mapped Locations";
while ((@row) = $sth->fetchrow_array) {
       $count0{$row[0]} = $row[1];
}

foreach $item (@library) {
print OUT "$item\t$mappable{$item}\t$name{$item}\t$count01{$item}\t$count001{$item}\t$count0001{$item}\t$count0{$item}\n";
}
close OUT;

print "9 complete - all results and summaries are placed in correct locations\n";

###############################################################################################
# calculate elapsed time
my $end_time = time;
@time=localtime($start_time);

$local_sec = $time[0];
$local_min = $time[1];
$local_hr = $time[2];
$local_day_of_month = $time[3];
$local_mon = $time[4]+1;
$local_yr = $time[5]+1900;
$local_weekday = $time[6];

# To print time as 7/15/2013 14:15:32
print "TAPDANCE v. TKS1 ended:\t$local_mon/$local_day_of_month/$local_yr\t$local_hr:$local_min:$local_sec\n";
print "Elapsed time in hours\t",($end_time - $start_time)/3600,"\n";
###############################################################################################

