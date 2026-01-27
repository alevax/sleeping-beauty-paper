# config.pl 
# Last updated 10/26/18 by Tim Starr

# This perl script is called by multiple scripts within the lib folder.
# This script will set up the database interface, set the table name suffixes,
#	set the percentages for CIS calculations, and set the annotation file
# The subroutine "resolve_barcodes" and "resolve_IRDR" will create tables that 
#	identify and remove barcodes and IRDR sequences from the input sequences

# 10/26/18: The database currently being used is TKS MySQL database at MSI

#######################################################################################
# Set up the database handler 
use DBI;
# Note: the DBI is a database access module for the Perl programming language. 
#	It defines a set of methods, variables, and conventions that provide a 
#	consistent database interface, independent of the actual database being used.
# 	I'm pretty sure the mysql driver DBD::MYSQL must also be installed based on

 $dbh;
 
 $db_name     = 'tapdance_database';
 $db_username = 'tapdance_username';
 $db_password = 'tapdance_password';
 $db_host     = 'localhost;mysql_local_infile=1';
 $db_type     = 'mysql';

 $data_source = "DBI:$db_type:database=$db_name;host=$db_host";

$dbh = DBI->connect($data_source, $db_username, $db_password,
{ RaiseError => 1, AutoCommit => 0 });

# I think the RaiseError and AutoCommit allows multiple active instances to be 
#	running in the db, but I'm not sure about this.

#######################################################################################
# Set the project name suffix and the annotation file
$proj = 'prostate';
# TAPDANCE will create ~18 tables with this suffix for the table name
$annotation_file = 'lib/mm10.bed';
# Set the bowtie directory path
$bowtie_path = '/Users/matteo/Desktop/bowtie-1.2.3/indexes/mm10/mm10';

#######################################################################################
# Set the stringency levels for including insertions and calculating CIS's

# The library_percent threshold will result in discarding mapped sequences where the
#	read count for that insertion is lower than the library_percent based on total
#	read counts assigned to that library.
$library_percent ='0.01';

# The CIS and cocis thresholds will require the pvalue to be below the pvalue listed
#	to be reported as a CIS.
$CIS_total_pvalue = '0.05';
$CIS_library_pvalue = '0.05';
$CIS_region_pvalue = '0.05';

$cocis_threshold ='0.001';

#################################################################################################

sub resolve_barcodes {
# This is used at "2.1" in the TAPDANCE.pl script
# Change the sequence length (sequence, #) to length of barcode +1. For example, if your barcode
#	is 12 bases, enter (sequence,13) in the line below
$sth = $dbh->prepare("create table illumina_decoded_$proj select library, id, substring(sequence,13) as decoded_sequence from barcode_$proj,illumina_raw_$proj where sequence like concat(seq,'%')");
$sth->execute;
  }

sub resolve_IRDR {
# This is used at "2.3" in the TAPDANCE.pl script
# Adjust the sequences in the statements below and run the program to determine the optimal search string
#	The optimal search string will find the most matches, without too many false matches.

# For example, the following two statements work for my myeloid project. Notice the first search string
#	only differs from the second search string by "...AAA..." vs. "...AA..."
#	Also, the underscores represent wild-cards.  In this instance there are 5 wildcards.
#	The (decoded_sequence,32) must match the length of the search string, which includes the wild-cards
#	and the % sign. (i.e. everything within the single quotes).

$sth = $dbh->prepare("create table illumina_without_IRDR_$proj select library,id,substring(decoded_sequence,30) as insertion_sequence, 'good' as type from illumina_decoded_$proj where decoded_sequence like '___TGTATGTAAACTTCCGACTTCAACTG%'");
$sth->execute;

}
