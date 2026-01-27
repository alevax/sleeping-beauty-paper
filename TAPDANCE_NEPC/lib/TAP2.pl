# Description: This is the TAp2.pl script. 
#	The goal of this is to count the number of CIS present in a NR set
# one file is required nr.txt sorted in chromosome order!

# Excluded chromosomes:
#	This script requires a file "chromo.tab" to be located in data/ directory. The chromosomes
#	should be in list form:
#		chr1
#		chr15
#	If no chromosomes are to be ommitted, use chrN

# Feature file:
#	The correct features file that matches the genome must be in the lib/ directory (e.g. mm9.bed).  
#	The name of the file should be entered in the config.pl file (e.g. @annotation_file = 'lib/mm9.bed').
#	This file is generally created from UCSC Table browser using the RefSeq gene names.
#	Note: You can manually add BADrepeat or BADEN2 regions to the features file
#		These will be considered in the analysis (not sure how this is done).
#	See section 10.2 below to change parameters for using the feature finder table.

# Version history:
# 	Aaron Lyman Sarver Oct 21, 2010
# 	This script takes over following the bowtie mapping and generates CIS.
# 	Update NOV15,2010 to pass library information through 
# 	$Fully modified december 1 2010i
# 	Fully modifed May 1 2010
# 	This version, TAP2_TKS1.pl was created 11/22/13 by Tim Starr

use File::Copy;
require 'config.pl';
require 'lib/feature_finder_methods.pl';

unless(-d "CIS") {
    mkdir("CIS");
}

open CIS, "> results/summary_CIS_$proj.txt";
close CIS;

print "Created CIS file and populated with summary_CIS_project.txt\n";

###################################################################################################
# 10.1 generate required tables
$sth = $dbh->prepare("drop table if exists cis_result_final_$proj");
$sth->execute;
$sth = $dbh->prepare("create table cis_result_final_$proj (chromo varchar(20), start int, stop int, name  varchar(20), number varchar(20), size  varchar(20), pvalue  double, total varchar(20), pvaluetotal double,region varchar(20), pvalueregion double, library_name varchar(5000), strand int, gene_name varchar(1000), type  varchar(20))");
$sth->execute;

$sth = $dbh->prepare("drop table if exists sort_$proj");
$sth->execute;
$sth = $dbh->prepare("create table sort_$proj (chromo varchar(20), start int, stop int, name  varchar(20), number varchar(20), size  varchar(20), pvalue  double, total varchar(20), pvaluetotal double,region varchar(20), pvalueregion double, library_name varchar(5000), strand int, gene_name varchar(1000), type  varchar(20))");
$sth->execute;
$sth = $dbh->prepare("drop table if exists sort2_$proj");
$sth->execute;
$sth = $dbh->prepare("create table sort2_$proj (chromo varchar(20), start int, stop int, name  varchar(20), number varchar(20), size  varchar(20), pvalue  double, total varchar(20), pvaluetotal double,region varchar(20), pvalueregion double, library_name varchar(5000), strand int, gene_name varchar(1000), type  varchar(20))");
$sth->execute;

# 10.2 load concatamer information from file
$sth = $dbh->prepare("drop table if exists chromo_$proj");
$sth->execute;
$sth = $dbh->prepare("create table chromo_$proj (chromo varchar(20) ) ");
$sth->execute;
$sth = $dbh->prepare("load DATA local INFILE 'data/chromo.tab' INTO TABLE chromo_$proj ");
$sth->execute;
$sth = $dbh->prepare("select distinct chromo from chromo_$proj");
$sth->execute;
print "10.2 complete - Loaded concatamer chromo\n";

# 10.3 generate donor exclusion portion of the SQL query
@exclude;
$query_chromo= '';
while ((@row) = $sth->fetchrow_array) {
       push (@exclude,$row[0]);
       print "\n$row[0]";
}
print @exclude, "\n";

foreach $desc(@exclude) {
    $query_chromo = $query_chromo."and chromo != '".$desc."' ";
}
print $query_chromo, "\n";
print "10.3 - loaded chromo.tab with omitted chromosomes\n";

# 10.4 load library metadata
$sth = $dbh->prepare("drop table if exists metadata_$proj");
$sth->execute;
$sth = $dbh->prepare("create table metadata_$proj (library varchar(20), descriptor varchar(20), type varchar (20) ) ");
$sth->execute;
$sth = $dbh->prepare("load DATA local INFILE 'data/metadata.tab' INTO TABLE metadata_$proj ");
$sth->execute;
$sth = $dbh->prepare("select distinct descriptor from metadata_$proj where type = 'cis'");
$sth->execute;
print "10.4 - Loaded metadata\n";

# 10.2 load the feature hash
print "10.2 - Loading features...\n";
my ($feature_file, $window, $fchrom, $fstart, $fend, $debug) = ($annotation_file, 20000, 0, 1, 2, 1);
my $feature_hash_ref = &feature_index(\$feature_file, \$window, \$fchrom, \$fstart, \$fend, \$debug);
print "10.2 - Features loaded.\n";

# 10.3 generate metadata based library portion of the SQL query
@descriptors;
$query_definition = '';
while ((@row) = $sth->fetchrow_array) {
       push (@descriptors,$row[0]);
       print "\n$row[0]";
}
print @descriptors;
foreach $desc(@descriptors) {
    $query_definition = '';
    $query = "select distinct library from metadata_$proj where descriptor = '$desc'";
	#print "\n$query\n";
    $sth = $dbh->prepare("$query");
    $sth->execute;
    
    while ((@row) = $sth->fetchrow_array) {
	$query_definition = $query_definition."or A.library = '$row[0]'  "; 
	#print "$row[0]";
    }
    print "\n$desc\n";
#print "$query_definition\n";
print "10.3 complete\n";


# 11.0 setup the looping feature
    @options = ('');
    foreach $library(@options) {
	@percent =("$library_percent");
	foreach $percent(@percent) {

# 11.1 query to select regions of interest
	    $sth = $dbh->prepare("select chromo, start, A.library,count, count/total,orient  from illumina_hist_$proj A, lib_mappable_$proj B where A.library = B.library and count/total > $percent $query_chromo and (1=2 $query_definition) order by chromo,start+0 ;");
	    $sth->execute;
	    
	    unless(-d "results") {
		mkdir("results");
	    }
	    unless(-d "results/$desc") {
		mkdir("results/$desc");
	    }
	    
	    open OUT, "> CIS/nr.txt";
	    $name = "$desc-nr-$proj-$percent";
	    open OUT2, "> results/$desc/$name.txt";
	    open OUT3, "> results/$desc/$name.BED";
	    print OUT3 "track type ='BED' name='$name' description='$name Insertions' visibility=2 itemRgb='On'";
	    while ((@row) = $sth->fetchrow_array) {
		print OUT "\n$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]";
		$start = $row[1]-50;
		$stop = $row[1]+50;
		print OUT2 "\n$row[0]\t$start\t$stop\t$row[2]\t$row[3]\t$row[4]\t$row[5]";
		print OUT3 "\n$row[0]\t$start\t$stop\t$row[2]\t$row[3]\t$row[5]\t$start\t$stop\t159,0,197";
	    }
	    close OUT;
	    
#11.2 Anotate the inserts from the .txt file
	    my ($interval_file, $out_file, $req_res_type, $direction, $interval_cols, $feature_cols, $display_distance) = ("results/$desc/$name.txt", "results/$desc/$name.ann.txt", "ALL", "BOTH", "1:2:3:4:10", "1:2:3:4:6", 0);
	    &process_intervals(\$interval_file, \$out_file, $feature_hash_ref, \$req_res_type, \$direction, \$window, \$html, \$cturl, \$db, \$interval_cols, \$feature_cols, \$debug, \$display_distance);
	    unlink("results/$desc/$name.txt");
            #generate a list of gense with high density inserts
            #if (($percent == '0.01') and ($library eq 'library')) {
            #system("perl lib/convert_nr_2_annot.pl results/$desc/$name.txt results/$desc/$name.ingenuity.txt");
            #}

	    $count = 0;
	    open SOURCE, "< CIS/nr.txt";
	    while (defined($line = <SOURCE>)) {
		$count++;
	    }
	    @libSizes =();
	    open CIS, ">> results/summary_CIS_$proj.txt";
	    print "\n\n11.3 - $name  Total number of NR regions $count\n";
	    print CIS "\n\n$name\nExcluded chromosomes $query_chromo\nTotal number of NR regions $count\n";
	    print CIS "CIS analyses calculate window sizes based on NR count\n";

#11.4 Ideal window size calculation
	    $n = 10000;
	    $f = 3;
	    while ($n < 300000) {
		$stat =poissonBonf($n,$count,$f);
		if ($stat > 0.05) {
		    $nf=$n-1000;
		    $stat2=poissonBonf($nf,$count,$f);
		    if ($stat2 < 0.05) {
			print "$stat2\t$nf\t$f\n";
			print CIS "$stat2\t$nf\t$f\n";
			unshift(@libSizes, $nf);
		    }
		    $f=$f+1;
		}
		$n = $n+1000;
	    }
	    $nf='301000';
	    unshift(@libSizes, $nf);
	    print @libSizes;
	    print CIS @libsizes;
	    close CIS;
	    if ($count < 2000) {
		@libSizes =('301000','200000','100000','50000','25000','12500');
	    }
	    if ($count > 200000) {
		@libSizes =('301000','200000','100000','50000','25000','12500');
	    }

#11.5 
	    foreach $lib(@libSizes) {
		
#11.5 call the cis_finder.pl script to calculate the CIS.
		system("perl lib/cis_finder.pl $lib");
		print "\n$lib\n";

		
#11.6 call the feature_finder.pl script to annotat the CIS regions
		my ($interval_file, $out_file, $req_res_type, $direction, $interval_cols, $feature_cols, $display_distance) = ("CIS/cis$lib.txt", "CIS/final_cis_named_$lib.txt", "ALL", "BOTH", "1:2:3:4:10", "1:2:3:4:6", 0); 
		&process_intervals(\$interval_file, \$out_file, $feature_hash_ref, \$req_res_type, \$direction, \$window, \$html, \$cturl, \$db, \$interval_cols, \$feature_cols, \$debug, \$display_distance);
	    }
	    
		#load the results into the database.and resolve overlap
	    $sth = $dbh->prepare("delete from sort_$proj");
	    $sth->execute;
	    $sth = $dbh->prepare("delete from sort2_$proj");
	    $sth->execute;
	    
#11.7
	    foreach $lib(@libSizes) {
		$sth = $dbh->prepare("load DATA local INFILE 'CIS/final_cis_named_$lib.txt' INTO TABLE sort2_$proj");
		$sth->execute;
		resolve_cis();
}

print "11 complete\n";

# Move and rename CIS to cis_result_final table 
	    $filename = $desc."_".$percent.$library;
	    $query2 = "delete from cis_result_final_$proj where type= '$filename'";
	    print "\n$query2\n";
	    $sth = $dbh->prepare("$query2");
	    $sth->execute;
	    $sth = $dbh->prepare("insert into cis_result_final_$proj select * from sort_$proj");
	    $sth->execute;
	    $query3 = "update cis_result_final_$proj set type = '$filename' where type is NULL";
	    print "$query3";
	    $sth = $dbh->prepare("$query3");
	    $sth->execute;
	    
#Export CIS list to the results subdirectory
	    $query4 = "select number,CONCAT(chromo,':',start,'-',stop)as pos, total, pvaluetotal,number, pvalue,region,pvalueregion,gene_name,library_name,strand from cis_result_final_$proj where type = '$filename'  and pvalueregion < '$CIS_region_pvalue' and pvaluetotal < '$CIS_total_pvalue' and pvalue < '$CIS_library_pvalue' and gene_name not like '%BAD%' order by pvalueregion;";
	    $sth = $dbh->prepare("$query4");
	    $sth->execute;
	    open CIS, "> results/$desc/cis_$name.xls";
	    print CIS "$filename\nnumber\tpos\t#inserts\tpvalueinsert#\t#library\tpvaluelibrary#\t#regions\tpvalueregion#\tgene_name\tlibrary_name\tnumber of inserts drive transcription on positive strand";
	    while ((@row) = $sth->fetchrow_array) {
		print CIS "\n$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\t";
	    }
	    close CIS;
	    
	    
	    open SOURCE, "< results/$desc/cis_$name.xls";
	    open CIS, "> results/$desc/plot_$name.wig";
	    print CIS "track type=wiggle_0 name=\"$proj.$filename\" description=\"$proj.$filename\" visibility=full autoScale=off";
	    while (defined($line = <SOURCE>)) {
		chomp $line;
		@field= split(/\t/, $line);
		@chromo=split(/:|-/, $field[1]);
		
		if ($field[1] =~ m/chr/) {
		    if ($field[7] < 1.0e-100) {
			$field[7] = 1.0e-100;}
		    
		    $log= -(log($field[7])/log(10));
		    print CIS "\n$chromo[0]\t$chromo[1]\t$chromo[2]\t$log";
		}
	    }
	    close CIS;

#11.11 Generate the cis summary file of top 25 library CIS
	    $query4 = "select number,CONCAT(chromo,':',start,'-',stop)as pos, total, pvaluetotal,number, pvalue,region,pvalueregion,gene_name,library_name,strand from cis_result_final_$proj where type = '$filename'  and pvalueregion < '$CIS_region_pvalue' and pvaluetotal < '$CIS_total_pvalue' and pvalue < '$CIS_library_pvalue' and gene_name not like '%BAD%' order by pvalueregion;"; 
	    print "$query4";
	    $sth = $dbh->prepare("$query4");
	    $sth->execute;
	    open CIS, ">> results/summary_CIS_$proj.txt";
	    print CIS "\n\n$filename\nnumber\tpos\t#inserts\tpvalueinsert#\t#library\tpvaluelibrary#\t#regions\tpvalueregion#\tgene_name\tlibrary_name\t+strand\t";
	    while ((@row) = $sth->fetchrow_array) {        
		print CIS "\n$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\t";
	    }
	    close CIS;

#generate matrix and carry out association analyses
	    if (($desc eq 'all') and ($percent == 0.0001) ) {
		print "\n\n Starting Association analyses for $name";
		mkdir "FISH";
		mkdir "results/Assoc";
		system ('perl lib/visualize.pl');
		copy("FISH/cis4cluster.txt", "results/Assoc/cis4cluster.txt") or die "Unable to copy file: $!\n";
		#system ('module load R');
		system ('perl lib/matrix2fisher.pl FISH/cis4cluster.txt');
		copy("FISH/R_result.txt", "results/Assoc/R_result.txt") or die "Unable to copy file: $!\n";
		copy("FISH/Fisher_pre_named.txt", "results/Assoc/Fisher_pre_named.txt") or die "Unable to copy file: $!\n";
		system ('perl lib/Fisher_clean.pl');
	    }

	}#close the threshhold loop
    }#close the libraries loop
}#close the metadata loop

print "\nTAP2.pl complete - Yeehaw!\n";

exit(0);

sub resolve_cis {
    $sth = $dbh->prepare("delete from A using sort_$proj A, sort2_$proj B where A.chromo = B.chromo and A.start <= B.start and A.stop >= B.start and A.pvalue >  B.pvalue");
    $sth->execute;
    
    $sth = $dbh->prepare("delete from A using sort_$proj A, sort2_$proj B where A.chromo = B.chromo and A.start <= B.stop and A.stop >= B.stop and A.pvalue >  B.pvalue");
    $sth->execute;
    $sth = $dbh->prepare("delete from B using sort_$proj A, sort2_$proj B where A.chromo = B.chromo and A.start <= B.start and A.stop >= B.start and A.pvalue <  B.pvalue; ");
    $sth->execute;
    $sth = $dbh->prepare("delete from B using sort_$proj A, sort2_$proj B where A.chromo = B.chromo and A.start <= B.stop and A.stop >= B.stop and A.pvalue <=
  B.pvalue; ");
    $sth->execute;
    
    $sth = $dbh->prepare("insert into sort_$proj select * from sort2_$proj");
    $sth->execute;
    
    $sth = $dbh->prepare("delete from sort2_$proj");
    $sth->execute;
}

sub poisson {
    return ((($_[1]/(2393726033/$_[0]))**($_[2]))*(2.71828**-($_[1]/(2393726033/$_[0
								     ])))/(fac($_[2])));
}

sub poissonBonf {
    return (((($_[1]/(2393726033/$_[0]))**($_[2]))*(2.71828**-($_[1]/(2393726033/$_[
									  0])))/(fac($_[2])))*$_[1]);
}

sub fac {
    my ($n) = @_;
    
    if ($n < 2) {
	return $n;
    }
    else {
	return $n * fac($n-1);
    }
}

