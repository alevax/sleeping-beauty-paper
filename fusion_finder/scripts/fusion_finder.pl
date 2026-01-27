#!/usr/bin/env perl
##################################################
#this is a script to identify improperly paired fusions given one end of the fusion
#Aaron L. Sarver NOV 21, 2013
#Sarver@umn.edu
#Version 3.1
#
# Modified for sleeping_beauty project structure (Dec 2025)
# - Uses relative paths from project root
# - Outputs to results/fusions/
#
###Requirements###################################
#input.txt
#loc_SB.txt
#samtools
#starts with a bam file indexed makes a list of RNA based insertions for each file and a bed file of all insertions
##################################################

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

# Get project root (parent of scripts/)
my $script_dir = dirname(abs_path($0));
my $project_dir = dirname($script_dir);
chdir($project_dir) or die "Cannot chdir to $project_dir: $!";

# Directories
my $results_dir = "results/fusions";
my $working_dir = "$results_dir/working";

# Create output directories
system("mkdir -p $results_dir");
system("mkdir -p $working_dir");

open INFO, "< input.txt" or die "Cannot open input.txt: $!";
while (defined(my $data = <INFO>)) {
    print $data;
    chomp $data;
    my @binfo = split(/\t/, $data);
    my $dir = $binfo[0];
    my $name = $binfo[1];
    my $sam = $binfo[2];

    open OUT2, "> $results_dir/RNA_insertion_$name.$sam.txt" or die "Cannot open output: $!";

    open INFO2, "< $sam.txt" or die "Cannot open $sam.txt: $!";
    while (defined(my $data2 = <INFO2>)) {
        print $data2;
        chomp $data2;
        my @loc = split(/\t/, $data2);
        my $loc1 = $loc[1];
        my $locname = $loc[0];

        system "samtools view $dir $loc1 -o $working_dir/location_raw";

        `cut -f7,8 $working_dir/location_raw|sort -k1,1 -k2n,2 > $working_dir/loc.txt`;

        open SOURCE, "< $working_dir/loc.txt";
        open OUT, "> $working_dir/nr.txt";
        my $count = 0;
        while (defined(my $line = <SOURCE>)) {
            chomp $line;
            my @field = split(/\t|\s+/, $line);
            $count++;
            print "$count\n";
            if ($field[0] ne "=") {
                if ($field[0] ne "*") {
                    my $round = 300*int(0.5+$field[1]/300);
                    print OUT "$field[0]\t$round\n";
                }
            }
        }
        close OUT;
        close SOURCE;

        `uniq -c $working_dir/nr.txt > $working_dir/cis1.txt`;

        open SOURCE, "< $working_dir/cis1.txt";
        while (defined(my $line = <SOURCE>)) {
            chomp $line;
            my @field = split(/\t|\s+/, $line);
            if (defined($field[1]) && $field[1] > 2) {
                my $end = $field[3]+150;
                my $start = $field[3]-150;
                my $val = `samtools view -c $dir $field[2]:$start-$end`;
                chomp $val;
                if ($val > 0 && $field[1]/$val > 0.01) {
                    $end = $field[3]+150;
                    $start = $field[3]-150;
                    print OUT2 "$field[2]\t$start\t$end\t$name\t$field[1]\t$field[1]\t$val\t$locname\n";
                }
            }
        }
        close SOURCE;
    }
    close INFO2;
    close OUT2;
}
close INFO;

`cat $results_dir/RNA* > $results_dir/fusions_all.txt`;
print "\nDone! Results in $results_dir/fusions_all.txt\n";
