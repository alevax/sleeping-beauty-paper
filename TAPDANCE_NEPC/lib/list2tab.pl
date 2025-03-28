#!/usr/bin/perl

# This script is called by visualize.pl and Fisher_clean.pl

# The Missing Textutils, Ondrej Bojar, obo@cuni.cz
# http://www.cuni.cz/~obo/textutils
#
# 'list2tab' builds a 2D table out of key-key-value triples.
#
# $Id: list2tab,v 1.4 2006/07/07 15:01:53 bojar Exp $
#

$rows = shift;
$cols = shift;
$data = shift;
$blankvalue = shift;
@rows = split /,/, $rows;
@cols = split /,/, $cols;
@data = split /,/, $data;

$blankvalue = "-" if $blankvalue eq "";

if (!$rows || !$cols || !$data) {
  my $help = <<EOH;
Sample usage:
  ./list2tab.pl 1,2 5,6 3,4 [default_value]  < datafile > tablefile
The output table will have lines labelled with values seen in columns 1 and 2,
columns labelled with values from columns 5,6 and the values in the interior of the table will come from columns 3,4.

Sample input:
GIN	Praha    	5
IOL	Praha    	20
GIN	Brno     	10
IOL	Nova Paka	2

Output produced by: "list2tab 2 1 3 none"

         	GIN 	IOL
Brno     	10  	none
Nova Paka	none	2
Praha    	5   	20
EOH
  print STDERR $help;
  exit $1;
}


while (<>) {
  chomp;
  @line = split /\t/;
  $key = $val = $datum = "";
  foreach $row (@rows) {
    $key .= $line[$row-1]."\t";
  }
  chop $key;
  foreach $col (@cols) {
    $val .= $line[$col-1]."\t";
  }
  chop $val;
  foreach $dat (@data) {
    $datum .= $line[$dat-1]."\t";
  }
  chop $datum;
  $datum =~ s/^ *//;
  $datum =~ s/ *$//;

  $table{"$key\t$val"} = $datum;
#print STDERR "Tabulka >>$key\t$val<<.....$datum\n";
  $keys{$key} = 1;
  push @svals, $val if (!$vals{$val});
  $vals{$val} = 1;
}

@skeys = sort {uc($a) cmp uc($b)} keys %keys;
#@svals = sort {uc($a) cmp uc($b)} keys %vals;

#print STDERR "Klice: ".join(",",@skeys)."\n";
#print STDERR "Hodnoty: ".join(",",@svals)."\n";


print "\t" x $#rows;

foreach $col (@svals) {
  print "\t$col";
}
print "\n";

foreach $row (@skeys) {
  print "$row";
  foreach $col (@svals) {
    $pos = "$row\t$col";
    print "\t$table{$pos}" if defined $table{$pos};
    print "\t$blankvalue" if ! defined $table{$pos};
  }
  print "\n";
}
