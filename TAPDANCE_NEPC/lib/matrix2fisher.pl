#!/usr/bin/perl

# This script is called by TAP2.pl

# This script takes a matrix as input and generates an R script to look for 
# simularity between rows.
#
# 'tab2list' is a complement to list2tab. It reads a 2-dimensional table and
# produces a list in the form: row label - column label - value.
#
use Getopt::Long;

$keycols = 1;
$usage = 0;
@cis;
%hash;
open OUT, "> FISH/Fisher_pre.txt";
open OUT2, "> FISH/Fisher_pre_named.txt";
GetOptions("keycols=i" => \$keycols, "usage" => \$usage, "help"=>\$usage);

if ($usage) {
  print STDERR "usage: tab2list < infile > outfile
  Options:
    --keycols=X  ... How many cols from left are to be kept (default: 1)\n";
  exit 0;
}

$keycols--;

$_ = <>;
chomp;
@headline = split /\t/, $_;
@anim = @headline;
shift(@anim);
while (<>) {
  chomp;
  @line = split /\t/, $_;
  $head = "";
  foreach my $i (0 .. $keycols) {
    $head .= "$line[$i]\t";
  }
  chop $head;
push (@cis,$head);

  foreach my $i ($keycols+1 .. $#line) {
 #   print "$head\t$headline[$i]\t$line[$i]\n";
 $value=$head."|".$headline[$i];
 $hash{$value}=$line[$i];

}
}
foreach (@cis){
#print "$_";
}

print "\n";
foreach (@anim){
#print "$_";
}

print "\n";
foreach (keys %hash) {
#        print "The key $_ $hash{$_}\n";
}

%seen;
$first=0;
$second=0;
$third=0;
$fourth=0;

foreach $x(@cis){
foreach $y(@cis){
$first=0;
$second=0;
$third=0;
$fourth=0;
if ($x eq $y){
next;}
$z = $x."|".$y;
if (exists $seen{$z}){
next;}
$seen{$z}= '1';
$z2 = $y."|".$x;
$seen{$z2}= '1';
foreach $a(@anim) {
$b =$x."|".$a;
$c =$y."|".$a;
if (($hash{$b}==0) and ($hash{$c}==0)){
$fourth = $fourth + 1;
}
elsif (($hash{$b}==0) and ($hash{$c}==1)){
$third = $third + 1;
}
elsif (($hash{$b}==1) and ($hash{$c}==0)){
$second = $second + 1;
}
elsif (($hash{$b}==1) and ($hash{$c}==1)){
$first = $first + 1;
}
#print "$a|$x|$y\n";
#print "$hash{$b}|$hash{$c}\n";
}
print OUT "object__object\t$first\t$second\t$third\t$fourth\n";
print OUT2 "$x"."__"."$y\t$first\t$second\t$third\t$fourth\n";
}
}
close OUT;
close OUT2;
open SOURCE, "< FISH/Fisher_pre.txt";
open OUT, "> FISH/Fisher.R";
while (defined($line = <SOURCE>)) {
chomp $line;
    @field= split(/\t/, $line);
print OUT "$field[0] <- matrix(c($field[1],$field[2],$field[3],$field[4]), nr=2, dimnames=list(c(\"cis\", \"nocis\"), c(\"pheno\", \"notpheno\")))\n$field[0]\nfisher.test($field[0]".",alternative = \"greater\")\n";
}
#system "module load R";
#system "R < FISH/Fisher.R --vanilla > FISH/fisher_dump.txt";
open (my $fisher_dump, ">", "FISH/fisher_dump.txt") || die "Unable to open FISH/fisher_dump.txt. $!\n";
open (my $out, "Rscript --vanilla FISH/Fisher.R | ");
while (<$out>) {
    print $fisher_dump $_;
}
close($out);
close($fisher_dump);
close OUT;

open SOURCE, "< FISH/fisher_dump.txt";
open OUT, "> FISH/R_result.txt";
while (defined($line = <SOURCE>)) {
   chomp ($line); 
   @field= split(/\t/, $line);
           if ($field[0] =~ m/data/){
           print OUT "$field[0]\t";
           }
           if ($field[0] =~ m/p-value/){
           print OUT "$field[0]\n";
           }
   }
close OUT;
print "cooperation analyses complete";

