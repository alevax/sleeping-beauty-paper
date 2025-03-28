# This script is called by TAP2.pl

# This script takes a set of chromosomes and positions from a text file and generates defined CIS'es
# version 5.0 ALS 
# Nov 15 2010

#  Open the presorted NR chromosome file and load it into a hash.
require 'config.pl';
$boxSize = $ARGV[0];

open OUT, "> CIS/cis$ARGV[0].txt";
open SOURCE, "< CIS/nr.txt";
while (defined($line = <SOURCE>)) {
$count++; 
chomp $line;
@field= split(/\t/, $line);
$hash[$count]= $field[1];
$chrom[$count]= $field[0];
$lib[$count]= $field[2];
$orient[$count]=$field[5];
}

#Count the number of nr positions within a box the size of which is defined below and add to a parallel hash.
foreach $n (1..$count) {
$cis_count = 0;
@lib_list = '';
@pos_list = '';
$lib_count = '';
$orient_count = '0';
$i =$n;
#print "$n\n";
$j = $hash["$n"]+$boxSize;
#print "pos=>".$hash[$n]."\n";
while (($hash["$n"] <= $hash["$i"]) and ($j > ($hash["$i"]+0)) and ($i <= $count)) {
push(@lib_list,$lib["$i"]);
push(@pos_list,$hash["$i"]);
if ($orient[$i] eq "+") {
$orient_count = $orient_count +1;
}

$i=$i+1;
#print "$i\n";
$cis_count=$cis_count+1;
#print "$cis_count\n";

}
$cis[$n]=$cis_count;
$strand[$n]=$orient_count;
#print "Cis_count=$cis_count\t";
#print @lib_list;
#This section rather hackily counts the number of libraries represented in the CIS window.
%seen=();
@uniq = ();
$libstring = '';
foreach $item (@lib_list) {
unless ($seen{$item}) {
$seen{$item}= 1;
push(@uniq,$item)
}
}
#print "\n@uniq\n";

$lib_count = -1;
foreach $i (@uniq) {
$lib_count++;
$libstring = $libstring.$i.'::';
}

#print "Library_count=$lib_count";

$lib[$n]=$lib_count;
$lib_name[$n]=$libstring;
#print "$cis[$n]\t$lib[$n]";

#This section rather hackily counts the the number of regions represented in the CIS window.
%seen=();
@uniq = ();
$posstring = '';
foreach $item (@pos_list) {
unless ($seen{$item}) {
$seen{$item}= 1;
push(@uniq,$item)
}
}
#print "\n@uniq\n";

$pos_count = -1;
foreach $i (@uniq) {
$pos_count++;
$posstring = $posstring.$i.'::';
}

#print "Library_count=$lib_count";

$pos[$n]=$pos_count;
$pos_name[$n]=$posstring;
#print "\t$pos[$n]\n";



}


$cis[0] = 0;

# select peaks and export them with complex logic scheme.
# update Sept22 seperate this into 4 cases.

foreach $n (1..$count) {
$prev = $n-1;
$this = $n;
$next = $n+1;
#print "\t$ARGV[1]\n";


if ($ARGV[1] eq library) {
#print "$ARGV[1]\t$cis[$this]\t$lib[$n]\n";
$cis[$this]=$lib[$this];
$cis[$next]=$lib[$next];
$cis[$prev]=$lib[$prev];
}
# outcome 1 no other cises are relevant i.e. nothing within + or - cis size.
if ((($hash[$this] - $hash[$prev] > $boxSize) or ($chrom[$this] ne $chrom[$prev]))  and (($hash[$next]-$hash[$this] > $boxSize) or ($chrom[$next] ne $chrom[$this])) and ($cis[$this] > 1)){ 
outputpart1();
}

# outcome 2 no other cis + cis size cis - cis exists and is smaller.
if (($hash[$this] - $hash[$prev] >0) and ($hash[$this] - $hash[$prev] <= $boxSize) and ($cis[$prev] < $cis[$this]) and ($hash[$next]-$hash[$this] > $boxSize) and ($cis[$this] > 1)){ 
outputpart1();
}

# outcome 3 no other cis - cis size + cis exists and is smaller are relevant nothing within + or - cis size.
if ((($hash[$this] - $hash[$prev] > $boxSize) or ($chrom[$this] ne $chrom[$prev]))  and ($hash[$next]-$hash[$this] >= 0)  and ($hash[$next]-$hash[$this] <= $boxSize) and ($cis[$this] >= $cis[$next]) and ($cis[$this] > 1)){ 
outputpart1();
}

# outcome 4 other cises on both sides.
if (($hash[$this] - $hash[$prev] >0) and ($hash[$this] - $hash[$prev] <= $boxSize)  and ($hash[$next]-$hash[$this] > 0) and ($hash[$next]-$hash[$this] <= $boxSize) and ($cis[$prev] < $cis[$this]) and ($cis[$this] >= $cis[$next]) and $cis[$this] > 1 ) {
outputpart1();
}
#print "$chrom[$n]\t$hash[$n]\t$cis[$n]\n";
}
#hack to prevent permanent looping...observedon one occasion unsure why...
$iter = 0;
$count3 = 1;
$count2 = 0;
while (($count3 ne $count2) and ($iter < 100)){
print "$iter";
$iter++;
close OUT;
#print "$count3\t\t\t\t$count2\n\n";
$count3 = $count2;
$count2 = 0;
%hash= ();
%chrom= ();
%cis=();
%libraries=();
%full_line=();
#open OUT, "> cis2.txt";
open SOURCE, "< CIS/cis$ARGV[0].txt";
while (defined($line = <SOURCE>)) {
$count2++; 
chomp $line;
$full_line[$count2]=$line;
@field= split(/\t/, $line);
$hash[$count2]= $field[1];
$chrom[$count2]= $field[0];
$cis[$count2]= $field[4];
#$libraries[$count2]=$field[11];
}
close SOURCE;
open OUT, "> CIS/cis$ARGV[0].txt";
# select peaks and export them with complex logic scheme.
# update Sept22 seperate this into 4 cases.

foreach $n (1..$count2) {
$prev = $n-1;
$this = $n;
$next = $n+1;

# outcome 1 no other cises are relevant i.e. nothing within + or - cis size.
if ((($hash[$this] - $hash[$prev] > $boxSize) or ($chrom[$this] ne $chrom[$prev]))  and (($hash[$this]-$hash[$next] < -$boxSize) or ($chrom[$next] ne $chrom[$this])) and ($cis[$this] > 1)){ 
output();
}

# outcome 2 no other cis + cis size cis - cis exists and is smaller.
if (($hash[$this] - $hash[$prev] >0) and ($hash[$this] - $hash[$prev] <= $boxSize) and ($cis[$prev] < $cis[$this]) and ($hash[$next]-$hash[$this] > $boxSize) and ($cis[$this] > 1)){ 
output();
}

# outcome 3 no other cis - cis size + cis exists and is smaller are relevant nothing within + or - cis size.
if ((($hash[$this] - $hash[$prev] > $boxSize) or ($chrom[$this] ne $chrom[$prev]))  and ($hash[$next]-$hash[$this] >= 0) and ($hash[$next]-$hash[$this] <= $boxSize) and ($cis[$this] >= $cis[$next]) and ($cis[$this] > 1)){ 
output();
}

# outcome 4 other cises on both sides.
if (($hash[$this] - $hash[$prev] >0) and ($hash[$this] - $hash[$prev] <= $boxSize) and ($hash[$next]-$hash[$this] >= 0) and ($hash[$next]-$hash[$this] <= $boxSize) and ($cis[$prev] < $cis[$this]) and ($cis[$this] >= $cis[$next]) and $cis[$this] > 1 ) {
output();
}
#print "$chrom[$n]\t$hash[$n]\t$cis[$n]\n";
}

}
close OUT;

#Subroutines
# provide data in the form of poisson('BoxSize','Total','count')
sub outputpart1 {
$pval = poissonBonf($boxSize,$count,$lib[$this]);
$pvalTotal = poissonBonf($boxSize,$count,$cis[$this]);
$pvalRegion = poissonBonf($boxSize,$count,$pos[$this]);
if ($pval < $CIS_total_pvalue) {
$end = $hash[$n]+$boxSize;
print OUT "$chrom[$n]\t$hash[$n]\t$end\t$chrom[$n]_$hash[$n]\t$lib[$this]\t$boxSize\t$pval\t$cis[$this]\t$pvalTotal\t$pos[$this]\t$pvalRegion\t$lib_name[$n]\t$strand[$this]\n";           
}
}

sub output {
print OUT "$full_line[$n]\n";
#if ($pval < 0.05) {
#$end = $hash[$n]+$boxSize;
#print OUT "$chrom[$n]\t$hash[$n]\t$end\t$chrom[$n]_$hash[$n]\t$cis[$this]\t$boxSize\t$pval\t$bonf\t$libraries[$n]\n";
#}
}
sub poisson {
return ((($_[1]/(2393726033/$_[0]))**($_[2]))*(2.71828**-($_[1]/(2393726033/$_[0])))/(fac($_[2])));
}

sub poissonBonf {
return (((($_[1]/(2393726033/$_[0]))**($_[2]))*(2.71828**-($_[1]/(2393726033/$_[0])))/(fac($_[2])))*$_[1]);
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


