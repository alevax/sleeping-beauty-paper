# This script is called by TAP2.pl and by TAP4.pl

use strict;
use List::Util qw[min max];

my $all = "ALL";
my $close = "CLOSEST";
my $first_idx = 0;
my $left = -1;
my $right = 1;
my $both = 0;
my $keys = "KEYS";
my $NOT_FOUND = -3;
my $NOT_FOUND_STR = "NOT FOUND";
my $ENCOMPASSING = -2;
#my $ENCOMPASSING_STR = "";
my $ENCOMPASSING_STR = "ENCOMPASSING FEATURE";
my $INTERNAL = -1;
#my $INTERNAL_STR = "";
my $INTERNAL_STR = "INTERNAL FEATURE";

return 1;

sub feature_index {
    my ($feature_file, $window, $fchrom_col, $fstart_col, $fend_col, $debug) = @_;
    open(FEATURES, "<${$feature_file}") || die "Unable to open features file, ${$feature_file}: $!\n";
    my %feature_hash = ();
    my $num_indexed = 0;
    while (<FEATURES>) {
	if (${$debug}) {
	    if (($num_indexed % 1000) == 0) {
		my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime(time);
		print sprintf("%02d-%02d %02d:%02d:%02d Features Indexed: %9s\n", ($mon+1), $day, $hour, $min, $sec, $num_indexed);
	    }
	}
	chomp;
	$num_indexed++;
	my @split_line = split("\t");
	if (!@split_line) { next; }
	my $interim_hash1 = int($split_line[${$fstart_col}]/${$window});
	my $interim_hash2 = int($split_line[${$fend_col}]/${$window});
	for (my $interim_hash = min($interim_hash1, $interim_hash2); $interim_hash <= max($interim_hash1, $interim_hash2); $interim_hash++) {
	    push @{$feature_hash{$split_line[${$fchrom_col}]}->{$interim_hash}->{$split_line[${$fstart_col}]}}, \@split_line;
	}
    }
    close(FEATURES);
    if (${$debug}) { print "$num_indexed indexed\n"; }

    my $min_keys_in_window;
    my $max_keys_in_window;
    my $num_idxs = 0;
    my $total_keys = 0;
    foreach my $chrom (keys %feature_hash) {
	foreach my $interim_hash (keys %{$feature_hash{$chrom}}) {
	    $num_idxs++;
	    my @sorted_keys = sort {$a <=> $b} keys %{$feature_hash{$chrom}->{$interim_hash}};
	    my $interim_key = join("_", ($interim_hash, $keys));
	    my $num_keys = $#sorted_keys+1;
	    #print OUT "interim_hash:$interim_hash interim_key:$interim_key num_keys:$num_keys @sorted_keys\n";
	    $feature_hash{$chrom}->{$interim_key} = \@sorted_keys;
	    $total_keys += $#sorted_keys + 1;
	    if (!defined($min_keys_in_window) || ($#sorted_keys+1) < $min_keys_in_window) { $min_keys_in_window = ($#sorted_keys+1); }
	    if (!defined($max_keys_in_window) || ($#sorted_keys+1) > $max_keys_in_window) { $max_keys_in_window = ($#sorted_keys+1); }
	}
    }
    if (${$debug}) {
	my $avg_keys = $total_keys/$num_idxs;
	print "Smallest # of keys in window: $min_keys_in_window, Largest: $max_keys_in_window $total_keys in $num_idxs windows for an average of $avg_keys\n"; 
    }
    
    return \%feature_hash;
}

sub process_intervals {
    my ($intervals_file, $out_file, $feature_hash_ref, $req_res_type_ref, $direction_ref, $window_ref, $html_ref, $cturl_ref, $db_ref, $icols_ref, $fcols_ref, $debug_ref, $dd_ref, $capture, $promoter) = @_;
    if (${$debug_ref}) { print sprintf("process_intervals(intervals:%s, out:%s, direction:%s, window:%s, icols_ref:%s)\n", ${$intervals_file}, ${$out_file}, ${$direction_ref}, ${$window_ref}, ${$icols_ref}); }
    open(INTERVALS, "<${$intervals_file}") || die "Unable to open intervals file, ${$intervals_file}: $!\n";
    open(my $out_h, ">${$out_file}") || die "Unable to open output file, ${$out_file}: $!\n";

    #if (${$debug_ref}) { print sprintf("%s features in index.\n", keys %{$feature_hash_ref}); }

    my @icols = split(":", ${$icols_ref});
    #pod2usage({ message => "Must specify 5 columns with -icol, chrom:start:end:label:strand.\n", exitval => 2}) if ($#icols != 4);  
    my $ichrom = $icols[0]-1;
    my $istart = $icols[1]-1;
    my $iend = $icols[2]-1;
    my $ilabel = $icols[3]-1;
    my $istrand = $icols[4]-1;

    my @fcols = split(":", ${$fcols_ref});
    #pod2usage({ message => "Must specify 5 columns with -fcol, chrom:start:end:label:strand.\n", exitval => 2}) if ($#fcols != 4);  
    my $fchrom = $fcols[0]-1;
    my $fstart = $fcols[1]-1;
    my $fend = $fcols[2]-1;
    my $flabel = $fcols[3]-1;
    my $fstrand = $fcols[4]-1;

    my $cturl;
    my $print = \&print_res;
    if (${$html_ref}) {
	$print = \&print_res_html;
	print $out_h "<HTML>\n<HEAD>\n\t<TITLE>Feature Finder Annotation</TITLE>\n</HEAD>\n<BODY>\n<TABLE>\n";
	if (defined(${$cturl_ref})) {
	    ${$cturl} =~ s/X26/\%26/g;
	}
    }

    my $count = 0;
    foreach (<INTERVALS>) {
	if (${$debug_ref}) {
	    if (($count % 1000) == 0) {
		my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime(time);
		print sprintf("%02d-%02d %02d:%02d:%02d Intervals Processed: %9s\n", ($mon+1), $day, $hour, $min, $sec, $count);
	    }
	    $count++;
	}
	chomp;
	my @split_line = split("\t");
	if (!@split_line) { next; }
	if (!defined(${$feature_hash_ref}{$split_line[$ichrom]})) { 
	    #if (${$debug_ref}) { print "."; } 
	    next; 
	}
	my $win_start = min($split_line[$istart], $split_line[$iend]) - ${$window_ref};
	if ($win_start < 0) { $win_start = 0; }
	my $start_idx = int($win_start/${$window_ref}) -1;
	if ($start_idx < 0) { $start_idx = 0; }
	my $win_end = max($split_line[$istart], $split_line[$iend]) + ${$window_ref};
	my $end_idx = int($win_end/${$window_ref}) +1;
        
	if (uc ${$req_res_type_ref} eq $all) {
	    &$print(\$out_h, \@split_line, \$istart, \$iend, &compare(\$out_h, $feature_hash_ref, \$all, \@split_line, \$split_line[$ichrom], \$win_start, \$win_end, \$start_idx, \$end_idx, \$both, \$ichrom, \$istart, \$iend, \$ilabel, \$istrand, \$fchrom, \$fstart, \$fend, \$flabel, $debug_ref), \$flabel, \$fstart, \$fend, \$fstrand, $window_ref, $promoter, $capture, $dd_ref, $db_ref, $cturl_ref, \$ichrom, \$istart, \$iend);
	} 
	else {
	    my %left_res = ();
	    my %right_res = ();
	    #if (${$debug_ref}) { print sprintf("STRAND:%s\tDIRECTION:%s\n", $split_line[$istrand], ${$direction_ref}); }
	    if ($split_line[$istrand] ne "+" && $split_line[$istrand] ne "-") { warn sprintf("Found \"%s\" as the strand of %s:%s-%s %s, please verify icols match your input data\n", $split_line[$istrand], $split_line[$ichrom], $split_line[$istart], $split_line[$iend], $split_line[$ilabel]); }
	    if (($split_line[$istrand] eq "+"  && uc ${$direction_ref} ne "DOWNSTREAM") || ($split_line[$istrand] eq "-" && uc ${$direction_ref} ne "UPSTREAM")) {
		%left_res = %{&compare(\$out_h, $feature_hash_ref, \$close, \@split_line, \$split_line[$ichrom], \$win_start, \$win_end, \$start_idx, \$end_idx, \$left, \$ichrom, \$istart, \$iend, \$ilabel, \$istrand, \$fchrom, \$fstart, \$fend, \$flabel, $debug_ref)};
	    }
	    if (($split_line[$istrand] eq "-" && uc ${$direction_ref} ne "DOWNSTREAM") || ($split_line[$istrand] eq "+" && uc ${$direction_ref} ne "UPSTREAM")) {
		%right_res = %{&compare(\$out_h, $feature_hash_ref, \$close, \@split_line, \$split_line[$ichrom], \$win_start, \$win_end, \$start_idx, \$end_idx, \$right, \$ichrom, \$istart, \$iend, \$ilabel, \$istrand, \$fchrom, \$fstart, \$fend, \$flabel, $debug_ref)};
	    }
	    if (keys %left_res == 0) { 
		&$print(\$out_h, \@split_line, \$istart, \$iend, \%right_res, \$flabel, \$fstart, \$fend, \$fstrand, $window_ref, $promoter, $capture, $dd_ref, $db_ref, $cturl_ref, \$ichrom);
	    } 
	    elsif (keys %right_res == 0) {
		&$print(\$out_h, \@split_line, \$istart, \$iend, \%left_res, \$flabel, \$fstart, \$fend, \$fstrand, $window_ref, $promoter, $capture, $dd_ref, $db_ref, $cturl_ref, \$ichrom);
	    }
	    else {
		my @left_keys = keys %left_res;
		my @right_keys = keys %right_res;
		if (uc ${$direction_ref} eq "BOTH") {
		    @left_res{keys %right_res} = values %right_res;
		    &$print(\$out_h, \@split_line, \$istart, \$iend, \%left_res, \$flabel, \$fstart, \$fend, \$fstrand, $window_ref, $promoter, $capture, $dd_ref, $db_ref, $cturl_ref, \$ichrom);
		}
		else {
		    if ($left_keys[0] <= $right_keys[0]) {
			&$print(\$out_h, \@split_line, \$istart, \$iend, \%left_res, \$flabel, \$fstart, \$fend, \$fstrand, $window_ref, $promoter, $capture, $dd_ref, $db_ref, $cturl_ref, \$ichrom);
		    }
		    elsif ($right_keys[0] < $left_keys[0]) {
			&$print(\$out_h, \@split_line, \$istart, \$iend, \%right_res, \$flabel, \$fstart, \$fend, \$fstrand, $window_ref, $promoter, $capture, $dd_ref, $db_ref, $cturl_ref, \$ichrom);
		    }
		}
	    }
	}
    }
    close(INTERVALS);
    if (${$html_ref}) { 
	print $out_h "</TABLE>\n</BODY>\n</HTML>\n";
    }
    close($out_h);
   
}

sub compare {
    my ($out_h_ref, $feature_hash_ref, $res_type, $interval_array_ref, $chrom, $window_start, $window_end, $start_idx, $end_idx, $direction, $ichrom, $istart, $iend, $ilabel, $istrand, $fchrom, $fstart, $fend, $flabel, $debug) = @_;
    #if (${$debug}) { print sprintf("compare(result type:%s, chrom:%s, win start:%s, win end:%s, direction:%s, ichrom_col:%s, fchrom_col:%s, ilabel_col:%s, flabel_col:%s, debug:%s)\n", ${$res_type}, ${$chrom}, ${$window_start}, ${$window_end}, ${$direction}, ${$ichrom}, ${$fchrom}, ${$ilabel}, ${$flabel}, ${$debug}); }
    my %result = ();
    my %reshash = ();
    my $feature_start;
    my $feature_end;
    my $interval_start;
    my $interval_end;
    my $interval_idx;
    my $best_distance = $NOT_FOUND;
    for ($interval_idx = ${$start_idx}; $interval_idx <= ${$end_idx}; $interval_idx++) {
	my $interval_key = join("_", ($interval_idx, $keys));
	if (!defined(${$feature_hash_ref}{${$chrom}}->{$interval_key})) { 
	    next; 
	}
	
	my @feature_keys = @{${$feature_hash_ref}{${$chrom}}->{$interval_key}};
	my $interval_midp = (($interval_array_ref->[${$istart}] + $interval_array_ref->[${$iend}])/2);
	foreach my $feature_idx (@feature_keys) {
	    foreach my $feature (@{${$feature_hash_ref}{${$chrom}}->{$interval_idx}->{$feature_idx}}) {
		$feature_start = min($feature->[${$fstart}], $feature->[${$fend}]);
		$feature_end = max($feature->[${$fstart}], $feature->[${$fend}]);
		$interval_start = min($interval_array_ref->[${$istart}], $interval_array_ref->[${$iend}]);
		$interval_end = max($interval_array_ref->[${$istart}], $interval_array_ref->[${$iend}]);
		my $reskey = join(":", ($interval_array_ref->[${$ichrom}], ${$window_start}, ${$window_end}, $interval_array_ref->[${$ilabel}], $feature->[${$fchrom}], $feature_start, $feature_end, $feature->[${$flabel}]));  
		my $distance;
		if ($feature_start <= $interval_start && $feature_end >= $interval_end) { $distance = $ENCOMPASSING; }
		elsif ($interval_start < $feature_start && $interval_end > $feature_end) { $distance = $INTERNAL; }
		else {
		  $distance = min(abs($interval_array_ref->[${$istart}] - $feature_start), abs($interval_array_ref->[${$istart}] - $feature_end), abs($interval_array_ref->[${$iend}] - $feature_start), abs($interval_array_ref->[${$iend}] - $feature_end));
		  #if (${$debug}) { print "FEATURE:$feature->[${$flabel}]\tDISTANCE:$distance\n"; }
		}	
		my $dir_str = "";
		my $feature_midp = ($feature_start + $feature_end)/2;
		if ($interval_midp >= $feature_midp) {
		    if ($interval_array_ref->[${$istrand}] eq "+") { $dir_str = "UPSTREAM"; }
		    else { $dir_str = "DOWNSTREAM" }
		}
		else {
		    if ($interval_array_ref->[${$istrand}] eq "-") { $dir_str = "UPSTREAM"; }
		    else { $dir_str = "DOWNSTREAM" }
		}

		if (${$res_type} eq $close) {
		    if (${$direction} == $left) {
			if ($interval_midp < $feature_midp) { next; }
		    }
		    elsif (${$direction} == $right) { 
			if ($interval_midp >= $feature_midp) { next; }
		    }

		    #print ${$out_h_ref} "distance:$distance best_distance:$best_distance\n";
		    if ($best_distance > $NOT_FOUND) {
			if ($distance < $best_distance) {
			    delete $result{$best_distance};
			    if (!defined($reshash{$reskey})) {
				$reshash{$reskey} = $feature;
				push @{$result{$distance}{$feature->[${$flabel}]}}, $feature, $dir_str;
				#$result{$distance}{$feature->[${$flabel}]} = $dir_str;
				$best_distance = $distance;
			    }
			} 
		    } 
		    else {
			if (!defined($reshash{$reskey})) {
			    $reshash{$reskey} = $feature;
			    push @{$result{$distance}{$feature->[${$flabel}]}}, $feature, $dir_str;
			    #$result{$distance}{$feature->[${$flabel}]} = $dir_str;
			    $best_distance= $distance;
			}
		    }
		}
		else {
		    # if the interval start - the window is less than the end of the feature and the interval end + the window is greater than the start of the feature it's included
		    if (${$window_start} <= $feature_end && ${$window_end} >= $feature_start) {
			if (!defined($reshash{$reskey})) {
			    $reshash{$reskey} = $feature;
			    push @{$result{$distance}{$feature->[${$flabel}]}}, $feature, $dir_str;
			    #$result{$distance}{$feature->[${$flabel}]} = $dir_str;
			} 
		    }
		}
	    }
	}
    }
    return \%result;
}

sub print_res { 
    my ($out_h_ref, $interval_array_ref, $istart, $iend, $result_hash_ref, $flabel, $fstart, $fend, $fstrand, $window, $promoter_len, $capture, $dd_ref) = @_;
    print ${$out_h_ref} join("\t", @{$interval_array_ref});
    my @sorted_keys = sort {$a <=> $b} keys %{$result_hash_ref};
    if (defined($sorted_keys[0])) {
	print ${$out_h_ref} sprintf("\t");
	my $left_end = $interval_array_ref->[${$istart}]; 
	my $right_end = $interval_array_ref->[${$iend}];
	my $fstart_pos;
	my $fend_pos;
	foreach my $next_key (@sorted_keys) {
	  foreach my $res (keys %{$result_hash_ref->{$next_key}}) {
	    #my @feature_keys = keys %{$res};
	    if (defined($capture) && ${$capture}) {
	      if (@{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstrand}] eq '+') {
		$fstart_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstart}] - ${$promoter_len};
	      } else {
		$fstart_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstart}];
	      }
	      if (@{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstrand}] eq '-') {
		$fend_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fend}] + ${$promoter_len};
	      } else {
		$fend_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fend}];
	      }
	      if ($fstart_pos < $left_end) { $left_end = $fstart_pos; }
	      if ($fend_pos > $right_end) { $right_end = $fend_pos; }
	    }
	    if (!${$dd_ref}) {
	      print ${$out_h_ref} sprintf("%s, ", $res);
	    } else {
	      if ($next_key == $ENCOMPASSING) {
		print ${$out_h_ref} sprintf("%s (%s), ", $res, $ENCOMPASSING_STR);
	      }
	      elsif ($next_key == $INTERNAL) {
		print ${$out_h_ref} sprintf("%s (%s), ", $res, $INTERNAL_STR);
	      }
	      else {
		print ${$out_h_ref} sprintf("%s (%s bp %s), ", $res, $next_key, @{$result_hash_ref->{$next_key}->{$res}}[1]);
	      }		    
	    }
	  }
	}
	if (defined($capture) && ${$capture}) {
	  print ${$out_h_ref} sprintf("\t%s\t%s\t%s", $left_end, $right_end, $right_end - $left_end);
	}
    }
    else {
	print ${$out_h_ref} sprintf("\tNo results within ${$window} bp window");
    }
    print ${$out_h_ref} sprintf("\n");
}

sub print_res_html {  
    my ($out_h_ref, $interval_array_ref, $istart, $iend, $result_hash_ref, $flabel, $fstart, $fend, $fstrand, $window, $promoter_len, $capture, $dd_ref, $db, $cturl, $ichrom) = @_;
    print ${$out_h_ref} sprintf("<TR>\n\t<TD>");
    my $win_start = min($interval_array_ref->[${$istart}], $interval_array_ref->[${$iend}]) - $window;
    my $win_end = max($interval_array_ref->[${$istart}], $interval_array_ref->[${$iend}]) + $window;
    print ${$out_h_ref} sprintf("<A HREF=\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s:%s-%s&display_app=ucsc&authz_method=display_at", $db, $interval_array_ref->[${$ichrom}], $win_start, $win_end);
    my $left_end = $interval_array_ref->[${$istart}];
    my $right_end = $interval_array_ref->[${$iend}];
    my $fstart_pos;
    my $fend_pos;
    if (defined($cturl)) {
	#print ${$out_h_ref} sprintf ("&hgt.customText=%s/display_as?id=%s&display_app=ucsc&authz_method=display_at", $server_url, $out_dbkey);
	print ${$out_h_ref} sprintf("&hgt.customText=%s", $cturl);
    }
    print ${$out_h_ref} sprintf("\">UCSC Visualization</A></TD><TD>");
    #custom track stuff something like : &hgt.customText=%s/display_as?id=%s&display_app=ucsc&authz_method=display_at, $server_url, $out_dbkey
    print ${$out_h_ref} join("</TD><TD>", @{$interval_array_ref});
    my @sorted_keys = sort {$a <=> $b} keys %{$result_hash_ref};
    if (defined($sorted_keys[0])) {
	foreach my $next_key (@sorted_keys) {
	    foreach my $res (keys %{$result_hash_ref->{$next_key}}) {
	      if (${$capture}) {
		if (@{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstrand}] eq '+') {
		  $fstart_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstart}] - ${$promoter_len};
		} else {
		  $fstart_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstart}];
		}
		if (@{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fstrand}] eq '-') {
		  $fend_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fend}] + ${$promoter_len};
		} else {
		  $fend_pos = @{$result_hash_ref->{$next_key}->{$res}}[0]->[${$fend}];
		}
		if ($fstart_pos < $left_end) { $left_end = $fstart_pos; }
		if ($fend_pos > $right_end) { $right_end = $fend_pos; }
	      }
	      
	      if (!${$dd_ref}) {
		print ${$out_h_ref} sprintf("<TD>%s</TD>", $res);
	      } else {
		if ($next_key == $ENCOMPASSING) {
		  print ${$out_h_ref} sprintf("<TD>%s (%s)</TD>", $res, $ENCOMPASSING_STR);
		}
		elsif ($next_key == $INTERNAL) {
		  print ${$out_h_ref} sprintf("<TD>%s (%s)</TD>", $res, $INTERNAL_STR);
		}
		else {
		  print ${$out_h_ref} sprintf("<TD>%s (%s bp %s)</TD>", $res, $next_key, @{$result_hash_ref->{$next_key}->{$res}}[0]);
		}
	      }
	    }
	}
	
	if (${$capture}) {
	  print ${$out_h_ref} sprintf("<TD>%s</TD><TD>%s</TD><TD>%s</TD>", $left_end, $right_end, $right_end - $left_end);
	}
    }
    else {
	print ${$out_h_ref} sprintf("<TD>No results within $window bp window</TD>");
    }
    print ${$out_h_ref} sprintf("</TD>\n</TR>\n");
}
