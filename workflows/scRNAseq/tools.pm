use File::Temp qw/ tempfile tempdir /;
use List::Util 'shuffle';

sub s_ascii2phred {
    # returns Phred quality of basecall in Solexa/Ilumina reads given ASCII symbol
    return 10 * log(1 + 10 ** ((ord($_[0]) - 64) / 10.0)) / log(10);
}

sub ascii2phred {
    # returns Phred quality of basecall given ASCII symbol
    return ord($_[0]) - 33;
}

sub phred2errprob {
    # returns error probability of a given Phred value
    return exp(-$_[0] * log(10)/10.0);
}

sub ascii2errprob {
    # returns error probability given ASCII symbol for basecall
    return phred2errprob(ascii2phred($_[0]));
}

sub s_ascii2errprob {
    # returns error probability given ASCII symbol for Solexa/Ilumina basecall
    return phred2errprob(s_ascii2phred($_[0]));
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
#    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}

sub max {
 my $a = shift;
 my $b = shift;
 return $a > $b ? $a : $b;
}

sub min {
 my $a = shift;
 my $b = shift;
 return $a < $b ? $a : $b;
}

sub hash_append {
    # receives a pointer to hash as first argument
    my $hash = shift;
    my $key = shift;
    my $el = shift;
    if (!exists($$hash{$key})){
	@{$$hash{$key}} = ($el);
    }else{
	push( @{$$hash{$key}},$el);
    }
    return;
}

sub array_append {
    my $array = shift;
    my $el = shift;
    if (scalar @$array == -1){
	@$array = ($el);
    }else{
	push(@$array, $el);
    }
} 

sub fill_histo {
    my ( $histo, $key ) = @_;
    if (!exists($$histo{$key})){
	$$histo{$key} = 1;
    }else{
	$$histo{$key} ++;
    }
}

sub print_num_histo {
    # receives a pointer to histo-hash as argument
    $histo = $_[0];
    $flag = 0;
    if ($#_ == 1){
	$out = $_[1];
	$flag = 1;
	open(OUT,">",$out);
    }
    foreach $k (sort {$a <=> $b} keys %$histo){
	if ($flag){
	    print OUT $k."\t".$$histo{$k}."\n";
	}else{
	    print $k."\t".$$histo{$k}."\n";
	}
    }
    if ($flag){
	close(OUT);
    }
    return;
}

sub print_hash {
    # receives a pointer to hash as argument
    $hash = $_[0];
    $flag = 0;
    if ($#_ == 1){
	$out = $_[1];
	$flag = 1;
	open(OUT,">",$out);
    }
    foreach $k (sort keys %$hash){
	if ($flag){
	    print OUT $k."\t".$$hash{$k}."\n";
	}else{
	    print $k."\t".$$hash{$k}."\n";
	}
    }
    if ($flag){
	close(OUT);
    }
    return;
}
sub fasta2hash {
    # receives a pointer to a hash as argument                                                                                                                                            
    $hash = $_[0];
    $file = $_[1];
    open(IN,"<",$file);
    $j = 0;
    print STDERR $file."\n";
    while(<IN>){
        chomp;
	next if $_ eq "";
	$j ++;
        print STDERR $j."\r";
        if (/^\>/){
            $key = $_;
            $key =~ s/^\>//g;
            $$hash{$key} = "";
        }else{
            $$hash{$key} = $$hash{$key}.$_;
        }
    }
    close(IN);
}

sub fasta2hash_old {
    # receives a pointer to a hash as argument
    $hash = $_[0];
    $file = $_[1];
    ($fh, $filename) = tempfile();
    open(IN,"<",$file);
    $j = 0;
    print STDERR $file."\n";
    open(TMP,">",$filename);
    while(<IN>){
        chomp;
	$j ++;
	print STDERR $j."\r";
	if (/^\>/){
	    print TMP $_."\n" if $j == 1;
	    print TMP "\n".$_."\n" if $j > 1;
	}else{
	    print TMP $_;
	}
    }
    print TMP "\n";
    close(TMP);
    close(IN);

    open(IN,"<",$filename);
    while(<IN>){
	chomp;
	if (/^\>/){
	    $_ =~ s/^\>+//;
	    $key = $_;
	    next;
	}
	if ( exists($$hash{$key}) ){
	    $$hash{$key} = $_ if length($_) > length($$hash{$key});
	}else{
	    $$hash{$key} = $_;
	}
    }
    close(IN);
    system("rm ".$filename)
}

sub mafref2hash {
    # receives a pointer to a hash as argument
    $hash = $_[0];
    $file = $_[1];

    open(IN,"<",$file);
    while(<IN>){
	chomp;
	if (/^\>\>/){
	    $key = $_;
	    $key =~ s/^\>\>//;

	}
	if ( !exists($$hash{$key}) ){
	    $$hash{$key} = $_."\n";
	}else{
	    $$hash{$key} = $$hash{$key}.$_."\n";
	}
    }
    close(IN);
}

sub revcompl {
    $seq = shift;
    $seq =~ tr/ACGTacgt/TGCAtgca/ ;
    return reverse($seq);
}

sub log2 {
    my $n = shift;
    return (log($n)/log(2));
}

sub calc_rpkm {
    my $count  = shift;
    my $length = shift;
    my $mapped_reads = shift;
    return "0" if $length * $mapped_reads == 0;
    return $count/$length/$mapped_reads * 1000 * 1000000;
}
sub clean_array {
    my $array = shift;
    my @clean_array = ();
    for $a (@{$array}){
	array_append(\@clean_array, $a) if $a ne "NA";
    }
    @{$array} = @clean_array;
}
sub average {
    @_ == 1 or die ('Sub usage: $average = average(\@array);');
    my ($array_ref) = @_;
    my @clean_array = @{$array_ref};
    clean_array(\@clean_array);
    my $sum;
    my $count = scalar @clean_array;
    foreach (@clean_array) { $sum += $_; }
    return $sum / $count;
}

sub median {
    @_ == 1 or die ('Sub usage: $median = median(\@array);');
    my ($array_ref) = @_;
    my @clean_array = @{$array_ref};
    clean_array(\@clean_array);
    my $count = scalar @clean_array;
    # Sort a COPY of the array, leaving the original untouched
    my @array = sort { $a <=> $b } @clean_array;
    if ($count % 2) {
	return $array[int($count/2)];
    } else {
	return ($array[$count/2] + $array[$count/2 - 1]) / 2;
    }
} 
sub average_trimmed {
    my $array_ref = shift;
    my $alpha = shift;
    my @clean_array = @{$array_ref};
    clean_array(\@clean_array);
    my $count = scalar @clean_array;
    my $sum;
    $cut = int($alpha * $count);
    my @array = sort { $a <=> $b } @clean_array;
    foreach (@array[$cut...$#array - $cut]){ $sum += $_; }
    return $sum/($count - 2*$alpha);
} 

sub get_name_gff {
    my @F     = split(/\t/,shift);
    my $check = shift;
# %check is initialized like that
#    $parent = "transcript"; parent record
#    $child  = "exon"; child record
#    $check{$parent} = "ID"; find name of parent
#    $check{$child} = "Parent"; find name of parent
    my @G     = split(/\;/,$F[8]);
    for $i (0..$#G){
#	print STDERR $G[$i]."\t".$F[2]."\t".$$check{$F[2]}."\n";
	if ($G[$i] =~ /^($$check{$F[2]})/){
	    my @g = split(/\=/,$G[$i]);
	    my @n = split(/\:/,$g[1]);
	    return $n[1] if $g[1] =~ /\:/;
	    return $n[0] if $g[1] !~ /\:/;
	}
    }
}
sub get_name_gff_SPALN {
    my @F     = split(/\t/,shift);
    my $check = shift;
# %check is initialized like that
#    $parent = "transcript"; parent record
#    $child  = "exon"; child record
#    $check{$parent} = "ID"; find name of parent
#    $check{$child} = "Parent"; find name of parent
    my @G     = split(/\;/,$F[8]);
    for $i (0..$#G){
	if ($G[$i] =~ /^($$check{$F[2]})/){
	    my @g = split(/\=/,$G[$i]);
	    my @n = split(/\:/,$g[1]);
	    return $n[1] if $g[1] =~ /\:/;
	    return $n[0] if $g[1] !~ /\:/;
	}
    }
}

sub get_entry_gtf {
    my @F     = split(/\t/,shift);
    my $entry = shift;
    my @G     = split(/\;/,$F[8]);
    foreach $g (@G){
	if ($g =~ /($entry)/){
	    #print STDERR "here 1:\t".$g."\n";
	    $g =~ s/($entry)//g;
	    $g =~ s/[\s\"\n]//g;
	    #print STDERR "here 2:\t".$g."\n";
	    #die $F[0]."\n";
	    return $g;
	}
    }

}

sub get_entry_gff {
    my @F     = split(/\t/,shift);
    my $entry = shift;
    my @G     = split(/\;/,$F[8]);
    foreach $g (@G){
	if ($g =~ /($entry)/){
	  ($dum,$name) = split(/\=/,$g);
	  $name =~ s/[\s\"\n]//g;
	  return $name;
	}
    }

}

sub hash_file_1 {
    my $hash = shift;
    my $file = shift;
    open(TMP,"<",$file);
    while(<TMP>){
	chomp;
	$$hash{$_} = 1;
    }
    close(TMP);
}

sub cut_string {
    my $str = shift;
    my $c   = shift;
    my $s   = 0;
    my @a   = ();
    my $flag = 1;
    my $pos = 0;
    while ($flag) {
	if (length(substr($str,$pos, length($str) - $pos)) > $c){
	    push(@a, substr($str,$pos,$c));
	    $pos = $pos + $c;
	}else{
	    push(@a, substr($str,$pos,length($str) - $pos));
	    $flag = 0;
	}
    }
    return @a 
}

sub round{
    my $n = shift;
    my $digits = shift;
    $rounded = sprintf("%.".$digits."f",$n );
    return $rounded;
}


sub seeded_shuffle {
    my @F = @_;
    my $old_seed = rand(2**32);
    srand($F[$#F]);
    my @shuffled = shuffle(@F[0..($#F - 1)]);
    srand($old_seed);
    return @shuffled;
}


sub random_array {
  my $n = shift;
  my $l = shift;
  my $u = shift;
  my %seen = ();
  my @a = ();
  my $i = 0;
  my $random = ();
  while( $i < $n ){
    $random = int( rand( $u-$l+1 ) ) + $l;
    if (!exists($seen{$random})){
      push(@a,$random);
      $seen{$random} = 1;
      $i ++
    }
  }
  return @a;
}


sub makedir {
    my $directory = shift;
    
    unless(-e $directory or mkdir $directory) {
        die "Unable to create $directoryn";
    }
}

1;


