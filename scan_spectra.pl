#!/usr/bin/perl
use Number::Closest;
use List::MoreUtils qw(first_index);
use List::Util qw(sum);

# Read clump l,b,v data
my @clump_names = read_col("Clumps_l_b_v.csv",1); 
my @clump_velo = read_col("Clumps_l_b_v.csv",4); 

my %velo; my $i = 0;
foreach my $name (@clump_names) {
	$velo{$name} = $clump_velo[$i];
	#print "velo{$name} = $clump_velo[$i]\n";
	$i++;
}

my $dir = "_FCRAO_12CO_spectra";

open(FH, '>', "$dir"."_Tpeak.csv") or die $!;
print FH "Clump,Tpeak\n";



my $clump = 0;
for ($clump=1; $clump<=145; $clump=$clump+1) {
	my $sn = sprintf("%03d",$clump);
	my $fn = "$dir/Clump$sn"."\.dat";
	my $name = "Clump$sn";
	if (not -f $fn) {
		print "$name,\n";
		print FH "$name,\n";
		next; 	
	};
	
	my $GC_velo = $velo{$name};
	
	# Read spectra
	my @vdata = read_col($fn,3); # читаем данные скорости
	my @idata = read_col($fn,4); # читаем данные интенсивности


	my $finder = Number::Closest->new(number => $GC_velo, numbers => \@vdata) ;
	

	# Get 3 closest points
	my $closest_three = $finder->find(3) ;
	my @closest = @$closest_three;
	my @peaks_closest = ();
	foreach my $clos (@closest) {
		my $ind = first_index { $_ == $clos } @vdata;
		my $peak_closest = sprintf("%5.2f",$idata[$ind]);
		#print "close value $clos\n";
		push(@peaks_closest,$peak_closest);
	}

	#print join(", ", @peaks_closest);
	#print "\n";
	my $peak_mean = sprintf("%5.2f",mean(@peaks_closest));
	#print "mean = $peak_mean\n";

	print "$name,$peak_mean\n";
	print FH "$name,$peak_mean\n";


}


sub mean {
    return sum(@_)/@_;
}


sub read_col {
	open F,$_[0];
	my $col_num = $_[1]; 
	<F>; # skip header
	my @fdata = <F>; chomp @fdata;
	my @res = ();
	foreach my $line (@fdata) {
		next if $line =~ /#/;
		next if length $line < 3;
		$line =~ s/\s+/ /g;
		my @vals = (); @vals = split(/[:,\s\/]+/,$line);
		#print $vals[$col_num-1]."\n";
		push @res, $vals[$col_num-1];
	}
	return @res;
}


close FH
