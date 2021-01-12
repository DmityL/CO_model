#!/usr/bin/perl
use Number::Closest;
use List::MoreUtils qw(first_index);
use List::Util qw(sum);

# Read clump l,b,v data
my @clump_names = read_col("Calc.csv",1);
my @clump_v = read_col("Calc.csv",4); 
my @clump_tex = read_col("Calc.csv",6); 
my @clump_tau12 = read_col("Calc.csv",9); 
my @clump_sigma = read_col("Calc.csv",11);

 my $i = 0;
foreach my $name (@clump_names) {
	if ($clump_tex[$i] and $clump_tau12[$i]) {
	
	#print "$name, Tex=$clump_tex[$i]\n";
	open(FH, '>', "_spect_comb/$name.dat.dat") or die $!;
print FH qq|tbg = 2.7
rat = 70
shift10 = 0

tx1 = $clump_tex[$i]
t1 = $clump_tau12[$i]
v1 = $clump_v[$i]
d1 = $clump_sigma[$i]
|;
	close FH;
	} else {print "In clump $name values are undefined\n"}
	$i++;
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
		my @vals = (); @vals = split(",",$line);
		#print $vals[$col_num-1]."\n";
		push @res, $vals[$col_num-1];
	}
	return @res;
}


close FH
