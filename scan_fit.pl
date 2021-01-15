#!/usr/bin/perl
use Number::Closest;
use List::MoreUtils qw(first_index);
use List::Util qw(sum);

# Read clump l,b,v data
my @clump_names = read_col("Clumps_l_b_v.csv",1); # читаем данные скорости
my @clump_velo = read_col("Clumps_l_b_v.csv",4); # читаем данные интенсивности

my %velo; my $i = 0;
foreach my $name (@clump_names) {
	$velo{$name} = $clump_velo[$i];
	#print "velo{$name} = $clump_velo[$i]\n";
	$i++;
}

my $dir = "_SMT_13CO_spectra";

open(FH, '>', "$dir"."_FWHM.csv") or die $!;
print FH "Clump,FWHM\n";
opendir(DIR, $dir);
my @clumps= readdir(DIR); 

my $clump = 0;
for ($clump=0; $clump<@clumps; $clump=$clump+1) {
	next if $clumps[$clump] eq '..';
	next if $clumps[$clump] eq '.';
	next if not $clumps[$clump] =~ m/\.fit/;
	my $name = $clumps[$clump]; $name =~ s/\.fit//g;
	my $GC_velo = $velo{$name};


	# Read fit
	open(F, '<', $dir."/".$clumps[$clump]) or die $!;
	while(my $line = <F>){
		# 
	   if ($line =~ m|\#MNW\s*Line width:\s*([0-9\.\-\+]+)\s*([0-9\.\-\+]+)|) {
		   print "$name,$1\n";
		    print FH "$name,$1\n";
	   }
	
	}

}



close FH;



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
