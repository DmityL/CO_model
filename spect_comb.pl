#!/usr/bin/perl
use strict;
use warnings;
# Only pull in PDL::Primitive:
use PDL::LiteF;        # loads less modules
use PDL::NiceSlice;    # preprocessor for easier pdl indexing syntax 
 
use PDL::Stats;

#my $dir = "+GaussClumps";
my @cubes = qw(
_FCRAO_12CO_spectra	0
_FCRAO_13CO_spectra	0
_SMT_12CO_spectra	0
_SMT_13CO_spectra	0
);

#_FCRAO_13CO_spectra_model	0

my $dv = 35;
my $cubes_num = @cubes; $cubes_num = $cubes_num / 2;
for (my $i = 0; $i<$cubes_num ; $i++) {
	$cubes[$i*2+1] = $dv*$i;
}


# read all files as clumps
opendir(DIR, $cubes[0]);
my @clumps= readdir(DIR); 

# Цикл по сгусткам (счётчик = cube)
my $clump = 0;
for ($clump=0; $clump<@clumps; $clump=$clump+1) {
	next if $clumps[$clump] eq '..';
	next if $clumps[$clump] eq '.';
	next if not $clumps[$clump] =~ m/\.dat/;
print "$clumps[$clump]\n";

my $i = 0;
my @xdata = (); # Скорость
my @ydata = (); # Интенсивность
my @zdata = (); # Сигма

# loop throught data cubes  (counter = i)
for ($i=0; $i<@cubes; $i=$i+2) {
	print "$cubes[$i]\n";
	my @vdata = read_col("./$cubes[$i]/$clumps[$clump]",2); # read velocity data
	my @idata = read_col("./$cubes[$i]/$clumps[$clump]",3); # read intensity data

	my @v_filter = ();
	my @i_filter = ();
	my @uncdata = ();	

	if ($idata[0]) {
	

	my $piddle = new PDL @idata;
	my $rms = 0;
	my ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($piddle->slice('0:70'));

	# add velocity shift for each line
	my $j = 0;
	foreach my $v (@vdata) {
		my $unc = 1;
		#$unc = 0.5 if ($v >= -18) and ($v <= -14);
		#$unc = 0.2 if ($v >= -20) and ($v <= -14) and ($i==0 or $i==1);
		if ($v > -35 and $v < 0) {
			push @v_filter,$v+ $cubes[$i+1];
			push @i_filter,$idata[$j];
			push @uncdata,$rms;

		}
		#print "i = $i, v = $v, idata=$idata[$j]\n";
	
		$j++;
	}

	}

#	Вычисляем стат. вес
#	foreach my $unc (@uncdata) {
#		$unc = 1 if $unc < 5; 
#		$unc = 0.7 if $unc >= 5;
#		$unc = 0.3 if $unc > 10;
#		$unc = 0.1 if $unc > 15;
#	}

	# combining data
	push @xdata, @v_filter;
	push @ydata, @i_filter;
	push @zdata, @uncdata;
}

# output resulted combined spectra
open FILE,">_spect_comb/$clumps[$clump].spec";
for ($i = 0; $i < @xdata; $i++) {
	print FILE "$xdata[$i]\t$ydata[$i]\t$zdata[$i]\n";	
}

}

sub read_col {
	open F,$_[0];
	my $col_num = $_[1];
	my @fdata = <F>; chomp @fdata;
	my @res = ();
	foreach my $line (@fdata) {
		next if $line =~ /#/;
		$line =~ s/    / /g;
		$line =~ s/   / /g;
		$line =~ s/  / /g;
		my ($xval,$yval,$zval) = split(" ",$line);
		push @res, $xval if $col_num == 1;
		push @res, $yval if $col_num == 2;
		push @res, $zval if $col_num == 3;
	}
	return @res;
}
