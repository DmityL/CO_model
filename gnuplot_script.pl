#!/usr/bin/perl
use strict;
use warnings;

open(FH, '>', "gnuplot_long_script.txt") or die $!;
print FH qq|
	set terminal pdfcairo portrait size 7, 5 enhanced color "Helvetica" 13
cd 'D:/YandexDisk/CO_approx'
set xlabel 'Velocity (km/s)'
set ylabel 'Ta (K)'
set key top noopaque
unset arrow
set grid lw 1 
show grid
|;
for (my $i=0; $i < 300; $i++ ) {

	my $is = sprintf("%03d",$i);
	my $ism = "$is"."m";
	my $ismb = "$is"."mb";
	if (-f "_FCRAO_13CO_spectra/Clump$is.dat") {
#set arrow from 0.75,0 to 0.75,1 nohead lc rgb 'red'

	print FH qq|
 set title 'Spectra of clump $is'
set output 'D:/YandexDisk/CO_approx/clump$is.pdf'
plot '_FCRAO_12CO_spectra/Clump$is.dat' using 2:3 with lines lw 2  title 'FCRAO 12CO(1-0)' ,\\
'_FCRAO_13CO_spectra/Clump$is.dat' using 2:3 with lines lw 2  title 'FCRAO 13CO(1-0)' ,\\
'_FCRAO_13CO_spectra_model/Clump$is.dat' using 2:3 with lines lw 2  title 'GAUSSCLUMP 13CO(1-0)' ,\\
'_SMT_12CO_spectra/Clump$is.dat' using 2:3 with lines lw 2  title 'SMT 12CO(2-1)' ,\\
'_SMT_13CO_spectra/Clump$is.dat' using 2:3 with lines lw 2  title 'SMT 13CO(2-1)'
set output

|;


# '_FCRAO_13CO_spectra/Clump$ism.dat' using 2:3 with lines lw 2  title 'FCRAO 13CO(1-0) model' ,\\
#'_FCRAO_13CO_spectra/Clump$ismb.dat' using 2:3 with lines lw 2  title 'FCRAO 13CO(1-0) mb06' ,\\
	}




}

close FH