#!/usr/bin/perl
use File::Path qw(make_path remove_tree);
use List::Util qw(max);
use File::chdir;
use File::Path; 
 # ======================= Input parameters ================================

my $input_dir = "input";
my $prepare_dir = "prepare";
my $work_dir = "calculations";

my @files = qw(12co_1to0.fits 13co_1to0.fits 12co_2to1.fits 13co_2to1.fits);
my @lines = qw(12CO10 13CO10 12CO21 13CO21);
my @beams = qw(45.0 46.0 32.0 32.0); # arcsec
my @main_beam_eff = qw(0.48 0.48 0.75 0.75);
my @clip_vals = qw(3 3 3 3); # in sigma level

my $transition_to_fit = "10";

my $ignore_convol = 1;
my $remove_before_start = 1;

my $base_coord_system_file = 0;
my $file_to_extract_clumps = 0;

my %restfreq = (
	'12CO10' => '115.2712018E+09',
	'13CO10' => '110.2013543E+09',
	'12CO21' => '230.5380000E+09',
	'13CO21' => '220.3986841E+09'
);




 # ======================== Actual code =================================



if ($remove_before_start) {
	rmtree $prepare_dir;
	rmtree $work_dir;
}

if (not -d $prepare_dir) {
	make_path("$prepare_dir");
}
if (not -d $work_dir) {
	make_path("$work_dir");
}

sub run {
	print "$_[0]\n";
	my $res = `$_[0]`;
	print $res."\n";
	return $res;
}

sub run_quiet {
	my $res = `$_[0]`;
	return $res;
}

my @last= ();

print "Converting from FITS to Miriad...\n";
my $i = 0;
foreach $file (@files) {
	my $name = $file; $name =~ s/\.fits//g;
	my $line = $lines[$i];
	my $beam = $beams[$i];
	my $freq = $restfreq{$line}; 
	run_quiet("fits in=$input_dir/$file out=$prepare_dir/$name op=xyin");
	run_quiet("puthd in=$prepare_dir/$name/restfreq value=$freq");
	run_quiet("puthd in=$prepare_dir/$name/bmaj value=$beam");
	run_quiet("puthd in=$prepare_dir/$name/bmin value=$beam");
	run_quiet("puthd in=$prepare_dir/$name/bpa value=0.0");
	run_quiet("puthd in=$prepare_dir/$name/bunit value='K'");
	push(@last,$name);
	$i++;
}
$CWD = "./$prepare_dir";

print "Convert from Ta to Tmb ...\n";
my $i = 0;
foreach $file (@last) {
	my $mbeff = $main_beam_eff[$i];
	run("maths exp=\"<$file>/$mbeff\" out=$file.mb");
	$file = $file.".mb";
	$i++;
}

my @rms = ();
print "Estimate RMS ...\n";
my $i = 0;
foreach $file (@last) {
	my $rms_text = run_quiet("sigest in=$file");
	my $rms_val = 0;
	if ($rms_text =~ m/Estimated rms is ([0-9\-\.E]+)/i) {
		$rms_val = $1;
	}
	push(@rms,$rms_val);
	print "For file $file RMS is $rms_val;\n";
	$i++;
}



if (not $ignore_convol) {
	my $max_beam = max(@beams);
	print "Convolving data cubes... (max. beam is $max_beam)\n";
	$i = 0;
	foreach $file (@last) {
		my $conv_size = sqrt($max_beam*$max_beam - $beams[$i]*$beams[$i]); 
		if ($conv_size > 1) {
			run("convol map=$file out=$file.convol fwhm=$conv_size");
			$file = $file.".convol";
		}
		$i++;
	}
}

my $cs0 = `gethd in=$last[$base_coord_system_file]/ctype1`;
print "Converting gal/eq ... (assuming base coordinate system of file $files[$base_coord_system_file]: $cs0)\n";
foreach $file (@last) {
	my $cs = `gethd in=$file/ctype1`;
	if ($cs ne $cs0 and ($cs0 =~ /ra/i and $cs =~ /glon/i or $cs0 =~ /glon/i and $cs =~ /ra/i)) {
		run("regrid in=$file out=$file.galeq options=galeqsw,offset");
		$file = $file.".galeq";

	}
}

print "Regridding all cubes to same grid (assuming base coordinate system of file $files[$base_coord_system_file])...\n";
$i = 0;
foreach $file (@last) {
	my $template = $last[$base_coord_system_file];
	if ($i != $base_coord_system_file) {
		run("regrid in=$file out=$file.regrid tin=$template axes=1,2"); # axes=1,2
		$file = $file.".regrid";
	}
	$i++;
}

$CWD = "..";

print "Converting from Miriad to FITS ... \n"; 
$i = 0;
foreach $file (@last) {
	my $name = $lines[$i];
	if (not -f "$work_dir/$name.fits") {
		run("fits in=$prepare_dir/$file out=$work_dir/$name.fits op=xyout");
	}
	if (not -d "$work_dir/$name") {
		run("cp -R $prepare_dir/$file $work_dir/$name");
	}
	$file = $name;
	$i++;
}

$CWD = "./$work_dir";

print "Calculating moments... \n";
$i = 0;
foreach $file (@last) {
	my $clip = $clip_vals[$i]*$rms[$i];
	run("moment in=$file out=$file.sum mom=0 clip=-10,$clip");
	run("fits in=$file.sum out=$file.sum.fits op=xyout");
	run("moment in=$file out=$file.peak mom=-2");
	run("fits in=$file.peak out=$file.peak.fits op=xyout");
	run("moment in=$file out=$file.vpeak mom=-3");
	run("fits in=$file.vpeak out=$file.vpeak.fits op=xyout");
	$i++;
}

my $optimize_var = 1;


print "Writing parameters for fit...\n";

if ($optimize_var == 2) {

my %nu = ("10" => "115.271E9", "21" => "230.538E9");
my %nu13 = ("10" => "110.201E9", "21" => "220.399E9");
my %T0 = ("10" => "5.532", "21" => "11.064");
my %T013 = ("10" => "5.288", "21" => "10.577");
my $Brot = 55.1E9; # 13CO
my %Eu13 = ("10" => "5.29", "21" => "15.87");
my %J = ("10" => "1", "21" => "2");



print "Fitting transition $t...\n";

open(F,">params_var2.dat");
print F qq|
{"data12_10": "12CO10.fits",
"data12_21": "12CO21.fits",
"data13_10": "13CO10.fits",
"data13_21": "12CO21.fits",
}
|;
#run("python ../optimize_var2.py");
close(F)

}

if ($optimize_var == 1) {

my %nu = ("10" => "115.271E9", "21" => "230.538E9");
my %nu13 = ("10" => "110.201E9", "21" => "220.399E9");
my %T0 = ("10" => "5.532", "21" => "11.064");
my %T013 = ("10" => "5.288", "21" => "10.577");
my $Brot = 55.1E9; # 13CO
my %Eu13 = ("10" => "5.29", "21" => "15.87");
my %J = ("10" => "1", "21" => "2");

foreach my $t (qw(10 21)) {

print "Fitting transition $t...\n";

open(F,">params.dat");
print F qq|
{"data12": "12CO$t.peak.fits",
"data13": "13CO$t.peak.fits",
"data12_int": "12CO$t.sum.fits",
"data13_int": "13CO$t.sum.fits",
"nu": $nu{$t},
"nu13": $nu13{$t},
"T0": $T0{$t},
"T013": $T013{$t},
"Brot": $Brot,
"Eu":$Eu13{$t},
"J":$J{$t} 
}
|;
run("python ../optimize.py");
close(F)
}

}

#print "Converting from FITS to NDF ... \n";

#my $name = $lines[$file_to_extract_clumps];
#run("\$STARLINK_DIR/bin/convert/fits2ndf $prepare_dir/$name.fits $prepare_dir/data");  

#print "Extracting clumps ... \n";
