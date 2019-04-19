#!/usr/bin/env perl

use strict;
use warnings;

my $file = $ARGV[0];
my $out_dir = $ARGV[1];



open(FH, $file);
$file =~ s/.*\///;
print "file: $file\n";
print "out_dir: $out_dir\n";
open(OUT, ">$out_dir/$file.filtered");


# 48119868        T       A,TAG   15118.34        AF=8.52273e-03,1.42045e-03
# 48119869        G       A,C     203198.33       AF=0.00000e+00,0.00000e+00
# 48119872        T       G       162596.74       AF=0.00000e+00
my $cnt = 0;
while(my $line = <FH>){
	
	if($cnt == 0){	
		if($line ne "POS\tREF\tALT\tQUAL\tINFO\tINFO\tINFO\tINFO\n"){
			die "[Error]: incorrect header in input file '$file'.\n    >> Please make sure the 1st line of your input file is:\n    POS     REF     ALT     QUAL    INFO	INFO	INFO	INFO\n    (tab-separated)\n";
		}
		print OUT $line;
		$cnt++;	
		next;
	}

	chomp($line);

	my @vals = split(' ', $line);
	my $pos = $vals[0];
	my $ref = $vals[1];
	my $alt = $vals[2];
	my $qual = $vals[3];
	my $ac = $vals[4];
	my $an = $vals[5];
	my $af = $vals[6];
	my $dp = $vals[7];
	
	my @alt_vals = split(',', $alt);
	#print join(", ", @alt_vals)."\n";
	
	$ac =~ s/AC=//;
	my @ac_vals = split(',', $ac);

	$af =~ s/AF=//;
	my @af_vals = split(',', $af);

	$an =~ s/AN=//;

	$dp =~ s/DP=//;
	my @dp_vals = split(',', $dp);


	for(my $i=0; $i<scalar @alt_vals; $i++){
		#print "$pos\t$ref\t$alt_vals[$i]\t$qual\t$ac_vals[$i]\t$af_vals[$i]\t$an\t$dp_vals[$i]\n";
		print OUT "$pos\t$ref\t$alt_vals[$i]\t$qual\t$ac_vals[$i]\t$af_vals[$i]\t$an\t$dp_vals[$i]\n";
	}
	

	if($cnt % 100000 == 0){
		print "$cnt\n";
	}
	$cnt++;
}

close(FH);
close(OUT);

print "<< $file expansion complete.\n";
