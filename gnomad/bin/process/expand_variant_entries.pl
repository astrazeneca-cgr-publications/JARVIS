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


my $cnt = 0;
while(my $line = <FH>){
	
	if($cnt == 0){	
		if($line ne "POS\tREF\tALT\tQUAL\tAC\tAF\tAN\tDP\n"){
			print "[Warning]: incorrect/non-existent header in input file '$file'.\n    >> Please make sure the 1st line of your input file is:\n    POS     REF     ALT     QUAL    AC	AF	AN	DP\n    (tab-separated) or the header line was knowingly ommited.\n";
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
	my $af = $vals[5];
	my $an = $vals[6];
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
