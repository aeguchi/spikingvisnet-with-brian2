#!/usr/bin/perl

	#
	#  ParamGenAndRunVisNet.pl
	#  VisBack
	#
	#  Created by Akihiro Eguchi on 15/12/15.
	#  Copyright 2015 OFTNAI. All rights reserved.
	#
	#  This is a perl script to be called by dakota optimizer

	use strict;
    use warnings;
    use POSIX;
	use File::Copy;
	use File::Copy 'cp';
	use Data::Dumper;
	#use Data::Compare;
	use Cwd 'abs_path';
	
	
	
	### LOADING INPUT ###
	my ($infi, $outdir) = @ARGV;
	
	open(my $fh, "<:encoding(UTF-8)", $infi)
		or die "Could not open";
	
	my @table;
	while (my $row = <$fh>){
		chomp $row;
		my $white = substr($row, 0, 1);
		my @rowInput;
		my $flag = 0;
		my $tmpStr = '';
		for (my $i=0;$i<length($row);$i++){
			my $char = substr($row, $i, 1);
			#print "$char\n";
			if ($char ne $white){
				$flag = 1;
				$tmpStr = $tmpStr . $char;
				#print($tmpStr);
				#print($char);
			}elsif($flag == 1){
				#print("$tmpStr\n");
				push(@rowInput,$tmpStr);
				$tmpStr = '';
				$flag = 0;
			}
			if($i==length($row)-1){
				#print("$tmpStr\n");
				push(@rowInput,$tmpStr);
				$tmpStr = '';
				$flag = 0;
			}
		}
	
	#	print("@rowInput\n");
		push(@table,[@rowInput]);	
	}
	close $fh;
	my $par1 = sprintf("%.3f", $table[1][0]);
	my $par2 = sprintf("%.3f", $table[2][0]);
	my $par3 = sprintf("%.3f", $table[3][0]);
	my $par4 = sprintf("%.3f", $table[4][0]);
	my $par5 = sprintf("%.3f", $table[5][0]);
#	my $ss1 = $table[3][0];
#	my $ss2 = $table[4][0];
#	my $ir1 = $table[5][0];
#	my $ir2 = $table[6][0];
	
	system("python ../../dakota_runMe.py $outdir $par1 $par2 $par3 $par4 $par5");
	