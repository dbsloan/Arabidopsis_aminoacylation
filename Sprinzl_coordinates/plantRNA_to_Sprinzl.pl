#!/usr/bin/env perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 input_file > output_file\n\n";

my $file = shift or die ($usage);

my @file_lines = file_to_array ($file);

shift @file_lines;

print "Gene\tPosition\tSprinzlPos\tNucleotide\tRegion\tNote\n";

while (my $line = shift (@file_lines)){
	
	my @sl = split (/\t/, $line);
	my $tRNA_pos = 1;
	my $Sprinzl_pos = 0;
	
	my $name = $sl[0];
	
	#process P -1 position (generally only in tRNA-His and only 1 nt)
	if ($sl[6]){
		if (length ($sl[6]) > 1){
			$Sprinzl_pos -= length ($sl[6]) - 1;
			for (my $i = 0; $i < length ($sl[6]); ++$i){
				print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[6], $i, 1), "\tP -1\tManual Check: more than one Pos -1 base\n";
				++$tRNA_pos;
				++$Sprinzl_pos;
			}
		}else{
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t$sl[6]\tP -1\t\n";
			++$tRNA_pos;
			++$Sprinzl_pos;
		}
	}else{
		++$Sprinzl_pos;
	}
	
	#process acceptor stem (5'). Should be 7 nt.
	if (length($sl[7]) > 7){
		for (my $i = 0; $i < length ($sl[7]); ++$i){
			if ($i >= 7){
				print "$name\t$tRNA_pos\t7x\t", substr($sl[7], $i, 1), "\t5\' acceptor stem\tManual Check: 5' acceptor stem does not contain 7 nt\n";			
			}else{
				print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[7], $i, 1), "\t5\' acceptor stem\tManual Check: 5' acceptor stem does not contain 7 nt\n";
			}
			++$tRNA_pos;
			++$Sprinzl_pos;
		}
	}else{
		for (my $i = 0; $i < length ($sl[7]); ++$i){
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[7], $i, 1), "\t5\' acceptor stem\t\n";
			++$tRNA_pos;
			++$Sprinzl_pos;		
		}
	}
	
	#set Sprinzl to 8 after acceptor stem regardless of number of nt so far.
	$Sprinzl_pos = 8;
	
	
	#process pos 8/9. Should be 2 nt
	if (length($sl[8]) != 2){
		for (my $i = 0; $i < length ($sl[8]); ++$i){
			if ($i >= 2){
				print "$name\t$tRNA_pos\t9x\t", substr($sl[8], $i, 1), "\tP 8\-9\tManual Check: Pos 8-9 is not 2 nt\n";		
			}else{
				print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[8], $i, 1), "\tP 8\-9\tManual Check: Pos 8-9 is not 2 nt\n";		
			}
			++$tRNA_pos;
			++$Sprinzl_pos;
		}
	}else{
		for (my $i = 0; $i < length ($sl[8]); ++$i){
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[8], $i, 1), "\tP 8\-9\t\n";
			++$tRNA_pos;
			++$Sprinzl_pos;
		}	
	}

	#set Sprinzl to 10 after pos 8/9 regardless of number of nt so far.
	$Sprinzl_pos = 10;
	
	#process D stem and d loop
	length($sl[9]) != length($sl[11]) and print STDERR "Warning: D stems are not of equal length for the $name";
	
	#first D stem
	for (my $i = 0; $i < length ($sl[9]); ++$i){
		print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[9], $i, 1), "\t5\' D stem\t\n";
		++$tRNA_pos;
		++$Sprinzl_pos;
	}	
	
	#d loop
	my @gg_dloop = ($sl[10] =~ /GG/g);
	my $gg_count = scalar(@gg_dloop);
	if ($gg_count){
		my $last_gg_pos = $-[0];
		my $dloop_warning = "";
		if ($gg_count > 1){
			$dloop_warning = "Manual Check: D loop contains multiple GG motifs";
		}
		for (my $i = 0; $i < $last_gg_pos; ++$i){
			if ($Sprinzl_pos > 17){
				print "$name\t$tRNA_pos\t17x\t", substr($sl[10], $i, 1), "\tD loop\t$dloop_warning\n";
			}else{
				print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[10], $i, 1), "\tD loop\t$dloop_warning\n";
			}
			++$tRNA_pos;
			++$Sprinzl_pos;
		}
		print "$name\t$tRNA_pos\t18\tG\tD loop\t$dloop_warning\n";
		++$tRNA_pos;
		print "$name\t$tRNA_pos\t19\tG\tD loop\t$dloop_warning\n";
		++$tRNA_pos;
		$Sprinzl_pos = 20;
			
		if ($last_gg_pos + 2 < length ($sl[10])){
				
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[10], $last_gg_pos + 2, 1), "\tD loop\t$dloop_warning\n";
			++$tRNA_pos;
			++$Sprinzl_pos;
				
			for (my $i = $last_gg_pos + 3; $i < length ($sl[10]); ++$i){
				
				if (length($sl[10]) - $i + length($sl[11]) > 5){
					print "$name\t$tRNA_pos\t20x\t", substr($sl[10], $i, 1), "\tD loop\t$dloop_warning\n";
					++$tRNA_pos
				}else{
					print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[10], $i, 1), "\tD loop\t$dloop_warning\n";
					++$tRNA_pos;
					++$Sprinzl_pos;
				}
			}
		}
					
	}else{
		for (my $i = 0; $i < length ($sl[10]); ++$i){
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[10], $i, 1), "\tD loop\tManual Check: D loop does not contain GG motif\n";
			++$tRNA_pos;
			++$Sprinzl_pos;
		}
	}
	
	#second D stem
	my $dstem_warning = "";
	
	if ($Sprinzl_pos + length ($sl[11]) - 1 != 25){
		$dstem_warning = "Manual Check: 3' D stem does not end at position 25";
	}
	
	for (my $i = 0; $i < length ($sl[11]); ++$i){
		print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[11], $i, 1), "\t3\' D stem\t$dstem_warning\n";
		++$tRNA_pos;
		++$Sprinzl_pos;
	}

	#set Sprinzl to 26 after second D stem regardless of number of nt so far.	
	$Sprinzl_pos = 26;
	
	#process P 26
	my $p26_warning = "";
	if (length ($sl[12]) > 1){
		$p26_warning = "Manual Check: Multiple nucleotides at P 26";
	}
	for (my $i = 0; $i < length ($sl[12]); ++$i){
		print "$name\t$tRNA_pos\t26\t", substr($sl[12], $i, 1), "\tP 26\t$p26_warning\n";
		++$tRNA_pos;
	}

	#set Sprinzl to 27 regardless of number of nt so far.	
	$Sprinzl_pos = 27;
	
	#process 5' anticodon stem
	my $antistem5_warning = "";
	if (length($sl[13]) != 5){
		$antistem5_warning = "Manual Check: 5\' anticodon stem is not 5 nt long";
	}
	
	for (my $i = 0; $i < length ($sl[13]); ++$i){
		if ($Sprinzl_pos > 31){
			print "$name\t$tRNA_pos\t31x\t", substr($sl[13], $i, 1), "\t5\' anticodon stem\tManual Check: 5\' anticodon stem is longer than 5 nt\n";	
		}else{
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[13], $i, 1), "\t5\' anticodon stem\t$antistem5_warning\n";
		}
		++$tRNA_pos;
		++$Sprinzl_pos;
	}

	#set Sprinzl to 32 regardless of number of nt so far.	
	$Sprinzl_pos = 32;

	#process anticodon loop
	my $antiloop_warning = "";
	if (length($sl[14]) != 7){
		$antiloop_warning = "Manual Check: Anticodon loop is not 7 nt long";
	}
	for (my $i = 0; $i < length ($sl[14]); ++$i){
		if ($Sprinzl_pos > 38){
			print "$name\t$tRNA_pos\t38x\t", substr($sl[14], $i, 1), "\tAnticodon loop\tManual Check: Anticodon loop is longer than 7 nt\n";	
		}else{
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[14], $i, 1), "\tAnticodon loop\t$antiloop_warning\n";
		}
		++$tRNA_pos;
		++$Sprinzl_pos;
	}
	
	#set Sprinzl to 39 regardless of number of nt so far.	
	$Sprinzl_pos = 39;
	
	#process 3' anticodon stem
	my $antistem3_warning = "";
	if (length($sl[15]) != 5){
		$antistem3_warning = "Manual Check: 3\' anticodon stem is not 5 nt long";
	}
	for (my $i = 0; $i < length ($sl[15]); ++$i){
		if ($Sprinzl_pos > 43){
			print "$name\t$tRNA_pos\t43x\t", substr($sl[15], $i, 1), "\t3\' anticodon stem\tManual Check: 3\' anticodon stem is longer than 5 nt\n";	
		}else{
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[15], $i, 1), "\t3\' anticodon stem\t$antistem3_warning\n";
		}
		++$tRNA_pos;
		++$Sprinzl_pos;
	}

	#set Sprinzl to 44 regardless of number of nt so far.	
	$Sprinzl_pos = 44;
	
	#process variable loop
	if (length($sl[16]) < 4){
		for (my $i = 0; $i < length ($sl[16]); ++$i){
			print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[16], $i, 1), "\tVariable loop\tManual Check: Variable loop is shorter than 4 nt\n";
			++$tRNA_pos;
			++$Sprinzl_pos;
		}
	}elsif(length($sl[16]) == 4){
		print "$name\t$tRNA_pos\t44\t", substr($sl[16], 0, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
		print "$name\t$tRNA_pos\t45\t", substr($sl[16], 1, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
		print "$name\t$tRNA_pos\t46\t", substr($sl[16], 2, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
		print "$name\t$tRNA_pos\t48\t", substr($sl[16], 3, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
	}else{
		#assuming that any V loop > 5 nt does have a position 47.
		print "$name\t$tRNA_pos\t44\t", substr($sl[16], 0, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
		print "$name\t$tRNA_pos\t45\t", substr($sl[16], 1, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
		
		for (my $i = 2; $i < length ($sl[16]) - 3; ++$i){
			print "$name\t$tRNA_pos\t45x\t", substr($sl[16], $i, 1), "\tVariable loop\t\n";
			++$tRNA_pos;
		}	
		
		print "$name\t$tRNA_pos\t46\t", substr($sl[16], -3, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
		print "$name\t$tRNA_pos\t47\t", substr($sl[16], -2, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
		print "$name\t$tRNA_pos\t48\t", substr($sl[16], -1, 1), "\tVariable loop\t\n";
		++$tRNA_pos;
	}
	
	#set Sprinzl to 49 regardless of number of nt so far.	
	$Sprinzl_pos = 49;
	
	#process first T stem
	my $tstem5_warning = "";
	if (length($sl[17]) != 5){
		$tstem5_warning = "Manual Check: 5\' T stem is not 5 nt long";
	}
	
	for (my $i = 0; $i < length ($sl[17]); ++$i){
		print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[17], $i, 1), "\t5\' T stem\t$tstem5_warning\n";
		++$tRNA_pos;
		++$Sprinzl_pos;
	}

	#set Sprinzl to 54 regardless of number of nt so far.	
	$Sprinzl_pos = 54;
	
	#process T loop
	my $tloop_warning = "";
	if (length($sl[18]) != 7){
		$tloop_warning = "Manual Check: T loop is not 7 nt long";
	}
	
	for (my $i = 0; $i < length ($sl[18]); ++$i){
		print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[18], $i, 1), "\tT loop\t$tloop_warning\n";
		++$tRNA_pos;
		++$Sprinzl_pos;
	}

	#set Sprinzl to 61 regardless of number of nt so far.	
	$Sprinzl_pos = 61;
	
	#process second T stem
	my $tstem3_warning = "";
	if (length($sl[19]) != 5){
		$tstem3_warning = "Manual Check: 3\' T stem is not 5 nt long";
	}
	
	for (my $i = 0; $i < length ($sl[19]); ++$i){
		print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[19], $i, 1), "\t3\' T stem\t$tstem3_warning\n";
		++$tRNA_pos;
		++$Sprinzl_pos;
	}

	#set Sprinzl to 66 regardless of number of nt so far.	
	$Sprinzl_pos = 66;
	
	#process second acceptor stem
	my $accstem3_warning = "";
	if (length($sl[20]) != 7){
		$accstem3_warning = "Manual Check: 3\' acceptor stem is not 7 nt long";
	}
	
	for (my $i = 0; $i < length ($sl[20]); ++$i){
		print "$name\t$tRNA_pos\t$Sprinzl_pos\t", substr($sl[20], $i, 1), "\t3\' acceptor stem\t$accstem3_warning\n";
		++$tRNA_pos;
		++$Sprinzl_pos;
	}

	#set Sprinzl to 73 regardless of number of nt so far.	
	$Sprinzl_pos = 73;
	
	#process discriminator base
	my $discrim_warning = "";
	if (length($sl[21]) != 1){
		$discrim_warning = "Manual Check: Disciminator base is not 1 nt long";
	}
	
	for (my $i = 0; $i < length ($sl[21]); ++$i){
		print "$name\t$tRNA_pos\t73\t", substr($sl[21], $i, 1), "\tDisciminator base\t$discrim_warning\n";
		++$tRNA_pos;
	}

	#add CCA
	print "$name\t$tRNA_pos\t74\tC\tCCA tail\t\n";
		++$tRNA_pos;
	print "$name\t$tRNA_pos\t75\tC\tCCA tail\t\n";
		++$tRNA_pos;
	print "$name\t$tRNA_pos\t76\tA\tCCA tail\t\n";
}

