#!/usr/bin/env perl;

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 sam_file sprinzl_file output_basename. Optionally, add --unique flag after the filename to restrict analysis to reads that had a single top mapping hit (no ties).\n\n";

my $sam_file = shift or die ($usage);
my $sprinzl_file = shift or die ($usage);
my $output_base = shift or die ($usage);
my $unique;
if ($unique = shift){
	$unique eq "--unique" or die ("\nERROR: the only option that can be provided after the same file name is --unique.\n\n")
}


my @sprinzl_lines = file_to_array ($sprinzl_file);
shift @sprinzl_lines;

my %sprinzl_HoH;

foreach (@sprinzl_lines){
	my @sl = split (/\t/, $_);
	$sprinzl_HoH{$sl[0]}->{$sl[1]} = $sl[2];
}



my $FH = open_file ($sam_file);

my %startPos_HoH;
my %endPos_HoH;

while (my $line = <$FH>){

	$line =~ /^\@/ and next;
	
	my @sl = split (/\t/, $line);

	my $seq = $sl[9];
	my $name = $sl[0];
	my $ref = $sl[2];
	my $start = $sl[3];
	my $cigar = $sl[5];
	my $AS = $sl[11];
	my $XS = $sl[12];
	
	if ($unique){
		my $top_score;
		my $second_score;
		
		if ($AS =~ /AS\:i\:([\-\d]+)$/){
			$top_score = $1;
		}else{
			die ("\nERROR: could not parse AS score value: $AS\n\n");
		}

		if ($XS =~ /XS\:i\:([\-\d]+)$/){
			$second_score = $1;
			$top_score > $second_score or next;
		}
	}

	++$startPos_HoH{$ref}->{$start};
	
	my $end_pos = $start - 1;
	
	while (length ($cigar) > 0){
		
		if ($cigar =~ /^(\d+[DIM])/){
			my $aln_piece = $1;
			my $aln_len = substr ($aln_piece, 0, -1);
			unless (substr ($aln_piece, -1) eq "I"){
				$end_pos += $aln_len;
			}
		
			if (length ($cigar) > length($aln_piece)){
				$cigar = substr ($cigar, length($aln_piece));
			}else{
				$cigar = "";
			}
		}else{
			die ("\nERROR: Could not parse CIGAR string $cigar\.\n\n")
		}		
	}
	
	++$endPos_HoH{$ref}->{$end_pos};
	

}

my $FHO1 = open_output("$output_base\.start_pos.txt");
my $FHO2 = open_output("$output_base\.end_pos.txt");

print $FHO1 "Gene\tPosition\tSprinzlPosition\tStartPosCount\n";

foreach my $gene (sort keys %startPos_HoH){
	my %pos_hash = %{$startPos_HoH{$gene}};
	foreach my $pos (sort { $a <=> $b } keys %pos_hash){
		print $FHO1 "$gene\t$pos\t", $sprinzl_HoH{$gene}->{$pos}, "\t", $pos_hash{$pos}, "\n";
	}
}

print $FHO2 "Gene\tPosition\tSprinzlPosition\tEndPosCount\n";

foreach my $gene (sort keys %endPos_HoH){
	my %pos_hash = %{$endPos_HoH{$gene}};
	foreach my $pos (sort { $a <=> $b } keys %pos_hash){
		print $FHO2 "$gene\t$pos\t", $sprinzl_HoH{$gene}->{$pos}, "\t", $pos_hash{$pos}, "\n";
	}
}

