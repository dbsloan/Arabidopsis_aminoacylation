#!/usr/bin/env perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 input_file\n\n";

my $file = shift or die ($usage);

my @lines = file_to_array($file);

shift @lines;

my %count_sum_HoHoH; #IsoDecoder -> Type -> Position

foreach (@lines){
	chomp $_;
	my @sl = split (/\t/, $_);

	$count_sum_HoHoH{$sl[2]}->{$sl[1]}->{$sl[3]} += abs($sl[4]);
}

my $print_string =  "\+ facetted_pos_scales( y = list (";

foreach my $isodec (sort keys %count_sum_HoHoH){

	my $max_count=0;
	
	my %type_HoH = %{$count_sum_HoHoH{$isodec}};
	
	foreach my $type (sort keys %type_HoH){
		my %pos_hash = %{$type_HoH{$type}};
		foreach my $pos (sort keys %pos_hash){
			abs($pos_hash{$pos}) > $max_count and $max_count = abs($pos_hash{$pos});
		}
	}

	$max_count = sprintf("%.1f", $max_count);

	$print_string .= "IsoDecoder == \"$isodec\" ~ scale_y_continuous(limits = c(" . -1 * $max_count / 1000 . "," . 1 * $max_count / 1000 . ")),";

}

$print_string = substr ($print_string, 0, -1);
$print_string .= "))\n";

print $print_string;