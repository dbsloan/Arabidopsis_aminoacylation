#!/usr/bin/env perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 group1_file_list group2_file_list group1_name group2_name genome_filter isodecoder_renaming_file exclusion_file\n\n";

my $group1_file_list = shift or die ($usage); #a test file with list of filenames included in treatment group 1 (positive on y axis)
my $group2_file_list = shift or die ($usage); #a test file with list of filenames included in treatment group 2 (negative on y axis)
my $group1_name = shift or die ($usage); #name for group1 to include in Type column in output
my $group2_name = shift or die ($usage); #name for group2 to include in Type column in output
my $genome_filter = shift or die ($usage); #name of genome for subsetting dataset. Only genes from this genome will be returned
my $isodecoder_renaming_file = shift or die ($usage); #file with scheme for renaming isodecoders (tab delimited text file with two columns: oldname/newname). e.g. eMet -> Met(e)
my $exclusion_file = shift or die ($usage); #file with list of gene names to exclude


my @group1_files = file_to_array($group1_file_list);
my @group2_files = file_to_array($group2_file_list);

my $group1_n = scalar @group1_files; #counts will be used to divide TPM values so that column plotting sums to *average* TPM.
my $group2_n = scalar @group2_files;

my @isodecoder_renaming_lines = file_to_array($isodecoder_renaming_file);
my %isodecoder_renaming_hash;
foreach (@isodecoder_renaming_lines){
	chomp $_;
	my @sl = split (/\t/, $_);
	$isodecoder_renaming_hash{$sl[0]} = $sl[1];
}

my @exclusion_lines = file_to_array($exclusion_file);
my %exclusion_hash;
foreach (@exclusion_lines){
	chomp $_;
	$exclusion_hash{$_} = 1;
}


print "File\tType\tIsoDecoder\tSprinzlPos\tTPM\n";

foreach my $file (@group1_files){
	
	chomp $file;
	my $total_read_count = 0;
	my %isodecoder_count_HoH;
	my @file_lines = file_to_array($file);
	shift @file_lines;
	foreach my $line (@file_lines){
		chomp $line;
		my @sl = split (/\t/, $line);
		exists ($exclusion_hash{$sl[0]}) and next;
		my @sn = split (/\-/, $sl[0]);
		$sn[0] eq $genome_filter or next;
		$total_read_count += $sl[3];
		
		my $SprinzlPos = $sl[2];
		if ($SprinzlPos =~ /^(\d+)x/){ # remove "x" from Sprinzl position strings
			$SprinzlPos = $1;
		}
		
		my $isodec = $sn[2];
		
		if (exists($isodecoder_renaming_hash{$isodec})){
			$isodec = $isodecoder_renaming_hash{$isodec};
		}
		
		$isodecoder_count_HoH{$isodec}->{$SprinzlPos} += $sl[3];
	}
	
	foreach my $isodecoder (sort keys %isodecoder_count_HoH){
		
		my %pos_count_hash = %{$isodecoder_count_HoH{$isodecoder}};
		
		foreach my $pos (sort {$a <=> $b} keys %pos_count_hash){
			print "$file\t$group1_name\t$isodecoder\t$pos\t", 1e6 * $pos_count_hash{$pos} / ($group1_n * $total_read_count), "\n";
		}
	}
}


foreach my $file (@group2_files){
	
	chomp $file;
	my $total_read_count = 0;
	my %isodecoder_count_HoH;
	my @file_lines = file_to_array($file);
	shift @file_lines;
	foreach my $line (@file_lines){
		chomp $line;
		my @sl = split (/\t/, $line);
		exists ($exclusion_hash{$sl[0]}) and next;
		my @sn = split (/\-/, $sl[0]);
		$sn[0] eq $genome_filter or next;
		$total_read_count += $sl[3];
		
		my $SprinzlPos = $sl[2];
		if ($SprinzlPos =~ /^(\d+)x/){ # remove "x" from Sprinzl position strings
			$SprinzlPos = $1;
		}
		
		my $isodec = $sn[2];
		
		if (exists($isodecoder_renaming_hash{$isodec})){
			$isodec = $isodecoder_renaming_hash{$isodec};
		}
		
		$isodecoder_count_HoH{$isodec}->{$SprinzlPos} += $sl[3];
	}
	
	foreach my $isodecoder (sort keys %isodecoder_count_HoH){
		
		my %pos_count_hash = %{$isodecoder_count_HoH{$isodecoder}};
		
		foreach my $pos (sort {$a <=> $b} keys %pos_count_hash){
			print "$file\t$group2_name\t$isodecoder\t$pos\t", -1e6 * $pos_count_hash{$pos} / ($group1_n * $total_read_count), "\n";
		}
	}
}