use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 file_list sprinzl_file min_read_count\n\n";


my $file_list = shift or die ($usage);
my $sprinzl_file = shift or die ($usage);
my $min_read_count = shift or die ($usage);


my @sprinzl_lines = file_to_array ($sprinzl_file);
shift @sprinzl_lines;

my %sprinzl_HoH;

foreach (@sprinzl_lines){
	my @sl = split (/\t/, $_);
	$sprinzl_HoH{$sl[0]}->{$sl[1]} = $sl[2];
}


my @files = file_to_array($file_list);

#print "@files";

my %Isodecoder_fam_total_CC_HoHoA; #1st key: isodecoder, 2nd key: genome, array positions file number, value total read count
my %Isodecoder_fam_total_CCA_HoHoA; #same structure as above

my $file_num = 0;

foreach (@files){
	chomp $_;
	my $FH = open_file ($_);
	my $first_line = <$FH>;
	
	while (my $line = <$FH>){
		my @sl = split (/\t/, $line);
		my @split_gene = split (/\-/, $sl[0]);
		$split_gene[2] eq "stemloop" and next;
		$Isodecoder_fam_total_CC_HoHoA{$split_gene[2]}->{$split_gene[0]}[$file_num] += $sl[5];
		$Isodecoder_fam_total_CCA_HoHoA{$split_gene[2]}->{$split_gene[0]}[$file_num] += $sl[6];
	}
	close $FH;
	++$file_num;
}


my %Isodecoder_fam_byPos_CC_HoHoHoA; #1st key: isodecoder, 2nd key: genome, 3rd key sprinzl pos, array positions file number, value total read count
my %Isodecoder_fam_byPos_CCA_HoHoHoA; #same structure as above
my %pos_tracker_HoHoH; #same as above, but no array. Tracking which positions have counts in any file (regardless of CC/CCA)

$file_num = 0;

foreach (@files){
	chomp $_;
	my $FH = open_file ($_);
	my $first_line = <$FH>;
	
	while (my $line = <$FH>){
		my @sl = split (/\t/, $line);
		my @split_gene = split (/\-/, $sl[0]);
		$split_gene[2] eq "stemloop" and next;
		my $sprinzl_pos = $sprinzl_HoH{$sl[0]}->{$sl[1]};
		if ($sprinzl_pos =~ /^(\d+)x/){
			$sprinzl_pos = $1;
		}
		$Isodecoder_fam_byPos_CC_HoHoHoA{$split_gene[2]}->{$split_gene[0]}->{$sprinzl_pos}[$file_num] += $sl[5];
		$Isodecoder_fam_byPos_CCA_HoHoHoA{$split_gene[2]}->{$split_gene[0]}->{$sprinzl_pos}[$file_num] += $sl[6];
		$pos_tracker_HoHoH{$split_gene[2]}->{$split_gene[0]}->{$sprinzl_pos} = 1;
	}
	close $FH;
	++$file_num;
}


#set all undefined values to 0.
foreach my $isodecoder (keys %pos_tracker_HoHoH){
	my %genome_HoH = %{$pos_tracker_HoHoH{$isodecoder}};
	foreach my $genome (sort keys %genome_HoH){
		
		my $fileCounter = 0;
		foreach my $file (@files){
			exists ($Isodecoder_fam_total_CC_HoHoA{$isodecoder}->{$genome}[$fileCounter]) or $Isodecoder_fam_total_CC_HoHoA{$isodecoder}->{$genome}[$fileCounter] = 0;
			exists ($Isodecoder_fam_total_CCA_HoHoA{$isodecoder}->{$genome}[$fileCounter]) or $Isodecoder_fam_total_CCA_HoHoA{$isodecoder}->{$genome}[$fileCounter] = 0;
			++$fileCounter;
		}
		
		my %pos_hash = %{$genome_HoH{$genome}};
		foreach my $pos (sort keys %pos_hash){
			my $fileCounter2 = 0;
			foreach my $file2 (@files){
				exists ($Isodecoder_fam_byPos_CC_HoHoHoA{$isodecoder}->{$genome}->{$pos}[$fileCounter2]) or $Isodecoder_fam_byPos_CC_HoHoHoA{$isodecoder}->{$genome}->{$pos}[$fileCounter2] = 0;
				exists ($Isodecoder_fam_byPos_CC_HoHoHoA{$isodecoder}->{$genome}->{$pos}[$fileCounter2]) or $Isodecoder_fam_byPos_CC_HoHoHoA{$isodecoder}->{$genome}->{$pos}[$fileCounter2] = 0;				
				++$fileCounter2;
			}
		}
	}
}


print "Isodecoder\tGenome\tPosition\tLibrary\tEndType\tAbundance\n";

foreach my $isodecoder (keys %pos_tracker_HoHoH){
	my %genome_HoH = %{$pos_tracker_HoHoH{$isodecoder}};
	foreach my $genome (sort keys %genome_HoH){
		
		my $meet_minimums = 1;
		
		my $fileCounter = 0;
		foreach my $file (@files){
			$Isodecoder_fam_total_CC_HoHoA{$isodecoder}->{$genome}[$fileCounter] >= $min_read_count or $meet_minimums = 10;
			$Isodecoder_fam_total_CCA_HoHoA{$isodecoder}->{$genome}[$fileCounter] >= $min_read_count or $meet_minimums = 10;
			++$fileCounter;
		}
		
		$meet_minimums or next;
		
		my %pos_hash = %{$genome_HoH{$genome}};
		foreach my $pos (sort {$a <=> $b} keys %pos_hash){
			my $fileCounter2 = 0;
			foreach my $file2 (@files){
				print "$isodecoder\t$genome\t$pos\t$file2\tCC\t";
				print -1 * $Isodecoder_fam_byPos_CC_HoHoHoA{$isodecoder}->{$genome}->{$pos}[$fileCounter2] / $Isodecoder_fam_total_CC_HoHoA{$isodecoder}->{$genome}[$fileCounter2];
				print "\n";
				print "$isodecoder\t$genome\t$pos\t$file2\tCCA\t";
				print $Isodecoder_fam_byPos_CCA_HoHoHoA{$isodecoder}->{$genome}->{$pos}[$fileCounter2] / $Isodecoder_fam_total_CCA_HoHoA{$isodecoder}->{$genome}[$fileCounter2];
				print "\n";
				++$fileCounter2;
			}
		}
	}
}