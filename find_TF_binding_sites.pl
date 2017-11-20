#!/usr/bin/perl
use warnings;
use strict;

my $count = {}; my $matrix = {}; my $TF;
#get base counts
(open MATRIX, "motif_matrix.txt") or die $!;
while (<MATRIX>)
{
	chomp;
	if ($_ =~ "Transcription Factor Name:")
	{
		my @line = split(": ", $_);
		$TF = $line[1];
	}
	elsif ($_ =~ /[a|c|g|t]/)
	{
		my @line = split("\t", $_);
		my $base = shift(@line);
		shift(@line); #remove |
		
		my $position = 1;
		foreach my $number (@line)
		{
			if ($base eq "a") {$base = "A"} elsif ($base eq "g") {$base = "G"} elsif ($base eq "c") {$base = "C"} elsif ($base eq "t") {$base = "T"}
			$count->{$TF}{$position}{$base}{$number} = "";
			$position++;
		}
	}
}

#make matrix
foreach my $factor (keys %{$count})
{
	foreach my $position (keys %{$count->{$factor}})
	{
		my $total = 0;
		foreach my $base (keys %{$count->{$factor}{$position}})
		{
			my @number = (keys %{$count->{$factor}{$position}{$base}}); my $number = "@number"; 
			$total += $number;
		}
		foreach my $base (keys %{$count->{$factor}{$position}})
		{
			my @number = (keys %{$count->{$factor}{$position}{$base}}); my $number = "@number";
			my $freq = ($number+($total*0.05))/($total+($total*0.2));
			$matrix->{$factor}{$position}{$base}{$freq} = "";
		}
	}
}

#load sequence
my $c = 0; my $sequence;
(open SEQ, "bglFpAST.fasta") or die $!;
while (<SEQ>)
{
	chomp;
	if ($c < 1)
	{
		$c++;
	}
	else
	{
		$sequence = $sequence.$_;
	}
}

my @split_seq = split("", $sequence);
my $max = scalar(@split_seq);

(open OUTPUT, ">QUERY_TFsites.txt") or die $!;
#find top TF binding sites
foreach my $trans (keys %{$matrix})
{
	my $max_score = 1;
	foreach my $pos (keys %{$matrix->{$trans}})
	{
		my @A = (keys %{$matrix->{$trans}{$pos}{"A"}}); my $A = "@A";
		my @T = (keys %{$matrix->{$trans}{$pos}{"T"}}); my $T = "@T";
		my @G = (keys %{$matrix->{$trans}{$pos}{"G"}}); my $G = "@G";
		my @C = (keys %{$matrix->{$trans}{$pos}{"C"}}); my $C = "@C";
		my @numbers = ("$A", "$T", "$C", "$G");
		my @sorted_numbers = sort{$a <=> $b}@numbers;
		my $top = pop(@sorted_numbers);
		$max_score = $max_score*$top;
	}
	my @positions = (keys %{$matrix->{$trans}});
	my $start = 0; my $length = scalar(@positions);
	while ($start+$length < $max)
	{
		my $query = substr($sequence, $start, $length);
		my @query = split("", $query);
		my $r = 1; my $prob = 1;
		foreach my $nuc (@query)
		{
			my @like = (keys %{$matrix->{$trans}{$r}{$nuc}}); my $like = "@like";
			$prob = $prob*$like;
			$r++;
		}
		my $score = $prob/$max_score;
		if ($score > 0.1)
		{
			print OUTPUT "$trans\t$start\t$length\t$score\t$prob\t$max_score\n";
		}
		$start++;
	}

}	




