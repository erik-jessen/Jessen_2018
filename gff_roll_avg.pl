#!/usr/bin/perl
# Reads from from a gff file and computes rolling average over a window 2*$len
# rollavgfiles.txt contains 2 columns: 1) name of file before rolling avg, 2) name of file wanted after rolling avg

use strict;
use warnings;

my $files = {};
open (FILE, '<rollavgfiles.txt') or die $!;
while (<FILE>)
{
	chomp;
	my ($input, $output) = split("\t", $_);
	$files->{$input}{$output} = "";
}

foreach my $in (keys %{$files})
{
	foreach my $out (keys %{$files->{$in}})
	{
		rollavg($in, $out);
	}
}

sub rollavg
{
my $input = shift;
my $output = shift;
open (FHI, '<', $input) or die $!;
open (FHO, '>', $output) or die $!;

my ($line, @bp,  $bpo, $y, $x, $w, $z, @pos, $inc);
my $line2="abc";
my @cols;
my @colsa;

my $count=0;

my $len=150;

my $bkgd=0; 

#-0.3795; ## specifies background correction for data

$bpo=0;


my $now = gmtime;

print $now,"\n";


#READ LINES FROM INPUT FILE AND CALCULATE


while ($line = <FHI>) 
{
	
	if (substr($line,0,1) eq "#") 
	{

			print FHO $line;

	}
	
	next if (substr($line,0,1) eq "#");

	chomp($line);

	@cols = split("\t",$line);
 
	#split on tabs
	
	$pos[$count]=$cols[3];
	
	$bp[$count]=$cols[5]-$bkgd;
	
	$count++;
	
}
$inc=@bp;
print "$count\t$inc\n";
close (FHI);
open (FHI, '<', $input) or die $!; 
## reopens file - be careful it's the right one

$count=0;

$y=0;


while ($line = <FHI>) 
{
	
	if (substr($line,0,1) eq "#") 
	{
		
		print FHO $line;
		
	}
	
	next if (substr($line,0,1) eq "#");	
	
	chomp($line);
	@cols = split("\t",$line); 

	#split on tabs
#	
	$cols[0]= "K12.fna";
	
	$cols[2]= $output;  

	## specifies new feature for file
	
	$cols[3]=int($cols[3]);
	
	$cols[4]=int($cols[4]);
	
	if (($y>=100) and ($y<=(378237-100))) 
	{
	
		for ($x=$y-100;$x<=$y+100;$x++) 
		{

	
			$w=$x;

			$z=$y;
	
			if ($x<0) 
			{
	
				$w=0;
	
			} 
			elsif ($x>378237) 
			{	
	
				$w=378237;
	
			} 
			else 
			{
	}

	
			if (($pos[$w]>=$pos[$y]-$len) and ($pos[$w]<=$pos[$y]+$len)) 
			{
	
				$bpo=$bpo+$bp[$x];
	
				$count++;
	
			}
	
		}
	
		$bpo=sprintf("%.2f", $bpo/$count);
	
		$count=0;
		
		print FHO "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$bpo\t$cols[6]\t$cols[7]\t$cols[8]\n";

	} 
	else 
	{

		print FHO "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\t$cols[7]\t$cols[8]\n";

		print "$y\t$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\t$cols[7]\t$cols[8]\n";
	
	}

	$bpo=0;
$y++;

}


$now = gmtime;
print $count, "\t",$now,"\n";
close (FHI);
close (FHO);
}