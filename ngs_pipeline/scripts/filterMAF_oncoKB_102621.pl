#!/usr/bin/perl
use strict;
use warnings;

open(IN, "$ARGV[0]") or die "cannot open file $ARGV[0]:$!\n";
my $title = <IN>;
my @title = split(/\t/,$title);
chomp $title;
print "$title\n";

my ($mutation_effect,$oncogenic,$LEVEL_1,$LEVEL_2,$LEVEL_3A,$LEVEL_3B,$LEVEL_4,$LEVEL_R1,$LEVEL_R2,$LEVEL_R3,$Highest_level,$citations);
for(0..$#title){
	$mutation_effect= $_ if $title[$_] =~ "MUTATION_EFFECT";
	$oncogenic	= $_ if $title[$_] =~ "ONCOGENIC";
	$LEVEL_1	= $_ if $title[$_] =~ "LEVEL_1";
	$LEVEL_2	= $_ if $title[$_] =~ "LEVEL_2";
	$LEVEL_3A	= $_ if $title[$_] =~ "LEVEL_3A";
	$LEVEL_3B	= $_ if $title[$_] =~ "LEVEL_3B";
	$LEVEL_4	= $_ if $title[$_] =~ "LEVEL_4";
	$LEVEL_R1	= $_ if $title[$_] =~ "LEVEL_R1";
	$LEVEL_R2	= $_ if $title[$_] =~ "LEVEL_R2";
	$LEVEL_R3	= $_ if $title[$_] =~ "LEVEL_R3";
	$Highest_level	= $_ if $title[$_] =~ "HIGHEST_LEVEL";
	$citations	= $_ if $title[$_] =~ "CITATIONS";

}
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	#if(defined $line[$mutation_effect] or defined $line[$oncogenic] or defined $line[$Highest_level] or defined $line[$citations]){
	#	print "$_\n";
	#}
	if($line[$mutation_effect] ne "Unknown"){
		print "$_\n";
	}elsif(defined $line[$oncogenic] or defined $line[$Highest_level] or defined $line[$citations]){
		print "$_\n";
	}
}
close IN;
