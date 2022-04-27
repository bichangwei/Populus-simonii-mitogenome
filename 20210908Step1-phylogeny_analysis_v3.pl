#!/usr/bin/perl -w
#########################################################################
# FileName: Gene_gain_loss.pl
# Version: dc7ec151-6ff8-4224-8965-faf94f2f5014
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Fri Sep  3 09:39:48 2021
#########################################################################

use strict;
use Getopt::Long;
use File::Basename;
my $red = "\033[0;31m";	my $end = "\033[0m";
my ($in,$out,$help);
GetOptions
(
	"in=s"=>\$in,	"out=s"=>\$out,	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file> <in.input file>
		-o <string> <out.output file>
INFO

die $usage if ($help || !$in || !$out);
open IN1,"$in" || die $!;
open GL,">>gene_list.txt";#gene_list.txt
$/="\n";
while(my $line1=<IN1>)
{
	chomp $line1;
	open IN2,"<$line1";
##output the transformated file to a file starting as 'trans'
	open OUT1,">trans$line1";
	my $gl_name=$line1;
	$gl_name=~ s/(.*)\.fa/$1/;
	print GL ">$gl_name\n";
	while(<IN2>)
	{
		chomp;
	    s/\borf25\b/atp4/ig;
		s/\borfB\b/atp8/ig;
		s/\borfX\b/mttB/ig;
		s/NAD/nad/g;
		s/ATP/atp/g;
		s/COX/cox/g;
		s/ccmFN\d+/ccmFn/ig;
		s/ccmFN/ccmFn/ig;
		s/ccmFC/ccmFc/ig;
		s/atp1-\d/atp1/g;
		s/atp6-\d/atp6/g;
		s/atp8-\d/atp8/g;
		s/ccmC-\d/ccmC/g;
		s/cob-\d/cob/g;
		s/cox2-\d/cox2/g;
		s/nad4L-\d/nad4L/g;
		s/rps19-\d/rps19/;
		s/rps7-\d/rps7/;
		s/sdh3-\d/sdh3/g;
#		s/mat-r/matR/g;
		s/mat-R/matR/g;
		s/\bcoxI\b/cox1/g;
		s/\bcoxIII\b/cox3/g;
		s/tatC/mttB/g;
#		s/nad4l/nad4L/g;
		if($_=~ /^>lcl.*gene=orf.*hypothetical protein.*/){
			next;
		}elsif($_ =~ /^>lcl.*gene=(.*?)\]/){
			print OUT1 ">$1\n";
			print GL "$1\n";
#			push(@gene_name,$1);
		}else{
			print OUT1 "$_\n";
			next;
#			push(@species_name,$_);
		}
	}
	print GL "\n";
	close IN2;
	close OUT1;
}
close IN1;
close GL;

open GE,"<gene_list.txt";
open SS,"> $out" || die $!;

my (@species_name,@gene_name);
while(<GE>)
{
	chomp;
	if($_ =~ />(.*)/)
	{
		push(@species_name,$1);
	}elsif($_ =~ /orf/ig){
		next;
	}elsif($_ =~ /^[a-zA-Z]/){
		push(@gene_name,$_);
	}
}
print SS "Gene\t";
foreach (@species_name){
	/([A-Z]).*_([a-z].*)/;
	print SS "$1_$2\t";
}
my %seen = ( ); 
foreach my $item (@gene_name) { 
#		print "$item\n";
		$seen{$item}++; 
} 
my @uniq = sort keys %seen; 
print SS "\n";
foreach my $gname (@uniq)	
{
	print SS "$gname";
	foreach my $sname (@species_name)
	{
			my $temp=0;
			open TS,"<trans$sname.fa";
			while(<TS>)
			{
				chomp;
				s/>(.*)/$1/;
				if($gname eq $_)
				{
					$temp++;
				}
			}
			if($temp > 0){
				print SS "\t$temp";  #输出具体的基因数量
			}else{
				print SS "\t0";
			}
#			close TS;
	}
	print SS "\n";
}
close TS;
close SS;
system("mkdir Trans && mv trans* Trans");
system("rm -rf gene_list.txt");




#############################################
##END
##程序结束
#############################################
