#!/usr/bin/perl -w
##########################################################################
## FileName: Step2_Generate_the_phylogeny_fa.pl
## Version: 3.0
## Author: Changwei Bi<bichwei@163.com>
##########################################################################
use strict;
use Getopt::Long;
use File::Basename;
my $red = "\033[0;31m";	my $end = "\033[0m";
my ($gene_name,$output,$help);

GetOptions( 
#	"input|i=s" => \$organelle_name,	#file
	"gene|g=s" => \$gene_name,		#file
	"out|o=s" => \$output,			#string
	"help|?" => \$help			#string
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-g <file> <conserved_gene.list file>
		-o <string> <prefix of your final result in fa format>
INFO

die $usage if ($help || !$gene_name || !$output);

##select_protein.pl
system('cp Trans/tran* .');
system('ls trans* >trans_species.txt');
open IN3,"<trans_species.txt";
my $name;
my $seq;
my @name;
my %name_seq;
while(my $tmp=<IN3>)
{
	chomp $tmp;
	open IN4,"<$tmp";
	open OUT2,">select$tmp";
	local $/=">";
	while(<IN4>)
	{
		chomp;
		my ($name,$seq)=split(/\n/,$_,2);
		if(defined($name)){
			$name_seq{$name}=$seq;
			push(@name,$name);
		}
	}
	local $/="\n";
	open (IN5,"<$gene_name") or die("ERROR,you should input a correct gene file\n");
	while(my $line=<IN5>)
	{
		chomp $line;
		if(grep {$_ =~ /$line/i} @name)
		{
			print OUT2 ">$line\n$name_seq{$line}\n";
		}
	}
	close IN4;
	close IN5;
	close OUT2;
	undef @name;
}
close IN3;

##merge_all_data.pl
system('ls select* >select_species.txt');
open IN6,"<select_species.txt";
while(my $line1=<IN6>)
{
	chomp $line1;
        open IN7,"<$line1";
	my $tj_name=$line1;
        $line1 =~ s/selecttrans(.*)\.fa/$1/;
	print  "$line1\t";
	system("grep '>' $tj_name | wc -l ");
	print  "\n";
	open OUT3,">final$line1";
	print OUT3 ">$line1\n";
	while(my $line2=<IN7>)
	{
		chomp $line2; 
		if($line2 =~ /^>/)
		{
			next;
		}else{
			print OUT3 "$line2"; 
		}
	}
	close IN7;
	print OUT3 "\n";
	close OUT3;
}
close IN6;

system("cat final* > $output.phylogeny.fa");	#The Final phylogeny file in fasta format

#delete all temp files
system('rm -rf trans*');
system('rm -rf select*');
system('rm -rf final*');
system('rm -rf Trans*');

#run.muscle_IQTREE.sh
system("muscle -align $output.phylogeny.fa -output $output.muscle.alignment.afa");
system("iqtree -s $output.muscle.alignment.afa -m MFP -B 1000 --bnni -T AUTO"); 
