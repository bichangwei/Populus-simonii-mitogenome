##Phylogeny Analysis Software of Plant Organelle
##Author: Changwei Bi, Nanjing Forestry University
##Version: 3.0
##If you have any queries, please feel free to contact us at the address below.
##E-mail:bichwei@163.com

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
my $red = "\033[0;31m";	my $end = "\033[0m";
my ($organelle_name,$gene_name,$output,$help);

GetOptions( 
	"input|i=s" => \$organelle_name,
	"gene|g=s" => \$gene_name,
	"out|o=s" => \$output,
	"help|?" => \$help);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file> <species.name file>
		-g <file> <conserved_gene.list file>
		-o <string> <fa.prefix file>
INFO

die $usage if ($help || !$organelle_name || !$gene_name || !$output);
##第一步：转换NCBI原蛋白质文件数据格式
##transformat_file.pl
######################################

open IN1,"<$organelle_name" || die $!;
##定义分隔符为换行符，因为之后会将分隔符定义为‘>’
$/="\n";
while(my $line1=<IN1>)
{
	chomp $line1;
	open IN2,"<$line1";
##将转换格式后的所有细胞器蛋白质文件放到trans开头的文件中
	open OUT1,">trans$line1";
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
##这一步只针对通过标准格式注释提交的NCBI细胞器基因组蛋白质文件
		if($_ =~ /^>lcl.*gene=(.*?)\]/){
			print OUT1 ">$1\n";
		}else{
			print OUT1 "$_\n";
		}
	}
	close IN2;
	close OUT1;
}
close IN1;

###########################################
##第二步：根据蛋白质文件筛选用于分析的基因
##select_protein.pl
###########################################

##将上一步生成的trans开头的文件名称放到trans_species.txt文件中，用于下一步分析
system('ls trans* >trans_species.txt');
open IN3,"<trans_species.txt";
##定义4个全局变量
my $name;
my $seq;
my @name;
my %name_seq;
while(my $tmp=<IN3>)
{
	chomp $tmp;
	open IN4,"<$tmp";
##第二步生成的文件放到select开头的文件中
	open OUT2,">select$tmp";
##重新定义分隔符为'>'，默认分隔符是换行符
	local $/=">";
	while(<IN4>)
	{
		chomp;
##对上一步格式转换后的文件进行处理，运用hash将名称与序列相连
		($name,$seq)=split(/\n/,$_,2);
		#print $seq,"\n";
		$name_seq{$name}=$seq;
##将所有的名称通过压栈的方式放到数组中，用于下面的提取分析
		push(@name,$name);
	}
##将分隔符还原为换行符
	local $/="\n";
##打开之前准备好的包含所有想要用来建树的基因文件
##线粒体一般不超过30个，叶绿体一般不低于50个
	open (IN5,"<$gene_name") or die("ERROR,you should input a correct gene file\n");
	while(my $line=<IN5>)
	{
		chomp $line;
##根据名字来获得基因	
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

###########################################
##第三步：将筛选出来的基因放到一个只包含物种名的文件中
##merge_all_data.pl
###########################################

##将上一步生成的select开头的文件名称放到select_species.txt文件中，用于下一步分析
system('ls select* >select_species.txt');
open IN6,"<select_species.txt";
while(my $line1=<IN6>)
{
	chomp $line1;
##统计所选择的基因是否在所有细胞器基因组中都有，如果数值都是相同的，则构建的进化树是可信的，否则是不可信的，需要重新选择基因
        open IN7,"<$line1";
	my $tj_name=$line1;
        $line1 =~ s/selecttrans(.*)\.fa/$1/;
	print  "$line1\t";
##统计基因数量，手工查看是否一致
	system("grep '>' $tj_name | wc -l ");
	print  "\n";
##用于检查运行过程中的一些bug
#	print "$line1\n";
#	print $species_name,"\n";
##最后将结果保存到以final开头的文件中
	open OUT3,">final$line1";
	print OUT3 ">$line1\n";
	while(my $line2=<IN7>)
	{
		chomp $line2; 
#		if($line2 =~ /^>\s+/)
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
#print OUT3 "\n";
close IN6;

##第四步：生成系统发育分析所需文件
system("cat final* > $output.phylogeny.fa");
###########################################

###########################################
##第五步：将程序运行过程中生成的中间文件分类
##删除一些无用的中间文件
##delete.pl
###########################################
system('mkdir Trans && mv trans* Trans');
system('mkdir Select && mv select* Select');
system('mkdir Final && mv final* Final');
#
#system('rm -rf trans*');
#system('rm -rf select*');
#system('rm -rf final*');

###########################################
##第六步：调用Muscle构建进化树
##use_muscle.pl
###########################################

#my $cmd2='/home/YTM/bcw/bin/muscle3.6_src/muscle -in phylogeny.fa -maxmb 5000 -clw -out tree.fa';
#system($cmd2);

#############################################
##END
##程序结束
#############################################
