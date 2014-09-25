#!/usr/bin/perl

use strict;
use warnings;
use 5.010;


$|=1; ## flush out the STDOUT buffer immediately

##############
## Get Opts ##
##############
use Getopt::Std;
use vars qw(
	$opt_h
	$opt_s
	$opt_m
	$opt_g
	$opt_d
	$opt_o
);


#########
# usage #
#########
my $usage = <<EOF;

Needleman-Wunsch Algorithm v.1.2

purpose: 	Global sequence alignment over the entire length of two sequences. Suitable when the two sequences
		are of similar length, with a significant degree of similiarity throughout.
	

bugs: 	If default values for -s and -m are used, please make sure only 1 *.fasta file and 1 *.matrix 
	file exist in the same directory.	
         
usage: NeedlemanWunsch_v1.2.pl [options]

options:

  -h	this help screen
  -s	<path to sequence filename>  default: any file with *.fasta extension 
	(only the 1st and 2nd sequence will be considered --> see README)
  -m	<path to matrix filename>  default: any file with *.matrix extension 
	(for FORMAT INFO of matrix file --> see README)
  -g	gap penalty [integer] default:-2
  -d	debug [0,1] default: 0 (no debug)
  -o	<name of output file> default: STDOUT (not implemented in this version)


     
      
Author: Fiorella Schischlik (fiorella.schischlik\@gmail.com)

For more information please see the README
 
EOF

###########################
# command line processing #
#####################################################
## check for valid arguments in commandline input  ##
#####################################################

getopts('hs:m:g:o:d:');

####### check if help was requested
die "$usage" if $opt_h;

####### check optional parameters and use defaults if not provided
$opt_h = 0					unless $opt_h;
$opt_m = `ls *.matrix`				unless $opt_m;   # default: any file with a *.matrix extension
$opt_s = `ls *.fasta`				unless $opt_s;   # default: any file with a *.fasta extension
$opt_g = -2                                 	unless $opt_g;   # default: -2
$opt_d = 0					unless $opt_d;   # default: 0=no debug 1=debug modus on
$opt_o = 0					unless $opt_o;   # default: print to STDOUT

print "\n";
print "++++++++++++++++\n";
print "+Options used: + \n";
print "++++++++++++++++\n";
print "h: $opt_h\n";
print "s: $opt_s";
print "m: $opt_m";
print "g: $opt_g\n";
print "d: $opt_d\n";
print "o: $opt_o\n";

##### debug modus recommended
my $debug = $opt_d;

##################################
# open and read in matrix file   #
##################################

my ( %score, @temp, $line, @header );

open (SCORE, $opt_m) || die "Could not open $opt_m: $!"; 

$line = <SCORE>;
chomp ($line);
@header = split (/\t/, $line);
while ( $line = <SCORE> ) {
    chomp ($line);
    @temp = split /\t/, $line;
    #$score{$temp[0]} = {};
    for ( my $i = 1 ; $i <= $#temp ; $i++ ) {
        $score{ $temp[0] }{ $header[$i] } = $temp[$i];
    }
}

close(SCORE);

##################################
# open and read in sequence file #
##################################

unless (open SEQFILE, "<$opt_s"){
	die "Could not open $opt_s: $!";}
####### read by FASTA record
local $/ = "\n>";  
####### only 1st and 2nd sequence are read into array @seqS and @seqT
my $switch = 1;	
my @seqS;
my @seqT;

while (my $seq = <SEQFILE>) {
	chomp $seq;
	####### remove FASTA header
    	$seq =~ s/^>*.+\n//;  
	####### remove endlines
    	$seq =~ s/\n//g;  
		if($switch){@seqS = split(//,$seq);}
		else {@seqT=split(//,$seq);last;}
		$switch = 0;
}

print "\n";
print "++++++++++++++++++\n";
print "+Sequences used: + \n";
print "++++++++++++++++++\n";
print "seqS: @seqS"."\n";
print "seqT: @seqT"."\n";

close(SEQFILE);

###############################
# Needleman-Wunsch Algorithm  #
#####################################################################
# Initialization of dynamic programming matrix and traceback matrix #
#####################################################################

my $gap = $opt_g;
my @matAlign;
my @traceBack;

for (my $i=0; $i<=@seqT; $i++){
	for (my $j=0; $j<=@seqS; $j++){ #####  matrix is filled per row
		$matAlign[$i][$j]=0;
		$traceBack[$i][$j]=0;
	}
}

if($debug){print "DIM seqT: ".@seqT." DIM seqS: ".@seqS."\n";}

##################################################################################################
# 1.) First row and first column of the score and traceback matrix are filled during initiatzion #
##################################################################################################

#### ROW
for(my $j=0;$j<=@seqS;$j++){

	$matAlign[0][$j]=$j*$gap;
	$traceBack[0][$j]=2;
}
##### COLUMN
for(my $i=0;$i<=@seqT;$i++){

	$matAlign[$i][0]=$i*$gap;
	$traceBack[$i][0]=1;
}

if($debug){print "DIM matAlign: ".@matAlign." DIM traceBack: ".@traceBack."\n";}

#### print to file if specificed in command line -o 
#### TO DO: this feature is not implemented, output to STDOUT

 			print "+++++++++++++++++++++++++++++++\n";
			print "+ DYNAMIC PROGRAMMING MATRIX:+ \n";
			print "+++++++++++++++++++++++++++++++\n";
			foreach my $ref_row (@matAlign){
				foreach my $column (@$ref_row) { print "$column"."\t"; }
				print "\n";
			}
			print "+++++++++++++++++++++++++++++++\n";
			print "+ TRACEBACK MATRIX:           + \n";
			print "+++++++++++++++++++++++++++++++\n";
			foreach my $ref_row (@traceBack){
				foreach my $column (@$ref_row) { print "$column"."\t"; }
				print "\n";
			}
		

###################################################################
# 2.) Calculation of the scores and filling the traceback matrix  #
###################################################################

##### @dir array stores the 0=from North [up]; 1=West [left]; 2=from NorthWest [diag]
my @dir=0;
my $pos=0;

##### declare subroutines
sub max(@);

##### subroutines
sub max(@) 
{
	#### Return the position in array of max number in array
	#### start with diag number as max number
	my @dir = @_;
	my $max = $dir[2];
	my $pos=2;
	
	for(my $i=1;$i>=0;$i--){
		
		if($dir[$i]>$max){
			$max=$dir[$i];
			$pos=$i;
		}
	}
	return $pos;
}

##### start loop
for(my $i=1;$i<=@seqT;$i++){
	
	for(my $j=1;$j<=@seqS;$j++){	

		if ($debug) {	print "seq: ".$seqS[$j-1]." <--> ".$seqT[$i-1];
			      	print " coord: i:".$i." j:".$j;}
		
		$dir[0]=$matAlign[$i-1][$j]+($gap); ## left
		$dir[1]=$matAlign[$i][$j-1]+($gap); ## up
		$dir[2]=$matAlign[$i-1][$j-1]+($score{$seqS[$j-1]}{$seqT[$i-1]}); ## diag

		$pos =max(@dir);

		if ($debug) {	print "\t[up]".$dir[0]."\t[left]".$dir[1]."\t[diag]".$dir[2];
				print " traceback matrix: ".($pos+1)."\n";}
		## matrix filling
		$matAlign[$i][$j]=$dir[$pos];
		$traceBack[$i][$j]=($pos+1);
	}
}

	print "+++++++++++++++++++++++++++++++\n";
	print "+ DYNAMIC PROGRAMMING MATRIX:+ \n";
	print "+++++++++++++++++++++++++++++++\n";
	foreach my $ref_row (@matAlign){
		foreach my $column (@$ref_row) { print "$column"."\t"; }
		print "\n";
	}

	$traceBack[0][0] = 0; ##### 0 for the stop sign (for traceback)
	
	print "+++++++++++++++++++++++++++++++\n";
	print "+ TRACEBACK MATRIX:           + \n";
	print "+++++++++++++++++++++++++++++++\n";
	##### the finished traceBack matrix has a 0=stop 1=TO North [up]; 2=TO West [left]; 3=TO NorthWest [diag]
	foreach my $ref_row (@traceBack){
		foreach my $column (@$ref_row) { print "$column "; }
		print "\n";
	}

##########################################################
## 3.) Deducing the alignment from the traceback matrix ##
##########################################################

##### start deducing alignment from southwest point of matrix
my @alignS;
my @alignT;
my $posx = @seqS; #### stores the current position in the matrix
my $posy = @seqT;

if($debug) {	my $dim_traceBack = @traceBack;
		print "DIM OF TRACEBACK MATRIX ".$dim_traceBack."\n";}

##### endless loop stops when northeast point of matrix is reached 
##### the traceback matrix has a 0 at $traceBack[0][0]
##### sequences are aligned backwards
while (1){
	
	##### if position $traceBack[0][0] is reached, programm finishes
	if ($traceBack[$posy][$posx] == 0){
		last;
	}
	
	else{
		
		if ($traceBack[$posy][$posx]==2){

				if($debug) {	print "TRACEBACK-MATRIX [2] and y: ".$posy." und x: ".$posx." value: ".$traceBack[$posy][$posx]."\n";
						print "seqS von y an pos: ".$seqS[$posy-1]."\n";
						print "seqT von x an pos: ".$seqT[$posx-1]."\n";}
			#### if max n is 2 [left] a gap is introduced in the $seqT
			unshift(@alignT,'-');
			unshift(@alignS,$seqS[$posx-1]);
		
				if($debug) {	print "posx: ".$posx."\n";
						print "posy: ".$posy."\n";
						say("S: @alignS");
						say("T: @alignT");
						say("");}
			$posx--;
			
		}
	
		if ($traceBack[$posy][$posx]==1){

				if($debug) {	print "TRACEBACK-MATRIX [1] and y: ".$posy." und x: ".$posx." value: ".$traceBack[$posy][$posx]."\n";
						print "seqS von y an pos: ".$seqS[$posx-1]."\n";
						print "seqT von x an pos: ".$seqT[$posy-1]."\n";}
			#### if max n is 1 [up] a gap is introduced in the $seqS
			unshift(@alignT,$seqT[$posy-1]);
			unshift(@alignS,'-');
		
				if($debug) {	print "posx: ".$posx."\n";
						print "posy: ".$posy."\n";
						say("S: @alignS");
						say("T: @alignT");
						say("");}
			$posy--;
				
		}
	
		if ($traceBack[$posy][$posx]==3){

				if($debug) {	print "TRACEBACK-MATRIX [3] and y: ".$posy." und y: ".$posx." value: ".$traceBack[$posy][$posx]."\n";
						print "seqS von y an pos: ".$seqS[$posy-1]."\n";
						print "seqT von x an pos: ".$seqT[$posx-1]."\n";}
			#### if max is 3 [diag] letters from two sequences are aligned
			unshift(@alignS,$seqS[$posx-1]);
			unshift(@alignT,$seqT[$posy-1]);
		
				if($debug) {	print "posx: ".$posx."\n";
						print "posy: ".$posy."\n";
						say("S: @alignS");
						say("T: @alignT");
						say("");}
		
			$posx--;
			$posy--;
		}
	}	
		
}

print "+++++++++++++++++++++++++++++++\n";
print "++ ORIGINAL SEQUENCE:        ++ \n";
print "+++++++++++++++++++++++++++++++\n";
say("@seqT");
say("@seqS");

print "+++++++++++++++++++++++++++++++\n";
print "++ ALIGNED SEQUENCE:         ++ \n";
print "+++++++++++++++++++++++++++++++\n";
say("@alignT");
say("@alignS");

exit;
