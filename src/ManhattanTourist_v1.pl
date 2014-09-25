#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
no strict 'refs';

$|=1; ## flush out the STDOUT buffer immediately

##############
## Get Opts ##
##############
use Getopt::Std;
use vars qw(
        $opt_h
        $opt_m
);


#########
# usage #
#########
my $usage = <<EOF;

Manhattan Tourist Problem Algorithm v.1

purpose: 	Find longest path between two vertices in a weighted DAG. Outputs longest path in G from source to sink.
		Input matrices can adopt any size but only certain format is accepted (see README)
         
usage: ManhattanTourist_v1.pl [options]

options:

  -h	this help screen
  -m	<path to matrix filename>  default: any file with *.asp extension 
	(for FORMAT INFO of matrix file --> see README)


Author: Fiorella Schischlik (fiorella.schischlik\@gmail.com)

For more information please see the README
 
EOF

###########################
# command line processing #
#####################################################
## check for valid arguments in commandline input  ##
#####################################################

getopts('hm:');

####### check if help was requested
die "$usage" if $opt_h;

####### check optional parameters and use defaults if not provided
$opt_h = 0                                 	 unless $opt_h;
$opt_m = `ls *.asp`                      		 unless $opt_m;   # default: any file with a *.matrix extension

print "\n";
print "++++++++++++++++\n";
print "+Options used: + \n";
print "++++++++++++++++\n";
print "h: $opt_h\n";
print "m: $opt_m";
print "\n";

####################################
# declare and specify subroutines  #
####################################

sub trim($);
sub max(@);

##### trim() function removes whitespace from the start and end of the string
##### trim() function was adapted from http://www.somacon.com/p114.php
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s*//; #changed, original s/^\s+//
	$string =~ s/\s*$//;
	return $string;
}

##### max() function, searches for the highest value in an array of numbers
sub max(@) 
{
	my @dir = @_;
	my $max = $dir[0];
	for(my $i=1;$i<@dir;$i++){
		if($dir[$i]>$max){
			$max=$dir[$i];
		}
		
	}
	return $max;
}


####################################
######## read in matrix file #######
####################################

##### declare global variables for symbolic references	 (see $matrix_name in push())
our @G_down;
our @G_right;
our @G_diag;

##### declare private variables
my $matrix_name;
	
open (MATRIX, $opt_m) || die "Could not open $opt_m: $!"; 
	
while (my $line = <MATRIX>) { 
		
	next if $line =~ /^-/; ##### skip "---" lines /^\s*$/

	if ($line =~ /^([A-Za-z]\w*)/) { 
			$matrix_name = $1;
	} 
		
	else { 
			my $trimLine=trim($line);
		
			my (@row) = split (/\s+/, $trimLine); 
			push(@{$matrix_name}, \@row); ##### $matrix_name is a symbolic reference
	} 
	
} 

close(MATRIX); 

############################################
######## Manhattan Tourist Algorithm #######
############################################

##### declare and initialize dynamic programming @matrix
my @matrix = ([0],[0]);

#### 
#my $dim_down = @G_down;
#my $dim_right = @G_right;

#print "dim G_down: $dim_down\n";
#print "dim G_right: $dim_right\n";

for (my $i=0; $i<=@G_diag; $i++){
	for (my $j=0; $j<=@G_diag; $j++){ ##### matrix is filled per row
		$matrix[$i][$j]=0;
	}
}

##### calculate all $matrix[$i][0] = first column
for(my $i=1;$i<=@G_diag;$i++){

	$matrix[$i][0]=$matrix[$i-1][0]+$G_down[$i-1][0];
}

##### calculate all $matrix[0][$i] = first row
for(my $j=1;$j<=@G_diag;$j++){

	$matrix[0][$j]=$matrix[0][$j-1]+$G_right[0][$j-1];
}


say("+++++++++++++++++++++++++++++++++++++++++++++++++++++++");
say("+Longest path matrix with first line and first column:+");
say("+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");

foreach my $ref_row (@matrix){
	foreach my $column (@$ref_row) { print "$column". "\t"; }
	print "\n";
}	

##### calculate s for all $matrix[$i][$j]
my @dir = 0; #### stores the values for path caluclated from north, west and diagonal
my $max  = 0; #### max value of path from north, west, diagonal

for(my $i=1;$i<@matrix;$i++){
	
	for(my $j=1;$j<@matrix;$j++){
		
		$dir[0]=$matrix[$j][$i-1]+$G_right[$j][$i-1];
		$dir[1]=$matrix[$j-1][$i]+$G_down[$j-1][$i];
		$dir[2]=$matrix[$j-1][$i-1]+$G_diag[$j-1][$i-1];
	
		##### $matrix[$i][$j] is the max number calculated from path north, west, diagonal
		$max = max(@dir);
		$matrix[$j][$i]=$max;

	}
}

say("\n+++++++++++++++++++++++");
say("+Longest path matrix: +");
say("+++++++++++++++++++++++\n");

foreach my $ref_row (@matrix){
	foreach my $column (@$ref_row) { print "$column"."\t"; }
	print "\n";
}

exit; 


