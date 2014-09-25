#!/usr/bin/perl

use strict;
use warnings;
use 5.010;

use Time::HiRes;

$|=1; ## flush out the STDOUT buffer immediately

##############
## Get Opts ##
##############
use Getopt::Std;
use vars qw(
        $opt_h
        $opt_n
	 $opt_s
);


#########
# usage #
#########
my $usage = <<EOF;

Towers of Hanoi Algorithm v.1

purpose:	Outputs a list of moves that solves the Towers of Hanoi. If -s is specified calculation time 
		for n to m disks are printed instead of moves.
         
usage: HanoiTowers_v1.pl [options]

options:

  -h	this help screen
  -n	number of disks [positive integer, >1] default: 4 <play only version>
  -s	number of max disk to use for time statistics [positive integer] default: 0 <statistics only version>


Author: Fiorella Schischlik (fiorella.schischlik\@gmail.com)
 
EOF

###########################
# command line processing #
#####################################################
## check for valid arguments in commandline input  ##
#####################################################

getopts('hs:n:');

####### check if help was requested
die "$usage" if $opt_h;

####### check optional parameters and use defaults if not provided
$opt_h = 0                                 	 unless $opt_h;
$opt_n = 4		                        	 unless $opt_n;   # default: 4
$opt_s = 0						 unless $opt_s;

print "\n";
print "++++++++++++++++\n";
print "+Options used: + \n";
print "++++++++++++++++\n";
print "h: $opt_h\n";
print "n: $opt_n\n";
print "s: $opt_s";
print "\n";

###########################
# start algorithm         #
###########################

###### number of disks
my $ndisk=$opt_n;

###### declare subroutines
sub HanoiTowers(@);

print "\n";
print "++++++++++++++++++\n";
print "+Required Moves: + \n";
print "++++++++++++++++++\n";

###### check for correct input, false input forces program to exit
if($ndisk == 0 || $ndisk ==1){print "\nAt least 2 disks are required. Rerun script.\n"; exit;}

sub HanoiTowers(@)
{
	my $ndisk=$_[0];
	my $fromPeg=$_[1];
	my $toPeg=$_[2];
	my $unusedPeg;
	
	if ($ndisk==1){
		if ($opt_s == 0){
			say("Move disk from peg $fromPeg to peg $toPeg");}
		return;
	}
	
	else{
		$unusedPeg = 6-$fromPeg-$toPeg;
		##### HanoiTowers function calls itself with some part of the task
		HanoiTowers($ndisk-1,$fromPeg, $unusedPeg);
		if ($opt_s == 0){
			say("Move disk from peg $fromPeg to peg $toPeg");}
		HanoiTowers($ndisk-1,$unusedPeg,$toPeg);
		return;
	}
	
	
		
}

###### start recursive HanoiTower subroutine
HanoiTowers($ndisk,1,3);

###### if -s was specified
###### moves are not printed, instead calculation time is for n to m disks is printed to STDOUT

if($opt_s){
	
	my @timestamps;
	my $timeStart;
	my $timEnd;
	
	##### this loop was gratefully provided by Marlies Dolezal 
	for (my $i = $opt_n; $i <= $opt_s; $i++){
            $timeStart = Time::HiRes::gettimeofday();
            HanoiTowers($i,1,3);                   
            $timEnd = Time::HiRes::gettimeofday();
            push(@timestamps, $timEnd - $timeStart);
	}
	
	print "\n++++++++++++++++++++++++++++++++++++++++++++\n";
	print "+Calculation time increases exponentially: +\n";
	print "++++++++++++++++++++++++++++++++++++++++++++\n";
	foreach my $time (@timestamps){
		print "$time\n";
	}
}

exit;


