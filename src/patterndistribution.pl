#!/usr/bin/perl -w

# Betrachte grafisch die Zuordnung bestimmter Muster zu ihren Sequenzen
# Ermoegliche auch die Ausgabe in einem Format, das durch die Statistik-
# Umgebung R gelesen werden kann.

use strict;
use Getopt::Long;

my ($patternfile,$sequencefile,$minRatio,$minSeqsHit,$filesuffixsep,$debug)
  =(undef,       undef,        undef,    undef,     ,"."	   , 0);
GetOptions(
	   "patterns=s"      => \$patternfile,
           "sequences=s"     => \$sequencefile,
	   "minRatio=i"	     => \$minRatio,
	   "minSeqsHit=i"    => \$minSeqsHit,
           "filesuffixsep=s" => \$filesuffixsep,
	   "debug"           => \$debug
	   )
or Getopt::Long::HelpMessage(2);

#filesuffixsep is unused for now


open(PATTERNS,"<",$patternfile) or die "Could not open list of patterns '$patternfile'\n";
open(SEQUENCES,"<",$sequencefile) or die "Could not open fasta file with sequences '$sequencefile'\n";

$minRatio=7 if !defined($minRatio);
$minSeqsHit=5 if !defined($minSeqsHit);

my %patterns;
my $npattern=0;
my $npatternAccepted=0;

while(<PATTERNS>) {
	chomp;
	next if /#/;
	next if /^\s/;
	next if /^$/;
	$npattern++;
	my $l=$_;
	my ($p,$sum,$ratio)=$l=~/(\S+)\s+(\d+).*R:([0-9.E+-]+)/i
		or die "Could not interpret line $l\n";
	#print STDERR "$p,$sum,$ratio\n";
	my ($seqsHit)=$sum   =~ /^(\d+)/ or die "Could not read sum ($sum) in $l\n";
	next unless  $seqsHit>= $minSeqsHit;
	next unless  $ratio  >= $minRatio;
	$p =~ s/\?/./g;
	next unless $p=~/[A-Z]/;
	if (exists($patterns{$p})) {
		$patterns{$p}++;
	}
	else {
		$patterns{$p}=1;
		$npatternAccepted++;
	}
}; 
close(PATTERNS);
print STDERR "Accepted $npatternAccepted of $npattern Patterns.\n";

$/=">";
my %sequence2no;
my %matches;
my %NumOfMatches;
my $num=0;
my $nmatches=0;
my @errors;
while(<SEQUENCES>) {
	chomp;
	next unless length($_)>4;
	$num++;

	my @lines=split(/[\n\r]+/,$_);
	my @l=split(/[ .]/,shift(@lines));
	my $name=shift(@l);
	$name =~ s/^>//;
	my $nameok=0;
	do {
		if(exists($sequence2no{$name})) {
			my $err="Attention: Name $name appears at least twice.\n";
			print STDERR $err;
			push @errors,$err;
			#$name="$name'";
			$nameok=-1;
		} else {
			$nameok=1;
		}
	} while (0==$nameok);

	next if -1 == $nameok;

	my $seq=join("",@lines);
	$sequence2no{$name}=$num;
# 	my $numMatchInSequence;
	foreach my $pat (keys %patterns) {
		if (!defined($matches{$pat})) {
			$matches{$pat}=[];
		}
# 		print STDERR "the current sequence is $seq \n";
		my @fragments = split /$pat/i, $seq;
		my $numMatchInSequence = $#fragments;
# 		print STDERR "the number of matches is $numMatchInSequence \n"; 
		if ($seq =~ /($pat)/i) {
			print $name,"\t",$pat,"\tmatched." if $debug;
# 			push (@{$NumOfMatches{$pat},$numMatchInSequence);
			my @tupel = ($name,$numMatchInSequence);
			push (@{$matches{$pat}},\@tupel);
			$nmatches++;
		}
		print "\n" if $debug;
	}
}
close(SEQUENCES);
print STDERR "Found ".($nmatches)." in $num sequences.\n";


# header only   -  lists all sequence IDs
#my $i=0;
foreach my $seqId (sort { $sequence2no{$a} <=> $sequence2no{$b} } keys %sequence2no) {
	#print "\t" if $i>0;
	#for(my $j=0; $j<$sequence2no{$seq}; $j++) {print " ";}
	print "\t\"$seqId\"";
	#print "\n";
	#$i++;
}
print "\n";

my @temparray=keys %sequence2no;
my $columlength = scalar(@temparray);

# contents
foreach my $m (sort keys %matches) {
	my @sequences=@{$matches{$m}};
	#my @positions;
	#foreach my $tmp_s (@sequences) {
#		# Converting tupel reference back to tupel
# 		my ($seqno,$nummatches)=@{$tmp_s} 
# 		push @positions,$sequence2no{$seqno};
# 	}

	my $lastpos=0;

	sub positionOfSequence($) {
		my $tmp_s=shift;
		my ($seqname,$bla)=@{$tmp_s};
		my $r=$sequence2no{$seqname};
		die "Could not retrieve '$seqname' from assignment to number.\n"
			unless defined($r);
		return($r);
	}

	sub numberOfMatchesInSequence($) {
		my $tmp_s=shift;
		my ($seqname,$numMatchInSequence)=@{$tmp_s};
		if ($numMatchInSequence>1900) {
			die "Number of matches in sequence $seqname is higher than 1900 --- error!!!\n";
		}
		return($numMatchInSequence);
	}

	my @sequencesSorted = sort {positionOfSequence($a)<=>positionOfSequence($b)} @sequences;

	print "\"$m\"";
	foreach my $seq (@sequencesSorted) {
		my $pos = positionOfSequence($seq);
		my $num = numberOfMatchesInSequence($seq);
		for(my $i=$lastpos; $i<$pos; $i++) {
			print "\t";
			if ($i<$pos-1) {
				print "0";
			}
		}
		#print "X";
# 		print $pos . ":";
# 		my $realnum = $num*0.5;
		print $num;
		$lastpos=$pos;
	}
	for(my $i=$lastpos; $i<$columlength; $i++)
	{
			print "\t";
			print "0";
	}	


	print "\n";
}

if ($#errors > -1) {
	foreach my $e (@errors) {
		print $e;
	}
	exit(-1);
}

__END__

=head1 NAME

patterndistribution -  Evaluation of patterns provided by SPEXS and creation a matrix indicating the hits.

=head1 SYNOPSIS

patterndistribution [options] --patterns <SPEXS output> --sequences <sequences to run on> 

Options:

   --help               brief help message
   --minRatio	        minimal ratio announced by SPEXS in comparson to second set
   --minSeqsHit         minimal number of sequences to match

=head1 OPTIONS

=over 8

=item B<--patterns>

The SPEXS output file.

The patterns are retrieved from the SPEXS output and over each single pattern it is decided if it should be used for the further analyses.

=item B<--sequences>

File with FASTA protein sequences on which the patterns from SPEXS are evaluated.

=item B<--minRatio>

The ratio of matches in the first set of sequences versus the second as announced from the SPEXS tool.

=item B<--minSeqsHit>

The minimal number of matches a pattern should have to be considered worthy for the analysis.

=item B<--help>

Print a brief help message to standard output and exits.

=back

=head1 DESCRIPTION

hmmtop2sequences interprets the output of the membrane topology prediction program HMMTOP in order to return the intracellular, extracellular and transmembraneous subsequences.

=head1 AUTHOR

Juliane Jäpel and Steffen Möller <moeller@inb.uni-luebeck.de>

=head1 COPYRIGHT

This programs is Copyright 2005, by Juliane and Steffen.

=cut

