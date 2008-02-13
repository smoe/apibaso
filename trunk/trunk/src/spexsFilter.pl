#!/usr/bin/perl -w

# Read in reference sequences

use strict;
use Getopt::Long;

my ($inputfile,$spexs_outputfile_a)
  =(undef,undef);
GetOptions(
	   "reference=s"      => \$inputfile,
	   "spexsfile=s"      => \$spexs_outputfile_a,
	)
or Getopt::Long::HelpMessage(2);


open(REFFILE,"<",$inputfile) or die "Could not open list of incoming patterns \n";
open(SPEXSFILE,"<",$spexs_outputfile_a) or die "Could not open .spexs file with spexs_output '$spexs_outputfile_a'\n";




$/=">";
my @seqs;
my @patterns;
my $debug=1;

while(<REFFILE>) {
	chomp;
	next if length($_)<4;
	my @lines=split(/\n/,$_);
	my $head=shift(@lines);
	my $seq=join("\n",@lines);
	$seq =~ s/\s*;[^\n]*//gm;
#	$seq =~ s/[ \0\n\r\t]+//gm;	
	$seq =~ s/>$//;	
	$seq =~ s/[\W]+//gm;	
	my @tmp=split(/ /,$head);
	my $prot=shift(@tmp);
	print STDERR $prot," ",$seq if $debug;
	if ($seq =~ /[A-Z]/) {
		push @seqs,$seq;
		print STDERR "$seq ok" if $debug;
	}
	else {
		print STDERR "ignored" if $debug;
	}
	print STDERR "\n" if $debug;
}
close(REFFILE);

# Read in SPEXS_OUTPUT file





$/="\n";




while(<SPEXSFILE>) {
	next if /^#/;
	my $line=$_;
	my @patterns=split(/ /,$line);
	my $found=0;
	my $pattern=shift(@patterns);
	next unless defined($pattern);
	next unless length($pattern)>1;
	# question mark bug of spexs?
	$pattern =~ s/\?//g;
	foreach my $seq (@seqs) {
		if ($seq =~ /$pattern/) {
			$found=1;
			print STDERR "found sequence '$seq' by '$pattern' \n";
			last:
		}
	}
	if ($found) {
		print $pattern." ".join(" ",@patterns);
	}
}
close(SPEXSFILE);