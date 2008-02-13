#!/usr/bin/perl -w

# Script to retrieve subsequences for the inner, outer and transmembrane moieties
# of a protein sequence annotated automatically by the transmembrane predictor HMMTOP.
#
# Copyright by Juliane J�pel and Steffen M�ller 2005

use strict;
use Getopt::Long;

my ($hmmtop,$sequences,$corename,$destdir,$debug,$membranesep,$filesuffixsep)=(undef,undef,undef,undef,undef,"@",".");

GetOptions("hmmtop=s"     => \$hmmtop,
	   "sequences=s"  => \$sequences,
	   "corename=s"   => \$corename,
	   "destdir=s"   => \$destdir,
	   "sep=s"        => \$membranesep,
	   "filesuffixsep=s"        => \$filesuffixsep,
	   "debug"        => \$debug)
or Getopt::Long::HelpMessage(2);
print STDOUT "  Get it started\n";
print STDERR "  Parameter of hmmtop2sequences:\n";
print STDERR "	hmmtop: $hmmtop\n" if defined($hmmtop);
print STDERR "	sequences: $sequences\n" if defined($sequences);
print STDERR "	corename: $corename\n" if defined($corename);
print STDERR "	membranese: $membranesep\n" if defined($membranesep);
print STDERR "	filesuffixsep: $filesuffixsep\n" if defined($filesuffixsep);
print STDERR "	debug: $debug\n" if defined($debug);




# first argument is the HMMTOP annotation
open(HMMTOP,"<",$hmmtop) or die "Could not open file $hmmtop.\n";
# second argument is the FASTA file that was also the source of the HMMTOP annotation
open(SEQ,"<",$sequences) or die "Could not open file $sequences.\n";

if (!defined($corename)) {
	($corename)=split(/.complete./,$sequences);
	print STDERR "Corename not defined - choosing '$corename'.\n";
}

# transmembrane regions
my %predictions;
# sidedness of insertion
my %sides;

# READ result of HMMTOP predictions and store results in hashes predictions and sides
while(<HMMTOP>) {
	my ($nterm,$numtm,$positions);
	my $name="";
	if (($name,$nterm,$numtm,$positions) = $_ =~ />HP:\s+\d+\s+(.*)\s+(OUT|IN)\s+(\d+)\s*(.*)/) {
		my ($id) = $name =~ /^(\S+)/m;
		$id =~ s/^>//;
		if (0==$numtm) {
			print STDERR "$id is soluble.\n";
		}
		else {
			print "$id,$nterm,$numtm,$positions\n";
			$predictions{$id}=$positions;
			$sides{$id}=$nterm;
			#print STDERR join(",",@pos);
		}
	}
	else {
		print "Line does not match regular expression: $_";
	}
}

$/=">";
my $sep=".";
open (OUT, ">$destdir/${corename}${filesuffixsep}out.fasta")  or die ("Could not open file $destdir/${corename}${filesuffixsep}out.fasta"." for writing.\n");
open (TM,  ">$destdir/${corename}${filesuffixsep}tm.fasta")   or die ("Could not open file $destdir/${corename}${filesuffixsep}tm.fasta"." for writing.\n");
open (TMIN,">$destdir/${corename}${filesuffixsep}tmin.fasta") or die ("Could not open file $destdir/${corename}${filesuffixsep}tmin.fasta"." for writing.\n");
open (IN,  ">$destdir/${corename}${filesuffixsep}in.fasta")   or die ("Could not open file $destdir/${corename}${filesuffixsep}in.fasta"." for writing.\n");

while (<SEQ>) {
	next if length($_) le 10;
	my @lines= split(/\n/,$_);
	my $tmp_id=shift @lines;
	my @ids= split(/[\t ]+/,$tmp_id);
	my $id = shift @ids;
	$id =~ s/^\>//;
	my $seq=join("",@lines);
	$seq =~ s/>$//;
	if (exists($predictions{$id})) {
		my $side=$sides{$id};
		my @pos=split(/\s+/,$predictions{$id});
		
		my $post;
		my $lastpos=0;
		my ($num_in,$num_out,$num_tmin,$num_tm)=(0,0,0,0);
		print OUT ">$id\nb";
		print TM  ">$id\nb";
		print IN  ">$id\nb";
		print TMIN  ">$id\nb";
		foreach my $b (@pos) {
			my $a=$lastpos;
			if (($post)=$side =~ /TM-(OUT|IN)/) {
				if ($a - $b > 25) {
					print STDERR " ! ! ! ! Programming error ! ! ! !\n";
					print STDERR " No such long membrane spanning regions are anticiapted.\n\n";
				}
				$side=$post;
				my $s=substr($seq,$a-1,$b-$a+1);
				print TM $membranesep if $num_tm;
				print TM $s;
				$num_tm++;
				if ($num_tmin>0 && $side =~ /TM-IN/) {
					print TMIN $membranesep;
				}
				print TMIN $s;
				if ($num_tmin>0 && $side =~ /TM-IN/) {
					$num_tmin++;
				}
			}
			elsif ($side =~ /IN/) {
				$a++;
				my $s=substr($seq,$a-1,$b-1-$a+1);
				$side="TM-OUT";
				print IN $membranesep if $num_in;
				print IN $s;
				print TMIN $s;
				$num_in++;
			}
			elsif ($side =~ /OUT/) {
				$a++;
				my $s=substr($seq,$a-1,$b-1-$a+1);
				$side="TM-IN";
				print OUT $membranesep if $num_out;
				print OUT $s;
				$num_out++;
			}
			else {
				print STDERR "Programming error: state $side does not "
				             . "exist in my plan.\n";
			}
			$lastpos=$b;
		}
		my $a=$lastpos;
		if ($side =~ /IN/) {
			$a++;
			my $s=substr($seq,$a-1);
			$side="TM-OUT";
			print IN $membranesep if $num_in;
			print IN $s;
		}
		elsif ($side =~ /OUT/) {
			$a++;
			my $s=substr($seq,$a-1);
			$side="TM-IN";
			print OUT $membranesep if $num_out;
			print OUT $s;
		}
		else {
			print STDERR "Programming error: state $side does not "
					. "exist in my plan.\n";
		}
		print OUT  "z\n";
		print IN   "z\n";
		print TM   "z\n";
		print TMIN "z\n";
	}
	else {
		print STDERR "Could not find predictions for id '$id'.\n";
	}
	$seq =~ tr/ \t\r\n//d;
	print "$id:$seq\n" if defined($debug);  
}

close OUT;
close TM;
close IN;
close HMMTOP;
close SEQ;

__END__

=head1 NAME

hmmtop2sequences -  Script to retrieve subsequences for the inner, outer and transmembrane moieties of a protein sequence annotated automatically by the transmembrane predictor HMMTOP.

=head1 SYNOPSIS

hmmtop2sequences [options] --hmmtop <hmmtop results file> --sequences <sequence file> [--corename <common name of output file>]

Options:

   --help		brief help message
   --sep <char>         separator 

=head1 OPTIONS

=over 8

=item B<--hmmtop>

Specifies the HMMTOP output file to be interpreted.

=item B<--sequences>

Specifies the FASTA sequences from which the HMMTOP annotation was derived. The order of files must be identical.

=item B<--corename>

Specifies the core filename of the output file to which "_in.fasta","_out.fasta","_tm.fasta" and "_tmin.fasta" are appended. The last stands for a combined presentation of transmembraneous and intracellular domains. The separator then only separates the outer moieties.

=item B<--sep>

The separator for membraneous and soluble fragments.

=item B<--help>

Print a brief help message to standard output and exits.

=back

=head1 DESCRIPTION

hmmtop2sequences interprets the output of the membrane topology prediction program HMMTOP in order to return the intracellular, extracellular and transmembraneous subsequences.

=head1 AUTHOR

Juliane J�pel and Steffen M�ller <moeller@inb.uni-luebeck.de>

=head1 COPYRIGHT

This programs is Copyright 2005, by Juliane and Steffen.

=cut
