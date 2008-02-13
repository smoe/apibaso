#!/usr/bin/perl -w
use strict;

my $output_rownames=0;
my $output_colnames=0;

my $skip_beginning=1;
my @exps;
if (0)  {
	while(<>) {
		if ($skip_beginning) {
			if (/Set_of_Sets/ ) {
				$skip_beginning=0;
			}
			next;
		}
		next if /^#/;
		next if /class.*fication/;

		my @tmp=split(/ /,$_);
		my $exp=$tmp[0];
		next unless $exp =~ /[A-Z.]/;

		#print "$exp\n";
		push @exps,$exp;
	}
}
else {
	my $first=0;
	while(<>) {
		if ($first) {
			$first=0;
			next;
		}
		next if /^#/;
		next if /class.*fication/;
		my @tmp=split(/[ \t]/,$_);
		my $exp=$tmp[0];
		print STDERR "$exp\n";
		next unless $exp =~ /[A-Z.]/;
		push @exps,$exp;
	}
}

sub generate_random_sequence 
{
   my $passwordsize = shift;
   my @alphanumeric = ('A','C'..'I','K'..'N','P'..'T','V','W','Y');
   my $randpassword = join('',map $alphanumeric[rand @alphanumeric], 0..$passwordsize);
   return $randpassword;
}


if ($output_colnames) {
	foreach my $e (@exps) {
		print "\t$e";
	}
	print "\n";
}

for(my $i=0; $i<2000; $i++) {
   print STDERR "$i\t" if 0==$i%40;
   my $s=generate_random_sequence(15);
   print $s if $output_rownames;
   foreach my $e (@exps) {
   	print "\t", (($s=~/$e/g)?1:0);
   }
   print "\n";
}

