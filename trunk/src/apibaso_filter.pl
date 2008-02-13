#!/usr/bin/perl -w

# be careful with what you do
use strict;

use Getopt::Long;
my ($xmlfile,$format,$db,$debug)
  =(undef,   undef,  undef, 0);
GetOptions("xmlfile=s"   => \$xmlfile,
           "format=s"    => \$format,
           "db=s"        => \$db,
	   "debug"       => \$debug)
or Getopt::Long::HelpMessage(2);



# BioPerl knows how to retrieve protein sequences from the Internet
use Bio::SeqIO;

# XML knows how to "understand" our data
use XML::Parser;
use XML::SimpleObject;

# the input file is for now set statically
$xmlfile = 'apibaso_test.xml' unless defined($xmlfile);

# a parser for the input file is created
my $parser = XML::Parser->new(ErrorContext => 2, Style => "Tree");
my $xso = XML::SimpleObject->new( $parser->parsefile($xmlfile) );

# initialise export format of protein sequences
$format="fasta" unless defined($format);
$db="swissprot" unless defined($db);
my %params = ( '-format' => $format);
my $remotedb;
# create communicator with remote protein database SWISS-PROT (-> UniProt)
eval {
    require Bio::DB::SwissProt;
    $remotedb = Bio::DB::SwissProt->new();
};
die($@) unless ! $@;

my $f=undef;
if ( -r "apical.complete.fasta") {
	unlink("apical.complete.fasta")      or die "Could not remove file apical.complete.fasta\n";
}
if ( -r "basolateral.complete.fasta") {
	unlink("basolateral.complete.fasta") or die "Could not remove file basolateral.complete.fasta\n";
}

# iterate over all entries in our data
foreach my $variant ($xso->child('api:variantlist')->children('api:variant')) {
    my $err=0;
    # comments to appear on the screen - progress indicators
    print $variant->attribute('gene').":";
    print $variant->attribute('location');
    if (!defined($variant->attribute('location'))) {
    	print "  Undefined location!\n";
	$err++;
    }
    if (!defined($variant->attribute('UniProt'))) {
    	print "  Undefined UniProt accession number!\n";
	$err++;
    }
    if (0<$err) {
    	print "  Skipping entry.\n\n";
    	next;
    }
    print ",".$variant->attribute('UniProt')."\n";
    
    # show retrieved protein sequence
    my $stream;
    #if( $remotedb->can('get_Stream_by_batch') ) {
    #    $stream = $remotedb->get_Stream_by_batch($variant->attribute('UniProt'));
    #} else {
        $stream = $remotedb->get_Stream_by_acc($variant->attribute('UniProt'));
    #}

    my $f=undef;
    # schaue was in location steht
    if ($variant->attribute('location') eq "apical") {
    # output file is dies
    	open($f, ">>apical.complete.fasta")      or die "Could not open file apical.complete.fasta for writing\n";
    }
    elsif ($variant->attribute('location') eq "basolateral") {
    # sonst jenes
	open($f,">>basolateral.complete.fasta") or die "Could not open file basolateral.complete.fasta for writing\n";
    }
    # und bei anderen Sonderfällen das
    else {
        # set output file
        $f = \*STDOUT;
    }
    $params{'-fh'} = $f;

    # empty stream of sequences (it is only a single one) to the prior determined file
    my $seqio = new Bio::SeqIO(%params);
    while( my $seq = $stream->next_seq ) {
        $seqio->write_seq($seq);
    }
    if ($variant->attribute('location') eq "apical" or $variant->attribute('location') eq "basolateral") {
	    close($f);
    }
}

__END__

=head1 NAME

apibaso_filter -  Use XML file to retrieve entries from UniProt database.

=head1 SYNOPSIS

to be written

=head1 OPTIONS

to be written

=head1 DESCRIPTION

The XML file describing protein location refers to the UniProt accession number and entry IDs. The protein sequence however is not explicitly stored. This script reads the entry and downloads the sequences online.

=head1 AUTHOR

Juliane Jäpel and Steffen Möller <moeller@inb.uni-luebeck.de>

=head1 COPYRIGHT

This programs is Copyright 2005-2008, by Juliane and Steffen.

=cut

