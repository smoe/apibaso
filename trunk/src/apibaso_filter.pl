#!/usr/bin/perl -w

use strict; # avoids errors that might otherwise be tolerated

use Getopt::Long;
my ($xmlfile,$format,$db,$debug,$destdir)
  =(undef,   undef,  undef, 0,  undef);
GetOptions("xmlfile=s"   => \$xmlfile,
           "format=s"    => \$format,
           "db=s"        => \$db,
	   "destdir=s"	 => \$destdir,
	   "debug"       => \$debug)
or Getopt::Long::HelpMessage(2);

# BioPerl knows how to retrieve protein sequences from the Internet
use Bio::SeqIO;

# XML knows how to "understand" our data
use XML::Parser;
use XML::SimpleObject;

# the input file is for now set statically
$xmlfile = 'apibaso_test.xml' unless defined($xmlfile);

$destdir = "." unless defined($destdir);

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

die "Destination directory '$destdir' does not exist.\n"
	unless -d "$destdir";

my $APICALFNAME="$destdir/apical.complete.fasta";
my $BASOFNAME="$destdir/basolateral.complete.fasta";

my $f=undef;
if ( -r "$APICALFNAME" ) {
	unlink($APICALFNAME) or die "Could not remove file $APICALFNAME\n";
	print STDERR "Cleared file '$APICALFNAME'.\n";
}
if ( -r "basolateral.complete.fasta") {
	unlink("$BASOFNAME") or die "Could not remove file $BASOFNAME";
	print STDERR "Cleared file '$BASOFNAME'.\n";
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
    	open($f, ">>$APICALFNAME")      or die "Could not open file $APICALFNAME for writing\n";
    }
    elsif ($variant->attribute('location') eq "basolateral") {
    # sonst jenes
	open($f,">>$BASOFNAME") or die "Could not open file $BASOFNAME for writing\n";
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

apibaso_filter --xmlfile apibaso.xml 

=head1 OPTIONS

=over 4

=item xmlfile <file> : specifies the file to be parsed

=format <val> : paramter expected by BioPerl, defaults to fasta

=item db <name> : parameter expected by BioPerl, defaults to swissprot

=item destdir <path> : specification of folder in which to store the result files

=item  debug : set for increased verbosity

=back


=head1 DESCRIPTION

The XML file describing protein location refers to the UniProt accession
number and entry IDs. The protein sequence however is not explicitly
stored. This script reads the entry with routines provided by the BioPerl
library and downloads the sequences online.

The script creates two files, apical.complete.fasta and basolateral.complete.fasta
in the directory that is specified by the destdir paramter.

=head1 AUTHOR

Juliane Jäpel and Steffen Möller <moeller@inb.uni-luebeck.de>

=head1 COPYRIGHT

This programs is Copyright 2005-2008, by Juliane and Steffen.

=cut

