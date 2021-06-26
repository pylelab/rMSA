#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $rootdir=dirname(abs_path(__FILE__));
my $bindir ="$rootdir/bin";
my $dbdir  ="$rootdir/database";
my $db1tsv ="$dbdir/rnacentral.tsv";
my $db2    ="$dbdir/nt";
my $max_split_seqs=5000;

my $docstring=<<EOF
getTaxonID.pl seq.afa seq.afa.tsv
    for rMSA format alignment seq.afa, map each hit to taxonID and output the
    result to seq.afa.tsv

getTaxonID.pl seq.afa seq.afa.tsv names.dmp
    in addition to taxonID, also map each hit to scientific name
    "names.dmp" can be downloaded by:
    \$ wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    \$ unzip taxdmp.zip names.dmp
EOF
;

if (@ARGV<2)
{
    print $docstring;
    exit(1);
}

my $infile =$ARGV[0];
my $outfile=$ARGV[1];
my $tmpfile=$outfile.".tmp";
my $namedmp="";
$namedmp   =$ARGV[2] if (@ARGV>2);

my @header_list=();
my @accession_list=();
my $accession="";
my @rc_list=();
my @nt_list=();
my $header;
my $cmd="grep -E '";
foreach $header(`grep '^>' $infile|sed 's/>//g'|cut -f1`)
{
    chomp($header);
    $accession="";
    if ($header=~/^(URS[A-Z0-9]+)/) # RNAcentral
    {
        $accession="$1";
        push(@rc_list,($accession)) if (! grep( /^$accession$/, @rc_list));
        $cmd.="$accession|";
    }
    elsif ($header=~/^([_A-Z0-9]+)[.]/)
    {
        $accession="$1";
        push(@nt_list,($accession)) if (! grep( /^$accession$/, @nt_list));
    }
    if (length $accession)
    {
        push(@header_list,($header));
        push(@accession_list,($accession));
    }
    else
    {
        print "skip unmappable entry >".$header."\n";
    }
}

printf "mapping %d RNAcentral accession(s)\n", scalar @rc_list;
my %taxon_dict;
%taxon_dict = map { $_ => "" } @accession_list;
my $taxonIDs;
if (scalar @rc_list)
{
    $cmd=substr($cmd,0,(length $cmd)-1)."' $db1tsv";
    my $success=0;
    #print "$cmd\n";
    foreach my $line(`$cmd`)
    {
        if ($line=~/^(URS[A-Z0-9]+)\s\d+\s([,\d]+)$/)
        {
            $accession="$1";
            $taxonIDs="$2";
            $taxon_dict{$accession}=$taxonIDs;
            $success++;
        }
    }
    if ($success==0)
    {
        printf "0 RNAcentral accession mapped. re-search full database.\n";
        foreach my $line(`cat $db1tsv`)
        {
            if ($line=~/^(URS[A-Z0-9]+)\s\d+\s([,\d]+)$/)
            {
                $accession="$1";
                $taxonIDs="$2";
                $taxon_dict{$accession}=$taxonIDs;
            }
        }
    }
}

printf "mapping %d NCBI nucleotide (NT) accession(s)\n", scalar @nt_list;
for (my $j=0;$j<scalar @nt_list;$j+=$max_split_seqs)
{
    #print "mapping entry $j to ".($j+$max_split_seqs)."\n";
    my $txt="";
    foreach (my $i=$j;$i<=$j+$max_split_seqs;$i++)
    {
        $txt.="$nt_list[$i]\n";
    }
    open(FP,">$tmpfile");
    print FP $txt;
    close(FP);
    foreach my $line(`$bindir/blastdbcmd -db $db2 -entry_batch $tmpfile -outfmt '%a %T'`)
    {
        if ($line=~/(\S+)\s(\S+)/)
        {
            $accession="$1";
            $taxonIDs="$2";
            $accession="$1" if ($accession=~/(\S+)[.]\d+/);
            $taxon_dict{$accession}=$taxonIDs;
        }
    }
}
system("rm $tmpfile");

my %name_dict;
my $name;
my $taxonID;
if (length $namedmp)
{
    printf "mapping scientific names\n";
    foreach my $line(`grep -P "\tscientific name\t" $namedmp`)
    {
        if ($line=~/^(\d+)\t\|\t([\S\s]+?)\t\|\t/)
        {
            $taxonID="$1";
            $name="$2";
            $name_dict{$taxonID}=$name;
        }
    }
    printf "read %d scientific names\n",scalar %name_dict;
}

printf "writing mapping file for %d hits\n", scalar @header_list;
my $txt="#accession\thit\ttaxonID\n";
$txt="#accession\thit\ttaxonID\tname\n" if (length $namedmp);
my $names;
for (my $i=0;$i<scalar @header_list;$i++)
{
    $header=$header_list[$i];
    $accession=$accession_list[$i];
    $taxonIDs=$taxon_dict{$accession};
    if (length $taxonIDs==0)
    {
        print "failed to map >$header\n";
    }
    if (length $namedmp)
    {
        $names="";
        foreach $taxonID(split(/,/,$taxonIDs))
        {
            $names.=",$name_dict{$taxonID}" if (exists($name_dict{$taxonID}));
        }
        $names=substr($names,1);
        $txt.="$accession\t$header\t$taxonIDs\t$names\n";
    }
    else
    {
        $txt.="$accession\t$header\t$taxonIDs\n";
    }
}
open(FP,">$outfile");
print FP $txt;
close(FP);

exit(0);
