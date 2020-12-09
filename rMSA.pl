#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $rootdir=dirname(abs_path(__FILE__));
my $bindir ="$rootdir/bin";
my $dbdir  ="$rootdir/database";
my $db0    ="$dbdir/Rfam.cm";
my $db1    ="$dbdir/rnacentral.fasta";
my $db2    ="$dbdir/nt";
my $db0to1 ="$dbdir/rfam_annotations.tsv.gz";
my $db0to2 ="$dbdir/Rfam.full_region.gz";

my $docstring=<<EOF
rMSA.pl seq.fasta \\
    -db0=$db0  \\
    -db1=$db1 \\
    -db2=$db2 \\
    -db0to1=$db0to1 \\
    -db0to2=$db0to2 \\
    -cpu=1 -tmpdir=/tmp/$ENV{USER}/rMSA_`date +%N`

    for query sequence seq.fasta, output MSA to seq.afa

Input:
    seq.fasta - single sequence query fasta
    db0       - infernal covariance models for all rfam families.
    db1       - (colon separated list of) blastn format sequence database(s)
                where only the watson strand will be searched
    db2       - (colon separated list of) blastn format sequence database(s)
                where both watson and crick strand will be searched
    db0to1    - mapping file from db0 to db1
    db0to2    - mapping file from db0 to db2
    cpu       - number of threads. default 1.
    tmpdir    - temporary folder
    ssfile    - optional dot brack secondary structure prediction file.
                default is to predict ss by RNAfold.
EOF
;

my $max_rfam_num     =300;    # maximum number of rfam families to parse
my $max_split_seqs   =5000;   # max number of blastn sequences to parse per batch
                              # control tmp file size. do not change final result
my $max_target_seqs  =20000;  # max number of blastn sequences to report
my $max_aln_seqs     =100000; # max number of blastn alignmnets to parse
my $max_hhfilter_seqs=5000;   # max number of hhfilter sequences to report
my $min_hhfilter_seqs=10;     # min number of hhfilter sequences to report
my $max_mafft_qinsi  =50000;  # max $Lch * $hitnum for mafft-qinsi
my $max_mafft_xinsi  =10000;  # max $Lch * $hitnum for mafft-xinsi
my $target_Nf        =128;
my $target_Nf_cov    =60;     # only sequences with cov>=$target_Nf_cov are
                              # included for Nf calculation

my $inputfasta="";
my $ssfile    ="";
my $cpu       =1;
my $tmpdir    =""; #"/tmp/$ENV{USER}/rMSA_$$\_".`date +%N`;

foreach (my $a=0;$a<@ARGV;$a++)
{
    if    ($ARGV[$a]=~/-db0=(\S+)/)    { $db0="$1"; }
    elsif ($ARGV[$a]=~/-db1=(\S+)/)    { $db1="$1"; }
    elsif ($ARGV[$a]=~/-db2=(\S+)/)    { $db2="$1"; }
    elsif ($ARGV[$a]=~/-db0to1=(\S+)/) { $db0to1="$1"; }
    elsif ($ARGV[$a]=~/-db0to2=(\S+)/) { $db0to2="$1"; }
    elsif ($ARGV[$a]=~/-ssfile=(\S+)/) { $ssfile="$1"; }
    elsif ($ARGV[$a]=~/-cpu=(\S+)/)    { $cpu="$1"; }
    elsif ($ARGV[$a]=~/-tmpdir=(\S+)/) { $tmpdir="$1"; }
    else                               { $inputfasta=$ARGV[$a]; }
}

my @db_list =();
my @db1_list=();
my @db2_list=();
push(@db1_list, split(/:/,"$db1"));
push(@db2_list, split(/:/,"$db2"));
push(@db_list,@db1_list);
push(@db_list,@db2_list);

#### check input ####
if (length $inputfasta == 0 || (scalar @db_list == 0))
{
    print "$docstring";
    exit(1);
}

if (!-s "$inputfasta")
{
    print "ERROR! No such file $inputfasta\n";
    exit(1);
}

$inputfasta=abs_path($inputfasta);
my $prefix=dirname($inputfasta);
foreach my $basename(split(/\./,basename($inputfasta)))
{
    $prefix.="/$basename";
    last;
}

foreach my $db(@db_list)
{
    if (!-s "$db")
    {
        print "ERROR! No such file $db\n";
        exit(1);
    }
    if (!-s "$db.nal")
    {
        print "ERROR! $db not in blastn format. Please run\n";
        print "$dbdir/script/makeblastdb -in $db -parse_seqids -hash_index -dbtype nucl\n";
        exit(1);
    }
}

if ($db0)
{
    if (!-s "$db0")
    {
        print "ERROR! No such file $db0\n";
        exit(1);
    }
    elsif (!-s "$db0.i1m")
    {
        print "ERROR! $db0 not formatted by cmpress. Please run:\n";
        print "$dbdir/script/cmpress $db0\n";
        exit(1);
    }
    elsif (length "$db0to1" + length "$db0to2"==0)
    {
        print "To use -db0, either -db0to1 or -db0to2 or both must be specified\n";
        exit(1);
    }
    elsif (length "$db0to1" >0 && ! -s "$db0to1")
    {
        print "No such file $db0to1\n";
        exit(1);
    }
    elsif (length "$db0to2" >0 && ! -s "$db0to2")
    {
        print "No such file $db0to2\n";
        exit(1);
    }
}

#### make tmp ####
if (length $tmpdir==0)
{
    $tmpdir="/tmp/$ENV{USER}/rMSA_$$\_".`date +%N`;
    chomp($tmpdir);
}
if (!-d "$tmpdir")
{
    print "creating $tmpdir\n";
    system("mkdir -p $tmpdir")
}
else
{
    print "$tmpdir already exists\n";
}
$tmpdir=abs_path($tmpdir);
if ($tmpdir eq dirname($inputfasta))
{
    print "ERROR! Query fasta cannot be put inside tmpdir\n";
    exit(1);
}

#### prepare query fasta ####
my $seqnum=`grep '^>' $inputfasta|wc -l`;
if ($seqnum>=2)
{
    print "ERROR! More than one sequence in $inputfasta.\n";
    exit(1)
}
my $sequence=`$bindir/fastaNA $inputfasta | $bindir/fastaOneLine -|tail -1`;
chomp($sequence);
my $Lch=length $sequence;
print ">".abs_path($inputfasta)."\tlength=$Lch\n$sequence\n";
if ($Lch==0)
{
    print "ERROR! Sequence length 0.\n";
    exit(1)
}
open(FP,">$tmpdir/seq.fasta");
print FP ">query\n$sequence\n";
close(FP);
my $task="blastn";
$task="blastn-short" if ($Lch<30);

#### prepare secondary structure ####
if (length $ssfile>0 && !-s "$prefix.cm")
{
    if (!-s "$ssfile")
    {
        print "ERROR! No such file $ssfile\n".
        exit(1);
    }
    system("cp $ssfile $tmpdir/input.dbn");
    my $ss="";
    foreach my $line(`grep -v '^>' $ssfile|grep -ohP "^\\S+"`)
    {
        chomp($line);
        $ss.=$line;
    }
    if ($ss!~/^[.()]+$/)
    {
        print "ERROR! Only the following 3 symbols are allowed:\n";
        print "( ) .\n";
        print "Pseudoknots are not allowed\n";
        exit(1);
    }
    my $Lss=length $ss;
    if ($Lss!=$Lch)
    {
        print "ERROR! $Lss!=$Lch. inconsistency in SS and sequence length\n";
        exit(1);
    }
    my $count_left =$ss=~ tr/(//;
    my $count_right=$ss=~ tr/)//;
    if ($count_left!=$count_right)
    {
        print "ERROR! $count_right!=$count_left. number of ) different from that of (\n";
        exit(1);
    }
    open(FP,">$tmpdir/RNAfold.dbn");
    print FP "$ss\n";
    close(FP);
}


#### perform rfam search ####
print "==== rfam pre-screening ====\n";
if (length $db0>0 && (!-s "$prefix.db.gz" || `zcat $prefix.db.gz|wc -l`+0==0)
                  && (!-s "$prefix.db0.gz"||`zcat $prefix.db0.gz|wc -l`+0==0))
{
    system("$bindir/cmscan --tblout $tmpdir/cmscan.tblout -o $tmpdir/cmscan.out --noali $db0 $tmpdir/seq.fasta");
    my @family_list=();
    my $pattern="";
    foreach my $line(`grep -v '#' $tmpdir/cmscan.tblout`)
    {
        if ($line=~/^\S+\s+(\S+)/)
        {
            my $family="$1";
            next if ( grep( /^$family$/, @family_list) );
            if (scalar @family_list==0) { $pattern.= "($family)"; }
            else                        { $pattern.="|($family)"; }
            push(@family_list, ("$family"));
            print "$family\n";
            last if (scalar @family_list>$max_rfam_num);
        }
    }

    ### map rfam match to db1 ###
    # URS-Id Rfam-Model-Id Score E-value Sequence-Start Sequence-Stop
    # Model-Start Model-Stop Rfam-Model-Description
    # (0-indexed)
    ### map rfam match to db2 ###
    # rfam_acc rfamseq_acc seq_start seq_end bit_score evalue_score
    # cm_start cm_end truncated type
    for (my $dd=1;$dd<=2;$dd++)
    {
        last if (scalar @family_list==0);
        my $db0tod="$db0to1";
        $db0tod="$db0to2" if ($dd==2);
        next if (length "$db0tod"==0 || !-s "$db0tod");
        my $cat="cat";
        $cat="zcat" if ("$db0tod"=~/.gz$/);
        &System("cp $db0tod $tmpdir/db0to$dd");
        my $tabfile="$tmpdir/rfam$dd.tab";
        open(FP,">$tabfile");
        foreach my $line(`$cat $tmpdir/db0to$dd|grep -P "$pattern"`)
        {
            if (($dd==1 && $line=~/^(\S+)\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)/)||
                ($dd==2 && $line=~/^\S+\s+(\S+)\s+(\d+)\s+(\d+)/))
            {
                my $saccver ="$1";
                my $sstart  ="$2";
                my $send    ="$3";
                if ($dd==1)
                {
                    $sstart++;
                    $send++;
                }
                print FP "$saccver\t$sstart\t$send\n";
            }
        }
        close(FP);
        &System("rm $tmpdir/db0to$dd");
        &retrieveSeq($tabfile, $db1, "rfam1") if ($dd==1);
        &retrieveSeq($tabfile, $db2, "rfam2") if ($dd==2);
    }

    if (-s "$tmpdir/rfam1.db" || -s "$tmpdir/rfam2.db")
    {
        &System("cat $tmpdir/rfam1.db $tmpdir/rfam1.db > $tmpdir/db0");
        &plain2gz("$tmpdir/db0", "$prefix.db0.gz");
    }
    system("wc -l $tmpdir/rfam*.db $tmpdir/db0");
}

#### perform blastn ####
print "==== $task pre-screening ====\n";

if (-s "$prefix.db.gz" && `zcat $prefix.db.gz|wc -l`+0>0)
{
    &gz2plain("$prefix.db.gz", "$tmpdir/db");
}
else
{
    foreach (my $d=0;$d<scalar @db_list;$d++)
    {
        my $db    =$db_list[$d];
        my $strand="plus";
        $strand   ="both" if ( grep( /^$db$/, @db2_list) );
        my $tabfile="$tmpdir/blastn$d.tab";
        &System("$bindir/blastn -num_threads $cpu -query $tmpdir/seq.fasta -strand $strand -db $db -out $tabfile -task $task -max_target_seqs $max_target_seqs -outfmt '6 saccver sstart send evalue bitscore nident staxids'");
        &retrieveSeq($tabfile, $db, "blastn$d");
    }
    
    if (length $db0>0 && -s "$prefix.db0.gz")
    {
        &gz2plain("$prefix.db0.gz", "$tmpdir/db0");
    }
    &System("cat $tmpdir/db0 $tmpdir/blastn*.db > $tmpdir/trim.db");
    &rmredundant_rawseq("$tmpdir/trim.db", "$tmpdir/db");
    &plain2gz("$tmpdir/db", "$prefix.db.gz");
    if (-s "$prefix.db.gz" && `zcat $prefix.db.gz|wc -l`+0>0)
    {
        &System("rm $prefix.db0.gz");
    }
}

#### nhmmer ####
my $hitnum=`grep '^>' $tmpdir/db|wc -l`+0;
print "==== local alignment of $hitnum sequences by nhmmer ====\n";

if (-s "$prefix.nhmmer.afa.gz" && `zcat $prefix.nhmmer.afa.gz |wc -l`+0>0)
{
    &gz2plain("$prefix.nhmmer.afa.gz", "$tmpdir/nhmmer.afa");
}
else
{
    &System("$bindir/qnhmmer --noali -A $tmpdir/nhmmer.a2m --cpu $cpu --watson $tmpdir/seq.fasta $tmpdir/db | grep 'no alignment saved'");
    &addQuery2a2m("$tmpdir/nhmmer.a2m","$tmpdir/nhmmer.unfilter.afa");
    &run_hhfilter($max_hhfilter_seqs,0,"$tmpdir/nhmmer.unfilter.afa","$tmpdir/nhmmer.afa");
    &plain2gz("$tmpdir/nhmmer.afa", "$prefix.nhmmer.afa.gz");
}

my $Nf=&run_calNf("$tmpdir/nhmmer.afa");
if ($Nf>=$target_Nf)
{
    print "output nhmmer.afa (Nf>=$Nf) as final MSA\n";
    &System("cp $tmpdir/nhmmer.afa $prefix.afa");
    exit(0);
}

#### cmbuild and cmcalibrate ####
print "==== making covariance model ====\n";
if (-s "$prefix.cm")
{
    &System("cp $prefix.cm $tmpdir/infernal.cm");
}
else
{
    &System("$bindir/RNAfold --noPS $tmpdir/seq.fasta | awk '{print \$1}' | tail -n +3 > $tmpdir/RNAfold.dbn") if (!-s "$tmpdir/RNAfold.dbn");
    &System("$bindir/reformat.pl fas sto $tmpdir/nhmmer.afa $tmpdir/nhmmer.sto");

    ## reformat ss with according to gaps in reference sequence of .sto file ##
    foreach my $i(`awk '{print \$2}' $tmpdir/nhmmer.sto | head -n5 | tail -n1 | grep -b -o - | sed 's/..\$//'`)
    {
        &System("sed -i \"s/./&-/$i\" $tmpdir/RNAfold.dbn");
    }
    ## add reformated ss from last step to .sto file ##
    my $txt="";
    foreach my $line(`head -n -1 $tmpdir/nhmmer.sto`)
    {
        chomp($line);
        $line.="E=0.0" if ($line=~/^#=GF DE/);
        $txt.="$line\n";
    }
    $txt.="#=GC SS_cons                     ";
    $txt.=`cat $tmpdir/RNAfold.dbn`;
    $txt.="//\n";
    open(FP,">$tmpdir/RNAfold.sto");
    print FP $txt;
    close(FP);

    &System("$bindir/cmbuild --hand -F $tmpdir/infernal.cm $tmpdir/RNAfold.sto");
    &System("$bindir/cmcalibrate --cpu $cpu $tmpdir/infernal.cm");
    &System("cp $tmpdir/infernal.cm $prefix.cm");
}

print "==== cmsearch hits from cmscan and blastn ====\n";
if (-s "$prefix.cmsearch.afa.gz" && `zcat $prefix.cmsearch.afa.gz|wc -l`+0>0)
{
    &gz2plain("$prefix.cmsearch.afa.gz", "$tmpdir/cmsearch.afa");
}
else
{
    &System("$bindir/qcmsearch --noali -A $tmpdir/cmsearch.a2m --cpu $cpu $tmpdir/infernal.cm $tmpdir/db|grep 'no alignment saved'");
    &addQuery2a2m("$tmpdir/cmsearch.a2m","$tmpdir/cmsearch.unfilter.afa");
    &run_hhfilter($max_hhfilter_seqs,$min_hhfilter_seqs,"$tmpdir/cmsearch.unfilter.afa","$tmpdir/cmsearch.afa");
    &plain2gz("$tmpdir/cmsearch.afa", "$prefix.cmsearch.afa.gz");
}

my $Nf=&run_calNf("$tmpdir/cmsearch.afa");
if ($Nf>=$target_Nf)
{
    print "output cmsearch.afa (Nf>=$Nf) as final MSA\n";
    &System("cp $tmpdir/cmsearch.afa $prefix.afa");
    exit(0);
}

#### cmsearch for db1 and db2 ####
print "==== cmsearch hits from cmsearch ====\n";
for (my $dd=1;$dd<=2;$dd++)
{
    if (-s "$prefix.db$dd.gz" && `zcat $prefix.db$dd.gz|wc -l`+0>0)
    {
        &gz2plain("$prefix.db$dd.gz", "$tmpdir/db$dd");
    }
    else
    {
        my @db_list =();
        push(@db_list,@db1_list) if ($dd==1);
        push(@db_list,@db2_list) if ($dd==2);
        foreach (my $d=0;$d<scalar @db_list;$d++)
        {
            my $db    =$db_list[$d];
            &System("$bindir/qcmsearch --noali -A $tmpdir/cmsearch$d.$dd.a2m --cpu $cpu --incE 10.0 $tmpdir/infernal.cm $db|grep 'no alignment saved'");
            my $tabfile="$tmpdir/cmsearch$d.$dd.tab";
            open(FP,">$tabfile");
            foreach my $line(`grep '^>' $tmpdir/cmsearch$d.$dd.a2m`)
            {
                if ($line=~/^>(\S+)\/(\d+)-(\d+)/)
                {
                    my $saccver ="$1";
                    my $sstart  ="$2";
                    my $send    ="$3";
                    print FP "$saccver\t$sstart\t$send\n";
                }
            }
            close(FP);
            &retrieveSeq($tabfile, $db, "cmsearch$d.$dd");
        }
        &System("cat $tmpdir/cmsearch*.$dd.db > $tmpdir/trim.$dd.db");
        &rmredundant_rawseq("$tmpdir/trim.$dd.db", "$tmpdir/db$dd");
        &plain2gz("$tmpdir/db$dd", "$prefix.db$dd.gz");
        system("wc -l $tmpdir/cmsearch*.$dd.db $tmpdir/db$dd");
    }

    if (-s "$prefix.cmsearch.$dd.afa.gz" && `zcat $prefix.cmsearch.$dd.afa.gz|wc -l`+0>0)
    {
        &gz2plain("$prefix.cmsearch.$dd.afa.gz", "$tmpdir/cmsearch.$dd.afa");
    }
    else
    {
        &System("cat $tmpdir/db $tmpdir/db1 > $tmpdir/trimall.db");
        if ($dd==2)
        {
            &System("cat $tmpdir/db $tmpdir/db1 $tmpdir/db2 > $tmpdir/trimall.db");
        }
        &rmredundant_rawseq("$tmpdir/trimall.db", "$tmpdir/dball");
        &System("$bindir/qcmsearch --noali -A $tmpdir/cmsearch.$dd.a2m --cpu $cpu $tmpdir/infernal.cm $tmpdir/dball|grep 'no alignment saved'");
        &addQuery2a2m("$tmpdir/cmsearch.$dd.a2m","$tmpdir/cmsearch.$dd.unfilter.afa");
        &run_hhfilter($max_hhfilter_seqs,$min_hhfilter_seqs,"$tmpdir/cmsearch.$dd.unfilter.afa","$tmpdir/cmsearch.$dd.afa");
        &plain2gz("$tmpdir/cmsearch.$dd.afa", "$prefix.cmsearch.$dd.afa.gz");
    }
    my $Nf=&run_calNf("$tmpdir/cmsearch.$dd.afa");
    my $hitnum=`grep '^>' $tmpdir/cmsearch.$dd.afa|wc -l`+0;
    if ($Nf>=$target_Nf || $hitnum>=$max_hhfilter_seqs)
    {
        print "output cmsearch.$dd.afa (Nf>=$Nf) as final MSA\n";
        &System("cp $tmpdir/cmsearch.$dd.afa $prefix.afa");
        exit(0);
    }
}

#### adjust --incE ####
print "==== adjust cmsearch evalue ====\n";
&System("cat $tmpdir/db $tmpdir/db1 $tmpdir/db2 > $tmpdir/trimall.db");
&rmredundant_rawseq("$tmpdir/trimall.db", "$tmpdir/dball");
foreach my $incE((0.1,1,10))
{
    if (-s "$prefix.cmsearch.2.$incE.afa.gz" && `zcat $prefix.cmsearch.2.$incE.afa.gz|wc -l`+0>0)
    {
        &gz2plain("$prefix.cmsearch.2.$incE.afa.gz", "$tmpdir/cmsearch.2.$incE.afa");
    }
    else
    {
        &System("$bindir/qcmsearch --noali -A $tmpdir/cmsearch.2.$incE.a2m --cpu $cpu --incE $incE $tmpdir/infernal.cm $tmpdir/dball|grep 'no alignment saved'");
        &addQuery2a2m("$tmpdir/cmsearch.2.$incE.a2m","$tmpdir/cmsearch.2.$incE.unfilter.afa");
        &run_hhfilter($max_hhfilter_seqs,$min_hhfilter_seqs,"$tmpdir/cmsearch.2.$incE.unfilter.afa","$tmpdir/cmsearch.2.$incE.afa");
        &plain2gz("$tmpdir/cmsearch.2.$incE.afa", "$prefix.cmsearch.2.$incE.afa.gz");
        my $Nf=&run_calNf("$tmpdir/cmsearch.2.$incE.afa");
        my $hitnum=`grep '^>' $tmpdir/cmsearch.2.$incE.afa|wc -l`+0;
        if ($incE==10 || $Nf>=$target_Nf || $hitnum>=$max_hhfilter_seqs)
        {
            print "output cmsearch.2.$incE.afa (Nf>=$Nf) as final MSA\n";
            &System("cp $tmpdir/cmsearch.2.$incE.afa $prefix.afa");
            exit(0);
        }
    }
}


#### clean up ####
&System("sync");
&System("rm -rf $tmpdir") if (-d "$tmpdir");
exit(0);

#### add query to nhmmer/cmsearch a2m alignment ####
sub addQuery2a2m
{
    my ($infile,$outfile)=@_;
    my $hitnum=`grep '^>' $infile|wc -l`+0;
    if ($hitnum==0)
    {
        &System("$bindir/fastaOneLine $tmpdir/seq.fasta $outfile");
    }
    else
    {
        #&System("$bindir/a3m2msa $infile | grep -ohP '^\\S+' > $outfile.tmp");
        &System("$bindir/a3m2msa $infile | grep -ohP '^\\S+' | $bindir/fastaNA - > $outfile.tmp");
        my $template_seq=`head -2 $outfile.tmp|tail -1`;
        chomp($template_seq);
        my $Lt=length $template_seq;
        if ($Lch != $Lt)
        {
            ## merge msa by clustalo ##
            print "echo $Lch != $Lt. realign by clustalo\n";
            &System("$bindir/clustalo --p1=$tmpdir/seq.fasta --p2=$outfile.tmp --is-profile --threads=$cpu | $bindir/RemoveNonQueryPosition - | $bindir/fixAlnX - N $outfile");
        }
        else
        {
            &System("cat $tmpdir/seq.fasta $outfile.tmp | $bindir/fastaOneLine - | $bindir/fixAlnX - N $outfile");
        }
        &System("rm $outfile.tmp");
    }
    return $hitnum;
}


#### retrive sequence from listed in tabfile from db, ####
#### output the sequences to $tmpdir/$tag.db          ####
sub retrieveSeq
{
    my ($tabfile, $db, $tag)=@_;
    
    &System("cut -f1 $tabfile|sort|uniq > $tmpdir/$tag.list");
    &System("split -l $max_split_seqs $tmpdir/$tag.list $tmpdir/$tag.list.split.");
    foreach my $suffix (`ls $tmpdir/|grep $tag.list.split.|sed 's/$tag.list.split.//g'`)
    {
        chomp($suffix);
        &System("$bindir/blastdbcmd -db $db -entry_batch $tmpdir/$tag.list.split.$suffix -out $tmpdir/$tag.db.$suffix");
        &System("sort -k4g $tabfile | head -$max_aln_seqs | $bindir/trimBlastN $tmpdir/$tag.db.$suffix - $Lch $tmpdir/$tag.trim.split.$suffix");
        &System("rm $tmpdir/$tag.db.$suffix $tmpdir/$tag.list.split.$suffix");
    }
    &System("cat $tmpdir/$tag.trim.split.* | $bindir/fasta2pfam - | sort -u -k2 | $bindir/pfam2fasta - > $tmpdir/$tag.db");
    &System("rm  $tmpdir/$tag.trim.split.*");
    return;
}

#### remove identical sequences from unaligned database ####
sub rmredundant_rawseq
{
    my ($infile,$outfile)=@_;
    my $throw_away_sequences=int(0.4*$Lch);
    &System("$bindir/cd-hit-est-2d -T $cpu -i $tmpdir/seq.fasta -i2 $infile -c 1.0 -o $tmpdir/cdhitest2d.db -l $throw_away_sequences -M 5000");
    &System("$bindir/cd-hit-est -T $cpu -i $tmpdir/cdhitest2d.db -c 1.0 -o $outfile -l $throw_away_sequences -M 5000");
    &System("rm $tmpdir/cdhitest2d.db");
    return;
}

#### remove redundant sequences from MSA ####
sub run_hhfilter
{
    my ($max_hhfilter_seqs,$min_hhfilter_seqs,$infile,$outfile)=@_;

    my $cov=50;
    my $id =99;
    &System("$bindir/hhfilter -i $infile -id $id -cov $cov -o $outfile");
    my $hitnum=`grep '^>' $outfile|wc -l`+0;
    if ($hitnum>$max_hhfilter_seqs)
    {
        print "too many sequences: $hitnum > $max_hhfilter_seqs.\n";
        my $Nf=&run_calNf("$outfile");
        if ($Nf>$target_Nf)
        {
            foreach $cov ((60,70))
            {
                &System("$bindir/hhfilter -i $infile -id $id -cov $cov -o $outfile.tmp");
                $Nf=&run_calNf("$outfile.tmp");
                last if ($Nf<=$target_Nf);
                &System("mv $outfile.tmp $outfile");
                #$hitnum=`grep '>' $outfile|wc -l`+0;
                #last if ($hitnum<=$max_hhfilter_seqs);
            }
        }
        if ($hitnum>$max_hhfilter_seqs)
        {
            foreach $id ((96,93,90))
            {
                &System("$bindir/hhfilter -i $infile -id $id -cov $cov -o $outfile");
                $hitnum=`grep '>' $outfile|wc -l`+0;
                last if ($hitnum<=$max_hhfilter_seqs);
            }
        }
    }
    elsif ($hitnum<$min_hhfilter_seqs)
    {
        print "too few sequences: $hitnum < $min_hhfilter_seqs.\n";
        foreach $cov ((40))
        {
            &System("$bindir/hhfilter -i $infile -id $id -cov $cov -o $outfile");
            $hitnum=`grep '>' $outfile|wc -l`+0;
            last if ($hitnum>=$max_hhfilter_seqs);
        }
    }
    return $outfile;
}

#### calculate Nf for MSA ####
sub run_calNf
{
    my ($infile)=@_;
    &System("$bindir/hhfilter -i $infile -id 99 -cov $target_Nf_cov -o $infile.$target_Nf_cov");
    my $target_Nf_tmp=$target_Nf+1;
    return `$bindir/fastNf $infile.$target_Nf_cov 0.8 0 $target_Nf_tmp`+0;
}

#### copy gzip compressed file to plain file ####
sub gz2plain
{
    my ($infile,$outfile)=@_;
    return if (-s "$outfile");
    if ($infile=~/.gz$/ && $outfile!~/.gz$/)
    {
        &System("cp $infile $outfile.gz");
        &System("zcat $outfile.gz > $outfile");
    }
    else
    {
        &System("cp $infile $outfile");
    }
    return;
}

#### copy plain file to gzip compressed file ####
sub plain2gz
{
    my ($infile,$outfile)=@_;
    return if (-s "$outfile");
    if ($infile!~/.gz$/ && $outfile=~/.gz$/)
    {
        &System("cat $infile | gzip - > $infile.gz");
        &System("cp $infile.gz $outfile");
    }
    else
    {
        &System("cp $infile $outfile");
    }
    return;
}

#### run system command ####
sub System
{
    my ($cmd)=@_;
    print "$cmd\n";
    return system("$cmd");
}
