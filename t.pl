use strict;
use Vcf;
use Getopt::Long;
use Data::Dumper;

open my $OD, "data/EIEEgenes_gencode19.gtf" or die $!;

open my $NW, "data/EIEEgenes_loutre.GRCh38-GRCh37-lift.mapped" or die $!;

open my $RF, "data/refseq_70_37backmap.gtf" or die $!;

my(%OD,%NW,%RF,%VC);
my $intron_offset = 8;
my $vcf_file;
my $help;

my $usage = qq{USAGE:
 perl t.pl --help
 perl t.pl --vcf_file path/to/vcf_file (required) --intron_offset an_integer (optional)
 eg. perl t.pl --vcf_file test.vcf.gz
 
 required perl modules : Vcf.pm
};

GetOptions(
    "help" => \$help,
    "vcf_file=s" => \$vcf_file,
    "intron_offset:i" => \$intron_offset,
  );

if (!$vcf_file) {
  print $usage;
  exit(1);
}

if($help) {
  print $usage;
  exit(0);
}

foreach my $fh([$OD,\%OD],[$NW,\%NW],[$RF,\%RF]) {
 fill($fh);

 foreach my $trans( keys %{ $fh->[1] } ) {
  my $region = $fh->[1]->{$trans}->{'trans_region'};
  my $biotype = $fh->[1]->{$trans}->{'trans_type'};
  my $vcf = Vcf->new(file => "$vcf_file", region=>"$region");
  $vcf->parse_header();
  while (my $x=$vcf->next_data_array()) {
   my $snp = "$$x[0]:$$x[1]:$$x[2]:$$x[3]:$$x[4]:$$x[5]::$$x[6]:$$x[7]:$$x[8]";
   my $exon_ct = scalar keys %{ $fh->[1]->{$trans}->{'exon'} };
   my $found = 0;
   my $var_size = length($$x[3]) - 1; # if the alt allele is a deletion
   foreach my $exon(keys %{ $fh->[1]->{$trans}->{'exon'} }) {
    last if $found;
    foreach my $exon_no(keys %{ $fh->[1]->{$trans}->{'exon'}->{$exon} }) {
     my $exon_start = $fh->[1]->{$trans}->{'exon'}->{$exon}->{$exon_no}->[1];
     my $exon_end = $fh->[1]->{$trans}->{'exon'}->{$exon}->{$exon_no}->[2];
     $VC{ $snp }{ $fh->[1]->{$trans}->{'gene_id'} }{ $trans }{'biotype'} = $biotype;
     if( ($$x[1] <= $exon_end) && ($$x[1] + $var_size >= $exon_start) ) { # var is in the exon
      $VC{ $snp }{ $fh->[1]->{$trans}->{'gene_id'} }{ $trans }{'exon'}{ $exon }++;
      $found = 1;
      last;
     }
     elsif($exon_no < $exon_ct && $exon_no > 1) { # var not in/near first or last exon   
      if( 
         ($$x[1] <= $exon_end + $intron_offset) && ($$x[1] > $exon_end)
               || 
         ($$x[1] + $var_size >= $exon_start - $intron_offset) && ($$x[1] + $var_size < $exon_start)
        ){
       $VC{ $snp }{ $fh->[1]->{$trans}->{'gene_id'} }{ $trans }{'splice'}{ $exon }++;
       $found = 1;
       last;
      }
     }
     # var is in/near first exon - only check downstream of exon end
     elsif($exon_no == 1 && ($$x[1] <= $exon_end + $intron_offset) && ($$x[1] > $exon_end)) { 
      $VC{ $snp }{ $fh->[1]->{$trans}->{'gene_id'} }{ $trans }{'splice'}{ $exon }++;
      $found = 1;
      last;
     }
     # var is in/near last exon - only check upstream of exon start
     elsif($exon_no == $exon_ct && ($$x[1] + $var_size >= $exon_start - $intron_offset) && ($$x[1] + $var_size < $exon_start)) {
      $VC{ $snp }{ $fh->[1]->{$trans}->{'gene_id'} }{ $trans }{'splice'}{ $exon }++;
      $found = 1;
      last;
     }
    }
   }
   if(!$found) { # it must be intronic
    $VC{ $snp }{ $fh->[1]->{$trans}->{'gene_id'} }{ $trans }{'intron'}++;
   }
   foreach my $cds(keys %{ $fh->[1]->{$trans}->{'CDS'} }) {
    foreach my $cds_no(keys %{ $fh->[1]->{$trans}->{'CDS'}->{$cds} }) {
     my($cds_fr,$cds_to) = @{ $fh->[1]->{$trans}->{'CDS'}->{$cds}->{$cds_no} }[1,2];
     if($$x[1] <= $cds_to && $$x[1] + $var_size >= $cds_fr) {
      $VC{ $snp }{ $fh->[1]->{$trans}->{'gene_id'} }{ $trans }{'CDS'}++;
     }
    }
   }
  }
 }
}

foreach my $snp(keys %VC) {
 my($sr,$pos,$snp_id,$ref,$alt,$qual,$filter,$info,$format) = split':', $snp;
 foreach my $gene(keys %{ $VC{ $snp } }) {
  foreach my $trans(keys %{ $VC{ $snp }{ $gene } }) {
   my $trans_biotype = $VC{ $snp }{ $gene }{ $trans }{ 'biotype' };
   my $cds = exists( $VC{ $snp }{ $gene }{ $trans }{ 'CDS' } ) ? 1 : 0;
   $trans_biotype = defined($trans_biotype) ? $trans_biotype : q{.};
   foreach my $int_ex(keys %{ $VC{ $snp }{ $gene }{ $trans } }) {
    next if ($int_ex eq 'biotype' || $int_ex eq 'CDS');
    if($int_ex eq 'intron') {
     print join("\t", $sr,$pos,$snp_id,$ref,$alt,$qual,$filter,$info,$format,$gene,$trans,$trans_biotype,$int_ex,q{.},q{.}), "\n";
    }
    else {
     foreach my $exon(keys %{ $VC{ $snp }{ $gene }{ $trans }{ $int_ex } }) {
      print join("\t", $sr,$pos,$snp_id,$ref,$alt,$qual,$filter,$info,$format,$gene,$trans,$trans_biotype,$int_ex,$exon,$cds), "\n";
     }
    }
   }
  }
 }
}

sub fill {
 my $fh = shift;
 my $H = $fh->[1];
 my $F = $fh->[0];
 while(<$F>) {
  chomp;
  my($sr,$gp,$feat,$fr,$to,$sc,$str,$frame,$attr)=split"\t",$_;
  $sr=~s/^chr//;
  my($gene_id,$trans_id,@rest) = split q{; }, $attr; 
  ($gene_id)=$gene_id=~/gene_id \"*([\w\.]+)/;
  ($trans_id)=$trans_id=~/transcript_id \"*([\w\.]+)/;
  next if($gene_id eq $trans_id && $gp ne 'RefSeq_transMap'); 
  if("$feat" eq 'transcript') {
   my $mfr = $fr - $intron_offset;
   my $mto = $to + $intron_offset;
   $H->{ $trans_id }{ 'trans_region' } = $sr . q{:} . $mfr . q{-} . $mto;
   $H->{ $trans_id }{ 'gene_id' } = $gene_id; 
  }
  my ($exon_numb,$gene_name,$trans_name,$trans_type);
  foreach my $ele(@rest) {
   if($ele=~/exon_number/) {
    ($exon_numb)=$ele=~/exon_number \"*(\d+)/;
   }
   elsif($ele=~/gene_name/) {
    ($gene_name)=$ele=~/gene_name \"*([\w\.]+)/;
    $H->{ $trans_id }{ 'gene_name' } = $gene_name;
   }
   elsif($ele=~/transcript_name/) {
    ($trans_name)=$ele=~/transcript_name \"*([\w\.]+)/;
    $H->{ $trans_id }{ 'trans_name' } = $trans_name;
   }
   elsif($ele=~/transcript_type/) {
    ($trans_type)=$ele=~/transcript_type \"*([\w\.]+)/;
    $H->{ $trans_id }{ 'trans_type' } = $trans_type;
   }
   elsif( my($ex_id)=$ele=~/exon_id \"*([\w\.]+)/ ) {
    push @{ $H->{ $trans_id }{ $feat }{ $ex_id }{ $exon_numb } }, $sr,$fr,$to,$str,$frame;
   }
  }
 }
}
