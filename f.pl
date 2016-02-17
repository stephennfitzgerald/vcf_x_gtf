use strict;
use warnings;
use Data::Dumper;


my $usage = qq{USAGE:
 perl f.pl in_file
};

open my $fh, '<', shift || help();

die $! unless( fileno($fh) );

my (%H,%K);

while(<$fh>) {
 my @cols=split"\t",$_;
 if( $cols[9] =~/^OTTHUMG/ ) { # it's current gencode
  $H{ $cols[0] }{ $cols[1] }{ 'gen_new' }{ $cols[12] }++;
  push @{ $K{ $cols[0] }{ $cols[1] }{ 'gen_new' } }, $_;
 }
 elsif( $cols[9] =~/^NM_/ ) { # it's a ref_seq
  $H{ $cols[0] }{ $cols[1] }{ 'ref_seq' }{ $cols[12] }++;
  push @{ $K{ $cols[0] }{ $cols[1] }{ 'ref_seq' } }, $_;
 }
 elsif( $cols[9] =~/^ENSG/ ) { # it's an old gencode
  $H{ $cols[0] }{ $cols[1] }{ 'gen_old' }{ $cols[12] }++;
  push @{ $K{ $cols[0] }{ $cols[1] }{ 'gen_old' } }, $_;
 }
}


foreach my $chr(sort keys %H) {
 foreach my $pos(keys %{ $H{$chr} }) {
  if(! exists($H{$chr}{$pos}{'gen_new'}{'exon'}) && ! exists($H{$chr}{$pos}{'gen_new'}{'splice'})) {
   delete($H{$chr}{$pos}); # if the new gencode annotation is neither an exon or a splice
   next;
  }
  if(1 < keys %{ $H{$chr}{$pos}{'ref_seq'} }) {
   delete($H{$chr}{$pos}{'ref_seq'}); # if there is > 1 ref_seq annotation
  }
  elsif(( keys %{ $H{$chr}{$pos}{'ref_seq'} })[0] ne 'intron') {
   delete($H{$chr}{$pos}{'ref_seq'}); # and if the single ref_seq annotation is not 'intron'
  }
  if(1 < keys %{ $H{$chr}{$pos}{'gen_old'} }) {
   delete($H{$chr}{$pos}{'gen_old'});
  }
  elsif(( keys %{ $H{$chr}{$pos}{'gen_old'} })[0] ne 'intron') {
   delete($H{$chr}{$pos}{'gen_old'});
  }
 }
}

foreach my $chr(sort keys %H) {
 foreach my $pos(keys %{ $H{$chr} }) {
  if(1 < keys %{ $H{$chr}{$pos} }) { # there has to be more that just the new gencode annotation 
   foreach my $type(keys %{ $K{$chr}{$pos} }) {
    foreach my $record(@{ $K{$chr}{$pos}{$type} }) {
     print $record;
    }
   }
  }
 }
}

sub help {
  print $usage;
  exit(1);
}
