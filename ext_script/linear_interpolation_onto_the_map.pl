use strict ; 
use warnings ; 

my %map ;
my %cm ; 
open IN, "</shared/projects/gametic_segregation_distortion/data/A.lyrata_map_clean_scaffold.txt" ;
<IN> ; ### trash the header line, format: Chr, bp, cm
while (<IN>) { 
	chomp ; 
	my @split = split ( /\t/, $_ ) ; 
	push @{ $map{$split[0]} }, $split[1] ;
	push @{ $cm{$split[0]} }, $split[2] ; 
}
close IN ; 

### this file is assumed to ahve format chr, bp, ... 
open IN, "<$ARGV[0]" ;
while (<IN>) { 
	chomp ; 
	my @split = split ( /\t/, $_ ) ; 
	if ( !exists( $map{$split[0]} ) ) { 
		print STDERR "Chromosome $split[0] not in map file... please check.\n" ;
		next ;
	}

	## skip those sites before the first marker and those after the last marker
	if ( $split[1] < ${$map{$split[0]}}[0] || $split[1] > ${$map{$split[0]}}[$#{$map{$split[0]}}] ) { 
		next ; 
	}	

	### otherwise figure out which markers we are stuck between
	foreach my $index ( 0..$#{$map{$split[0]}} - 1 ) {
		if ( $split[1] >= ${$map{$split[0]}}[$index] && $split[1] < ${$map{$split[0]}}[$index+1] ) { 
			my $scalar = ( $split[1] - ${$map{$split[0]}}[$index] ) / ( ${$map{$split[0]}}[$index+1] - ${$map{$split[0]}}[$index] ) ; 
			my $cm_pos = $scalar * ( ${ $cm{$split[0]} }[$index+1] - ${ $cm{$split[0]} }[$index] ) + ${ $cm{$split[0]} }[$index] ; 
			print $split[0], "\t", $split[1], "\t", $cm_pos ; 
			foreach ( 2..$#split ) { 
				print "\t", $split[$_] ; 
			}
			print "\n" ; 
		}
	}
}
