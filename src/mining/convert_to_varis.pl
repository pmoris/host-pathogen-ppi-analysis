use strict;
use warnings;

open(IN,'ppi_go2name_clustering_filtered.txt');
open(OUT,'>varis_ppi_go2name_clustering_filtered.txt');

print OUT "lhs,rhs,confidence\n";

while(<IN>){
	#('b_A_1,b_F_6>a_A_-5,a_N_-1,a_V_-4', 221, (0.7189189189189189, 266))
	if(m/\('(.+)>(.+)',\s-?\d+,\s\((\d\.\d+),/){
		my$lhs = $1;
		my$rhs = $2;
		my$conf = $3;
		$lhs =~ s/,/ /gi;
		$rhs =~ s/,/ /gi;
		print OUT $lhs.','.$rhs.','.$conf."\n";
	}
} 
