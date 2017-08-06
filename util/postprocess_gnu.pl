#!/usr/bin/perl
#
#

$unfold_inp_file=$ARGV[0];
$unfold_log_file=$ARGV[1];
print "The unfold_inp_file is $unfold_inp_file\n";
print "The unfold_log_file is $unfold_log_file\n";

print "processing...\n";
print "\n";
print "OUTPUT:\n";
print "1. bands3d.dat	# which contains kpts,en,weitht every line.\n";

open(FILE_in, $unfold_inp_file);

$_=<FILE_in>;
$seedname=$_;
$_=<FILE_in>;
@tmp=split;
$e_min = $tmp[0];
$e_max = $tmp[1];
$de    = $tmp[2];
$nepts_total = ($e_max-$e_min)/$de + 1;

for($i=0;$i<2;$i++) {	#skip  lines
  $_=<FILE_in>;
}

$_=<FILE_in>;
@tmp=split;
$nkpts_path_border = $tmp[0];		#number of k-points
$nk_per_path       = $tmp[1];		#nubmer of kpts per k-path
$nkpts_total = ($nkpts_path_border-1)*$nk_per_path+1;

#warn "$nkpts_path_border $nk_per_path\n";

$_=<FILE_in>;
@kpt_path_border_old = split;
$k_path_length[0]=0;
for($i=0;$i<$nkpts_path_border-1;$i++){
	$_=<FILE_in>;
	@kpt_path_border_new = split;
	#increment of every two kpts in k path.
	for($xi=0;$xi<3;$xi++){
		$increment[$xi] = ( $kpt_path_border_new[$xi] - $kpt_path_border_old[$xi] ) / $nk_per_path;
	}
	@kpt_old=@kpt_path_border_old;
	for($j=0;$j<$nk_per_path;$j++){
		#calculate each kpoints in this current path, and the k_path_length from the original kpt.
		for($xi=0;$xi<3;$xi++){
			$kpt_new[$xi] = $kpt_old[$xi] + $increment[$xi];
		}
		$k_path_length[$i*$nk_per_path+$j+1] = $k_path_length[$i*$nk_per_path+$j] + sqrt( ($kpt_new[0]-$kpt_old[0])**2 + ($kpt_new[1]-$kpt_old[1])**2 + ($kpt_new[2]-$kpt_old[2])**2 );
		@kpt_old=@kpt_new;
	}
	@kpt_path_border_old = @kpt_path_border_new;
}
close(FILE_in);


$en[0]=$e_min;
for($j=0;$j<$nepts_total-1;$j++){
	$en[$j+1]=$en[$j]+$de;
}




open(FILE_in, $unfold_log_file);

while($_=<FILE_in>){
	if(/=========================================================/){
		$_=<FILE_in>;
		@tmp=split;
		#$nkpts_total=$tmp[0];
		#$nepts_total=$tmp[1]; # number of energy points.
		if($nkpts_total ne $tmp[0] or $nepts_total ne $tmp[1]){
			warn "Input files do not match each other.";
			print "$tmp[0] $tmp[1]\n";
			print "$nkpts_total $nepts_total\n";
			exit 1;
		}
		$nlines_per_kpt = int($nepts_total / 20);
		if($nepts_total%20 ne 0){
			$nlines_per_kpt += 1;
		}
		for($i=0;$i<$nkpts_total;$i++){
			$_=<FILE_in>;     #skip the line indicate the kpt.
			for($j=0;$j<$nlines_per_kpt;$j++){
				$_=<FILE_in>;
				@tmp=split;
				for($k=0;$k<@tmp;$k++){
					$weight[$i][$j*20+$k]=$tmp[$k];
				}
			}
		}
	}
}
close(FILE_in);

open(FILE_out, ">bands3d.dat");
for($i=0;$i<$nkpts_total;$i++){
	for($j=0;$j<$nepts_total;$j++){
		printf(FILE_out "%12.9f %7.3f %12.9f\n",$k_path_length[$i],$en[$j],$weight[$i][$j]);
	}
}
close(FILE_out);
