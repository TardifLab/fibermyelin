#!/usr/bin/perl -w

### Preprocess ir-diff data
#
##########################################################
##########################################################
## Preprocess data (i.e. conversions to various formats)
##
###########################################
#
###########################################
# Created by Ilana Leppert Feb 2019
###########################################
require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
use List::MoreUtils qw(zip);

if($0 =~ /[\/A-Za-z_0-9-]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Prepares the data for the ir-diff pipeline


Usage: $program -d <dcmdir>
-help for options

USAGE

my @args_table = (["-d","string",1,\$dcmdir,"name of dicom directory"]);
#["-clobber","boolean",undef,\$clobber,"clobber all output files (currently not implemented, sorry)"],#HERE not implemented);


Getopt::Tabular::SetHelp ($Usage, '');
my @args;
GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
die $Usage unless $#ARGV >=0;

print "---Get the list of sequences from the minc files---\n";
`lsminc *.mnc.gz > mnclist` unless -e "mnclist";

print "---Sort the dicom by sequence name and then sequence number-----\n";
print "sort_dicom.pl $dcmdir/MR* '0018,1030' $dcmdir//\n";
`sort_dicom.pl $dcmdir/MR* '0018,1030' $dcmdir/`;
`sort_dicom.pl $dcmdir/MR* '0020,0011' $dcmdir/`;


print "\n-----Convert the ir-diffusion and hardi to nii----------\n";
`mkdir nii`;
$done=`\\ls nii/*midiff*fov1*.nii*`; chomp($done); 
if ($done eq "") {
	print "dcm2nii -o nii $dcmdir/*fov1/*\n";
    `dcm2nii -o nii $dcmdir/*fov1/*` ;
}
$done=`\\ls nii/*midiff*fov2*.nii*`; chomp($done); 
if ($done eq "") {
    print "dcm2nii -o nii $dcmdir/*fov2/*\n";
    `dcm2nii -o nii $dcmdir/*fov2/*`;
}
$done=`\\ls nii/*cmrrmbep2d*.nii*`; chomp($done); 
if ($done eq "") {
    print "dcm2nii -o nii $dcmdir/cmrr_mbep2d_diff_multib*/*\n";
    `dcm2nii -o nii $dcmdir/cmrr_mbep2d_diff_multib*/*`;
}

#$hardi_nii = `\\ls nii/*cmrrmbep2ddiff*.nii*`; chomp($hardi_nii);
$hardi_bvec = `\\ls nii/*cmrrmbep2ddiff*.bvec*`; chomp($hardi_bvec);
$hardi_bval = `\\ls nii/*cmrrmbep2ddiff*.bval*`; chomp($hardi_bval);

print "\n----Grab a dicom file to be able to read the timing from it------\n";
$dicom = `\\ls $dcmdir/*_fov1/MR* | grep MR -m 1`;

print "\n-----Convert the HARDI diffusion to mif----------\n";
$hardid = "hardi-analysis/";
`mkdir $hardid` unless -e $hardid;
$done=`\\ls hardi-analysis/dwi.mif`; chomp($done); 
unless (-e $done) {
    print "mrconvert $dcmdir/cmrr_mbep2d_diff_multib*/ hardi-analysis/dwi.mif\n";
    `mrconvert $dcmdir/cmrr_mbep2d_diff_multib*/ hardi-analysis/dwi.mif` ;
}

chomp($dicom);
print "------\n";
#print "irdiff_pipeline.pl -irdiff nii/*{fov1,fov2}*.nii.gz -bvecs nii/*{fov1,fov2}*.bvec -bvals nii/*{fov1,fov2}*.bval -hardi hardi-analysis/dwi.mif -dcm $dicom > log\n\n";
print "irdiff_pipeline.pl -irdiff nii/*{fov1,fov2}*.nii.gz -hardi hardi-analysis/dwi.mif -bvecs $hardi_bvec -bvals $hardi_bval -dcm $dicom > log\n\n";

print "------\n";