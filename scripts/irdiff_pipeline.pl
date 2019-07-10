#!/usr/bin/perl -w

### Process ir-diff data
#
##########################################################
##########################################################
##
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

Runs the ir-diff pipeline.


Usage: $program
-help for options

USAGE
#if($ARGV[0]=~/help/i){print $Usage; exit(1);}

#presets: DO: must create these from the inputs. these are currently hardcoded; links are for last thing processed
$acq="/data_/tardiflab/ilana/midiff_mgh/acqparams.txt"; #didn't AP-PA so no distortion correction for now...

# TODO: this assume IR-idff always acquired the same way, (2 b=0, 30 Siemens directions, 4 averages)
# the bvec here has that 1st b=0 removed, have to remove it from the volume in the code, this could be done in the code
$bvecs_ir="/data_/tardiflab/ilana/midiff_mgh/30bvecs-4avg-rmb0.txt";# 1 b=0, 30 directions, from Siemens built-in
$bvals_ir="/data_/tardiflab/ilana/midiff_mgh/30bvals-4avg-rmb0.txt";# 1 b=0, 30 directions, b=1000
#$bvecs_ir="/data_/tardiflab/ilana/midiff_mgh/30bvecs-8avg-rmb0.txt";# 1 b=0, 30 directions, from Siemens built-in
#$bvals_ir="/data_/tardiflab/ilana/midiff_mgh/30bvals-8avg-rmb0.txt";# 1 b=0, 30 directions, b=1000

# TODO: right now assume the HARDI is always done the same way, the 108 direction protocol from UNF
#$scheme="/data_/tardiflab/ilana/midiff_mgh/NODDI.scheme";
# Create custom one based on bvals and results from eddy

# index file that tells eddy which line/of the lines in the acqparams.txt file that are relevant for the data passed into eddy.
$index="/data_/tardiflab/ilana/midiff_mgh/index.txt"; # 

#DO : niak (matlab interface for nii) needs unzippped nii files, so have to unzip before any calls, should find a way around this...

$afd_thresh = 0.1;

my @args_table = (#["-clobber","boolean",undef,\$clobber,"clobber all output files (currently not implemented, sorry)"],#HERE not implemented
		  ["-irdiff","string",2,\@irdiff,"ir-diffusion data, 2 files:fov1 fov2"],
    #[
#["-anat","string",2,\@anatfiles,"anatomical (T1w pre-contrast) nii file and minc file."],
    ["-hardi","string",1,\$hardi,"HARDI diffusion images in .mif format (gradient info included)"],
    ["-bvecs","string",1,\$bvecs,"bvec files for hardi"],
     ["-bvals","string",1,\$bvals,"bval files for hardi"],
    ["-dcm","string",1,\$dicom,"name of a single dicom file from the ir-diff epxeriment (to read timing info)"],
   #-acqparams","string",1,\$acq,"acq params for topup"],
);



Getopt::Tabular::SetHelp ($Usage, '');
my @args;
GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
die $Usage unless $#ARGV >=0;

#my $subdir=$args[0]; #the input dir
chomp($unique=`date +%y%m%d-%H:%M`); 
print "--> Date is $unique\n";
print "**** $program @ARGV ****\n";

##IR-dff
print "==============IR-diff processing====================\n\n";
`mkdir ir-analysis` unless -e 'ir-analysis';

print "----- unshuffle the data (this might take a minute or 2...)-----\n";
#unzip the .nii!
# unshuffle_diff(img,dcm)
$fov1 = $irdiff[0];
$fov1 =~ s/\.gz$//g; #remove .gz at end 
`gunzip $irdiff[0]` unless -e $fov1;

$ufov1 = "ir-analysis/unshuffle_fov1.nii";
$ufov2 = "ir-analysis/unshuffle_fov2.nii";
print "\nunshuffle_diff('$fov1','$dicom','$ufov1')\n";
run_matlab("unshuffle_diff('$fov1','$dicom','$ufov1')") unless -e $ufov1;

$fov2 = $irdiff[1];
$fov2 =~ s/\.gz$//g; #remove .gz at end 
`gunzip $irdiff[1]` unless -e $fov2;
print "\nunshuffle_diff('$fov2','$dicom','$ufov2')\n";
run_matlab("unshuffle_diff('$fov2','$dicom','$ufov2')") unless -e $ufov2;

# this create the TI.txt file, which should be the same for both runs as well as an 

print "\n---Now remove 1st b=0 (which is just there for steady-state)----\n";
# we need to know some params
$numavg = `more info.txt | grep numavg | grep -Eo '[0-9]+'`;chomp($numavg);
$numdiff = `more info.txt | grep directions | grep -Eo '[0-9]+'`;chomp($numdiff);
$totframes = $numavg*$numdiff;
print "fslroi $ufov1 ir-analysis/unshuffle_fov1-remb0.nii $numavg $totframes-$numavg\n";
`fslroi $ufov1 ir-analysis/unshuffle_fov1-remb0.nii $numavg $totframes-$numavg` unless -e "ir-analysis/unshuffle_fov1-remb0.nii";
print "fslroi $ufov2 ir-analysis/unshuffle_fov2-remb0.nii $numavg $totframes-$numavg\n";
`fslroi $ufov2 ir-analysis/unshuffle_fov2-remb0.nii $numavg $totframes-$numavg` unless -e "ir-analysis/unshuffle_fov2-remb0.nii";

print "-------Now combine the 2 datasets-------\n";
# we need to muck around or else it does not combine properly
`gunzip ir-analysis/*.gz` unless -e "ir-analysis/unshuffle_fov1-remb0.nii";
print "\ncomb_fov('ir-analysis/unshuffle_fov1-remb0.nii','ir-analysis/unshuffle_fov2-remb0.nii','TI.txt','ir-analysis/ir-diff.nii')\n";
run_matlab("comb_fov('ir-analysis/unshuffle_fov1-remb0.nii','ir-analysis/unshuffle_fov2-remb0.nii','TI.txt','ir-analysis/ir-diff.nii')") unless -e "ir-analysis/ir-diff.nii";

print "---Make sure the spacing is correct---\n";
print "mrconvert tmp-ir-diff.nii  -vox ,,2 ir-analysis/ir-diff.nii\n";
`mrconvert tmp-ir-diff.nii  -vox ,,2 ir-analysis/ir-diff.nii` unless -e "ir-analysis/ir-diff.nii";
#--> fit b=0
#fslroi ir-diff-reformat.nii b0-ir-diff.nii 0 8


## TODO, right now i'm just using the hard-coded ones but that won't work if we don't have the siemens 30directions
#print "-----Remove the extra b=0s in the bvec bval of ir-diff-----\n";
## open vector and bvalue files
#open(BVEC,"< $bvec") or die "can't open $bvec: $!";
#@bvecs = <BVEC>; #slurp all files, each line at an index
#open(BVAL,"< $bval") or die "can't open $bval: $!";
#@bval = <BVAL>; #slurp all files, each line at an index
#we want to get rid of the 1st b=0 at the beginning of each average


print "\n=============HARDI processing===========\n";
chdir("hardi-analysis");

 # extract the b=0s
print "dwiextract -bzero ../$hardi b0s_AP.mif\n";
`dwiextract -bzero ../$hardi b0s_AP.mif`;
## Preprocessing
### denoising
print "dwidenoise ../$hardi dwi_dn.mif -noise noise.mif\n";
`dwidenoise ../$hardi dwi_dn.mif -noise noise.mif` unless -e "dwi_dn.mif";
print "mrcalc ../$hardi dwi_dn.mif -subtract res.mif\n";
`mrcalc ../$hardi dwi_dn.mif -subtract res.mif`;
### Mask the data
print "dwi2mask dwi_dn.mif mask.mif\nmrview  dwi_dn.mif -roi.load mask.mif\n";
`dwi2mask dwi_dn.mif mask.mif` unless -e "mask.mif";
`mrview  dwi_dn.mif -roi.load mask.mif`;
### MEan B0 for visualization
print "dwiextract dwi_dn.mif - -shell 15 | mrmath - mean meanb0.mif -axis 3\n"; 
`dwiextract dwi_dn.mif - -shell 15 | mrmath - mean meanb0.mif -axis 3`;


###### NODDI processing
print "\n-----NODDI processing-------\n";
`mkdir noddi` unless -e 'noddi';
print "mrconvert mask.mif noddi/brain_mask.nii\n";
`mrconvert mask.mif noddi/brain_mask.nii`;
print "mrconvert  dwi_dn.mif noddi/dwi_dn.nii\n";
`mrconvert  dwi_dn.mif noddi/dwi_dn.nii`;

$eddy =0;
if ($eddy) {
  print "eddy_openmp --data_is_shelled --imain=noddi/dwi_dn.nii --mask=noddi/brain_mask.nii --acqp=$acq --index=$index --bvecs=../$bvecs --bvals=../$bvals --out=noddi/diff_eddy-corr\n";
  `eddy_openmp --data_is_shelled --imain=noddi/dwi_dn.nii --mask=noddi/brain_mask.nii --acqp=$acq --index=$index --bvecs=../$bvecs --bvals=../$bvals --out=noddi/diff_eddy-corr` unless -e "noddi/diff_eddy-corr.eddy_parameters";


  #---ugggh the mrconvert flipped it upside down... grrrr
  print "fslswapdim noddi/brain_mask.nii  x -y z noddi/brain_mask-swap.nii\n";
  `fslswapdim noddi/brain_mask.nii  x -y z noddi/brain_mask-swap.nii`;
  `gunzip noddi/brain_mask-swap.nii.gz` unless -e "noddi/brain_mask-swap.nii";
  print "fslswapdim noddi/diff_eddy-corr.nii  x -y z noddi/diff_eddy-corr-swap.nii\n";
  `fslswapdim noddi/diff_eddy-corr.nii  x -y z noddi/diff_eddy-corr-swap.nii` unless -e "noddi/diff_eddy-corr-swap.nii.gz";
  `gunzip noddi/diff_eddy-corr-swap.nii.gz` unless -e "noddi/diff_eddy-corr-swap.nii";

  $eddy_nogz = "noddi/diff_eddy-corr-swap.nii";

  # redo the scheme file based on modified bvecs
  print "bvecs_bvals2camino.pl -vec noddi/diff_eddy-corr.eddy_rotated_bvecs -val  ../$bvals -o NODDI.scheme\n";
  `bvecs_bvals2camino.pl -vec noddi/diff_eddy-corr.eddy_rotated_bvecs -val  ../$bvals -o NODDI.scheme`;

  print "AMICO_process('./','','$eddy_nogz','noddi/brain_mask-swap.nii','NODDI.scheme')\n";
  run_matlab("AMICO_process('./','','$eddy_nogz','noddi/brain_mask-swap.nii','NODDI.scheme')") unless -e "AMICO/NODDI/FIT_ICVF.nii";

  ###### mrtrix processing
  print "\n-----Fiber orientation/AFD processing (mrtrix)-------\n";
  # use the denoised, eddy-corrected output as input to this step
  print "mrconvert -fslgrad noddi/diff_eddy-corr.eddy_rotated_bvecs ../$bvals  $eddy_nogz  dwi_dn_ed.mif\n";
  `mrconvert -fslgrad noddi/diff_eddy-corr.eddy_rotated_bvecs ../$bvals $eddy_nogz dwi_dn_ed.mif`;

  $dwi = "dwi_dn_ed.mif";

}else{
  #---ugggh the mrconvert flipped it upside down... grrrr
  print "fslswapdim noddi/brain_mask.nii  x -y z noddi/brain_mask-swap.nii\n";
  `fslswapdim noddi/brain_mask.nii  x -y z noddi/brain_mask-swap.nii`;
  `gunzip noddi/brain_mask-swap.nii.gz` unless -e "noddi/brain_mask-swap.nii";
  print "fslswapdim noddi/dwi_dn.nii  x -y z noddi/dwi_dn-swap.nii\n";
  `fslswapdim noddi/dwi_dn.nii  x -y z noddi/dwi_dn-swap.nii` unless -e "noddi/dwi_dn-swap.nii.gz";
  `gunzip noddi/dwi_dn-swap.nii.gz` unless -e "noddi/dwi_dn-swap.nii";
  
  print "bvecs_bvals2camino.pl -vec ../$bvecs -val  ../$bvals -o NODDI.scheme\n";
  `bvecs_bvals2camino.pl -vec ../$bvecs -val  ../$bvals -o NODDI.scheme`;

  print "AMICO_process('./','','noddi/dwi_dn-swap.nii','noddi/brain_mask-swap.nii','NODDI.scheme')\n";
  run_matlab("AMICO_process('./','','noddi/dwi_dn-swap.nii','noddi/brain_mask-swap.nii','NODDI.scheme')") unless -e "AMICO/NODDI/FIT_ICVF.nii";

  ###### mrtrix processing
  print "\n-----Fiber orientation/AFD processing (mrtrix)-------\n";
  $dwi = "dwi_dn.mif";
}


#if this file exists don't run through everyting
$done=0;
if (-e "fixel_dir/afd.mif"){print "HARDI processing seems to be done...\n"; $done=1;}
if ($done==0){

### Fiber response estimation  
print "dwi2response tournier $dwi resp.txt\nshview resp.txt\n";
`dwi2response tournier $dwi resp.txt` unless -e "resp.txt";
`shview resp.txt`;

### Constrained spherical deconvolution  
#$ dwi2fod csd Input DWI Input response text file Output FOD image -mask Input DWI mask  
#$ mrview Input DWI -odf.load_sh Output FOD image
#*will only use outer shell*  
print "dwi2fod csd $dwi resp.txt fod.mif -mask mask.mif\nmrview meanb0.mif -odf.load_sh fod.mif\n";
`dwi2fod csd $dwi resp.txt fod.mif -mask mask.mif` unless -e "fod.mif";
`mrview meanb0.mif -odf.load_sh fod.mif`;

### Segment FOD images to estimate fixels and their apparent fibre density
print "fod2fixel -mask mask.mif fod.mif -afd afd.mif fixel_dir\n";
`fod2fixel -mask mask.mif fod.mif -afd afd.mif fixel_dir` unless -e "fixel_dir/afd.mif";

#Integral of the apparent fiber density (AFD) over a FOD 'lobe', i.e. a fiber density per population in the voxel  

#Colored by Apparent Fiber Density  
`mrview meanb0.mif -fixel.load fixel_dir/afd.mif`  ;
#Colored by direction  
`mrview meanb0.mif -fixel.load fixel_dir/directions.mif`;

}

## Now prepare everything for the fibermyelin processing
print "\n=====Pre-process for fiber-myelin=====\n";
## this comes from Jen's fibermyelin.py instructions

#convert the AFD and directions files to non-sparse format:
#fixel2voxel <afd.mif> split_data <afd_voxel.mif> ...
print "fixel2voxel fixel_dir/afd.mif split_data afd_voxel.mif\n";
`fixel2voxel fixel_dir/afd.mif split_data afd_voxel.mif` unless -e "split_data afd_voxel.mif";
#fixel2voxel <directions.mif> split_dir <directions_voxel.mif> ...
print "fixel2voxel fixel_dir/directions.mif split_dir directions_voxel.mif\n";
`fixel2voxel fixel_dir/directions.mif split_dir directions_voxel.mif` unless -e "directions_voxel.mif";
#-crop if necessary, mrcrop...0 offset

#for axial:#mrconvert -stride 1,2,3,4 <afd> ...
print "mrconvert -stride 1,2,3,4 afd_voxel.mif afd_voxel_strides.nii\n";
`mrconvert -stride 1,2,3,4 afd_voxel.mif afd_voxel_strides.nii` unless -e "afd_voxel_strides.nii";
print "mrconvert -stride 1,2,3,4 directions_voxel.mif directions_voxel_strides.nii\n";
`mrconvert -stride 1,2,3,4 directions_voxel.mif directions_voxel_strides.nii` unless -e "directions_voxel_strides.nii";


#for sagittal: 
#mrconvert -stride 3,1,2 <file.mif> <file.nii>  

#steps: ir-diff
#make sure it came directly from dcm
#axial:
#need to mrconvert -stride 1,2,3,4 
print "mrconvert -stride 1,2,3,4 ../ir-analysis/ir-diff.nii ../ir-analysis/ir-diff-strides.nii\n";
`mrconvert -stride 1,2,3,4 ../ir-analysis/ir-diff.nii ../ir-analysis/ir-diff-strides.nii` unless -e "../ir-analysis/ir-diff-strides.nii";

#if sagittal, mrconvert -stride 3,1,2,4 
#then, **flip in x, using fix_IRdiff_sag.py** OR just set strides to -3,1,2,4
#I think this is because of something in Ilana's unshuffling script, and we could change that and not flip here 

#do same to mask as to files that match it (first 3 dims)
#e.g., for axial, I needed to set strides to 1,2,3 
#for sagittal, 3,1,2

#NB!! NODDI processing erases the strides.  So, what I did was resample the file, reversing the negative stride direction (x)
#it might work to use strides -1,2,3, but I haven't tested it.
#I used script fix_NODDI_trans
print "fslswapdim AMICO/NODDI/FIT_ICVF.nii -x y z AMICO/NODDI/FIT_ICVF-strides.nii\n";
`fslswapdim AMICO/NODDI/FIT_ICVF.nii -x y z AMICO/NODDI/FIT_ICVF-strides.nii` unless -e "AMICO/NODDI/FIT_ICVF-strides.nii.gz";
print "fslswapdim noddi/brain_mask-swap.nii  -x y z noddi/brain_mask-swap-strides.nii\n";
`fslswapdim noddi/brain_mask-swap.nii -x y z  noddi/brain_mask-swap-strides.nii` unless -e "noddi/brain_mask-swap-strides.nii.gz";

print "\n---------Check everything is lined up!-------\n";
print "fsleyes ../ir-analysis/ir-diff-strides.nii* afd_voxel_strides.nii* AMICO/NODDI/FIT_ICVF-strides.nii* noddi/brain_mask-swap-strides.nii*\n";
`fsleyes ../ir-analysis/ir-diff-strides.nii* afd_voxel_strides.nii* AMICO/NODDI/FIT_ICVF-strides.nii* noddi/brain_mask-swap-strides.nii*`;



print "\n================Now draw the WM region you want to look at (call it roi.mif)=============\n";
print "dwi2tensor -mask mask.mif $dwi dt.mif \n";
`dwi2tensor -mask mask.mif $dwi dt.mif `;
print "tensor2metric dt.mif -vector V1.mif -modulate FA\n";
`tensor2metric dt.mif -vector V1.mif -modulate FA`;
`mrview V1.mif -fixel.load fixel_dir/afd.mif`;
#fod2dec [ options ]  input output
#print "mrview meanb0.mif -fixel.load fixel_dir/directions.mif";
#`mrview meanb0.mif -fixel.load fixel_dir/directions.mif`;

# remember to convert it properly!
print "mrconvert -stride 1,2,3 roi.mif roi_strides.nii\n";
`mrconvert -stride 1,2,3 roi.mif roi_strides.nii`;
$mask = "roi_strides.nii";

$TR = `more ../info.txt | grep TR | grep -Eo '[0-9]+'`;chomp($TR); #in ms

#print "---------Please enter the afd threshold (default=0.1)   : ";
#$afd_thresh = <STDIN>;
#chomp($afd_thresh);
$afd_thresh=0.2;
print "fibermyelin_pipeline.py -t1 T1.nii -Dpar Dpar.nii -mask $mask -vic AMICO/NODDI/FIT_ICVF-strides.nii.gz -afd afd_voxel_strides.nii -afdthresh $afd_thresh -dirs directions_voxel_strides.nii -IRdiff ../ir-analysis/ir-diff-strides.nii -TIs ../TIcomb.txt -bvals $bvals_ir -bvecs $bvecs_ir -TR $TR  -fixel fixel_dir-output\n";

`fibermyelin_pipeline.py -t1 T1.nii -Dpar Dpar.nii -mask $mask -vic AMICO/NODDI/FIT_ICVF-strides.nii.gz -afd afd_voxel_strides.nii -afdthresh $afd_thresh -dirs directions_voxel_strides.nii -IRdiff ../ir-analysis/ir-diff-strides.nii -TIs ../TIcomb.txt -bvals $bvals_ir -bvecs $bvecs_ir -TR $TR  -fixel fixel_dir-output`;

# mrview V1.mif -fixel.load fixel_dir-output/t1_fixel.nii

## Now generate a file with the sorted T1 by direction
# right now it's hard-coded in fibermyelin_pipeline.py but it could be an input
print "fibermyelin_pipeline.py -t1 T1.nii -Dpar Dpar.nii -mask $mask -vic AMICO/NODDI/FIT_ICVF-strides.nii.gz -afd afd_voxel_strides.nii -afdthresh $afd_thresh -dirs directions_voxel_strides.nii -IRdiff ../ir-analysis/ir-diff-strides.nii -TIs ../TIcomb.txt -bvals $bvals_ir -bvecs $bvecs_ir -TR $TR  -fixel fixel_dir-output -sort > t1sort.txt\n";

`fibermyelin_pipeline.py -t1 T1.nii -Dpar Dpar.nii -mask $mask -vic AMICO/NODDI/FIT_ICVF-strides.nii.gz -afd afd_voxel_strides.nii -afdthresh $afd_thresh -dirs directions_voxel_strides.nii -IRdiff ../ir-analysis/ir-diff-strides.nii -TIs ../TIcomb.txt -bvals $bvals_ir -bvecs $bvecs_ir -TR $TR  -fixel fixel_dir-output -sort > t1sort.txt`;


sub run_matlab {
  my ($command)=@_;
  #system("(echo \"$command\"| matlab -nosplash $logging)");
  open MATLAB,"|matlab -nosplash -nojvm -nodisplay " or die "Can't start matlab!\n";
  #print MATLAB "addpath('$mydir')\n";
  #print MATLAB "addpath(genpath('$niak_dir'))\n" if $niak_dir;
  print MATLAB $command,"\n";
  close MATLAB;
  die "Error in MATLAB!\n" if $?;
}