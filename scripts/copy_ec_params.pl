#!/usr/bin/perl -w

### Copy eddy current parameters from the b1000 shell of the HARDI
## so that they can be used by IR-diff
#
##########################################################
##########################################################
# Note that this only works for the MWC dataset, it expects the data a certain way
# meaning we have the exact same directions in both datasets
#NB: the data are assumed to have been acquired with diffusion fastest varying, and unshuffled to TI fastest varying.\n\
#    i.e., bvals/bvecs have diffusion fastest varying, and input IRdiff images do not.\n\

# We have to:
#   - get the b1000 shell
#   - remove all the extra b=0
#   - set all registration parameters to 0 (1st 6 columns), we only want EC params
#   - repeat each entry 4 times because IR-diff has TI fastest varying
#   - make an appropriate bvals and bvecs with TI fastest varying

###########################################
# Created by Ilana Leppert Feb 2022
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

Get the eddy parameter corrections from the b1000 shell to pass to eddy for IR-diff


Usage: $program -ec_in output.eddy_parameters -bvals_in bvals-diff-fast -bvecs_in bvecs-diff-fast -ec_out only-b1000-ec-params.txt -bvals new-bvals -bvecs new-bvecs
-help for options

USAGE

my @args_table = (["-ec_in","string",1,\$ec_in,"name of input ec parameters"],
                ["-bvals_in","string",1,\$bvals_in,"name of input bvals(of full HARDI)"],
                ["-bvecs_in","string",1,\$bvecs_in,"name of input bvecs(of full HARDI)"],
                ["-ec_out","string",1,\$ec_out,"name of output EC parameters"],
                ["-bvals","string",1,\$bvals,"name of output bvals for IR-diff (TI fastest varying)"],
                ["-bvecs","string",1,\$bvecs,"name of output bvecs for IR-diff (TI fastest varying)"]
            );
#["-clobber","boolean",undef,\$clobber,"clobber all output files (currently not implemented, sorry)"],#HERE not implemented);


Getopt::Tabular::SetHelp ($Usage, '');
my @args;
GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
die $Usage unless $#ARGV >=0;

open(FILE,"< $ec_in") or die "can't open $ec_in: $!\n";
@EC_in = <FILE>; #slurp
close(FILE);
open(FILE,"< $bvals_in") or die "can't open $bvals_in: $!\n";
$B = <FILE>; #slurp
@Bvals_in = split(/ /,$B);
close(FILE);
open(FILE,"< $bvecs_in") or die "can't open $bvecs_in: $!\n";
@dirs = <FILE>; #slurp this will have 3 elements
@x = split(/ /,$dirs[0]);
@y = split(/ /,$dirs[1]);
@z = split(/ /,$dirs[2]);
close(FILE);

#create output files
open(ECOUT,"> $ec_out") or die "can't create $ec_out: $!\n";
open(BVALS,"> $bvals") or die "can't create $bvals: $!\n";
open(BVECS,"> $bvecs") or die "can't create $bvecs: $!\n";

#add an input for the b=0, 4 times for each each TI
for ($j=0; $j<4; $j++){
    print   ECOUT ("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n");
    push(@bval,0);
    push(@x_irdiff,0);
    push(@y_irdiff,0);
    push(@z_irdiff,0);
}
#start at 1 because the 1st frame it's actually a b=0 from IR-diff
for ($i=1; $i<scalar(@Bvals_in); $i++){
    if (($Bvals_in[$i] > 900) & ($Bvals_in[$i] < 1500)){
        #print "--$i bval = $Bvals_in[$i] \n $EC_in[$i]\n";
        @EC_params = split(/  /,$EC_in[$i]);
        @only_EC = @EC_params[6..$#EC_params]; #the 1st 6 are the movt params which we don't want to keep
        #print "only EC params @only_EC\n";
        #print "xdir is: $x[$i]\n";
        for ($j=0; $j<4; $j++){ #repeat 4timeson once for each TI (TI fastest varying)
            # print 6 0s for motion and then append these EC params
            print ECOUT ("0 0 0 0 0 0 @only_EC");
            # get the diff direction that goes with it (4 times)
            push(@bval,1000);
            push(@x_irdiff,$x[$i]);
            push(@y_irdiff,$y[$i]);
            push(@z_irdiff,$z[$i]);
        }
    }
}
close ECOUT;

#now create that bvals and bvecs we need
print BVALS ("@bval");
print BVECS ("@x_irdiff\n@y_irdiff\n@z_irdiff");

