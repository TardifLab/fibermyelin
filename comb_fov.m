function comb_fov(vol1,vol2,TI, output)
% Want to combine 2 volumes acquired with an FOV shift

% vol1 and vol2: 2 interleaved FOVs, 
%  TI : txt file with list of TIs for each slice (will also be interleaved)
% output: name of output nii file with interleaved volumes and corrected spacing
%
% will also output TIcomb.txt, the combined TI text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create text file with combined list of TIs for each slice
TIs = importdata(TI);
[r,c] = size(TIs);
combTI= zeros(2*r,c);
combTI(1:2:end,:)=TIs;
combTI(2:2:end,:)=TIs;
save('TIcomb.txt', 'combTI', '-ASCII');

[HDR,VOL1] = niak_read_nifti(vol1);
[HDR,VOL2] = niak_read_nifti(vol2);
[x,y,z,t]=size(VOL1);

% make a volume double the number of slices
vol_tmp=zeros(x,y,2*z,t);
vol_tmp(:,:,2:2:end,:) = VOL1; % if volume 2 was shifted down from volume 1
vol_tmp(:,:,1:2:end,:) = VOL2;

hdr=HDR;
tmp='tmp-ir-diff.nii';
name=strcat(tmp);
hdr.file_name=name;
hdr.info.dimensions=[x y 2*z t];
niak_write_nifti(hdr,vol_tmp);

%in the terminal, change the slice spacing to 1/2 what the original was
vsize = hdr.info.voxel_size(3)/2;
cmd=sprintf('mrconvert %s -vox ,,%f %s',name,vsize,output);
cmd
system(cmd);


% the other file from mrtrix is meed up, read it in here and write out
% again, can't find another way to modify the starts!
%[HDR,VOL3] = niak_read_nifti('ir-diff-reformat.nii');
%[HDR2,VOL4] = niak_read_nifti('meanb0-24slices.nii');
%HDR.file_name='meanb0-24slices-mat.nii';
%niak_write_nifti(HDR,VOL4);
