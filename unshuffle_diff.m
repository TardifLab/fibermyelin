function unshuffle_diff(file,shift,spacing,TR,numavg,diffdirs,minTI,mode,mb)
[HDR,VOL] = niak_read_nifti(file);
 
%mode='interleave';
% shift=11;
% spacing=99.98; %[ms] TI spacing as output by the sequence
% TR=5700;
% diffdirs=31; %include b=0
% numavg=5;
% mb=3; multi-band factor
%
%dim='sag';
dim='tra';

[x,y,z,t]=size(VOL);

switch mode
    case{'interleaved'}
        if mod(z,2)==0 %even # of slices
          sidx = [2:2:z,1:2:z];%slice index
        else
           sidx = [1:2:z,2:2:z];%slice index
        end
    case{'interleaved-mb'}
        if mod(z/mb,2)~=0 %odd #slices within a block
            %sidx = [1:2:z,2:2:z];%slice index
        else %% even #slices
            sidx_gr1 = [2:2:z/2, 1:2:z/2];
            %sidx_gr2 = [z/2+2:2:z, z/2+1:2:z];
            %sidx = [sidx_gr1;sidx_gr2];
            if (rem(mb, 2) == 0) && (z >= 6)
                z_ = z/mb ;% this shift happens within a group
                % even slices are in the first half of SliceOrder, odd slices in the last half SliceOrder
                % With slices in this order, even SMS factors will generate cross talk if there is more than 1 repetition
                % This is how siemens reorders the lines
%                 % last line of even number slices moves to middle of group
                sidx_gr1(1:(z_/2)) =      sidx_gr1([1:floor(z_/4)    (z_/2) (floor(z_/4) + 1):(z_/2 - 1)]            );
                % first line of odd number slices moves to middle of the group
                sidx_gr1((z_/2 + 1):end) = sidx_gr1([2:(floor(z_/4)+1)     1      (floor(z_/4) + 2):(z_/2)    ] + z_/2);
                sidx = [sidx_gr1;sidx_gr1+z/2];
            end
        end
    case{'ascending'}
        sidx = 1:z;%slice index
end

%TI times
ti=minTI:spacing:TR;
TI=ti(1:z/mb);

% need to flip in z dimension if sagittal
if strcmp(dim,'sag')
    vol_tmp=zeros(x,y,z,t);
    for v=0:(z-1)
        vol_tmp(:,:,v+1,:)=VOL(:,:,end-v,:);
    end
    VOL=vol_tmp;
end

%% Un-shuffle the data by TI time
% The final volume will have b0s all TIs, d1 all TIs, etc
v=zeros(x,y,z,t); %unshuffled volume
% Shift at each average (shift after having acquired all the directions)
if mb>1
    %indr=reshape(sidx,z/mb,mb)'; %reshape so that each column is a MB stack
    indr=sidx;
    ind=indr;
    for i=1:numavg-1
        ind=cat(3,ind,circshift(indr,-shift*i,2)); %index of slices for each average, TR increases along the 3rd dimension
    end
else
    ind=sidx;
    %ind1=sidx1;
    %ind2=sidx2;
    for i=1:numavg-1
        ind=[ind; circshift(sidx,-shift*i)]; %index of slices for each average
        %ind1=[ind1; circshift(sidx1,-shift*i)]; %index of slices for each average
        %ind2=[ind2; circshift(sidx2,-shift*i)]; %index of slices for each average       
    end
end
%ind=cat(2,ind1,ind2);
% Sort by diffusion direction i.e. we'll have b0 at all TIs, d1 at all TIs
for d=1:diffdirs
    didx=d:diffdirs:numavg*diffdirs; %index of a particualr direction e.g. all the b=0s
    dir=VOL(:,:,:,didx); %get all time points of 1 direction
    
    TIm=zeros(z,numavg); %matrix of TIs, different set of TIs for each slice
    % Each slice will have different set ot TIs (if shift by more than 1)
    for s=1:z
        if mb>1
            [r,t_s_idx,t] = ind2sub(size(ind),find(ind == s));%idx of TI for slice s
        else
            [t_s_idx,c] = find(ind'==s); %idx of TI for slice s
        end
        t_s=TI(t_s_idx); %TIs for slice s
        [t_s_idx_sort,I]=sort(t_s_idx); %sort chronologically in time
        v(:,:,s,(d-1)*numavg+1:(d-1)*numavg+numavg) = dir(:,:,s,I); %reorder
        TIm(s,:) = sort(t_s);
        for i=1:numavg-1
            TRa(s,i)=TR-t_s(i)+t_s(i+1);% each slice's TR (it's the same for every average)
        end
    end
end

hdr=HDR;
name=strcat('unshuffle_diff_',file);
hdr.file_name=name;
niak_write_nifti(hdr,v);
save('TI.txt', 'TIm', '-ASCII');
save('TRa.txt', 'TRa', '-ASCII');

 