function pct_diff_nakul(dcmfolder1, dcmfolder2, outfolder)

%%% SUMMARY
%%% Convert a pair of dicom images and compute percent difference in signal
%%% intensities between them.

%%% CONTACT INFO
%%% Author:     Patrick Fisher, 20211003
%%% Email:      patrick@nru.dk

%%% INPUTS
%%% dcmfolder1: folder containing dicom files for precontrast image
%%% dcmfolder2: folder containing dicom files for postcontrast image
%%% outfolder:  where to save output data (e.g, batch script, nifti images,
%%%             pct diff images)

%%% EXAMPLES

%%% pct_diff_nakul()
%%% pct_diff_nakul(dcmfolder1, dcmfolder2, outfolder)

%% Process input

if nargin == 0
    dcmfolder1 = uigetdir(pwd, 'Select folder containing pre-contrast dicom files');
    dcmfolder2 = uigetdir(pwd, 'Select folder containing post-contrast dicom files');
    outfolder = uigetdir(pwd, 'Select output directory');
elseif nargin == 3
    if any([~exist(dcmfolder1,'dir') ~exist(dcmfolder2,'dir') ~exist(outfolder,'dir')])
        error('Not all specified folder exist')
    end
else
    error(['Unexpected number of inputs (' num2str(nargin) ')'])
end

% convert precontrast
unix(['/usr/local/fsl/fslpython/envs/fslpython/bin/dcm2niix' ...
    ' -o ' outfolder ...
    ' -f precontrast' ...
    ' -w 1 ' ...
    dcmfolder1]);

% convert postcontrast
unix(['/usr/local/fsl/fslpython/envs/fslpython/bin/dcm2niix' ...
    ' -o ' outfolder ...
    ' -f postcontrast ' ...
    ' -w 1 ' ...
    dcmfolder2]);

% image file names
i1 = fullfile(outfolder, 'precontrast.nii');
i2 = fullfile(outfolder, 'postcontrast.nii');

% coregister/reslice spm12 batch script
run_spm12
matlabbatch{1} = struct();
matlabbatch{1}.spm.spatial.coreg.estwrite.ref{1} =              i1; % precontrast is target image
matlabbatch{1}.spm.spatial.coreg.estwrite.source{1} =           i2; % postcontrast is move
matlabbatch{1}.spm.spatial.coreg.estwrite.other{1} =            '';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun =   'nmi'; % norm mut info
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep =        [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol =        [repmat(0.02,1,3) repmat(0.001,1,3) repmat(0.01,1,3) repmat(0.001,1,3)]; % spm12 default
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm =       [7 7]; % spm12 default
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp =     5; % 5th degree spline
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap =       [0 0 0]; % no wrap
matlabbatch{1}.spm.spatial.coreg.setwrite.roptions.mask =       0; % no mask
matlabbatch{1}.spm.spatial.coreg.setwrite.roptions.prefix =     'r'; % filename prefix

batchfile = fullfile(outfolder, 'batch_reslice.mat');
save(batchfile,'matlabbatch')
spm_jobman('run',batchfile) % run batch job

[pn,fn,ext] = fileparts(i2);
ri2 = fullfile(pn, ['r' fn ext]); % resliced postcontrast image

hdr1 = spm_vol(i1); % precontrast hdr file
hdr2 = spm_vol(ri2); % postcontrast hdr file
v1 = spm_read_vols(hdr1); % precontrast voxel values
v2 = spm_read_vols(hdr2); % postcontrast voxel values

% hdr and voxel values for pct diff img
hdrnew = hdr1; % borrow hdr info from precontrast
hdrnew.fname = fullfile(outfolder, 'pctdiff.nii'); % change filename
hdrnew.private.dat.fname = hdrnew.fname; % change filename
vnew = zeros(hdr1.dim); % voxel values

% limit pct diff computatino to subset of voxels
v1m = mean(mean(mean(v1)));
v2m = mean(mean(mean(v2)));
elem = find(v1>(v1m*0.5) & v2>(v2m*0.5)); % compute pct diff only at voxels above 50% of mean image

% compute pct diff
for i = 1:numel(elem)
    vnew(elem(i)) = ((v2(elem(i))-v1(elem(i)))/v1(elem(i)))*100; % ((post - pre)/pre)*100
end

spm_write_vol(hdrnew,vnew); % write new image file

fprintf('Done!\n')
