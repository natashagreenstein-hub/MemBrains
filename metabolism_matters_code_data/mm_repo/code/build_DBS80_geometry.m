%% build_DBS80_geometry.m — region centroids + distance matrix for DBS80
%
%  Produces DBS80_geometry.mat (used by NEMO_oxglc's spatial-autocorrelation
%  null) directly from the DBS80 atlas volume, so the geometry is generated in
%  MATLAB and matches the Methods text:
%
%     "Region centroids were computed as the centre of mass of each parcel in
%      the DBS80 atlas and expressed in MNI millimetre coordinates, and D_ij
%      denotes the Euclidean distance between centroids i and j."
%
%  For each label L = 1..N it takes the centre of mass of the parcel in
%  (0-based) voxel coordinates, maps it to world/MNI millimetres through the
%  NIfTI sform affine, and assembles the N×N Euclidean distance matrix D.
%
%  Output (saved to projDir):
%     coords : N×3  centroid coordinates in MNI mm
%     D      : N×N  Euclidean distance between centroids (mm)
%
%  Requires the Image Processing Toolbox (niftiread / niftiinfo), already used
%  elsewhere in the pipeline (ssim).

clear; close all;
projDir = '/Users/natashag/Documents/Claude/Projects/MemBrain';
N = 80;                         % DBS80 parcels (labels 1..80)

% --- locate the atlas volume (first existing candidate wins) ---
cands = { ...
    fullfile(projDir, 'dbs80symm_2mm.nii.gz'), ...
    '/Users/natashag/Downloads/MemBrain/Maps/derivatives/dbs80symm_2mm.nii.gz', ...
    '/Users/natashag/Downloads/Quantum Consciousness/Tommaso/David/dbs80symm_2mm.nii.gz'};
atlasPath = '';
for c = 1:numel(cands)
    if isfile(cands{c}), atlasPath = cands{c}; break; end
end
assert(~isempty(atlasPath), ...
    ['dbs80symm_2mm.nii.gz not found. Edit the `cands` list to point at the ' ...
     'DBS80 atlas volume on your machine.']);
fprintf('atlas: %s\n', atlasPath);

% --- read volume + sform affine (0-based voxel -> MNI mm) ---
info = niftiinfo(atlasPath);
V    = round(double(niftiread(atlasPath)));     % integer parcel labels
A    = double([info.raw.srow_x; info.raw.srow_y; info.raw.srow_z; 0 0 0 1]);
assert(abs(det(A(1:3,1:3))) > 0, 'sform affine is singular/absent in the header.');

% --- centre of mass per label -> MNI mm ---
coords = nan(N, 3);
nvox   = zeros(N, 1);
for L = 1:N
    idx = find(V == L);
    nvox(L) = numel(idx);
    if isempty(idx), continue; end
    [ii, jj, kk] = ind2sub(size(V), idx);       % 1-based subscripts
    com0  = [mean(ii) - 1, mean(jj) - 1, mean(kk) - 1];   % 0-based centre of mass
    world = A * [com0, 1].';                     % -> homogeneous world coords
    coords(L, :) = world(1:3).';
end
missing = find(nvox == 0);
assert(isempty(missing), 'labels with no voxels: %s', mat2str(missing(:).'));

% --- Euclidean distance matrix (pure MATLAB, no Statistics Toolbox) ---
D = sqrt(sum((reshape(coords, [N 1 3]) - reshape(coords, [1 N 3])).^2, 3));
D = 0.5 * (D + D.');                             % enforce exact symmetry

% --- sanity checks ---
offdiag = D(~eye(N, 'logical'));
fprintf('parcels found: %d/%d\n', sum(nvox > 0), N);
fprintf('distances (mm): min off-diag %.1f | median %.1f | max %.1f\n', ...
        min(offdiag), median(offdiag), max(offdiag));
% region 1 (left caudal ant. cingulate) and 80 (right homologue) should mirror in x:
fprintf('homotopy check  reg 1 = [%.1f %.1f %.1f]   reg 80 = [%.1f %.1f %.1f]\n', ...
        coords(1,:), coords(80,:));

% --- save ---
atlas_source = atlasPath; %#ok<NASGU>
outFile = fullfile(projDir, 'DBS80_geometry.mat');
save(outFile, 'coords', 'D', 'atlas_source', '-v7');
fprintf('wrote %s  (coords %dx3, D %dx%d)\n', outFile, N, N, N);
