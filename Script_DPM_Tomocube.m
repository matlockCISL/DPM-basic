%{
    Main script for Extracting phase, absorption from DPM measurements

    Script Initial Date: 210922, revised 220722
    Script Author: Alex Matlock
    Purpose: This script extracts the phase and absorption from a set of
    raw interferograms. The script assumes the dataset has a reference
    image for removing the background and can process multiple
    fields-of-view (FOV) as well as multiple measurements of the same FOV
    within the same folder. 
%}
clear all;
close all;
clc;

normImg = @(img) (img - min(min(img)))./(max(max(img)) - min(min(img)));
%% Add folders to path
addpath('Functions\');  % Adds folder containing relevant processing functions

%% Folder Name Declarations

% Set data, reference, and save directories
genfol = ['G:\My Drive\Data\Grating Measurements\250902\Raw Data\Wvgd_Lg_1']; % Generic folder for raw data location
svdir = 'G:\My Drive\Data\Grating Measurements\250902\Proc\Wvgd_Lg_'; % save folder location
svlbl = 'tomo_test';  % User-defined label to add to save folder name

lbl_data = '001';  % Name for raw data measurements (assumed to be identical or the same name with different numbers for all folders used in processing)
lFOV = [1,1];  % Set FOV range to process for folder (use [1 1] for processing only first FOV)
nmeas = 1;  % Number of images to process within each FOV folder
nref = 1;  % Number of reference images to process within the reffol location

%% Variable Declarations

% Microscope parameters
pm.dx = 5.5385;  % x Pixel size at camera plane (um)
pm.dy = 5.5385;  % y Pixel size at camera plane (um)
pm.Mo = 58.3;  % Microscope objective magnification
pm.Mf = 1;  % 4F System magnification
pm.Mtot = pm.Mo * pm.Mf;  % Total system magnification
pm.lmd = 0.532;  % System imaging wavelength (um)
pm.NA = 1.2;  % System collection NA
pm.grt = pm.Mo * 16.28e-2;  %7.96e-2;  %7.94e-2 * Mo;  %7.91e-2 * Mo;  % Grating line pairs per um after projection to image plane
pm.yshift = 5;  % y-axis translation correction based on exp. observation (pixel)

k0 = 2*pi/pm.lmd;  % Generate wavenumber

% Toggle image saving
tog.saveImg = 1;
tog.auto = 1;  % Toggles automatic center finding for Hilbert transform
tog.order = 1;  % Sets whether to use +1 (Right) or -1 (left) Fourier spectra
tog.halo = 0;  % Toggle halo artifact removal for data
tog.patch2D = 0; % Toggle patchwise 2D filtering of image
tog.selfRef = 2; % Toggle where nonzero value defines level of polynomial order fitting to implement (max 5)

%% Update save label
svlbl = ['Proc_' svlbl];
if(tog.halo)
    svlbl = [svlbl '_Halo'];
end
if(tog.patch2D)
    svlbl = [svlbl '_patch2DFit'];
elseif(tog.selfRef > 0)
    svlbl = [svlbl '_selfFit_Poly' num2str(tog.selfRef)];
else
    svlbl = [svlbl '_Global2DFit_nm'];
end
%% Load Background image and extract reference phase, amplitude

for nt = 1:nref
    disp(['Processing Reference Frame ' num2str(nt) '...']);
    if(tog.selfRef > 0)
        bkgnd = double(imread([genfol '\Interferogram\' lbl_data '.png']));
    else
        bkgnd = double(imread([genfol '\Background\' lbl_data '.png']));
    end
    bkgnd = normImg(bkgnd);
    % Allocate reference image stacks on first reference iteration
    if (nt == 1)
       A_ref = zeros(size(bkgnd));
       P_ref = zeros(size(bkgnd));
       recon_ref = zeros(size(bkgnd));
    end

    % Extract background field from interferogram
    [recon, pm.grt, pm.yshift]  = extractPA(bkgnd, k0, pm, tog);
    
    %Separate out magnitude
    A = abs(recon);

    % Separate out and unwrap phase
    Phi = unwrap2(angle((recon)));
    
    % Accumulate absorption, phase to obtain average
    A_ref = A_ref + A;
    P_ref = P_ref + Phi;
end




% Obtain average reference image, image size
P_ref = P_ref/nref;
A_ref = A_ref/nref;
sz = size(P_ref);

if(tog.selfRef > 0)
% % Lowpass filter the reference phase, ax = [-size(A_ref,1)/2:size(A_ref,2)/2-1]';
x = [-size(A_ref,1)/2:size(A_ref,2)/2-1]';
[xA,yA,zA] = prepareSurfaceData(x,x,A_ref);
fA = fit([xA,yA],zA,['poly' num2str(tog.selfRef * 10 + tog.selfRef)]);
A_ref = reshape(feval(fA, [xA, yA]), [size(A_ref,1), size(A_ref,2)]);

% Filter phase
[xP,yP,zP] = prepareSurfaceData(x,x,P_ref);
fP = fit([xP,yP],zP,['poly' num2str(tog.selfRef * 10 + tog.selfRef)]);
P_ref = reshape(feval(fP, [xP,yP]), [size(P_ref,1), size(P_ref,2)]);
% % 
end

clear Phi_f A_f P_nm A_nm fA fP xP yP zP Phi A recon x xA yA zA

%% Load Images and extract phase, amplitude
for nf = lFOV(1):lFOV(2)
    for nm = 1:nmeas
        disp(['FOV ' num2str(nf) ' of ' num2str(lFOV(2)) ', Image ' num2str(nm)]);
        
        % Load image (switch commented out code if you have multiple images in same FOV)
%         img = double(imread([genfol num2str(nf) '\' lbl_data sprintf('%08d', nm-1) '.tiff']));
        img = double(imread([genfol '\Interferogram\' lbl_data '.png']));

        img = normImg(img);
        %% Extract phase and amplitude information from images

        % Recover complex field
        [recon, pm.grt, pm.yshift] = extractPA(img, k0, pm, tog);
   
        %Separate out magnitude
        A = abs(recon);

        % Separate out and unwrap phase
        Phi = unwrap2(angle((recon)));

        % Remove reference information
        Phi = real(Phi - P_ref);
        A =-log(A./A_ref);      
        
    % Toggle global 2D fit or patchwise 2D fitting for background removal
    if(tog.patch2D)
        
        %         % Perform patch-based correction
        [Phi_stack, M, sep] = segmentObj(Phi, 1024, 512);
        [A_stack, M, sep] = segmentObj(A, 1024, 512);
        for k = 1:size(Phi_stack, 3)
            A_tmp = exp(-A_stack(:,:,k));
            x = [-size(A_tmp,2)/2:size(A_tmp,2)/2-1]';
            y = [-size(A_tmp,1)/2:size(A_tmp,1)/2-1]';
            [xA,yA,zA] = prepareSurfaceData(x,y,A_tmp);
            fA = fit([xA,yA],zA,'poly22');
            A_nm = reshape(feval(fA, [xA, yA]), [size(A_tmp,1), size(A_tmp,2)]);
            A_stack(:,:,k) = real(-log(A_tmp./A_nm));

            [xP,yP,zP] = prepareSurfaceData(x,y,Phi_stack(:,:,k));
            fP = fit([xP,yP],zP,'poly22');
            P_nm = reshape(feval(fP, [xP,yP]), [size(Phi_stack(:,:,k),1), size(Phi_stack(:,:,k),2)]);
            Phi_stack(:,:,k) = Phi_stack(:,:,k) - P_nm;

        end
        A = blendObj(A_stack, [M, M], [sep, sep]);
        Phi = blendObj(Phi_stack, [M, M], [sep, sep]);
    else
        % Apply 2D polynomial filter to remove low-frequency features
        [Phi, A] = polyfit_2D(Phi, A);
    end

        Phi_h = removeHalo(Phi, pm, 0.005);
        A_h = removeHalo(A, pm, 0.005);

        % If toggled, uses halo removal in data. If not, keeps halo removal
        % images for easier segmentation
        if(tog.halo)
            % Median filter images to PSF level
            A_f = (medfilt2(real(A_h), [3, 3]));
            Phi_f = (medfilt2(Phi_h, [3, 3]));
        else
            % Median filter images to PSF level
            A_f = (medfilt2(real(A), [3, 3]));
            Phi_f = (medfilt2(Phi, [3, 3]));
        end

        %% Save data into .mat, tiff images

        % Make save folder containing processing parameters, grating, etc.
        svfol = [svlbl '_NA' num2str(pm.NA) '_grt' ...
                  num2str(pm.grt/pm.Mo * 1e3) '_ys' num2str(pm.yshift)];


        mkdir([svdir num2str(nf) '\' svfol]);
        % save([svdir num2str(nf) '\' svfol '\Reconstruction_' num2str(nm-1) '.mat'],'Phi_f', 'A_f', 'pm', '-v7.3');
        save([svdir num2str(nf) '\' svfol '\Reconstruction_' num2str(nm-1) '.mat'],'Phi_f','A_f', 'Phi_h','A_h','pm', '-v7.3');
        if(tog.saveImg)
            mkdir([svdir num2str(nf) '\' svfol '\Phase']);
            mkdir([svdir num2str(nf) '\' svfol '\Abs']);

            phifile = Tiff([svdir num2str(nf) '\' svfol '\Phase\Phase_' num2str(nm-1) '.tiff'],'w');
            writeTiff_32bGray(phifile, Phi);

            ampfile = Tiff([svdir num2str(nf) '\' svfol '\Abs\Abs_' num2str(nm-1) '.tiff'],'w');
            writeTiff_32bGray(ampfile, A);
        end
        clear raw_img img Phi_f A_f P_nm A_nm
    end  % End of measurement within FOV loop
    
    
end  % end of FOV loop
