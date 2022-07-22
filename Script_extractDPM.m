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
%% Add folders to path
addpath('Functions\');  % Adds folder containing relevant processing functions

%% Folder Name Declarations

% Set data, reference, and save directories
genfol = ['G:\My Drive\Data\RBC_SCD\220720\25x_430nm_5um_30pBSA_p90_FOV_']; % Generic folder for raw data location
reffol = 'G:\My Drive\Data\RBC_SCD\220720\25x_430nm_5um_30pBSA_p90_Blank_1'; % Folder location for blank/reference image
svdir = 'G:\My Drive\Data\RBC_SCD\220720\25x_430nm_5um_30pBSA_p90_Processed\FOV_'; % save folder location

lbl_data = 'Raw';  % Name for raw data measurements (assumed to be identical or the same name with different numbers for all folders used in processing)
lFOV = [1, 20];  % Set FOV range to process for folder (use [1 1] for processing only first FOV)
nmeas = 1;  % Number of images to process within each FOV folder
nref = 1;  % Number of reference images to process within the reffol location

%% Variable Declarations

% Microscope parameters
pm.dx = 4.5;  % x Pixel size at camera plane (um)
pm.dy = 4.5;  % y Pixel size at camera plane (um)
pm.Mo = 25;  % Microscope objective magnification
pm.Mf = 300/75;  % 4F System magnification
pm.Mtot = pm.Mo * pm.Mf;  % Total system magnification
pm.lmd = 0.43;  % System imaging wavelength (um)
pm.NA = 1;  % System collection NA
pm.grt = pm.Mo * 16.3e-2;  %7.96e-2;  %7.94e-2 * Mo;  %7.91e-2 * Mo;  % Grating line pairs per um after projection to image plane
pm.yshift = 4;  % y-axis translation correction based on exp. observation (pixel)

k0 = 2*pi/pm.lmd;  % Generate wavenumber

% Toggle image saving
tog.saveImg = 1;

%% Make save folder containing processing parameters, grating, etc.
svfol = ['Proc_2DPolyFit_NA' num2str(pm.NA) '_grt' num2str(pm.grt/pm.Mo * 1e3) '_ys' num2str(pm.yshift)];

%% Load Background image and extract reference phase, amplitude

for nt = 1:nref
    disp(['Processing Reference Frame ' num2str(nt) '...']);
    bkgnd = double(imread([reffol '\' lbl_data '.tiff']));
    
    % Allocate reference image stacks on first reference iteration
    if (nt == 1)
       A_ref = zeros(size(bkgnd));
       P_ref = zeros(size(bkgnd));
       recon_ref = zeros(size(bkgnd));
    end

    % Extract background field from interferogram
    recon = extractPA(bkgnd, k0, pm);
    
    %Separate out magnitude
    A = abs(recon);

    % Separate out and unwrap phase
    Phi = unwrap2(angle(gather(recon)));
    
    % Accumulate absorption, phase to obtain average
    A_ref = A_ref + A;
    P_ref = P_ref + Phi;
end


% Obtain average reference image, image size
P_ref = P_ref/nref;
A_ref = A_ref/nref;
sz = size(P_ref);
clear Phi_f A_f P_nm A_nm fA fP xP yP zP Phi A recon bkgnd x xA yA zA

%% Load Images and extract phase, amplitude
for nf = lFOV(1):lFOV(2)
    for nm = 1:nmeas
        disp(['FOV ' num2str(nf) ' of ' num2str(lFOV(2)) ', Image ' num2str(nm)]);
        
        % Load image (switch commented out code if you have multiple images in same FOV)
%         img = double(imread([genfol num2str(nf) '\' lbl_data sprintf('%08d', nm-1) '.tiff']));
        img = double(imread([genfol num2str(nf) '\' lbl_data '.tiff']));


        %% Extract phase and amplitude information from images

        % Recover complex field
        recon = extractPA(img, k0, pm);
   
        %Separate out magnitude
        A = abs(recon);

        % Separate out and unwrap phase
        Phi = unwrap2(angle(gather(recon)));

        % Remove reference information
        Phi = real(Phi - P_ref);
        A =-log(A./A_ref);      
        
        % Apply 2D polynomial filter to remove low-frequency features
        [Phi, A] = polyfit_2D(Phi, A);

        % Median filter images to PSF level
        A_f = (medfilt2(real(A), [3, 3]));
        Phi_f = (medfilt2(Phi, [3, 3]));
         
        %% Save data into .mat, tiff images
        mkdir([svdir num2str(nf) '\' svfol]);
        save([svdir num2str(nf) '\' svfol '\Reconstruction_' num2str(nm-1) '.mat'],'Phi_f', 'A_f', 'pm', '-v7.3');

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
