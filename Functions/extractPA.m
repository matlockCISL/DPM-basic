function [recon, grt, yshift]  = extractPA(img, k0, pm, tog)
%% phase_recover function
%
%   Purpose: This function recovers the complex field from the DPM image.
%   This code currently assumes a square image, processing rectangular
%   images will require modification. 
%
%       Inputs:
%               img - MxM matrix containing normalized DPM image.
%               k0 - scalar value defining wavenumber
%               pm - structure containing the following parameters
%                   NA - scalar value defining collection objective NA
%                   dx - scalar value containing x-axis pixel size at image
%                        plane (um)
%                   grt - scalar value for grating frequency positioned at the
%                         imaging plane (um^-1)
%                   y_shift - scalar value denoting number of pixels to shift
%                             the spectra in the y direction.
%       Outputs:
%               recon - MxM matrix containing recovered complex field.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Handle Declaration
FT = @(x) ifftshift(fft2(fftshift(x)));
iFT = @(x) ifftshift(ifft2(fftshift(x)));


%% Variable declaration
[M,N]=size(img); % the image size of the image

dkx = 2*pi/(pm.dx/pm.Mtot .* M);  % Obtain k-vector sampling
pk_max = round(pm.NA*k0 / dkx);  % maximum pixel value for NA passband


fimg=FT(img); % Calculate the Frequency domain of the image

[NX,NY] = meshgrid(-M/2:M/2-1,-N/2:N/2-1); % Generate pixel grid for making spatial filter


%% Extract cross-correlation term, shift to center



% Automatically find field term or do it manually
if(tog.auto)
    if(tog.order == 1)
        % Perform peak search only in right half of image
        pmask = NX > 0;
    else
        pmask = NX < 0;
    end
    tmp = log10(abs(fimg)) .* pmask;
    shift = zeros(1,2);
    [~, shift(1)] = max(max(tmp, [], 1), [], 2);
    [~, shift(2)] = max(max(tmp, [], 2), [], 1);
    
    shift = shift - [M/2+1, N/2+1];

    % pm.grt = (shift(1) * dkx) / (2 * pi);
    % pm.yshift = shift(2);
else
    % grating shift position
    grt_max = round((2 * pi * pm.grt)/dkx);  

    % Find center of translated signal
    shift = [grt_max, pm.yshift];
end    
    % Filter out cross-correlation term
    mask_1 = sqrt(((NX-shift(1))./pk_max).^2 + ((NY+shift(2))./pk_max).^2) <= 1;
    
    % Apply filter to extract cross-correlation signal, shift to image
    % center
    fimg = circshift(fimg .* mask_1, [shift(2), -tog.order * shift(1)]);
    % 
    % figure(11);
    % imagesc(log10(abs(fimg)));axis image off;
    
    % Apply iFT to obtain complex field
    recon=iFT(fimg);
    grt = pm.grt;
    yshift = pm.yshift;
end

