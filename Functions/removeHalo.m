function [img] = removeHalo(img, pm, fwhm)
%% removeHalo Function
%
% Purpose: This function removes halo artifacts from a phase image using
%          the approach described in the 2018 work by Kandel et al. 
%          "Real-time halo correction in phase contrast imaging". The basic
%          design of this approach is to generate 3 filters in the Fourier 
%          space removing low frequency features in the image via a Hilbert
%          transform and extract only the maximum values from each filtered
%          image. This approach assumes a separation of the object's low 
%          and high spatial frequency features that allows for the image to
%          be filtered based on its directional spatial derivative 
%          information.
%
% Inputs: img - MxM image to be processed for halo removal.
%         pm - structure containing microscope parameters
%         fwhm - laser wavelength full width half maximum (um)     
%           

%% Handle Declarations
FT = @(x) fftshift(fft2(ifftshift(x)));
iFT = @(x) fftshift(ifft2(ifftshift(x)));

%% Variable initialization
sz = size(img);  % Obtain image pixel dimensions
du = 2 * pi * [1/(sz(1) * (pm.dx/pm.Mtot)), 1/(sz(2) * (pm.dy/pm.Mtot))];  % Determine K-space sampling
Lc = pm.lmd^2/(2*fwhm);  % Calculate laser coherence length


% Generate k-space grids
[U,V] = meshgrid(du(2) .* [-sz(2)/2:sz(2)/2-1], du(1) .* [-sz(1)/2:sz(1)/2-1]);
diagUV = abs(U + V);

%% Generate |k| filters
if(length(sz) > 2)
    k_0 = repmat((abs(V) .* (abs(V) < 1/Lc) + 1/Lc .* (abs(V) >= 1/Lc)), [1 1 sz(3)]);
    k_90 = repmat((abs(U) .* (abs(U) < 1/Lc) + 1/Lc .* (abs(U) >= 1/Lc)), [1 1 sz(3)]);
    k_45 = repmat(((diagUV .* (diagUV < 1/Lc)) + 1/Lc .* (diagUV >= 1/Lc)), [1 1 sz(3)]);
else
    k_0 = (abs(V) .* (abs(V) < 1/Lc) + 1/Lc .* (abs(V) >= 1/Lc));
    k_90 = (abs(U) .* (abs(U) < 1/Lc) + 1/Lc .* (abs(U) >= 1/Lc));
    k_45 = ((diagUV .* (diagUV < 1/Lc)) + 1/Lc .* (diagUV >= 1/Lc));
end

%% Apply filters
% Obtain 2D Fourier transforms of reconstructions
fimg = FT(img);

% Apply filters to each object slice
if(length(sz) > 2)
    tmpimg(:,:,:,1) = real(iFT(k_0 .* fimg));
    tmpimg(:,:,:,2) = real(iFT(k_45 .* fimg));
    tmpimg(:,:,:,3) = real(iFT(k_90 .* fimg));
    tmpimg(:,:,:,4) = img;
    img = max(tmpimg, [], 4);
else
    tmpimg(:,:,1) = real(iFT(k_0 .* fimg));
    tmpimg(:,:,2) = real(iFT(k_45 .* fimg));
    tmpimg(:,:,3) = real(iFT(k_90 .* fimg));
    tmpimg(:,:,4) = img;
    img = max(tmpimg, [], 3);
end

end