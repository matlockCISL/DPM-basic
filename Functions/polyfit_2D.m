function [P,A] = polyfit_2D(P, A)
%% polyfit_2D Function
%
% Purpose: This function is a wrapper providing a 2D polynomial fit to the
% phase and absorption of a measurement. 
%
% Inputs: P - MxM image containing the measured phase
%         A - MxM image containing the measured amplitude
%
% Outputs: P - MxM image of the filtered phase
%          A - MxM image of the filtered absorption
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain image size, make pixel index over image FOV
sz = size(P);
idx = [-sz(1)/2:sz(2)/2-1]';

% Convert absorption back to amplitude
A = exp(-A); 

% Filter absorption
[xA,yA,zA] = prepareSurfaceData(idx,idx,A);
fA = fit([xA,yA],zA,'poly22');
A_nm = reshape(feval(fA, [xA, yA]), [sz(1), sz(2)]);
A = -log(A./A_nm);

% Filter phase
[xP,yP,zP] = prepareSurfaceData(idx,idx,P);
fP = fit([xP,yP],zP,'poly22');
P_nm = reshape(feval(fP, [xP,yP]), [sz(1), sz(2)]);
P = P - P_nm;

end  % EOF