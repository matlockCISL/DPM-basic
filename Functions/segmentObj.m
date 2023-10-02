function [stack, M, sep] = segmentObj(obj, subsz, sep)
%% segmentObj Function
%
%   Purpose: This function segments the image of interest into a 3D image
%   stack based on the image subset size (subsz) they are using and the
%   pixel overlap (sep) between each image segment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate M
M = floor((size(obj,1) + sep - 2 * subsz)/(subsz - sep));

% Calculate starting pixel for subset generation
st_px = 0.5 * ( size(obj, 1) - ((M+2)*subsz - (M+1)*sep)) + 1;

%% Generate image stack
cnt = 1;
for k = 1:(M+2)
    for j = 1:(M+2)
        stack(:,:,cnt) = obj(st_px+(k-1)*(subsz-sep):(st_px+subsz-1)+(subsz-sep)*(k-1), ...
                                    st_px+(j-1)*(subsz-sep):(st_px+subsz-1)+(subsz-sep)*(j-1), :);
        cnt = cnt + 1;
    end
end

end