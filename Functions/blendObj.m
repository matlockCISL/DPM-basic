function out = blendObj(obj, M, sep)
%%%
%
%Function: blendObj
%Purpose: This function implements alpha blending on the 4D input stack to
%         recover a full FOV image.
%
%Inputs: obj: 4D image stack to be blended into 3D image
%        M: 2x1 vector denoting number of middle tiles for blended object
%               M = [Vert, Horz]. If value of M = -1, means no stitching occurs along that
%               direction
%        sep: 2x1 vector denoting amount of pixel overlap between image tiles
%               sep = [Vert Sep, Horz Sep]
% 
%%%
% Check if 2D or 3D stitching is required
if(ndims(obj) == 3)
    tog2D = 1;
else
    tog2D = 0;
end
imsz = [size(obj,1), size(obj,2), size(obj,3)];
finalsz = [(M(1)+2) * imsz(1) - (M(1)+1) * sep(1), ...
           (M(2)+2) * imsz(2) - (M(2)+1) * sep(2)];

%%  Generate blending filters

blender = ones(imsz(1), imsz(2), 4);

% Generate decreasing filter
dropx = sep(2):-1:1;
dropy = sep(1):-1:1;
blender(:, end-sep(2)+1:end, 1) = 1/sep(2) * repmat(dropx,imsz(1),1);
blender(end-sep(1)+1:end, :, 3) = 1/sep(1) * repmat(dropy',[1, imsz(2)]);

%Generate increasing filters
dropx = 1:sep(2);
dropy = 1:sep(1);
blender(:, 1:sep(2), 2) = 1/sep(2) * repmat(dropx,imsz(1),1);
blender(1:sep(1), :, 4) = 1/sep(1) * repmat(dropy',[1, imsz(2)]);

clear dropx dropy

%% Stitch image together
cnt = 1;

% Check for 2D or 3D blending
if(tog2D)
    out = zeros(finalsz(1), finalsz(2));
    for k = 1:M(1)+2
        for q = 1:M(2)+2
            if(q > 1)
                obj(:,:,cnt) = blender(:,:,2) .* obj(:,:,cnt);
            end
            if (q < (M(2)+2))
                obj(:,:,cnt) = blender(:,:,1) .* obj(:,:,cnt);
            end
            if (k < (M(1)+2))
                obj(:,:,cnt) = blender(:,:,3) .* obj(:,:,cnt);
            end
            if (k > 1)
                obj(:,:,cnt) = blender(:,:,4) .* obj(:,:,cnt);
            end

            out(1+(k-1)*(imsz(1)-sep(1)):imsz(1)+(k-1)*(imsz(1)-sep(1)), ...
                       1+(q-1)*(imsz(2) - sep(2)):imsz(2) + (q-1)*(imsz(2)-sep(2))) ...
                     = out(1+(k-1)*(imsz(1)-sep(1)):imsz(1)+(k-1)*(imsz(1)-sep(1)), ...
                       1+(q-1)*(imsz(2)-sep(2)):imsz(2)+(q-1)*(imsz(2)-sep(2))) ... 
                       +obj(:,:,cnt);
            cnt = cnt+1;
        end
    end

else
    out = zeros(finalsz(1), finalsz(2), imsz(3));
    
    for k = 1:M(1)+2
        for q = 1:M(2)+2
            if(q > 1)
                tmp = repmat(blender(:,:,2),[1, 1, imsz(3)]);
                obj(:,:,:,cnt) = tmp .* obj(:,:,:,cnt);
            end
            if (q < (M(2)+2))
                tmp = repmat(blender(:,:,1),[1, 1, imsz(3)]);
                obj(:,:,:,cnt) = tmp .* obj(:,:,:,cnt);
            end
            if (k < (M(1)+2))
                tmp = repmat(blender(:,:,3),[1, 1, imsz(3)]);
                obj(:,:,:,cnt) = tmp .* obj(:,:,:,cnt);
            end
            if (k > 1)
                tmp = repmat(blender(:,:,4),[1, 1, imsz(3)]);
                obj(:,:,:,cnt) = tmp .* obj(:,:,:,cnt);
            end

            out(1+(k-1)*(imsz(1)-sep(1)):imsz(1)+(k-1)*(imsz(1)-sep(1)), ...
                       1+(q-1)*(imsz(2) - sep(2)):imsz(2) + (q-1)*(imsz(2)-sep(2)), :) ...
                     = out(1+(k-1)*(imsz(1)-sep(1)):imsz(1)+(k-1)*(imsz(1)-sep(1)), ...
                                1+(q-1)*(imsz(2)-sep(2)):imsz(2)+(q-1)*(imsz(2)-sep(2)),:) ... 
                       +obj(:, :, :, cnt);
            cnt = cnt+1;
        end
    end
end  % end of 2D vs 3D toggle

end