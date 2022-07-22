%
% make complex matrix(2byte integer) for unwrapping from 4 phase images
% 2byte-integer(short) real-part, and 2byte-interger(short) imaginary-part
% 29 Aug. 2002  Hidena Iwai
% wph :   wrapping phase, Output file from Labview (int16)
% uph : unwrapping phase
%
function Y = unwrap2(X)
MX = size(X); Y=[];
fname1 = tempname;
fname2 = tempname;
fid = fopen(fname1, 'wb');
fwrite(fid, X, 'float');
fclose(fid);
disp(cmd);
%cmd = sprintf('goldstein -input %s -format float -output %s -xsize %d -ysize %d', fname1, fname2, MX(1), MX(2))
[s, w] = dos(cmd);
fid = fopen(fname2, 'rb');Y = fread(fid, MX, 'float');fclose(fid);
delete fname1;delete fname2;

