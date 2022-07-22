function writeTiff_32bGray(file, img)

    setTag(file,'Photometric',Tiff.Photometric.MinIsBlack);
    setTag(file,'Compression', Tiff.Compression.None);
    setTag(file,'ImageWidth', size(img,2));
    setTag(file,'ImageLength',size(img,1));
    setTag(file,'BitsPerSample', 32);
    setTag(file,'SamplesPerPixel', 1);
    setTag(file,'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    setTag(file,'SampleFormat',Tiff.SampleFormat.IEEEFP);
    write(file, single(img));
    close(file);
end  % EOF