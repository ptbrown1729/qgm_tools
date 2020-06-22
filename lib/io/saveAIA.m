function [] = saveAIA(Img, SavePath, Vals, Keys)
%[] = saveAIA(Img, SavePath, Vals, Keys)
%Saves an image as .aia file.
%should finish this...
lay = size(Img, 3);
row = size(Img, 1);
col = size(Img, 2);

try
    fid = fopen(SavePath,'w');
    %unused lines.
    fwrite(fid, 'AIA');
    fwrite(fid, 2, 'uint16'); %unsure what this is for
    %write image dimensions
    fwrite(fid, row, 'uint16');
    fwrite(fid, col, 'uint16');
    fwrite(fid, lay, 'uint16');
    
    %write image
    for ii=1:lay
        for jj=1:row
            fwrite(fid, Img(jj,:,ii), 'uint16');
        end
    end
    
    %write metadata
    if length(Vals) == length(Keys)
        for ii = 1:length(Vals)
            fwrite(fid, Keys{ii});
            fwrite(fid, char(0), '*char'); %strings should be null terminated
            fwrite(fid, Vals(ii), 'double');
        end
    else
        disp('Lengths of Keys and Vals did not match.')
    end
    
    fclose(fid);
    
catch
    fclose(fid);
end

end

