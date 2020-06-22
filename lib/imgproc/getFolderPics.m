function [images] = getFolderPics(folder_path, mode, nsequence)
% Load all images in a given folder, and average sequences of images if
% desired. If image files contain stacks of 3 images, automatically take OD
% image.
%
% FolderPath: path to folder storing image data
%
% Mode: allowed values are
%   'Single', which returns averages of individual pictures,
%   
%   'Normal', which returns averages of OD pictures if number of layers is
%   compatible with OD
%   
%   'Sequential':  assumes you've taken repeats of a sequences of nsequence
%   pictures, and averages pictures {1, NSequence+1,2*NSequence+1,...} then
%   separately {2,NSequence+2,2*NSequence+2,...}
%
%   'Raw': return all pictures, no averaging or od taking.
%
% NSequence: period of repetition of the sequence taken
%
% If an image file contains multiple images, these will be returned as an
% extra dimension in the array. e.g. if there are two images per file, Avg
% will be an array with size [Ny,Nx,Nimg] and Images will be an array with
% size [Ny,Nx,Nimg,Nfiles].
%
% TODO: average before taking OD image, instead of after.
allowed_modes = {'Normal', 'Sequential', 'Single', 'Raw'}; 

if ~exist('mode', 'var')
    mode = 'Normal';
end

% check that mode is an allowed mode
fn = @(c1, c2) strcmp(c1, c2);
if ~any(cellfun(@(c) fn(mode, c), allowed_modes))
    error('Mode was %s, which is not an allowed value', mode);
end

% warn about incompatible argument
if ~strcmp(mode, 'Sequential') && exist('nsequence', 'var')
    warning('Mode is not sequential, but NSequence argument is present. NSequence argument will be ignored');
end

start_path = pwd;
cd(folder_path)
files = [dir('*.aia'); dir('*.fits')];

nfiles = length(files);
first_path = fullfile(folder_path, files(1).name);
first_img = readimg(first_path);
image_err_positions = [];

if size(first_img, 3) == 1 && strcmp(mode, 'Normal')
    %If first image has three pictures, assume all are single OD pics.
    [vsize, hsize] = size(first_img);
    images = zeros(vsize, hsize, nfiles);
    images(:, :, 1) = first_img;
    for ii = 2 : 1 : length(files)
        try
            temp_path = fullfile(folder_path, files(ii).name);
        	temp_img = readimg(temp_path);
            images(:, :, ii) = temp_img;
        catch
            images(:, :, ii) = zeros(vsize, hsize);
            cat(2, image_err_positions, ii);
        end
    end
    
    [Ny, Nx, Nimg] = size(images);

elseif size(first_img, 3) == 3 && strcmp(mode, 'Normal')
    %If first image has three pictures, assume all are single OD pics.
    FirstOD = getOD(first_img);
    [vsize,hsize] = size(FirstOD);
    images = zeros(vsize, hsize, nfiles);
    images(:, :, 1) = FirstOD;
    for ii = 2 : 1 : length(files)
        try
            temp_path=fullfile(folder_path,files(ii).name);
        	temp_img = readimg(temp_path);
            images(:, :, ii) = getOD(temp_img);
        catch
            images(:, :, ii) = zeros(vsize,hsize);
            cat(2,image_err_positions,ii);
        end
    end
    
    avg = sum(images, 3) / size(images, 3);
    [Ny, Nx, Nimg] = size(images);
    images = avg;
    
elseif size(first_img, 3) == 4 && strcmp(mode, 'Normal')
    %If first image has three pictures, assume all are single OD pics.
    FirstOD = getOD(first_img);
    [vsize,hsize] = size(FirstOD);
    images = zeros(vsize,hsize,nfiles);
    images(:,:,1) = FirstOD;
    for ii=2:1:length(files)
        try
            temp_path=fullfile(folder_path,files(ii).name);
        	temp_img = readimg(temp_path);
            images(:,:,ii)=getOD(temp_img(:,:,[2,3,4]));
        catch
            images(:,:,ii) = zeros(vsize,hsize);
            cat(2,image_err_positions,ii);
        end
    end
    
    avg = sum(images,3)/size(images,3);
    [Ny,Nx,Nimg] = size(images);
    imavg = avg;
%     CascadedImage = reshape(Images,[Ny,Nx*Nimg]);
elseif size(first_img,3) == 5 && strcmp(mode,'Normal')
    %assume double shots. Returns average array of size [Ny,Nx,2] and all
    %images as array of size [Ny,Nx,2,Nimgs]
    FirstOD1 = getOD(first_img(:,:,[2,4,5]));
    FirstOD2 = getOD(first_img(:,:,[3,4,5]));
    [vsize,hsize] = size(FirstOD1);
    Images1 = zeros(vsize,hsize,nfiles);
    Images2 = zeros(vsize,hsize,nfiles);
    Images1(:,:,1) = FirstOD1;
    Images2(:,:,1) = FirstOD2;
    for ii=2:1:length(files)
        try
        temp_path=fullfile(folder_path,files(ii).name);
        temp_img = readimg(temp_path);
        Images1(:,:,ii)=getOD(temp_img(:,:,[2,4,5]));
        Images2(:,:,ii) = getOD(temp_img(:,:,[3,4,5]));
        catch
            Images1(:,:,ii) = zeros(vsize,hsize);
            Images2(:,:,ii) = zeros(vsize,hsize);
            cat(2,image_err_positions,ii);
        end
    end
    
    images = cat(4,Images1,Images2);
    Avg1 = sum(Images1,3)/size(Images1,3);
    Avg2 = sum(Images2,3)/size(Images2,3);
    avg = cat(3,Avg1,Avg2);
    [Ny,Nx,Nimg] = size(Images1);
    images = avg;
%     CascadedImage1 = reshape(Images1,[Ny,Nx*Nimg]);
%     CascadedImage2 = reshape(Images2,[Ny,Nx*Nimg]);
%     CascadedImage = cat(3,CascadedImage1, CascadedImage2);
    
elseif size(first_img, 3) == 6 && strcmp(mode, 'Normal')
    %assume double shots. Returns average array of size [Ny,Nx,2] and all
    %images as array of size [Ny,Nx,2,Nimgs]
    FirstOD1 = getOD(first_img(:, :, [2,4,6]));
    FirstOD2 = getOD(first_img(:, :, [3,5,6]));
    [vsize,hsize] = size(FirstOD1);
    Images1 = zeros(vsize,hsize,nfiles);
    Images2 = zeros(vsize,hsize,nfiles);
    Images1(:,:,1) = FirstOD1;
    Images2(:,:,1) = FirstOD2;
    %for ii=4:2:length(Files)
    for ii = 2:1:length(files)
        try
                temp_path=fullfile(folder_path,files(ii).name);
                temp_img = readimg(temp_path);
                Images1(:,:,ii) = getOD(temp_img(:,:,[2,4,6]));
                Images2(:,:,ii) = getOD(temp_img(:,:,[3,5,6]));
        catch
            Images1(:,:,ii) = zeros(vsize,hsize);
            Images2(:,:,ii) = zeros(vsize,hsize);
            cat(2,image_err_positions,ii);
        end
    end
    
    images = cat(4, Images1, Images2);
    Avg1 = sum(Images1,3)/size(Images1,3);
    Avg2 = sum(Images2,3)/size(Images2,3);
    avg = cat(3,Avg1,Avg2);
    [Ny,Nx,Nimg] = size(Images1);
    images = avg;

elseif size(first_img,3) == 7 && strcmp(mode,'Normal')
        %assume double shots. Returns average array of size [Ny,Nx,2] and all
    %images as array of size [Ny,Nx,2,Nimgs]
    FirstOD1 = getOD(first_img(:,:,[2,5,7]));
    FirstOD2 = getOD(first_img(:,:,[3,6,7]));
    [vsize,hsize] = size(FirstOD1);
    Images1 = zeros(vsize,hsize,nfiles);
    Images2 = zeros(vsize,hsize,nfiles);
    Images1(:,:,1) = FirstOD1;
    Images2(:,:,1) = FirstOD2;
    %for ii=4:2:length(Files)
    for ii = 2:1:length(files)
        try
                temp_path=fullfile(folder_path,files(ii).name);
                temp_img = readimg(temp_path);
                Images1(:,:,ii)=getOD(temp_img(:,:,[2,5,7]));
                Images2(:,:,ii) = getOD(temp_img(:,:,[3,6,7]));
        catch
            Images1(:,:,ii) = zeros(vsize,hsize);
            Images2(:,:,ii) = zeros(vsize,hsize);
            cat(2,image_err_positions,ii);
        end
    end
    
    images = cat(4,Images1,Images2);
    Avg1 = sum(Images1,3)/size(Images1,3);
    Avg2 = sum(Images2,3)/size(Images2,3);
    avg = cat(3, Avg1, Avg2);
    images = avg;
    
elseif strcmp(mode, 'Single')%||size(FirstImgs,3)~=3||size(FirstImgs,3)~=5||size(FirstImgs,3)~=6
    %returns average of each individual frame as an array of size
    %[Ny,Nx,Nimg]. Returns all images as array of size [Ny,Nx,Nimg,Nfiles]
    [vsize,hsize,layers] = size(first_img);
    images = zeros(vsize,hsize,layers,nfiles);
    images(:,:,:,1) = first_img;
    %for ii = 4:2:length(Files)
    for ii = 2:1:length(files)
        try
            temp_path=fullfile(folder_path,files(ii).name);
            temp_img = readimg(temp_path);
            images(:,:,:,ii) = temp_img;
        catch
            images(:,:,:,ii) = zeros(vsize,hsize,layers);
            cat(2,image_err_positions,ii);
        end
    end
    avg = sum(images, 4) / size(images,4);
    %CascadedImage = 0;
    images = avg;
    
elseif strcmp(mode, 'Sequential')
    % 
    [vsize, hsize, layers] = size(first_img);
    images = zeros(vsize, hsize, layers, nfiles);
    images(:, :, :, 1) = first_img;
    
    current_seq = 1;
    for ii = 2 : 1 : length(files)
        current_seq = mod(current_seq + 1, nsequence);
        try
            temp_path = fullfile(folder_path, files(ii).name);
            temp_img = readimg(temp_path);
            images(:, :, :, ii) = temp_img;
        catch
            images(:, :, :, ii) = zeros(vsize, hsize, layers);
            cat(2, image_err_positions, ii);
        end
    end
    
    n_partial_sequences = ceil(length(files) / nsequence);
    seq_counter = kron(ones(1, n_partial_sequences), 1:nsequence);
    seq_counter = seq_counter(1:length(files));
    avg = zeros(vsize, hsize, layers, nsequence);
    
    for jj = 1 : nsequence
        current_seq_images = images(:, :, :, seq_counter == jj);
        avg(:, :, :, jj) = sum(current_seq_images, 4) / size(current_seq_images, 4);
    end
    
    images = avg;
       
elseif strcmp(mode, 'Raw')
    % no averaging or od image taking. Just raw images.
    
    [vsize, hsize, layers] = size(first_img);
    images = zeros(vsize, hsize, layers, nfiles);
    
    for ii = 2 : 1 : length(files)
        try
            temp_path = fullfile(folder_path, files(ii).name);
            temp_img = readimg(temp_path);
            images(:, :, :, ii) = temp_img;
        catch
            images(:, :, :, ii) = zeros(vsize, hsize, layers);
            cat(2, image_err_positions, ii);
        end
    end
    
else
    error('incorrect mode argument.');
        
end
cd(start_path);

end

