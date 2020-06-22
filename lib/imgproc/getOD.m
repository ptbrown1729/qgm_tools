function [OD] = getOD(images) 
% Get absorption image from array of images for several different formats
% for stacks of images.
% 
% Most common format: first layer is atoms, second layer is beam, third layer is dark.
%
% Other common sizes return a stack of OD images.
%Change function so that 1. can accept a stack of [Ny,Nx,3,NImgs] and
%process correctly.
%For different cases, instead of rewriting the code, reformulate so that
%array size is [Ny,Nx,3,NImgs]

%elseif size(Images,3)==5 && size(Images,4)==1
%Images = cat(4,Images(:,:,[2,4,5]),Images(:,:,[3,4,5]));
%elseif size(Images,3)==6 && size(Images,4)==1
%Images = cat(4,Images(:,:,[2,4,6]),Images(:,:,[3,5,6]));
%elseif size(Images,3)==7 && size(Images,4)==1
%Images = cat(4,Images(:,:,[2,5,7]),Images(:,:,[3,6,7]));

if ndims(images) == 4
    % loop over image
    od_first = getOD(images(:, :, :, 1));
    OD = zeros([size(od_first, 1), size(od_first, 2), size(images, 4)]);
    OD(:, :, 1) = od_first;
    
    for ii = 2 : size(images, 4)
        OD(:, :, ii) = getOD(images(:, :, :, ii));
    end
    
else

    % single image
    if (size(images, 3) == 3)
        atoms = images(:,:,1);
        beam = images(:,:,2);
        dark = images(:,:,3);

        divided = (atoms-dark)./(beam-dark);

        divided(isnan(divided))=1;
        divided(divided<=0)=1;
        divided(isinf(divided))=1;

        OD = -log(divided);

    elseif (size(images, 3) == 4)
        atoms = images(:,:,2);
        beam = images(:,:,3);
        dark = images(:,:,4);

        divided = (atoms-dark)./(beam-dark);

        divided(isnan(divided))=1;
        divided(divided<=0)=1;
        divided(isinf(divided))=1;

        OD = -log(divided);

    elseif (size(images,3) == 5)
        atoms1 = images(:,:,2);
        atoms2 = images(:,:,3);
        beam1 = images(:,:,4);
        dark = images(:,:,5);

        divided1 = (atoms1-dark)./(beam1-dark);
        divided2 = (atoms2-dark)./(beam1-dark);

        divided1(isnan(divided1)|divided1<=0|isinf(divided1))=1;
        divided2(isnan(divided2)|divided2<=0|isinf(divided2))=1;


        OD = cat(3,-log(divided1),-log(divided2));

    elseif (size(images,3) == 6)
        atoms1 = images(:,:,2);
        atoms2 = images(:,:,3);
        beam1 = images(:,:,4);
        beam2 = images(:,:,5);
        dark = images(:,:,6);

        divided1 = (atoms1-dark)./(beam1-dark);
        divided2 = (atoms2-dark)./(beam2-dark);

        divided1(isnan(divided1)|divided1<=0|isinf(divided1))=1;
        divided2(isnan(divided2)|divided2<=0|isinf(divided2))=1;

        OD = cat(3,-log(divided1),-log(divided2));

    elseif (size(images,3) == 7)
        atoms1 = images(:,:,2);
        atoms2 = images(:,:,3);
        beam1 = images(:,:,5);
        beam2 = images(:,:,6);
        dark = images(:,:,7);

        divided1 = (atoms1-dark)./(beam1-dark);
        divided2 = (atoms2-dark)./(beam2-dark);

        divided1(isnan(divided1)|divided1<=0|isinf(divided1))=1;
        divided2(isnan(divided2)|divided2<=0|isinf(divided2))=1;

        OD = cat(3,-log(divided1),-log(divided2));
        %OD = -log(divided1);
        %OD2 = -log(divided2);
    else
        error('error in getOD: not implemented for %d images',size(images,3));
    end

end

end

