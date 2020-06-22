function [ImageData] = reconstructionExpt(ImageDataIn,Settings,Constants,FHandle)
%ImageData = reconstructionExpt(ImageDataIn,Settings,Constants)
%Wrapper function. Run Peter Schauss' reconstruction function from my
%real-time analysis code.


%convert from DatasetClass members used to identify datasets to cell arrays
%{Year,Month,Day,DataSetIndex,FileIndex,PictureIndex}
Dset = ImageDataIn.Dataset;
% try
%     Dset.PictureIndex = Settings.ActivePictures;
% catch
%     Dset.PictureIndex = 2;
% end
if isempty(Settings.ActivePictures)
    Settings.ActivePictures = 2;
end
dataset = {Dset.Year, Dset.Month, Dset.Day, Dset.DataSetIndex, Dset.FileIndex, Settings.ActivePictures};
imageDataArray = [];

%%%Peter Schauss code
userpath('clear') % reset user path
addpath('\\128.112.86.75\lithium\Data analysis\PeterS\reconstruction');
addpath('\\128.112.86.75\lithium\Data analysis\PeterS\reconstruction\misc');
addpath('\\128.112.86.75\lithium\Data analysis\PeterS\reconstruction\find_radii');
% if (~isempty(mfilename('fullpath')))
%     cd(fileparts(mfilename('fullpath')));
% end
% addpath(genpath('directory')) % add subfolders to user path

% addpath 'jlab';
%#ok<*NBRAK>

only_fast_preview = 0;
show_rounded_occ = 0;
showImageWithLattice = 0;
debugPhase = 0;


%dataset = {2016 8 22 3 1 2};


psfType = 10; % type 2: double gaussian; type 10: interpolated
maxIterations = 20;
minStepsize = 0.03;
borderPx = 5;


% normedData=load('Z:\singleatoms\data\_2010\normedEtalonData.dat');
% normedData=[]; % suppress normedData to reconstruct images with sizes other than 512x512 !
%k = 2
resultsFilename = 'results_manual.txt';
fid = fopen(resultsFilename, 'a');
fprintf (fid, imageData.getResultHeader());
% fprintf (fid, 'nr datasetName nAtoms nCounts centerX centerY phase1 phase1Err phase2 phase2Err ');
% fprintf (fid, 'ampl1 ampl1Err ampl2 ampl2Err sigma1 sigma1Err sigma2 sigma2Err\n');
datasetNo = 1;

diary('reconstr_manual.log');
% diary on;


[filenames,picsel] = getFilenames(dataset);
%filenames
%pause(0.5)


numPics = 0;
for l = 1:numel(filenames)
    images = readimg(filenames{l});
    %     if l == 1
    %         picsel = {size(images,3)};
    %     end
    directoryName = [fileparts(filenames{l}) '/'];
    fprintf('read %d pictures %d x %d\n',size(images,3),size(images,1),size(images,2))
    
    
    
    %close all
    
    angles = zeros(1,2);
    [angles(1), angles(2), latticeConstant, date] = lattice_angle_history(dataset{l,1},dataset{l,2},dataset{l,3});
    fprintf('using lattice angles from %s\n',date);
    
    for n=1:numel(picsel{l})
        %close all;
        name = sprintf ('%04d_%02d_%02d_%03d_%03d_%03d', dataset{l,1}, dataset{l,2},dataset{l,3},dataset{l,4},dataset{l,5}, picsel{l}(n));
        if (n > length(images))
            fprintf(['image ' name ' not found\n']);
            break;
        end
        
        %         if (all(images{n} == 0))
        %             fprintf('image %s is empty\n',name);
        %             continue;
        %         end
        
        fprintf ('%s (picture %d)\n', name, picsel{l}(n));
        
        
        
        %         centerF = imageData.findCenterOfMass (images{n}, 200);
        %         center = round(centerF);
        %         center = [220 280];
        im = images(1+borderPx:end-borderPx,1+borderPx:end-borderPx,picsel{l}(n)); % TODO kick out weird border data
        %         figure
        %         imagesc(im);
        %         figure
        %         imagesc(imfilter(im,fspecial('average',[2 2])));
        %         im = imfilter(im,fspecial('average',[2 2]));
        %         return;
        if (false) % find noise region
            figure
            test1 = circshift(im,[1 0]);
            test2 = circshift(im,[0 1]);
            test3 = medfilt2(im-test1,[2,2]).^2+medfilt2(im-test2,[2,2]).^2;
            test4 = test3-median(median(test3));
            imagesc(test4);
            subplot(2,2,4);
            test5 = 1.0.*(medfilt2(test4,[2 2])>mean(test4(:))+2*std(test4(:)));
            noiseregion = imfilter(test5,fspecial('disk',15))>0;
            im(noiseregion) = 0;
            %         binsize = 2;
            %         im = imresize(im, binsize, 'box');
        end
        occ = [];
        
        %         load('dummy_theta_-44_44_phi_0p7_0p2_1_new.mat');
        %         angles = [-44 44];
        
        %         load('Z:/singleatoms/Worksheets_Scripts/matlab/peter/imageReconstruction/tests/phasefittest/dummy_theta_-45_45_phi_0p3_0p3.mat');
        %         angles = [-45 45];
        %         load('reconstructed.mat','imgrec'); % use fake image
        
        %         load('generated_simulated_images/dummy1.mat','imgrec','occ'); % use fake image
        %         load('generated_simulated_images/dummy_nice.mat','imgrec','occ'); % use fake image
        %         load('generated_simulated_images/dummy_good1.mat','imgrec','occ'); % use fake image
        %         load('generated_simulated_images/dummy_bad2.mat','imgrec','occ'); % use fake image
        %         im = imgrec;
        
        %         figure
        %         imagesc(occ);
        %         return;
        %         myImage = im(center(2)-h:center(2)+h-1, center(1)-h:center(1)+h-1);
        
        %angles = 90 - angles;
        
        
        img = imageData ();
        img.showFinalRadiusFit = false;
        img.debugPhase = debugPhase;
        binning = 1;
        img.assignImage (im, [], [directoryName name],binning);
        
        addrCoordsFile = sprintf ('%s%s_4656_subtr_%03d_addrPoints.txt',directoryName,dataset{l,2},picsel{l}(n));
        if (exist(addrCoordsFile,'file') == 2)
            addrCoords = load(addrCoordsFile);
            img.setAddressedCoords(addrCoords);
        end
        
        
        % find single atoms
        img.findSingleAtomsInitial();
        
        % use found single atoms for PSF fitting
        [amplF, amplErr, sigmaF, sigmaErr, exitflag] = img.fitPSF();
        if (exitflag <= 0)
            if (exitflag == 0)
                fprintf('no isolated atoms found for psf fitting\n');
                %                 figure(1001);
                %                 imagesc(im);
                %                 axis square;
            else
                fprintf('psf fitting failed\n');
            end
        end
        fprintf ('ampl1 = %1.3f +- %1.3f   (%1.3f) \nampl2 = %1.3f +- %1.3f  (%1.3f)\nfullamp: %.1f+-%.1f\n\nsigma1 = %1.3f +- %1.3f\nsigma2 = %1.3f +- %1.3f\n',...
            amplF(1), amplErr(1), amplF(1) / sum(amplF),  amplF(2), amplErr(2), amplF(2) / sum(amplF),sum(amplF),sum(abs(amplErr)),  sigmaF(1), sigmaErr(1), sigmaF(2), sigmaErr(2));
        %         sigma = [2.02 6];
        %         sigma = [3.0 10];
        sigma = sigmaF;
        %         counts = amplF;
        ampl = [95 5];
        ampl = [90 10];
        %         ampl = [92 8];
        %         ampl = [87.4 12.6];
        ampl2 = ampl /  sum (ampl);
        ampl2 = amplF/sum(amplF);
        
        singleAtomCount = sum(amplF);
        
        
        img.setGridParameters (latticeConstant, angles);
        
        maxSigmaDeviationForPhaseFit = 0.4;
        % fit phase to found single atoms.
        [phases, stdDev] = img.phaseDifferenceFromSingleAtoms  (singleAtomCount/2, sigmaF(1)*(1-maxSigmaDeviationForPhaseFit),sigmaF(1)*(1+maxSigmaDeviationForPhaseFit));%1.5, 2.6);
        if (any(isnan(phases)))
            %             [phases, stdDev, guessedSingleAtomCount] = find_phase_from_cuts(im,angles(1),angles(2),mean(latticeConstant),'showCutFits',false,'showFinalResults',false,'debug',false);
            %             fprintf('use phase from cuts algorithm to determine phases without single atoms\n');
            %             singleAtomCount = guessedSingleAtomCount;
            fprintf('###### find_phase_from_cuts not working yet - phases random ######\n');
            phases = [0 0];
            stdDev = [3 3];
            %             singleAtomCount = 50;
        end
        if (any(isnan(sigma)))
            sigma = [1.71 1];%[2.26 1]; %default
        end
        if (any(isnan(ampl2)))
            ampl2 = [1 0]; % default
        end
        %         sigma = [1.71 1];%[1.65 1]; %default
        %         ampl2 = [1 0]; % default
        %psfType = 10;2; % dbl gauss
        
        img.setLatticePhase(phases,stdDev,1,showImageWithLattice);
        img.setPSFparameters (sigma, ampl2, psfType, 14);   % type 2 double gaussian
        
        img.reconstructImage (minStepsize, singleAtomCount, maxIterations, 0);
        %         img.reconstructImage (minStepsize, singleAtomCount, maxIterations, 0);
        %         if (showImageWithLattice)
        %         	img.showImageWithLattice();
        %
        %             figure
        %             imagesc(img.occupations);
        %             title('occupation grid');
        %             axis square tight;
        %             colorbar;
        %             imgocc = img.occupations;
        %             save('occupations_raw.mat','imgocc');
        %             imgrec = img.imageReconstructed;
        %             save('reconstructed.mat','imgrec');
        %         end
        
        
        
        
        if (~isempty(occ))
            if (any(size(img.occupations)~=size(occ)))
                fprintf('cannot compare to real occupation of simulated image: size mismatch\n');
            else
                
                
                if (size(img.occupations) ~=  size(occ))
                    fprintf('error: cannot compare reconstructed occupation to orignal one: size mismatch\n');
                    size(img.occupations)
                    size(occ)
                end
                %             occSmall = occ(1:(size(occ,1)-1),1:(size(occ,2)-1));
                
                nAtoms = sum(sum(img.occupationsRounded));
                
                subplot(2,2,1);
                imagesc(img.occupations-occ);
                colorbar
                title('difference between original and reconstructed');
                axis square;
                subplot(2,2,2);
                imagesc(img.occupations);
                colorbar
                title('reconstructed occupations');
                axis square;
                subplot(2,2,3);
                imagesc(occ);
                colorbar
                title('real occupations');
                axis square;
                % clamp original occ for comparison to min/max
                occClamped = occ;
                occClamped(occClamped>0&occClamped<img.amplHistCutoff) = 0;
                
                accError = sum(sum(abs(img.occupations-occ)));
                fprintf('accumulated error: %.1f = %.2f %%\n',accError,100*accError/nAtoms);
                accError = sum(sum(abs(img.occupations-occClamped)));
                fprintf('accumulated error (original clamped): %.1f = %.2f %%\n',accError,100*accError/nAtoms);
                
                atomError = sum(sum(abs(round(img.occupations)-round(occ))));
                fprintf('atom error (round): %.1f = %.2f %%\n',atomError,100*atomError/nAtoms);
                atomError = sum(sum(abs(round(img.occupations)-round(occClamped))));
                fprintf('atom error (round, original clamped): %.1f = %.2f %%\n',atomError,100*atomError/nAtoms);
                
                atomError = sum(sum(abs((round(img.occupations)>0)-(round(occClamped)>0))));
                fprintf('atom error (round, original clamped, 0/1): %.1f = %.2f %%\n',atomError,100*atomError/nAtoms);
                
                %             atomError = sum(sum(abs(~~ceil(img.occupations)-~~ceil(occ))));
                %             fprintf('atom error (ceil mod 2): %.1f = %.2f %%\n',atomError,100*atomError/nAtoms);
                
            end
        end
        
        % with all atoms (<minOcc only shown, but not included)
        saveImageBoolean = 0;
        
        % with all atoms
        %         img.showResults (img.occupations,300, true, saveImageBoolean);
        
        if (~only_fast_preview)
            saveImageBoolean=1;
        end
        % show with small atoms removed
        h = img.showResults (img.occupationsCutoff,300, true, saveImageBoolean);% default 300
        %some inconvenient way to put the results into my figure, instead
        %of the figure produced by the imageData object.
        warning('off','all');
        copyobj(allchild(h),FHandle);
        warning('on','all');
        close(h);
        
        %         fprintf('corrected meanCounts: %.0f (mean relative fitted signal was: %.2f)\n',img.meanCounts,img.meanRelativeSignal);
        
        saveImageBoolean=0;
        % show rounded occupation
        if (show_rounded_occ)
            img.showResults (img.occupationsRounded,300, false, saveImageBoolean);
        end
        
        if (~only_fast_preview)
            saveImageBoolean=1;
            img.showAddressingResults(200, saveImageBoolean);
        end
        
        % save raw image to summary folder
        rawPictureDirectoryName = 'raw_data_summary';
        dayDirname = fileparts(fileparts(directoryName));
        destinationFolder = [dayDirname '/' rawPictureDirectoryName '/'];
        if (exist(destinationFolder,'dir') ~= 7)
            mkdir(dayDirname,rawPictureDirectoryName)
        end
        img.saveRawImage(destinationFolder,300);
        
        
        resultLine = img.getResultLine(datasetNo, dataset{l,1}, dataset{l,2}, n);
        fprintf (fid, strrep(resultLine, '%', '%%')); % comment % signs.
        
        for i=1:length(img.messages)
            messageStr = fprintf('#### %s\n',img.messages{i});
        end
        
        datasetNo = datasetNo + 1;
        imageDataArray = cat(2,imageDataArray,img);
    end
    
end

fclose(fid);

diary off;

ImageData = ImageDataIn;
ImageData.Settings = Settings;
ImageData.ExtraDataStructure = imageDataArray;

if length(imageDataArray)>1
    Anum = zeros(1,length(imageDataArray));
    PSFSigmaPix = zeros(1,length(imageDataArray));
    Diff = zeros(1,length(imageDataArray)-1);
    Lost = zeros(1,length(imageDataArray)-1);
    Hopping = zeros(1,length(imageDataArray)-1);
    
    for ii = 1:length(imageDataArray)
        ImgCurrent = imageDataArray(ii);
        OccCurrent = ImgCurrent.occupationsRounded;
        OccCurrent(OccCurrent>1)=1;
        
        Anum(ii) = sum(sum(OccCurrent));
        if ~isempty(ImgCurrent.PSF)
            PSFSigmaPix(ii) = ImgCurrent.PSF.sigmaEstimate;
        else
            PSFSigmaPix(ii) = 0;
        end
        if ii~=1
            
            
            ImgPrevious = imageDataArray(ii-1);
            OccPrevious = ImgPrevious.occupationsRounded;
            OccPrevious(OccPrevious>1)=1;
            
            try
                DiffArray = xor(OccCurrent,OccPrevious);
            catch
                DiffArray = 0;
            end
            Diff(ii-1) = sum(sum(DiffArray))/sum(sum(OccPrevious));
            Lost(ii-1) = 1 - sum(sum(OccCurrent))/sum(sum(OccPrevious));
            Hopping(ii-1) = Diff(ii-1)-Lost(ii-1);
        end
    end
    
    DiffP = 100*mean(Diff);
    LostP = 100*mean(Lost);
    HoppingP = 100*mean(Hopping);
    
    AnumAvg = mean(Anum);
    PSFSigmaPixAvg = mean(PSFSigmaPix(PSFSigmaPix~=0));
    
    fprintf('Occupation change = %0.2f%%\n',DiffP);
    fprintf('Loss rate = %0.2f%%\n',LostP);
    fprintf('Inferred hopping rate = %0.2f%%\n',HoppingP);
else
        AnumAvg = sum(sum(imageDataArray.occupationsRounded));
        PSFSigmaPixAvg = imageDataArray.PSF.sigmaEstimate;
        DiffP = 0;
        LostP = 0;
        HoppingP = 0;
    
end

fprintf('Atom Number = %0.2f\n\n',AnumAvg);
% fprintf('PSF Sigma = %0.2f pix\n',PSFSigmaPixAvg); %seems like this is
% not right...
% PSFSigmaSIAvg = PSFSigmaPixAvg*Settings.PixelSize*Settings.HardwareBinSizeH/Settings.Magnification;
% fprintf('PSF Sigma = %0.2f um\n',PSFSigmaSIAvg*1e6);

ImageData.Vals = cat(2,ImageData.Vals,[DiffP,LostP,HoppingP,AnumAvg]);
ImageData.Keys = cat(2,ImageData.Keys,{'Difference in occupation rate','Loss rate','Hopping rate','AnumAvg'});

end

