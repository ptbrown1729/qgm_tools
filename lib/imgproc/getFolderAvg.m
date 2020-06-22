function [Avgs, Nfiles] = getFolderAvg(FolderPath, Mode, NSequence)
%[Avgs,Nfiles] = getFolderAvg(FolderPath,Mode,NSequence)
%Sequential mode assumes you've taken repeats of a sequences of NSequence
%pictures, and averages pictures {1, NSequence+1,2*NSequence+1,...} then
%separately {2,NSequence+2,2*NSequence+2,...}

%Ensure sensible arguments.
if ~exist('Mode','var')
    Mode = 'Sequential';
    NSequence = 1;
end

if strcmp(Mode,'Single')
    Mode = 'Sequential';
end

%Get list of files in directory.
OriginalPath = pwd;
cd(FolderPath)
Files=[dir('*.aia'); dir('*.fits')];
NTotalFiles = length(Files);


%Stores number of files currently in each average.
AveragedImages = zeros(1,NSequence);
%Get first image so know size.
FirstPath = fullfile(FolderPath,Files(1).name);
FirstImg = readimg(FirstPath);
Nfiles = 0;

[Vsize,Hsize,Lay] = size(FirstImg);
if strcmp(Mode,'Sequential')
    Avgs = zeros(Vsize,Hsize,Lay,NSequence);
    NCurrSeq = 1;
    for ii = 1:length(Files)
        %try
            TmpPath=fullfile(FolderPath,Files(ii).name);
            lastwarn(''); 
            IgnoreImg = 0;
            try
                TmpImgs = readimg(TmpPath);                
                %Unfortunately, when fitsread fails, it does not throw and
                %exception, but gives a warning. This code is a workaround
                %causing that warning to give an exception so I can catch
                %it and deal with it.
                if ~isempty(lastwarn)
                    error(warn);
                end
            catch
                IgnoreImg = 1;
                fprintf('Ignored Image %i \n',ii);
            end
            if ~IgnoreImg
                N = AveragedImages(NCurrSeq);
                Avgs(:,:,:,NCurrSeq) = Avgs(:,:,:,NCurrSeq)*N/(N+1)+TmpImgs/(N+1);
                AveragedImages(NCurrSeq) = AveragedImages(NCurrSeq)+1;
                Nfiles = Nfiles+1;
            end
%         catch Exception
%             disp('Caught Exception in getFolderAvg');
%         end
        NCurrSeq = mod(NCurrSeq,NSequence)+1;
    end
    
elseif strcmp(Mode,'OD')
    
end

if Nfiles~=NTotalFiles
    disp('Some files were ignored using getFolderAvg')
end
cd(OriginalPath);

end