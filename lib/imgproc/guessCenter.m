function [Cx,Cy,Sx,Sy] = guessCenter(M)
%[Cx,Cy,Sx,Sy] = guessCenter(M)

NImgs = size(M,3);

Cx = zeros(1,NImgs);
Cy = zeros(1,NImgs);
Sx = zeros(1,NImgs);
Sy = zeros(1,NImgs);

NBinLarge = 10;

for ii = 1:NImgs
    BinnedM = binImg(M(:,:,ii),NBinLarge,NBinLarge);

    [Max,Index] = max(BinnedM(:));
    [CyMax,CxMax] = ind2sub(size(M),Index);
    
    Fractions = 0:0.05:0.9;
    CxTmp = zeros(1,length(Fractions));
    StdsCxTmp = zeros(1,length(Fractions));
    
    CyTmp = zeros(1,length(Fractions));
    StdsCyTmp = zeros(1,length(Fractions));
    
    SxTmp = zeros(1,length(Fractions));
    SyTmp = zeros(1,length(Fractions));
    
    BinnedMTmp = BinnedM;
    
    %exclude increasing amounts of the picture and get center
    for jj = 1:length(Fractions)
        BinnedMTmp(BinnedM<Fractions(jj)*Max) = 0;
        [CxCurrent,CyCurrent,SxCurrent,SyCurrent] = getCenterMass(BinnedMTmp);
        CxTmp(jj) = CxCurrent; 
        CyTmp(jj) = CyCurrent;
        SxTmp(jj) = SxCurrent;
        SyTmp(jj) = SyCurrent;
    end
    
    %Trying to find point where Cx/Cy guesses become similar.
    %hopefully here have excluded all noise.
    for jj = 1:length(Fractions)
        StdsCxTmp(jj) = std(CxTmp(jj:end));
        StdsCyTmp(jj) = std(CyTmp(jj:end));
    end
    
    CounterArray = 1:length(Fractions);
    GoodIndices = CounterArray(StdsCxTmp<2 & StdsCyTmp<2);
    if ~isempty(GoodIndices)
        I = GoodIndices(1);
    else
        I = ceil(length(Fractions)/2);
    end
    
    Cx(ii) = CxTmp(I); Cy(ii) = CyTmp(I);
    Sx(ii) = SxTmp(I); Sy(ii) = SyTmp(I);
    
%     figure
%     plot(CxTmp)
%     figure
%     plot(StdsCxTmp);
%     figure
%     plot(CyTmp)
%     figure
%     plot(StdsCyTmp);
%     
%     disp(CyGoodGuesses)
%     disp(CxGoodGuesses)
    
    

end

end