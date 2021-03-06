% MI_Dataset_List = {{2016 10 28 28 1 2};{2016 11 1 31 1 2};{2016 11 1 37 1 2};{2016 10 31 34 1 2};{2016 10 31 23 1 2};{2016 11 1 42 1 2};{2016 10 31 53 1 2}};%;{}};
% SzSz_Blow1_List = {{2016 10 28 24 1 2};{2016 11 1 27 1 2};{2016 11 1 33 1 2};{2016 10 31 30 1 2};{2016 10 31 19 1 2};{2016 11 1 38 1 2};{2016 10 31 49 1 2}};%;{2016 11 11 24 1 2}};
% SzSz_Blow2_List = {{2016 10 28 26 1 2};{2016 11 1 29 1 2};{2016 11 1 35 1 2};{2016 10 31 32 1 2};{2016 10 31 21 1 2};{2016 11 1 40 1 2};{2016 10 31 51 1 2}};%;{2016 11 11 23 1 2}};
% SxSx_Blow1_List = {{2016 10 28 27 1 2};{2016 11 1 30 1 2};{2016 11 1 36 1 2};{2016 10 31 33 1 2};{2016 10 31 22 1 2};{2016 11 1 41 1 2};{2016 10 31 52 1 2}};%;{}};
% SxSx_Blow2_List = {{2016 10 28 25 1 2};{2016 11 1 28 1 2};{2016 11 1 34 1 2};{2016 10 31 31 1 2};{2016 10 31 20 1 2};{2016 11 1 39 1 2};{2016 10 31 50 1 2}};%;{}};

Save = 1;

DateList = {{2016 10 28},{2016 11 1}, {2016 11 1}, {2016 10 31}, {2016 10 31},{2016 11 1},{2016 10 31}};
DSets = {24:28,27:31,33:37,30:34,19:23,38:42,49:53};

BinEdges = sqrt(linspace(0,30^2,20));

if ~exist('FG','var')||~FG.IsInitialized
    FG = NonIntFG(1);
end

for ii = 1:length(DSets)
   set =  CorrDataSetRep(DateList{ii},DSets{ii},[],BinEdges);
   set.NonIntFG = FG;
   set.InitFG = 1;
   if Save
       set.saveToFile();
       set.savesummaryGraphs('./Figs');
   end
end