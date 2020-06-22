classdef ImageDataClass
    %Settings class which defines relevant experimental parameters.
    
    properties
        
        Dataset = DatasetClass;
        Settings = SettingsClass;
        
        %other information
        FilePath
        TimeStamp
        Vals
        Keys        
        
        %coordinate information
        ROIMask
        XcoordsOriginal
        YcoordsOriginal
        XcoordsProc
        YcoordsProc
        
        %images and fits
        OriginalImages
        ProcessedImages
        Fits
        
        %for other uses...
        ExtraDataStructure
        ExtraDataStructure2
        ExtraDataStructure3      
       
    end
    
    methods    
    end

end
