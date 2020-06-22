classdef SettingsClass
    %Settings class which defines relevant experimental parameters. These
    %are generated by the settings.m file stored in each data set. The
    %constructor to this class takes the path to such a file.
    
    properties
        %%%Folder Properties
        CurrentFolder % maybe these don't really belong in settings class?
        Experiment
        %%%Hardware Properties.
        Camera
        Manufacturer
        PictureType
        ActivePictures
        QuantumEfficiency
        CountsPerPhoton
        ImagingDuration
        PixelSize %TODO: deprecate this in favor of PixelSizeH and PixelSizeV.
        PixelSizeH
        PixelSizeV
        ImagingAxis
        Magnification
        HardwareBinSizeH
        HardwareBinSizeV
        
        %%%Software Properties.
        %binning
        SoftwareBinSizeH
        SoftwareBinSizeV
        %stripes due to readout noise
        RemoveStripes = 0
        StripesOrientation
        StripesRegionLimits
        %OD normalization
        NormalizeOD = 0
        NormRegionLimits
        WindowLimits
        %%%Fitting Properties.
        Fitting = 0
        FitInitialParameters %[135,46,2462,...]. e.g.
        FitFixedParameters %[0,1,0,...]. Boolean for each fit parameter. If 1, will fix to initial parameter.
        FitUpperLims = []
        FitLowerLims = []
        %these ones deprecated...
        UseFitAutoGuess = 0
        GridParams
        GridFixedParams
        AzimuthalAverageInitialParameters
        AzimuthalAverageFixedParameters
        %Image Display Properties.
        %Possible inconsistency...the actual images send to showPlots come
        %from the Experiment analysis file. So have to know something about
        %this.
        
        ColorMap
        %Assume that number of fits is same as number of images.
        % TODO: get rid of this image positioning information. I think that
        % this is too specific to be defined in the settings class. In the 
        % future I would like to change how the figures are handled by 
        % removing the "showPlots" function and generating the desired
        % figure from the experiment function
        ImageMultiplicity %Form [0,1,2,1,...]. Number of times to show image.
        ImageLocations % = [[2,2]; [3,3]; [3,2]; ...]
        ImageLimits % = [[0,1.5]; [0, 0.5]; ...]In order, same as number of images with multiplicity.
        FitMultiplicity %Form [0,1,0,3,...]. Number of times to show image.
        FitLocations % similar to images
        FitLimits
        ResidualBool %similar to images, except no multiplicity option
        ResidualLocations
        ResidualLimits
        
        %Debug mode
        MarkerLocations
        DebugMode = 0
        
    end
    
    methods
        function obj = SettingsClass(FilePath)
            %initialization function. Same functionality as
            %classFromSettingsFile function. Takes a settings.m file as an
            %input argument. Variables in the settings.m file which have
            %the same name as class properties are added to the class.
            if exist('FilePath','var')
                if exist(FilePath, 'file')
                    run(FilePath);
                    obj = SettingsClass;
                    
                    Fields = fieldnames(obj);
                    for ii = 1:length(Fields)
                        if exist(char(Fields(ii)), 'var')
                            FieldStr = char(Fields(ii));
                            eval(cat(2, 'obj.', FieldStr, '=', FieldStr, ';'));
                        end
                    end
                
                end
            end
        end
        
        %These ignore software binning, since our fn preserves the number
        %and hence size of these pixels.
        function EffPixelSizeH = getEffPixelSizeH(obj)
            EffPixelSizeH = [obj.PixelSize] * [obj.HardwareBinSizeH] / [obj.Magnification];
        end
        
        function EffPixelSizeV = getEffPixelSizeV(obj)
            EffPixelSizeV = [obj.PixelSize] * [obj.HardwareBinSizeV] / [obj.Magnification];
        end
        
        function EffPixelArea = getEffPixelArea(obj)
            EffPixelArea = ([obj.PixelSize] / [obj.Magnification])^2 * [obj.HardwareBinSizeV] * [obj.HardwareBinSizeH];
        end
        
    end
    
end

