classdef DataSummaryV2 < handle
    
    
    properties
        
        SavePath = 'Data3'
        
        Polarization = []
        PolarizationUnc = []
        
        NumNeighbors = 4;
        
        %central densities
        MI_CentralDensity = []
        MI_CentralDensityUnc = []
        Sz_Up_CentralDensity = []
        Sz_Up_CentralDensityUnc = []
        Sz_Down_CentralDensity = []
        Sz_Down_CentralDensityUnc = []
        Sx_Up_CentralDensity = []
        Sx_Up_CentralDensityUnc = []
        Sx_Down_CentralDensity = []
        Sx_Down_CentralDensityUnc = []
        
        %moment correlators
        MI_MaxMomentCorrNN = []
        MI_MaxMomentCorrNNUnc = []
        
        %Spin correlation matrices
        %first azimuthal average bin only...
        Sz_SpinCorrMat = []
        Sz_SpinCorrMatUnc = []
        Sx_SpinCorrMat = []
        Sx_SpinCorrMatUnc = []
        SpinCorrMat_Anisotropy = []
        SpinCorrMat_AnisotropyUnc = []
        
        %NN
        Sz_SpinCorrNN = [];
        Sz_SpinCorrNNUnc = [];
        Sx_SpinCorrNN = [];
        Sx_SpinCorrNNUnc = [];
        AnisotropyNN = [];
        AnisotropyNNUnc = [];
        
        %NNN
        Sz_SpinCorrNNN = [];
        Sz_SpinCorrNNNUnc = [];
        Sx_SpinCorrNNN = [];
        Sx_SpinCorrNNNUnc = [];
        AnisotropyNNN = [];
        AnisotropyNNNUnc = [];
        
        %NNNN
        Sz_SpinCorrNNNN = [];
        Sz_SpinCorrNNNNUnc = [];
        Sx_SpinCorrNNNN = [];
        Sx_SpinCorrNNNNUnc = [];
        AnisotropyNNNN = [];
        AnisotropyNNNNUnc = [];
        
        %Charge/density/moment correlators
        
        
        %azimuthal average parameters
        BinEdges = []
        
        
    end
    
    methods
        function obj = DataSummaryV2(CorrDataSetStack)
            
            if exist('CorrDataSetStack','var')
                NDSets = length(CorrDataSetStack);
                obj.Polarization = zeros(length(NDSets),1);
                
                for ii = 1:NDSets
                    CurrentDataSet = CorrDataSetStack(ii);
                    if ii == 1
                        CorrDims =ndims(CurrentDataSet.SzSz_Corr);
                        obj.BinEdges = CurrentDataSet.BinEdges;
                    else
                        if ~isequal(CurrentDataSet.BinEdges,obj.BinEdges)
                            error('BinEdges changed while looping over dataset stack');
                        end
                    end
                    
                    if ~isempty(CurrentDataSet.GlobalPolarization)
                        obj.Polarization(ii) = CurrentDataSet.GlobalPolarization;
                    else
                        obj.Polarization(ii) = 0;
                    end
                    obj.MI_CentralDensity = cat(1,obj.MI_CentralDensity,CurrentDataSet.MI.Occs_AzAvg(1));
                    obj.MI_CentralDensityUnc = cat(1,obj.MI_CentralDensityUnc,CurrentDataSet.MI.Occs_AzAvgUnc(1));
                    obj.Sz_Up_CentralDensity = cat(1,obj.Sz_Up_CentralDensity,CurrentDataSet.Up_Sz.Occs_AzAvg(1));
                    obj.Sz_Up_CentralDensityUnc = cat(1,obj.Sz_Up_CentralDensityUnc,CurrentDataSet.Up_Sz.Occs_AzAvgUnc(1));
                    obj.Sz_Down_CentralDensity = cat(1,obj.Sz_Down_CentralDensity,CurrentDataSet.Down_Sz.Occs_AzAvg(1));
                    obj.Sz_Down_CentralDensityUnc = cat(1,obj.Sz_Down_CentralDensityUnc,CurrentDataSet.Down_Sz.Occs_AzAvgUnc(1));
                    obj.Sx_Up_CentralDensity = cat(1,obj.Sx_Up_CentralDensity,CurrentDataSet.Up_Sx.Occs_AzAvg(1));
                    obj.Sx_Up_CentralDensityUnc = cat(1,obj.Sx_Up_CentralDensityUnc,CurrentDataSet.Up_Sx.Occs_AzAvgUnc(1));
                    obj.Sx_Down_CentralDensity = cat(1,obj.Sx_Down_CentralDensity,CurrentDataSet.Down_Sx.Occs_AzAvg(1));
                    obj.Sx_Down_CentralDensityUnc = cat(1,obj.Sx_Down_CentralDensityUnc,CurrentDataSet.Down_Sx.Occs_AzAvgUnc(1));
                    
                    [MaxMomentCorrNN,I] = max(squeeze(CurrentDataSet.MI.Density_Corr_AzAvg(5,6,:)));
                    MaxMomentCorrNNUnc = squeeze(CurrentDataSet.MI.Density_Corr_AzAvg(5,6,I));
                    obj.MI_MaxMomentCorrNN = cat(1,obj.MI_MaxMomentCorrNN,MaxMomentCorrNN);
                    obj.MI_MaxMomentCorrNNUnc = cat(1,obj.MI_MaxMomentCorrNNUnc,MaxMomentCorrNNUnc);
                    
                    obj.Sz_SpinCorrMat = cat(CorrDims,obj.Sz_SpinCorrMat,CurrentDataSet.SzSz_Corr(:,:,1));
                    obj.Sz_SpinCorrMatUnc = cat(CorrDims,obj.Sz_SpinCorrMatUnc,CurrentDataSet.SzSz_CorrUnc(:,:,1));
                    obj.Sx_SpinCorrMat = cat(CorrDims,obj.Sx_SpinCorrMat,CurrentDataSet.SxSx_Corr(:,:,1));
                    obj.Sx_SpinCorrMatUnc = cat(CorrDims,obj.Sx_SpinCorrMatUnc,CurrentDataSet.SxSx_CorrUnc(:,:,1));
                    obj.SpinCorrMat_Anisotropy = cat(CorrDims,obj.SpinCorrMat_Anisotropy,CurrentDataSet.Anisotropy(:,:,1));
                    obj.SpinCorrMat_AnisotropyUnc = cat(CorrDims,obj.SpinCorrMat_AnisotropyUnc,CurrentDataSet.AnisotropyUnc(:,:,1));
                    
                    obj.getOtherCorrelators();
                end
            end
            
        end
        
        function append(obj,CurrentDataSet)
            if ~isempty(CurrentDataSet.GlobalPolarization)
                obj.Polarization = cat(2,obj.Polarization,CurrentDataSet.GlobalPolarization);
            else
                obj.Polarization = cat(2,obj.Polarization,0);
            end
            
            obj.MI_CentralDensity = cat(1,obj.MI_CentralDensity,CurrentDataSet.MI.Occs_AzAvg(1));
            obj.MI_CentralDensityUnc = cat(1,obj.MI_CentralDensity,CurrentDataSet.MI.Occs_AzAvgUnc(1));
            obj.Sz_Up_CentralDensity = cat(1,obj.Sz_Up_CentralDensity,CurrentDataSet.Up_Sz.Occs_AzAvg(1));
            obj.Sz_Up_CentralDensityUnc = cat(1,obj.Sz_Up_CentralDensityUnc,CurrentDataSet.Up_Sz.Occs_AzAvgUnc(1));
            obj.Sz_Down_CentralDensity = cat(1,obj.Sz_Down_CentralDensity,CurrentDataSet.Down_Sz.Occs_AzAvg(1));
            obj.Sz_Down_CentralDensityUnc = cat(1,obj.Sz_Down_CentralDensityUnc,CurrentDataSet.Down_Sz.Occs_AzAvgUnc(1));
            obj.Sx_Up_CentralDensity = cat(1,obj.Sx_Up_CentralDensity,CurrentDataSet.Up_Sx.Occs_AzAvg(1));
            obj.Sx_Up_CentralDensityUnc = cat(1,obj.Sx_Up_CentralDensityUnc,CurrentDataSet.Up_Sx.Occs_AzAvgUnc(1));
            obj.Sx_Down_CentralDensity = cat(1,obj.Sx_Down_CentralDensity,CurrentDataSet.Down_Sx.Occs_AzAvg(1));
            obj.Sx_Down_CentralDensityUnc = cat(1,obj.Sx_Down_CentralDensityUnc,CurrentDataSet.Down_Sx.Occs_AzAvgUnc(1));
            
            [MaxMomentCorrNN,I] = max(squeeze(CurrentDataSet.MI.Density_Corr_AzAvg(5,6,:)));
            MaxMomentCorrNNUnc = squeeze(CurrentDataSet.MI.Density_Corr_AzAvgUnc(5,6,I));
            obj.MI_MaxMomentCorrNN = cat(1,obj.MI_MaxMomentCorrNN,MaxMomentCorrNN);
            obj.MI_MaxMomentCorrNNUnc = cat(1,obj.MI_MaxMomentCorrNNUnc,MaxMomentCorrNNUnc);
            
            obj.Sz_SpinCorrMat = cat(3,obj.Sz_SpinCorrMat,CurrentDataSet.SzSz_Corr(:,:,1));
            obj.Sz_SpinCorrMatUnc = cat(3,obj.Sz_SpinCorrMatUnc,CurrentDataSet.SzSz_CorrUnc(:,:,1));
            obj.Sx_SpinCorrMat = cat(3,obj.Sx_SpinCorrMat,CurrentDataSet.SxSx_Corr(:,:,1));
            obj.Sx_SpinCorrMatUnc = cat(3,obj.Sx_SpinCorrMatUnc,CurrentDataSet.SxSx_CorrUnc(:,:,1));
            obj.SpinCorrMat_Anisotropy = cat(3,obj.SpinCorrMat_Anisotropy,CurrentDataSet.Anisotropy(:,:,1));
            obj.SpinCorrMat_AnisotropyUnc = cat(3,obj.SpinCorrMat_AnisotropyUnc,CurrentDataSet.AnisotropyUnc(:,:,1));
            obj.getOtherCorrelators()
            
        end
        
        
        function getOtherCorrelators(obj)
            %NN
            obj.Sz_SpinCorrNN = squeeze(obj.Sz_SpinCorrMat(obj.NumNeighbors+1,obj.NumNeighbors+2,:));
            obj.Sz_SpinCorrNNUnc = squeeze(obj.Sz_SpinCorrMatUnc(obj.NumNeighbors+1,obj.NumNeighbors+2,:));
            obj.Sx_SpinCorrNN = squeeze(obj.Sx_SpinCorrMat(obj.NumNeighbors+1,obj.NumNeighbors+2,:));
            obj.Sx_SpinCorrNNUnc = squeeze(obj.Sx_SpinCorrMatUnc(obj.NumNeighbors+1,obj.NumNeighbors+2,:));
            obj.AnisotropyNN = squeeze(obj.SpinCorrMat_Anisotropy(obj.NumNeighbors+1,obj.NumNeighbors+2,:));
            obj.AnisotropyNNUnc = squeeze(obj.SpinCorrMat_AnisotropyUnc(obj.NumNeighbors+1,obj.NumNeighbors+2,:));
            %NNN
            obj.Sz_SpinCorrNNN = squeeze(obj.Sz_SpinCorrMat(obj.NumNeighbors+2,obj.NumNeighbors+2,:));
            obj.Sz_SpinCorrNNNUnc = squeeze(obj.Sz_SpinCorrMatUnc(obj.NumNeighbors+2,obj.NumNeighbors+2,:));
            obj.Sx_SpinCorrNNN = squeeze(obj.Sx_SpinCorrMat(obj.NumNeighbors+2,obj.NumNeighbors+2,:));
            obj.Sx_SpinCorrNNNUnc = squeeze(obj.Sx_SpinCorrMatUnc(obj.NumNeighbors+2,obj.NumNeighbors+2,:));
            obj.AnisotropyNNN = squeeze(obj.SpinCorrMat_Anisotropy(obj.NumNeighbors+2,obj.NumNeighbors+2,:));
            obj.AnisotropyNNNUnc = squeeze(obj.SpinCorrMat_AnisotropyUnc(obj.NumNeighbors+2,obj.NumNeighbors+2,:));
            %NNNN
            obj.Sz_SpinCorrNNNN = squeeze(obj.Sz_SpinCorrMat(obj.NumNeighbors+1,obj.NumNeighbors+3,:));
            obj.Sz_SpinCorrNNNNUnc = squeeze(obj.Sz_SpinCorrMatUnc(obj.NumNeighbors+1,obj.NumNeighbors+3,:));
            obj.Sx_SpinCorrNNNN = squeeze(obj.Sx_SpinCorrMat(obj.NumNeighbors+1,obj.NumNeighbors+3,:));
            obj.Sx_SpinCorrNNNNUnc = squeeze(obj.Sx_SpinCorrMatUnc(obj.NumNeighbors+1,obj.NumNeighbors+3,:));
            obj.AnisotropyNNNN = squeeze(obj.SpinCorrMat_Anisotropy(obj.NumNeighbors+1,obj.NumNeighbors+3,:));
            obj.AnisotropyNNNNUnc = squeeze(obj.SpinCorrMat_AnisotropyUnc(obj.NumNeighbors+1,obj.NumNeighbors+3,:));
            
            
        end
        
        function showMaxMomentCorr(obj)
            FigName = 'Maximum Moment Correlator Vs. Moment';
            h = figure('name',FigName);
            errorbar(obj.Polarization,obj.MI_MaxMomentCorrNN,obj.MI_MaxMomentCorrNNUnc,'r-o');
            grid on;
        end
        
        function saveMaxMomentCorr(obj)
            %make sure these are the correct shape...
            Pol = obj.Polarization;
            PolUnc = obj.PolarizationUnc;
            if length(Pol)>1 && size(Pol,1)==1
                Pol = transpose(Pol);
            end
            
            if length(PolUnc)>1 && size(PolUnc,1)==1
                PolUnc = transpose(PolUnc);
            end
            
            FName = 'MomentCorrelatorsData.txt';
            FPath = fullfile(obj.SavePath,FName);
            Titles = {'Pol','PolUnc','M(1,0)','M(1,0)Unc'};
            try
                FileID = fopen(FPath,'w');
                for ii = 1:length(Titles)-1
                    fprintf(FileID,'%s ',Titles{ii});
                end
                fprintf(FileID,'%s\n',Titles{end});
                fclose(FileID);
            catch
                try
                    fclose(FileID);
                catch
                end
            end
            FullMat = cat(2,Pol,PolUnc,obj.MI_MaxMomentCorrNN,obj.MI_MaxMomentCorrNNUnc);
            dlmwrite(FPath,FullMat,'-append','delimiter',' ');
        end
        
        %routines for displaying and saving data
        
        function showCorrVsP(obj)
            FigName = 'Correlators Vs. P';
            h = figure('name',FigName);
            
            if length(obj.Polarization)>length(obj.Sx_SpinCorrNN)
                Pol = obj.Polarization(1:length(obj.Sx_SpinCorrNN));
            else
                Pol = obj.Polarization;
            end
            errorbar(Pol,obj.Sz_SpinCorrNN,obj.Sz_SpinCorrNNUnc,'r-o')
            hold on;
            errorbar(Pol,obj.Sx_SpinCorrNN,obj.Sx_SpinCorrNNUnc,'b-o')
            errorbar(Pol,obj.Sz_SpinCorrNNN,obj.Sz_SpinCorrNNNUnc,'r-.')
            errorbar(Pol,obj.Sx_SpinCorrNNN,obj.Sx_SpinCorrNNNUnc,'b-.')
            xlim([0,1])
            grid on;
            hold off;
            legend({'C_z(1,0)','C_x(1,0)','C_z(1,1)','C_x(1,1)'})
            FigName2 = 'Anisotropy Vs. P';
            h = figure('name',FigName2);
            errorbar(Pol,obj.AnisotropyNN,obj.AnisotropyNNUnc,'g-o')
            hold on;
            errorbar(Pol,obj.AnisotropyNNN,obj.AnisotropyNNUnc,'r-o')
            hold off;
            xlim([0,1])
            ylim([-0.05,1.5])
            grid on;
            
        end
        
        function saveNNCorrs(obj)
            Pol = obj.Polarization;
            PolUnc = obj.PolarizationUnc;
            if length(obj.Polarization)>length(obj.Sx_SpinCorrNN)
                Pol = obj.Polarization(1:length(obj.Sx_SpinCorrNN));
                PolUnc = obj.PolarizationUnc(1:length(obj.Sx_SpinCorrNN));
            end
            
            %make sure these are the correct shape...
            if length(Pol)>1 && size(Pol,1)==1
                Pol = transpose(Pol);
            end
            
            if length(PolUnc)>1 && size(PolUnc,1)==1
                PolUnc = transpose(PolUnc);
            end
            
            FName = 'SpinCorrelatorsData.txt';
            FPath = fullfile(obj.SavePath,FName);
            
            Titles = {'Pol','PolUnc','Cz(1,0)','Cz(1,0)Unc','Cx(1,0)','Cx(1,0)Unc','A(1,0)','A(1,0)Unc','Cz(1,1)','Cz(1,1)Unc','Cx(1,1)','Cx(1,1)Unc','A(1,1)','A(1,1)Unc','Cz(2,0)','Cz(2,0)Unc','Cx(2,0)','Cx(2,0)Unc','A(2,0)','A(2,0)Unc'};
            try
                FileID = fopen(FPath,'w');
                for ii = 1:length(Titles)-1
                    fprintf(FileID,'%s ',Titles{ii});
                end
                fprintf(FileID,'%s\n',Titles{end});
                fclose(FileID);
            catch
                try
                    fclose(FileID);
                catch
                end
            end
            FullMat = cat(2,Pol,PolUnc,obj.Sz_SpinCorrNN,obj.Sz_SpinCorrNNUnc,obj.Sx_SpinCorrNN,obj.Sx_SpinCorrNNUnc,obj.AnisotropyNN,obj.AnisotropyNNUnc,obj.Sz_SpinCorrNNN,obj.Sz_SpinCorrNNNUnc,obj.Sx_SpinCorrNNN,obj.Sx_SpinCorrNNNUnc,obj.AnisotropyNNN,obj.AnisotropyNNNUnc,obj.Sz_SpinCorrNNNN,obj.Sz_SpinCorrNNNNUnc,obj.Sx_SpinCorrNNNN,obj.Sx_SpinCorrNNNNUnc,obj.AnisotropyNNNN,obj.AnisotropyNNNNUnc);
            dlmwrite(FPath,FullMat,'-append','delimiter',' ');
        end
        
        function showCorrMat(obj)
            if ~exist('MySpecialColorMap','var')
                run('load_MySpecialColorMap.m');
            end
            Lims = [-0.2,0.2];
            UncLims = [0,0.02];
            for ii = 1:size(obj.Sz_SpinCorrMat,3)
                NCols = 3;
                NRows = 3;
                FigName = sprintf('%d',ii);
                h = figure('name',FigName);
                subplot(NRows,NCols,1)
                imagesc(obj.Sz_SpinCorrMat(:,:,ii),Lims); axis equal; axis image;
                title('SzSz')
                subplot(NRows,NCols,2)
                imagesc(obj.Sz_SpinCorrMatUnc(:,:,ii),UncLims); axis equal; axis image;
                title('SzSz Unc')
                subplot(NRows,NCols,3)
                imagesc(abs(obj.Sz_SpinCorrMat(:,:,ii)./obj.Sz_SpinCorrMatUnc(:,:,ii)),[0,1]); axis equal; axis image;
                title('SzSz/SzSz Unc')
                subplot(NRows,NCols,4)
                imagesc(obj.Sx_SpinCorrMat(:,:,ii),Lims); axis equal; axis image;
                title('SxSx')
                subplot(NRows,NCols,5)
                imagesc(obj.Sx_SpinCorrMatUnc(:,:,ii),UncLims); axis equal; axis image;
                subplot(NRows,NCols,6)
                imagesc(abs(obj.Sx_SpinCorrMat(:,:,ii)./obj.Sx_SpinCorrMatUnc(:,:,ii)),[0,1]); axis equal; axis image;
                subplot(NRows,NCols,7)
                imagesc(abs(obj.SpinCorrMat_Anisotropy(:,:,ii)),[0,1]); axis equal; axis image;
                title('Anisotropy')
                subplot(NRows,NCols,8)
                imagesc(obj.SpinCorrMat_AnisotropyUnc(:,:,ii),[0,1]); axis equal; axis image;
                subplot(NRows,NCols,9)
                imagesc(abs(obj.SpinCorrMat_Anisotropy(:,:,ii)./obj.SpinCorrMat_AnisotropyUnc(:,:,ii)),[0,1]); axis equal; axis image;
                colormap(MySpecialColorMap);
                
            end
        end
        
        function generateCorrMatTextFiles(obj)
            for ii = 1:length(obj.Polarization)
                %Sz correlators
                Mat = obj.Sz_SpinCorrMat(:,:,ii);
                FileDir = fullfile(obj.SavePath,'Sz_Correlators');
                if ~exist(FileDir,'dir')
                    mkdir(FileDir);
                end
                FileName = sprintf('SzSz_CorrMat_P=%0.2f.txt',obj.Polarization(ii));
                FullPath = fullfile(FileDir,FileName);
                dlmwrite(FullPath,Mat,'delimiter',' ')
                %Sz uncertainties
                MatUnc = obj.Sz_SpinCorrMatUnc(:,:,ii);
                FileName = sprintf('SzSz_CorrMat_Unc_P=%0.2f.txt',obj.Polarization(ii));
                FullPath = fullfile(FileDir,FileName);
                dlmwrite(FullPath,MatUnc,'delimiter',' ')
                %Sx correlators
                Mat = obj.Sx_SpinCorrMat(:,:,ii);
                FileDir = fullfile(obj.SavePath,'Sx_Correlators');
                if ~exist(FileDir,'dir')
                    mkdir(FileDir);
                end
                FileName = sprintf('SxSx_CorrMat_P=%0.2f.txt',obj.Polarization(ii));
                FullPath = fullfile(FileDir,FileName);
                dlmwrite(FullPath,Mat,'delimiter',' ')
                %Sx uncertainties
                MatUnc = obj.Sx_SpinCorrMatUnc(:,:,ii);
                FileName = sprintf('SxSx_CorrMat_Unc_P=%0.2f.txt',obj.Polarization(ii));
                FullPath = fullfile(FileDir,FileName);
                dlmwrite(FullPath,MatUnc,'delimiter',' ')
                %Anisotropy
                Mat = obj.SpinCorrMat_Anisotropy(:,:,ii);
                FileDir = fullfile(obj.SavePath,'Anisotropy');
                if ~exist(FileDir,'dir')
                    mkdir(FileDir);
                end
                FileName = sprintf('Anisotropy_Mat_P=%0.2f.txt',obj.Polarization(ii));
                FullPath = fullfile(FileDir,FileName);
                dlmwrite(FullPath,Mat,'delimiter',' ')
                %Anisotropy uncertainties
                MatUnc = obj.SpinCorrMat_AnisotropyUnc(:,:,ii);
                FileName = sprintf('Anisotropy_Mat_Unc_P=%0.2f.txt',obj.Polarization(ii));
                FullPath = fullfile(FileDir,FileName);
                dlmwrite(FullPath,MatUnc,'delimiter',' ')
            end
        end
    end
    
end

