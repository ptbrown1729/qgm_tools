classdef CorrDataSetRep < handle
    %Class to store set of five folders need to compute spin correlators.
    %These individual data sets are MI, Up_Sz, Down_Sz, Up_Sx, and Down_Sx
    %These are stored in DataFolder class objects. This class is a
    %convenient container for those...Previously had only a single class,
    %but it was too much work to write out five separate properties for
    %e.g. MI, Up_Sz,... correlators.
    
    properties
        %Dataset and data usage properties.
        Date = {}
        Folders = []
        ExcludePictures = {}
        
        %Analysis settings
        UseEfficiencies = 1
        DetectionEff = 0.96
        
        %
        InitFG = 0
        NonIntFG
        
        %physical constants
        mLi = 9.9883e-27; %kg
        h = 6.6261e-34; %J.s
        hbar = 1.0546e-34; %J.s
        kb =  1.3806e-23;
        a = 1064e-9/sqrt(2); %m
        
        %Azimuthal average params
        BinEdges
        BinCenters
        MeanBinDist
        MeanBinDistUnc  
        
        %polarizations
        GPol
        GPolUnc
        Sz_LocalPolarization
        Sz_LocalPolarizationUnc
        Sx_LocalPolarization
        Sx_LocalPolarizationUnc
        
        %Datasets
        MI = DataFolder();
        Up_Sz = DataFolder();
        Down_Sz = DataFolder();
        Up_Sx = DataFolder();
        Down_Sx = DataFolder();
        
        %Densities
        Sz_Up_Dens
        Sz_Up_DensUnc
        Sz_Down_Dens
        Sz_Down_DensUnc
        Sx_Up_Dens
        Sx_Up_DensUnc
        Sx_Down_Dens
        Sx_Down_DensUnc
        MI_Dens
        MI_DensUnc
        %Derived densities
        Sz_Dens
        Sz_DensUnc
        Sx_Dens
        Sx_DensUnc  
        
        %individual component correlators
        SzSz_Up_Corr
        SzSz_Up_CorrUnc
        SzSz_Down_Corr
        SzSz_Down_CorrUnc
        SxSx_Up_Corr
        SxSx_Up_CorrUnc
        SxSx_Down_Corr
        SxSx_Down_CorrUnc
        MI_Corr
        MI_CorrUnc
        
        %Derived Correlators...combining all pictures
        SzSz_Corr
        SzSz_CorrUnc
        SxSx_Corr
        SxSx_CorrUnc
        Anisotropy
        AnisotropyUnc
        
        %Nearest neighbor derived correlators, for convenience
        SzSz_NNCorr
        SzSz_NNCorrUnc
        SxSx_NNCorr
        SxSx_NNCorrUnc
        NNAnisotropy
        NNAnisotropyUnc
        
        %Also NN corr
        SzSz_NNNCorr
        SzSz_NNNCorrUnc
        SxSx_NNNCorr
        SxSx_NNNCorrUnc
        NNNAnisotropy
        NNNAnisotropyUnc
        
        %also interested in momentum space.
        Static_Structure_FactorSz
        Static_Structure_FactorSzUnc
        Static_Structure_FactorSx
        Static_Structure_FactorSxUnc
        Qxs
        Qys
        
    end
    
    methods (Access = public)
        
        function obj = CorrDataSetRep(DateCell,DataSets,ExcludePictures,BinEdges)
            %         function obj = CorrDataSetV2(DataSetList,BinEdges)
            
            if exist('DateCell','var')
                if ~exist('BinEdges','var')
                    BinEdges = sqrt(linspace(0,30^2,30));
                end
                
                if ~exist('ExcludePictures','var')||(length(ExcludePictures)~=5)
                    ExcludePictures = {[],[],[],[],[]};
                end
                
                obj.BinEdges = BinEdges;
                obj.BinCenters = 0.5*(BinEdges(2:end)+BinEdges(1:end-1));
                
                obj.Date = DateCell;
                obj.Folders = DataSets;
                DataSetList = {DateCell,DateCell,DateCell,DateCell,DateCell};
                for ii = 1:5
                    DataSetList{ii}{4} = DataSets(ii);
                    DataSetList{ii}{5} = 1;
                    DataSetList{ii}{6} = 1;
                end
                
                

                obj.MI = DataFolder(DataSetList{5},ExcludePictures{5},BinEdges);
                Grid = obj.MI.DistGrid;
                obj.MeanBinDist = obj.MI.BinAvg;
                obj.MeanBinDistUnc = obj.MI.BinUnc;
                
                
                obj.Up_Sz = DataFolder(DataSetList{3},ExcludePictures{3},BinEdges,'external',Grid);
                obj.Down_Sz = DataFolder(DataSetList{1},ExcludePictures{1},BinEdges,'external',Grid);
                obj.Up_Sx = DataFolder(DataSetList{2},ExcludePictures{2},BinEdges,'external',Grid);
                obj.Down_Sx = DataFolder(DataSetList{4},ExcludePictures{4},BinEdges,'external',Grid);
                
                
                obj.doAnalysis();
            end
            
        end
        
        function InitializeFG(obj)
            tic;
            obj.InitFG = 1;
            %instantiate class.
            obj.NonIntFG = NonIntFG(1); 
            InitT = toc;
            fprintf('Took %0.2f s to initialize NonIntFG class \n',InitT);
        end
        
        function doAnalysis(obj)
            
            if ~obj.UseEfficiencies
                obj.DetectionEff = 1;
            end
            
            %Get densities and correct for detection efficiency.
            obj.MI_Dens = obj.MI.Occs_AzAvg/obj.DetectionEff;
            obj.MI_DensUnc = obj.MI.Occs_AzAvgUnc/obj.DetectionEff;
            obj.Sz_Up_Dens = obj.Up_Sz.Occs_AzAvg/obj.DetectionEff;
            obj.Sz_Up_DensUnc = obj.Up_Sz.Occs_AzAvgUnc/obj.DetectionEff;
            obj.Sz_Down_Dens = obj.Down_Sz.Occs_AzAvg/obj.DetectionEff;
            obj.Sz_Down_DensUnc = obj.Down_Sz.Occs_AzAvgUnc/obj.DetectionEff;
            obj.Sx_Up_Dens = obj.Up_Sx.Occs_AzAvg/obj.DetectionEff;
            obj.Sx_Up_DensUnc = obj.Up_Sx.Occs_AzAvgUnc/obj.DetectionEff;
            obj.Sx_Down_Dens = obj.Down_Sx.Occs_AzAvg/obj.DetectionEff;
            obj.Sx_Down_DensUnc = obj.Down_Sx.Occs_AzAvgUnc/obj.DetectionEff;
            
            obj.Sz_Dens = obj.Sz_Up_Dens+obj.Sz_Down_Dens;
            obj.Sz_DensUnc = sqrt(obj.Sz_Up_DensUnc.^2+obj.Sz_Down_DensUnc.^2);
            obj.Sx_Dens = obj.Sx_Up_Dens+obj.Sx_Down_Dens;
            obj.Sx_DensUnc = sqrt(obj.Sx_Up_DensUnc.^2+obj.Sx_Down_DensUnc.^2);
            
            %get local polarization, Sx and Sz
            %Uncertainty_C = [(2*B/(A+B)^2 * Unc_A)^2 + (2*A/(A+B)^2 * Unc_B)^2]^0.5
            obj.Sz_LocalPolarization = (obj.Sz_Up_Dens-obj.Sz_Down_Dens)./(obj.Sz_Up_Dens+obj.Sz_Down_Dens);
            obj.Sz_LocalPolarizationUnc = sqrt((2*obj.Sz_Down_Dens./(obj.Sz_Up_Dens+obj.Sz_Down_Dens).^2.*obj.Sz_Up_DensUnc).^2+(2*obj.Sz_Up_Dens./(obj.Sz_Up_Dens+obj.Sz_Down_Dens).^2.*obj.Sz_Down_DensUnc).^2);
            obj.Sx_LocalPolarization = (obj.Sx_Up_Dens-obj.Sx_Down_Dens)./(obj.Sx_Up_Dens+obj.Sx_Down_Dens);
            obj.Sx_LocalPolarizationUnc = sqrt((2*obj.Sx_Down_Dens./(obj.Sx_Up_Dens+obj.Sx_Down_Dens).^2.*obj.Sx_Up_DensUnc).^2+(2*obj.Sx_Up_Dens./(obj.Sx_Up_Dens+obj.Sx_Down_Dens).^2.*obj.Sx_Down_DensUnc).^2);
            
            %global polarization
            NupSingles = obj.Up_Sz.MeanAtomNum/obj.DetectionEff;
            NupSingles_Unc = obj.Up_Sz.AtomNumSD/obj.DetectionEff;
            NdownSingles = obj.Down_Sz.MeanAtomNum/obj.DetectionEff;
            NdownSingles_Unc = obj.Down_Sz.AtomNumSD/obj.DetectionEff;
            obj.GPol = (NupSingles-NdownSingles)./(NupSingles+NdownSingles);
            obj.GPolUnc = sqrt((2*NupSingles./(NdownSingles+NupSingles).^2.*NdownSingles_Unc).^2+(2*NdownSingles./(NupSingles+NdownSingles).^2.*NupSingles_Unc).^2);
            
            %get correlators and correct for efficiencies
            obj.SzSz_Up_Corr = obj.Up_Sz.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.SzSz_Up_CorrUnc = obj.Up_Sz.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            obj.SzSz_Down_Corr = obj.Down_Sz.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.SzSz_Down_CorrUnc = obj.Down_Sz.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            obj.SxSx_Up_Corr = obj.Up_Sx.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.SxSx_Up_CorrUnc = obj.Up_Sx.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            obj.SxSx_Down_Corr = obj.Down_Sx.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.SxSx_Down_CorrUnc = obj.Down_Sx.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            obj.MI_Corr = obj.MI.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.MI_CorrUnc = obj.MI.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            
            %Derived correlators
            obj.SzSz_Corr = 2*(obj.SzSz_Up_Corr+obj.SzSz_Down_Corr) - obj.MI_Corr;
            obj.SzSz_CorrUnc = sqrt((2*obj.SzSz_Up_CorrUnc).^2+(2*obj.SzSz_Down_CorrUnc).^2+obj.MI_CorrUnc.^2);
            
            obj.SxSx_Corr = 2*(obj.SxSx_Up_Corr+obj.SxSx_Down_Corr) - obj.MI_Corr;
            obj.SxSx_CorrUnc = sqrt((2*obj.SxSx_Up_CorrUnc).^2+(2*obj.SxSx_Down_CorrUnc).^2+obj.MI_CorrUnc.^2);
            
            obj.Anisotropy = 1-obj.SzSz_Corr./obj.SxSx_Corr;
            obj.AnisotropyUnc = abs(obj.SzSz_Corr./obj.SxSx_Corr).*sqrt((obj.SzSz_CorrUnc./obj.SzSz_Corr).^2+(obj.SxSx_CorrUnc./obj.SxSx_Corr).^2);
                       
            %NNcorr
            obj.SzSz_NNCorr = squeeze(obj.SzSz_Corr(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.SzSz_NNCorrUnc = squeeze(obj.SzSz_CorrUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.SxSx_NNCorr = squeeze(obj.SxSx_Corr(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.SxSx_NNCorrUnc = squeeze(obj.SxSx_CorrUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.NNAnisotropy = squeeze(obj.Anisotropy(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.NNAnisotropyUnc = squeeze(obj.AnisotropyUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
        
            %NNNCorr
            obj.SzSz_NNNCorr = squeeze(obj.SzSz_Corr(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.SzSz_NNNCorrUnc = squeeze(obj.SzSz_CorrUnc(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.SxSx_NNNCorr = squeeze(obj.SxSx_Corr(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.SxSx_NNNCorrUnc = squeeze(obj.SxSx_CorrUnc(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.NNNAnisotropy = squeeze(obj.Anisotropy(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
            obj.NNNAnisotropyUnc = squeeze(obj.AnisotropyUnc(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
    
        
        end
        
        function FigHandle = fitFGToEdge(obj,nCutoff,InitParams,FixedParams)
            %FigHandle = fitFGToEdge(obj,nCutoff,InitParams,FixedParams)
            %Params = [DeltaMu,T]
            
            if ~obj.InitFG
                fprintf('Non-interacting Fermi gas class not initialized \n')
                fprintf('Initializing...\n');
                obj.InitializeFG();
            end
            
            if ~exist('Cutoff','var')
                nCutoff = 0.25;
            end
            
            if ~exist('InitParams','var')
                InitParams = [0.3,0.5];
            end
            
            if ~exist('FixedParams','var')
                FixedParams = [0,0];
            end
            
            nstart = 0.05;

                %experimental quantities to fit.
                %to estimate real density, can invert densities.
                %gives some negative nums in the square root?
%                 a = obj.Sz_Up_Dens;
%                 b = obj.Sz_Down_Dens;
%                 n1 = 0.5*((1+a-b) - sqrt((a-b).^2 + (1-2*a-2*b)));
%                 n3 = 0.5*((1+b-a) - sqrt((a-b).^2 + (1-2*a-2*b)));
%                 n = n1+n3;
                n = obj.MI_Dens;
                nTofit = n(n<nCutoff & n>nstart);
                
                pTofit = obj.Sz_LocalPolarization(n<nCutoff & n>nstart);
                puncTofit = obj.Sz_LocalPolarizationUnc(n<nCutoff & n>nstart);
                puncTofit(isnan(puncTofit)) = 1;
                puncTofit(puncTofit==0) = 1;
                
                Corr1 = squeeze(obj.SzSz_Up_Corr(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
                Corr1Unc = squeeze(obj.SzSz_Up_CorrUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
                Corr2 = squeeze(obj.SzSz_Down_Corr(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
                Corr2Unc = squeeze(obj.SzSz_Down_CorrUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
                DensCorr = squeeze(obj.MI_Corr(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
                DensCorrUnc = squeeze(obj.MI_CorrUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
                
                
                Corr1ToFit = Corr1(n<nCutoff & n>nstart);
                Corr1UncToFit = Corr1Unc(n<nCutoff & n>nstart);
                Corr2ToFit = Corr2(n<nCutoff & n>nstart);
                Corr2UncToFit = Corr2Unc(n<nCutoff & n>nstart);
                DensCorrToFit = DensCorr(n<nCutoff & n>nstart);
                DensCorrUncToFit = DensCorrUnc(n<nCutoff & n>nstart);

                %<d_i d_j>_c = n_up^2*<n_up_i n_up_j>_c + n_down^2*<n_down_i n_down_j>_c + n_down^2*n_up^2
%                 DoublesCorrInferred = OnesCorr(FitParams,RInterp).*ThreesDensity(FitParams,RInterp).^2 + ThreesCorr(FitParams,RInterp).*OnesDensity(FitParams,RInterp).^2 + OnesCorr(FitParams,RInterp).*ThreesCorr(FitParams,RInterp);
                %This correlator is even worse combinatorally...
                %<(n_i_up + n_i_down - 2*n_i_up*n_i_down)(n_j_up +n_j_down- 2*n_j_up*n_j_down>_c
                % = <n_i_up n_j_up>_c + <n_i_down n_j_down>_c + 4<d_i d_j> ...
                % + <n_i_up n_j_down>_c + <n_i_down n_j_up>_c ...
                % - 2*<(n_i_up + n_i_down) d_j>_c - 2*<d_i (n_j_up + n_j_down)>_c
                %We already have the first line...our measured correlators
                %and the doubles. The second line vanishes for
                %non-interacting particles. 
                %We can evaluate the <n_i_up d_j> term using Wick's theorem.
                %Conveniently, only one term contributes ...
                %<n_i_up d_j> = - <c^dag_i_up c_j_up><c^dag_j_upc_i_up><c^dag_j_down c_j_down> ...
                % = <n_i_up n_j_up>_c <n_down>
                %where this last identity comes from using Wick's theorem
                %in reverse...i.e. we already found this term when we were
                %calculating the density correlations for a single
                %species.
                %%%Ones singles...similar to above
                %<(n_i_up-n_i_up*n_i_down)(n_j_up-n_j_up*n_j_down)>_c
                % = <n_i_up n_j_up>_c + <d_i d_j>_c - <n_i_up d_j> - <d_i n_j_up>
                %= <n_i_up n_j_up>_c + <d_i d_j>_c - 2*<n_i_up n_j_up>_c*n_down
                
%                 SinglesCorrInferred = OnesCorr(FitParams,RInterp)+ThreesCorr(FitParams,RInterp) + 4*DoublesCorrInferred -4*ThreesCorr(FitParams,RInterp).*OnesDensity(FitParams,RInterp) - 4*OnesCorr(FitParams,RInterp).*ThreesDensity(FitParams,RInterp);

%                 n1fn = @(singles,dmu,T) obj.NonIntFG.n1fn_singles_dmu_T(singles,dmu,T); %@(n,dmu,T) 0.5*(pfnReg(n,dmu,T).*n + n);
%                 n3fn = @(singles,dmu,T) obj.NonIntFG.n2fn_singles_dmu_T(singles,dmu,T); %n - n1fn(n,dmu,T);
%                 DoublesCorrfn = @(singles,dmu,T) obj.NonIntFG.c1fn_singles_dmu_T(singles,dmu,T).*n3fn(singles,dmu,T).^2 + obj.NonIntFG.c2fn_singles_dmu_T(singles,dmu,T).*n1fn(singles,dmu,T).^2 + obj.NonIntFG.c1fn_singles_dmu_T(singles,dmu,T).*obj.NonIntFG.c2fn_singles_dmu_T(singles,dmu,T); 
%                 SinglesCorrfn = @(singles,dmu,T) obj.NonIntFG.c1fn_singles_dmu_T(singles,dmu,T)+obj.NonIntFG.c2fn_singles_dmu_T(singles,dmu,T) + 4*DoublesCorrfn(singles,dmu,T) - 4*obj.NonIntFG.c2fn_singles_dmu_T(singles,dmu,T).*n1fn(singles,dmu,T) - 4*obj.NonIntFG.c1fn_singles_dmu_T(singles,dmu,T).*n3fn(singles,dmu,T);

                
                n1fn = @(n,dmu,T) obj.NonIntFG.n1fn_n_dmu_T(n,dmu,T); %@(n,dmu,T) 0.5*(pfnReg(n,dmu,T).*n + n);
                n3fn = @(n,dmu,T) obj.NonIntFG.n2fn_n_dmu_T(n,dmu,T); %n - n1fn(n,dmu,T);
                DoublesCorrfn = @(n,dmu,T) obj.NonIntFG.c1fn_n_dmu_T(n,dmu,T).*n3fn(n,dmu,T).^2 + obj.NonIntFG.c2fn_n_dmu_T(n,dmu,T).*n1fn(n,dmu,T).^2 + obj.NonIntFG.c1fn_n_dmu_T(n,dmu,T).*obj.NonIntFG.c2fn_n_dmu_T(n,dmu,T); 
                SinglesCorrfn = @(n,dmu,T) obj.NonIntFG.c1fn_n_dmu_T(n,dmu,T)+obj.NonIntFG.c2fn_n_dmu_T(n,dmu,T) + 4*DoublesCorrfn(n,dmu,T) - 4*obj.NonIntFG.c2fn_n_dmu_T(n,dmu,T).*n1fn(n,dmu,T) - 4*obj.NonIntFG.c1fn_n_dmu_T(n,dmu,T).*n3fn(n,dmu,T);
                n1SinglesCorrfn = @(n,dmu,T) obj.NonIntFG.c1fn_n_dmu_T(n,dmu,T).*(1-2*n3fn(n,dmu,T)) + DoublesCorrfn(n,dmu,T);
                n3SinglesCorrfn = @(n,dmu,T) obj.NonIntFG.c2fn_n_dmu_T(n,dmu,T).*(1-2*n1fn(n,dmu,T)) + DoublesCorrfn(n,dmu,T);
                %actually want these all in terms of nsingles...
                
                FitP = @(P)(obj.NonIntFG.pfn_n_dmu_T(nTofit,P(1),P(2))-pTofit)./puncTofit;
                FitC1 = @(P)(n1SinglesCorrfn(nTofit,P(1),P(2))-Corr1ToFit)./Corr1UncToFit;
                FitC2 = @(P)(n3SinglesCorrfn(nTofit,P(1),P(2))-Corr2ToFit)./Corr2UncToFit;
                FitSinglesCorr = @(P) (SinglesCorrfn(nTofit,P(1),P(2))-DensCorrToFit)./DensCorrUncToFit;
                
                
                Fitfn = @(P) [FitP(P.*(1-FixedParams)+InitParams.*FixedParams),FitC1(P.*(1-FixedParams)+InitParams.*FixedParams),FitC2(P.*(1-FixedParams)+InitParams.*FixedParams)];
                FitP = lsqnonlin(Fitfn,InitParams);
                T = FitP(2);
                DelMu = FitP(1);
                
                ninterp = transpose(linspace(min(obj.MI_Dens),max(obj.MI_Dens),100));
                
                FigHandle = figure;
                NRows = 2;
                NCols = 2;
                subplot(NRows,NCols,1)
                errorbar(obj.MI_Dens,obj.Sz_LocalPolarization,obj.Sz_LocalPolarizationUnc)
                hold on;
                plot(ninterp,obj.NonIntFG.pfn_n_dmu_T(ninterp,DelMu,T));
                grid on;
                xlabel('Singles Density n^s')
                ylabel('Singles Polarization, p^s')
                legend({'Experiment','Non-Interacting FG'})
                
                subplot(NRows,NCols,2)
                errorbar(obj.MI_Dens,DensCorr,DensCorrUnc);
                hold on;
                plot(ninterp,SinglesCorrfn(ninterp,DelMu,T));
                grid on;
                xlabel('Singles Density n^s')
                ylabel('Singles Correlator')
                
                subplot(NRows,NCols,3)
                errorbar(obj.MI_Dens,Corr1,Corr1Unc);
                hold on;
                plot(ninterp,n1SinglesCorrfn(ninterp,DelMu,T));
                grid on;
                xlabel('Singles Density n^s')
                ylabel('Ones Correlator')
                
                subplot(NRows,NCols,4)
                errorbar(obj.MI_Dens,Corr2,Corr2Unc);
                hold on;
                plot(ninterp,n3SinglesCorrfn(ninterp,DelMu,T));
                grid on;
                xlabel('Singles Density n^s')
                ylabel('Ones Correlator')
                
                suptitle(sprintf('Non-Interacting FG fit to edge, using region from n^s = %0.1f to %0.1f \nT = %0.2f t, DelMu %0.2f t',nstart,nCutoff,T,DelMu));
        end
        
        function getStructureFactor(obj)
            %first allowed q values.
            %start with only q = 0.
            NumQs = 10;
            [Qxs,Qys] = meshgrid(linspace(0,pi,NumQs),linspace(0,pi,NumQs));
%             Qs = cat(2,transpose(linspace(0,pi,10)),transpose(linspace(0,pi,10)));

            StructFactSzTemp = zeros(NumQs,NumQs,obj.MI.NBins);
            StructFactSzUncTemp = zeros(NumQs,NumQs,obj.MI.NBins);
            StructFactSxTemp = zeros(NumQs,NumQs,obj.MI.NBins);
            StructFactSxUncTemp = zeros(NumQs,NumQs,obj.MI.NBins);
            
            Indices = (1:(2*obj.MI.NumNeighbors+1))-obj.MI.CenterIndex_CorrMatrix;
            [IndicesX,IndicesY] = meshgrid(Indices,Indices);
            
            for ii = 1:NumQs
                for jj = 1:NumQs
                    ExpSingle = exp(1i*Qxs(ii,jj)*IndicesX+1i*Qys(ii,jj)*IndicesY);
                    Exp = repmat(ExpSingle,[1,1,obj.MI.NBins]);
                    StructFactSzTemp(ii,jj,:) = squeeze(sum(sum(obj.SzSz_Corr.*Exp,2),1));
                    StructFactSzUncTemp(ii,jj,:) = sqrt(squeeze(sum(sum(obj.SzSz_CorrUnc.^2,2),1)));
                    StructFactSxTemp(ii,jj,:) = squeeze(sum(sum(obj.SzSz_Corr.*Exp,2),1));
                    StructFactSxUncTemp(ii,jj,:) = sqrt(squeeze(sum(sum(obj.SzSz_CorrUnc.^2,2),1)));
                end
            end
            
            obj.Static_Structure_FactorSz = StructFactSzTemp; %squeeze(sum(sum(obj.Density_Corr_AzAvg,2),1));
            obj.Static_Structure_FactorSzUnc = StructFactSzUncTemp; %sqrt(squeeze(sum(sum(obj.Density_Corr_AzAvgUnc.^2,2),1)));
            obj.Static_Structure_FactorSx = StructFactSxTemp;
            obj.Static_Structure_FactorSxUnc = StructFactSxUncTemp;
            obj.Qxs = Qxs;
            obj.Qys = Qys;
        end
        
        function showStructureFactor(obj,Bin)
            if ~exist('Bin','var')
                Bin = 1;
            end
            
            FigHandle = figure('name','Static Structure Factor');
            NRows = 2;
            NCols = 2;
            
            subplot(NRows,NCols,1)
            imagesc(real(obj.Static_Structure_FactorSz(:,:,Bin)));
%             colormap(obj.MI.ColorMap)
            ax = gca;
            ax.XTickLabel = round(100*obj.Qxs(1,:)/pi)/100; %hack-y way to get two decimal places?
            ax.YTickLabel = round(100*obj.Qys(:,1)/pi)/100;
            xlabel('qx (pi)')
            ylabel('qy (pi)')
            axis equal;
            axis image;
            caxis([0,5]);
            colorbar;
            title('S_z(q)');
            
            subplot(NRows,NCols,3)
            imagesc(real(obj.Static_Structure_FactorSz(:,:,Bin))./real(obj.Static_Structure_FactorSzUnc(:,:,Bin)));
            caxis([0,5]);
            title('S_z Unc');
            axis equal; axis image;
            colorbar;
            
            subplot(NRows,NCols,2)
            imagesc(real(obj.Static_Structure_FactorSx(:,:,Bin)));
%             colormap(obj.MI.ColorMap)
            ax = gca;
            ax.XTickLabel = round(100*obj.Qxs(1,:)/pi)/100;
            ax.YTickLabel = round(100*obj.Qys(:,1)/pi)/100;
            xlabel('qx (pi)')
            ylabel('qy (pi)')
            axis equal;
            axis image;
            caxis([0,5]);
            colorbar;
            title('S_x(q)');
            
            
            subplot(NRows,NCols,4)
            imagesc(real(obj.Static_Structure_FactorSx(:,:,Bin))./real(obj.Static_Structure_FactorSxUnc(:,:,Bin)));
            caxis([0,5]);
            title('S_z Unc');
            axis equal; axis image;
            colorbar;
            
            
            
        end
        
        function FigHandle = showCorr(obj)
            Title = sprintf('%d/%d/%d, Folders %d, %d, %d, %d, and %d',obj.Date{1},obj.Date{2},obj.Date{3},obj.Folders(1),obj.Folders(2),obj.Folders(3),obj.Folders(4),obj.Folders(5));
            FigHandle = figure('name',Title);
            NRows = 2;
            NCols = 2;
            subplot(NRows,NCols,1)
            errorbar(obj.MeanBinDist,obj.SzSz_NNCorr,obj.SzSz_NNCorrUnc);
            hold on;
            errorbar(obj.MeanBinDist,obj.SxSx_NNCorr,obj.SxSx_NNCorrUnc);
            SinglesCorr = squeeze(obj.MI_Corr(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            SinglesCorrUnc = squeeze(obj.MI_CorrUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            errorbar(obj.MeanBinDist,SinglesCorr,SinglesCorrUnc);
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('Spin Correlator')
            legend({'Sz Corr','Sx Corr','Singles Density Corr'})
            
            subplot(NRows,NCols,2)
            errorbar(obj.MI_Dens,obj.SzSz_NNCorr,obj.SzSz_NNCorrUnc);
            hold on;
            errorbar(obj.MI_Dens,obj.SxSx_NNCorr,obj.SxSx_NNCorrUnc);
            errorbar(obj.MI_Dens,SinglesCorr,SinglesCorrUnc);
            grid on;
            xlabel('Singles Density')
            ylabel('Spin Correlator')
            legend({'Sz Corr','Sx Corr','Singles Density Corr'})
            
            subplot(NRows,NCols,3)
            errorbar(obj.MeanBinDist,obj.NNAnisotropy,obj.NNAnisotropyUnc);
            hold on;
            grid on;
            ylim([0,1])
            xlabel('Singles Density')
            ylabel('Anisotropy = 1-C^z/C^x')
            
            subplot(NRows,NCols,4)
            errorbar(obj.MI_Dens,obj.NNAnisotropy,obj.NNAnisotropyUnc);
            hold on;
            grid on;
            ylim([0,1])
            xlabel('Singles Density')
            ylabel('Anisotropy = 1-C^z/C^x')
            
            suptitle(Title);
            
            
        end
        
        function FigHandle = showCorrMat(obj,Bin)
            if ~exist('Bin','var')
                Bin = 1;
            end
            
            Title = sprintf('%d/%d/%d, Folders %d, %d, %d, %d, and %d',obj.Date{1},obj.Date{2},obj.Date{3},obj.Folders(1),obj.Folders(2),obj.Folders(3),obj.Folders(4),obj.Folders(5));
            FigHandle = figure('name',Title);
            imagesc(squeeze(obj.SzSz_Corr(:,:,Bin)),[-0.2,0.2]);
            axis equal; axis image;
            obj.MI.getColorMap()
            colormap(obj.MI.ColorMap)
        end
        
        function FigHandle = showAzAvg(obj)
            %FigHandle = showAzAvg(obj)
            %Display all data sets azimuthal averages.
            
            Title = sprintf('%d/%d/%d, Folders %d, %d, %d, %d, and %d',obj.Date{1},obj.Date{2},obj.Date{3},obj.Folders(1),obj.Folders(2),obj.Folders(3),obj.Folders(4),obj.Folders(5));
            FigHandle = figure('name',Title);
            NRows = 2;
            NCols = 2;
            
            subplot(NRows,NCols,1)
                errorbar(obj.MeanBinDist,obj.Sz_Up_Dens,obj.Sz_Up_DensUnc);
                hold on;
                errorbar(obj.MeanBinDist,obj.Sz_Down_Dens,obj.Sz_Down_DensUnc);
                errorbar(obj.MeanBinDist,obj.Sx_Up_Dens,obj.Sx_Up_DensUnc);
                errorbar(obj.MeanBinDist,obj.Sx_Down_Dens,obj.Sx_Down_DensUnc);
                errorbar(obj.MeanBinDist,obj.MI_Dens,obj.MI_DensUnc);
                ylim([0,1])
                grid on;
                legend({'One Singles,Sz','Three Singles, Sz','Ones Singles, Sx','Three Singles, Sx','Singles, Sz'});
                xlabel('Lattice Sites')
                ylabel('Density')
                title('Az Avg');
            subplot(NRows,NCols,2);
                errorbar(obj.MeanBinDist,obj.Sz_LocalPolarization,obj.Sz_LocalPolarizationUnc,'k');
                grid on;
                hold on;
                plot([obj.MeanBinDist(1),obj.MeanBinDist(end)],[obj.GPol,obj.GPol],'b-')
                plot([obj.MeanBinDist(1),obj.MeanBinDist(end)],[obj.GPol+obj.GPolUnc,obj.GPol+obj.GPolUnc],'b--')
                plot([obj.MeanBinDist(1),obj.MeanBinDist(end)],[obj.GPol-obj.GPolUnc,obj.GPol-obj.GPolUnc],'b--')
                ylim([0,1])
                xlabel('Lattice Sites')
                ylabel('Polarization')
                title(sprintf('P_g = %0.2f +/- %0.2f',obj.GPol,obj.GPolUnc));
            subplot(NRows,NCols,3);
                plot(obj.Up_Sz.PictureNumberInFolder,obj.Up_Sz.AtomNumbers,'o')
                hold on;
                plot(obj.Down_Sz.PictureNumberInFolder,obj.Down_Sz.AtomNumbers,'o');
                plot(obj.Up_Sx.PictureNumberInFolder,obj.Up_Sx.AtomNumbers,'o');
                plot(obj.Down_Sx.PictureNumberInFolder,obj.Down_Sx.AtomNumbers,'o');
                plot(obj.MI.PictureNumberInFolder,obj.MI.AtomNumbers,'o');
                ax = gca;
                ax.YLim(1) = 0;
                grid on;
                xlabel('Picture # in Folder')
                ylabel('Atom Number')
                legend({'One Singles,Sz','Three Singles, Sz','Ones Singles, Sx','Three Singles, Sx','Singles, Sz'});
                title('Atom Number Stability')
                
            subplot(NRows,NCols,4)
                errorbar(obj.MeanBinDist,obj.MI_Dens,obj.MI_DensUnc);
                hold on;
                errorbar(obj.MeanBinDist,obj.Sz_Dens,obj.Sz_DensUnc);
                errorbar(obj.MeanBinDist,obj.Sx_Dens,obj.Sx_DensUnc);
                errorbar(obj.MeanBinDist,obj.Sx_LocalPolarization,obj.Sx_LocalPolarization);
                grid on;
                xlabel('Lattice Sites')
                ylabel('Density')
                title(sprintf('DetEff = %0.2f',obj.DetectionEff));
                legend({'Ns, Sz','N1s+N3s, Sz','N1s+N3s, Sx','Pol, Sx'})
                
            suptitle(Title);
        end
        
        function shrinkToSave(obj)
            %shrink all data sets to a reasonable size so they can be
            %saved...
            obj.MI.shrinkToSave();
            obj.Up_Sz.shrinkToSave();
            obj.Down_Sz.shrinkToSave();
            obj.Up_Sx.shrinkToSave();
            obj.Down_Sx.shrinkToSave();
        end
        
        function reconstituteAllData(obj)
            %for use reconstituting all data sets, bc correlator full
            %matrices are too large to save (~300mb each), but might need
            %these to do more processing. They can easily be obtained
            %again.
            obj.MI.processData();
            obj.Up_Sz.processData();
            obj.Down_Sz.processData();
            obj.Up_Sx.processData();
            obj.Down_Sx.processData();
        end
        
        function savesummaryGraphs(obj,SavePath)
             if ~exist('SavePath','var')
                SavePath = '.';
            end
            
            if ~exist('SavePath','dir') && ~strcmp(SavePath,'.') && ~isempty(SavePath);
                mkdir(SavePath);
            end
            
            if ~exist('FName','var')
                FName = sprintf('%d_%d_%d_Folders=%d-%d_PGlobal=%0.2f+-%0.2f.fig',obj.Date{1},obj.Date{2},obj.Date{3},min(obj.Folders),max(obj.Folders),obj.GPol,obj.GPolUnc);
            end
            
            %for density plots
            FPath = fullfile(SavePath,strcat('Density_',FName));
            FigHandle = obj.showAzAvg();
            savefig(FigHandle,FPath);
            
            
            %for correlator plots
            FPath = fullfile(SavePath,strcat('Correlators_',FName));
            FigHandle = obj.showCorr();
            savefig(FigHandle,FPath);
            
            %for fg fit
            FPath = fullfile(SavePath,strcat('NonIntFGFit_',FName));
            FigHandle = obj.fitFGToEdge();
            savefig(FigHandle,FPath);
        end
        
        function saveToFile(obj,SavePath,FName)
            %save each radial profile to single text file.
            if ~exist('SavePath','var')
                SavePath = '';
            end
            
            if ~exist('SavePath','dir') && ~isempty(SavePath);
                mkdir(SavePath);
            end
            
            if ~exist('FName','var')
                FName = sprintf('%d_%d_%d_Folders=%d-%d_PGlobal=%0.2f+-%0.2f.txt',obj.Date{1},obj.Date{2},obj.Date{3},min(obj.Folders),max(obj.Folders),obj.GPol,obj.GPolUnc);
            end
            FPath = fullfile(SavePath,FName);
            
            
            %write text file.
            Labels = {'RadPosition','RadPositionUnc',...
                'LocalPolarization','LocalPolarizationUnc',...
                'SzUpSinglesDensity','SzUpSinglesDensityUnc',...
                'SzDownSinglesDensity','SzDownSinglesDensityUnc',...
                'SzSinglesDensity','SzSinglesDensityUnc',...
                'SxUpSinglesDensity','SxUpSinglesDensityUnc',...
                'SxDownSinglesDensity','SxDownSinglesDensityUnc',...
                'SxSinglesDensity','SxSinglesDensityUnc',...
                'MISinglesDensity','MISinglesDensityUnc',...
                'SzSzNNCorr','SzSzNNCorrUnc',...
                'SxSxNNCorr','SxSxNNCorrUnc',...
                'AnisotropyNN','AnisotropyNNUnc',...
                'SinglesNNCorr','SinglesNNCorrUnc',...
                'SzSzNNNCorr','SzSzNNNCorrUnc',...
                'SxSxNNNCorr','SxSxNNNCorrUnc',...
                'AnisotropyNNN','AnisotropyNNNUnc',...
                'SinglesNNNCorr','SinglesNNNCorrUnc',};
            
            SzSinglesDensity = obj.Sz_Up_Dens + obj.Sz_Down_Dens;
            SzSinglesDensityUnc = sqrt(obj.Sz_Up_DensUnc.^2+obj.Sz_Down_DensUnc.^2);
            SxSinglesDensity = obj.Sx_Up_Dens + obj.Sx_Down_Dens;
            SxSinglesDensityUnc = sqrt(obj.Sx_Up_DensUnc.^2+obj.Sx_Down_DensUnc.^2);
            
            MI_NNCorr = squeeze(obj.MI_Corr(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            MI_NNCorrUnc = squeeze(obj.MI_CorrUnc(obj.MI.CenterIndex_CorrMatrix,obj.MI.CenterIndex_CorrMatrix+1,:));
            MI_NNNCorr = squeeze(obj.MI_Corr(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
            MI_NNNCorrUnc = squeeze(obj.MI_CorrUnc(obj.MI.CenterIndex_CorrMatrix+1,obj.MI.CenterIndex_CorrMatrix+1,:));
            
            Mat = cat(2,...
            obj.MeanBinDist,obj.MeanBinDistUnc,...
            obj.Sz_LocalPolarization,obj.Sz_LocalPolarizationUnc,...
            obj.Sz_Up_Dens,obj.Sz_Up_DensUnc,...
            obj.Sz_Down_Dens,obj.Sz_Down_Dens,...
            SzSinglesDensity,SzSinglesDensityUnc,...
            obj.Sx_Up_Dens,obj.Sx_Up_DensUnc,...
            obj.Sx_Down_Dens,obj.Sx_Down_Dens,...
            SxSinglesDensity,SxSinglesDensityUnc,...
            obj.MI_Dens,obj.MI_DensUnc,...
            obj.SzSz_NNCorr,obj.SzSz_NNCorrUnc,...
            obj.SxSx_NNCorr,obj.SxSx_NNCorrUnc,...
            obj.NNAnisotropy,obj.NNAnisotropyUnc,...
            MI_NNCorr,MI_NNCorrUnc,...
            obj.SzSz_NNNCorr,obj.SzSz_NNNCorrUnc,...
            obj.SxSx_NNNCorr,obj.SxSx_NNNCorrUnc,...
            obj.NNNAnisotropy,obj.NNNAnisotropyUnc,...
            MI_NNNCorr,MI_NNNCorrUnc);
            
            
            %write to file
            try
                FileID = fopen(FPath,'w');
                for ii = 1:length(Labels)-1
                    fprintf(FileID,'%s ',Labels{ii});
                end
                fprintf(FileID,'%s\n',Labels{end});
                fclose(FileID);
            catch
                try
                    fclose(FileID);
                catch
                end
            end
            
            dlmwrite(FPath,Mat,'-append','delimiter',' ');
        end
 
    end
    
    %methods that we don't need to expose to be called externally.
    methods(Access = private)
    end
    
end