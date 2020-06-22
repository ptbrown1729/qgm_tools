classdef CorrDataSetRep_v2 < handle
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
        Desc_Str = ''
        NumNeighbors = 4
        AzAvgType = 'spatial' %spatial or density or external
        
        %Analysis settings
        UseEfficiencies = 0
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
        
        %Datasets
        n = DataFolder();
        nup = DataFolder();
        ndn = DataFolder();
        
        %Densities
        Sz_Up_Dens
        Sz_Up_DensUnc
        Sz_Down_Dens
        Sz_Down_DensUnc
        ns_Dens
        ns_DensUnc
        %Derived densities
        Sz_Dens
        Sz_DensUnc
        
        %individual component correlators
        SzSz_Up_Corr
        SzSz_Up_CorrUnc
        SzSz_Down_Corr
        SzSz_Down_CorrUnc
        ns_Corr
        ns_CorrUnc
        
        %Derived Correlators...combining all pictures
        SzSz_CorrUnSymm
        SzSz_CorrUnSymmUnc
        SzSz_Corr
        SzSz_CorrUnc
        
        %Nearest neighbor derived correlators, for convenience
        SzSz_NNCorr
        SzSz_NNCorrUnc
        
        %Also NN corr
        SzSz_NNNCorr
        SzSz_NNNCorrUnc
        
        %also interested in momentum space.
        Static_Structure_FactorSz
        Static_Structure_FactorSzUnc
        Qxs
        Qys
        
    end
    
    methods (Access = public)
        
        function obj = CorrDataSetRep_v2(DateCell,DataSets,ExcludePictures,BinEdges,AzAvgType)
            %         function obj = CorrDataSetV2(DataSetList,BinEdges)
            
            if exist('DateCell','var')
                if ~exist('BinEdges','var')
                    BinEdges = sqrt(linspace(0,30^2,30));
                end
                
                if ~exist('ExcludePictures','var')||(length(ExcludePictures)~=length(DataSets))
                    ExcludePictures = {[],[],[]};
                end
                
                if ~exist('AzAvgType','var')
                    AzAvgType = 'spatial';
                end
                obj.AzAvgType = 'spatial';
                
                obj.BinEdges = BinEdges;
                obj.BinCenters = 0.5*(BinEdges(2:end)+BinEdges(1:end-1));
                
                obj.Date = DateCell;
                obj.Folders = DataSets;
                DataSetList = {DateCell,DateCell,DateCell};
                for ii = 1:3
                    DataSetList{ii}{4} = DataSets(ii);
                    DataSetList{ii}{5} = 1;
                    DataSetList{ii}{6} = 1;
                end
                
                
                obj.n = DataFolder();
                obj.n.NumNeighbors = obj.NumNeighbors;
                obj.n.initialize(DataSetList{3},ExcludePictures{3},BinEdges,AzAvgType);
%                 obj.n = DataFolder(DataSetList{3},ExcludePictures{3},BinEdges);
                
                
                Grid = obj.n.DistGrid;
                obj.MeanBinDist = obj.n.BinAvg;
                obj.MeanBinDistUnc = obj.n.BinUnc;
                
                obj.nup = DataFolder();
                obj.nup.NumNeighbors = obj.NumNeighbors;
                obj.nup.initialize(DataSetList{2},ExcludePictures{2},BinEdges,'external',Grid);
%                 obj.nup = DataFolder(DataSetList{2},ExcludePictures{2},BinEdges,'external',Grid);
                
                obj.ndn = DataFolder();
                obj.ndn.NumNeighbors = obj.NumNeighbors;
                obj.ndn.initialize(DataSetList{1},ExcludePictures{1},BinEdges,'external',Grid);
%                 obj.ndn = DataFolder(DataSetList{1},ExcludePictures{1},BinEdges,'external',Grid);
                
%                 obj.Desc_Str = sprintf('%d/%d/%d, Folders %d, %d, and %d',obj.Date{1},obj.Date{2},obj.Date{3},obj.Folders(1),obj.Folders(2),obj.Folders(3));
                obj.Desc_Str = sprintf('%04d_%02d_%02d_%03d-%03d', obj.Date{1}, obj.Date{2}, obj.Date{3},...
                                                                   min(obj.Folders), max(obj.Folders));
                
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
            obj.ns_Dens = obj.n.Occs_AzAvg/obj.DetectionEff;
            obj.ns_DensUnc = obj.n.Occs_AzAvgUnc/obj.DetectionEff;
            obj.Sz_Up_Dens = obj.nup.Occs_AzAvg/obj.DetectionEff;
            obj.Sz_Up_DensUnc = obj.nup.Occs_AzAvgUnc/obj.DetectionEff;
            obj.Sz_Down_Dens = obj.ndn.Occs_AzAvg/obj.DetectionEff;
            obj.Sz_Down_DensUnc = obj.ndn.Occs_AzAvgUnc/obj.DetectionEff;
            
            
            obj.Sz_Dens = obj.Sz_Up_Dens+obj.Sz_Down_Dens;
            obj.Sz_DensUnc = sqrt(obj.Sz_Up_DensUnc.^2+obj.Sz_Down_DensUnc.^2);
            
            %get local polarization, Sx and Sz
            %Uncertainty_C = [(2*B/(A+B)^2 * Unc_A)^2 + (2*A/(A+B)^2 * Unc_B)^2]^0.5
            obj.Sz_LocalPolarization = (obj.Sz_Up_Dens-obj.Sz_Down_Dens)./(obj.Sz_Up_Dens+obj.Sz_Down_Dens);
            obj.Sz_LocalPolarizationUnc = sqrt((2*obj.Sz_Down_Dens./(obj.Sz_Up_Dens+obj.Sz_Down_Dens).^2.*obj.Sz_Up_DensUnc).^2+(2*obj.Sz_Up_Dens./(obj.Sz_Up_Dens+obj.Sz_Down_Dens).^2.*obj.Sz_Down_DensUnc).^2);
            
            %global polarization
            NupSingles = obj.nup.MeanAtomNum/obj.DetectionEff;
            NupSingles_Unc = obj.nup.AtomNumSD/obj.DetectionEff;
            NdownSingles = obj.ndn.MeanAtomNum/obj.DetectionEff;
            NdownSingles_Unc = obj.ndn.AtomNumSD/obj.DetectionEff;
            obj.GPol = (NupSingles-NdownSingles)./(NupSingles+NdownSingles);
            obj.GPolUnc = sqrt((2*NupSingles./(NdownSingles+NupSingles).^2.*NdownSingles_Unc).^2+(2*NdownSingles./(NupSingles+NdownSingles).^2.*NupSingles_Unc).^2);
            
            %get correlators and correct for efficiencies
            obj.SzSz_Up_Corr = obj.nup.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.SzSz_Up_CorrUnc = obj.nup.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            obj.SzSz_Down_Corr = obj.ndn.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.SzSz_Down_CorrUnc = obj.ndn.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            obj.ns_Corr = obj.n.Density_Corr_AzAvg/obj.DetectionEff.^2;
            obj.ns_CorrUnc = obj.n.Density_Corr_AzAvgUnc/obj.DetectionEff.^2;
            
            %Derived correlators
            obj.SzSz_CorrUnSymm = 2*(obj.nup.Density_Corr_AzAvg_NotSymmetrized + obj.ndn.Density_Corr_AzAvg_NotSymmetrized)/obj.DetectionEff.^2 - obj.n.Density_Corr_AzAvg_NotSymmetrized/obj.DetectionEff.^2;
            obj.SzSz_CorrUnSymmUnc = sqrt((2*obj.nup.Density_Corr_AzAvgUnc_NotSymmetrized).^2+(2*obj.ndn.Density_Corr_AzAvgUnc_NotSymmetrized).^2+(obj.n.Density_Corr_AzAvgUnc_NotSymmetrized).^2)./obj.DetectionEff.^2;
            obj.SzSz_Corr = 2*(obj.SzSz_Up_Corr+obj.SzSz_Down_Corr) - obj.ns_Corr;
            obj.SzSz_CorrUnc = sqrt((2*obj.SzSz_Up_CorrUnc).^2+(2*obj.SzSz_Down_CorrUnc).^2+obj.ns_CorrUnc.^2);
            
            %NNcorr
            obj.SzSz_NNCorr = squeeze(obj.SzSz_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            obj.SzSz_NNCorrUnc = squeeze(obj.SzSz_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            
            
            %NNNCorr
            obj.SzSz_NNNCorr = squeeze(obj.SzSz_Corr(obj.n.CenterIndex_CorrMatrix+1,obj.n.CenterIndex_CorrMatrix+1,:));
            obj.SzSz_NNNCorrUnc = squeeze(obj.SzSz_CorrUnc(obj.n.CenterIndex_CorrMatrix+1,obj.n.CenterIndex_CorrMatrix+1,:));
            
            
            
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
            n = obj.ns_Dens;
            nTofit = n(n<nCutoff & n>nstart);
            
            pTofit = obj.Sz_LocalPolarization(n<nCutoff & n>nstart);
            puncTofit = obj.Sz_LocalPolarizationUnc(n<nCutoff & n>nstart);
            puncTofit(isnan(puncTofit)) = 1;
            puncTofit(puncTofit==0) = 1;
            
            Corr1 = squeeze(obj.SzSz_Up_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            Corr1Unc = squeeze(obj.SzSz_Up_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            Corr2 = squeeze(obj.SzSz_Down_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            Corr2Unc = squeeze(obj.SzSz_Down_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            DensCorr = squeeze(obj.ns_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            DensCorrUnc = squeeze(obj.ns_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            
            
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
            
            ninterp = transpose(linspace(min(obj.ns_Dens),max(obj.ns_Dens),100));
            
            FigHandle = figure;
            NRows = 2;
            NCols = 2;
            subplot(NRows,NCols,1)
            errorbar(obj.ns_Dens,obj.Sz_LocalPolarization,obj.Sz_LocalPolarizationUnc)
            hold on;
            plot(ninterp,obj.NonIntFG.pfn_n_dmu_T(ninterp,DelMu,T));
            grid on;
            xlabel('Singles Density n^s')
            ylabel('Singles Polarization, p^s')
            legend({'Experiment','Non-Interacting FG'})
            
            subplot(NRows,NCols,2)
            errorbar(obj.ns_Dens,DensCorr,DensCorrUnc);
            hold on;
            plot(ninterp,SinglesCorrfn(ninterp,DelMu,T));
            grid on;
            xlabel('Singles Density n^s')
            ylabel('Singles Correlator')
            
            subplot(NRows,NCols,3)
            errorbar(obj.ns_Dens,Corr1,Corr1Unc);
            hold on;
            plot(ninterp,n1SinglesCorrfn(ninterp,DelMu,T));
            grid on;
            xlabel('Singles Density n^s')
            ylabel('Ones Correlator')
            
            subplot(NRows,NCols,4)
            errorbar(obj.ns_Dens,Corr2,Corr2Unc);
            hold on;
            plot(ninterp,n3SinglesCorrfn(ninterp,DelMu,T));
            grid on;
            xlabel('Singles Density n^s')
            ylabel('Ones Correlator')
            
            suptitle(sprintf('Non-Interacting FG fit to edge, using region from n^s = %0.1f to %0.1f \nT = %0.2f t, DelMu %0.2f t',nstart,nCutoff,T,DelMu));
        end
        
        function [Qxs,Qys,StructFactSzTemp] = getStructureFactor(obj)
            %first allowed q values.
            %start with only q = 0.
         
            [obj.Qxs, obj.Qys, obj.Static_Structure_FactorSz, obj.Static_Structure_FactorSzUnc] = get_structure_fact(obj.SzSz_Corr, obj.SzSz_CorrUnc);
            Qxs = obj.Qxs;
            Qys = obj.Qys;
            StructFactSzTemp = obj.Static_Structure_FactorSz;
        end
        
        function [Qxs,Qys,StructFactSz,StructFactSzUnc] = getStructureFactor2(obj)
            %Get structure factor from individual structure factors...

            %can get structure factors from adding up...
            obj.n.getStructureFactor();
            obj.nup.getStructureFactor();
            obj.ndn.getStructureFactor();
            
            StructFactSz = 2*(obj.nup.Static_Structure_Factor+obj.ndn.Static_Structure_Factor)-obj.n.Static_Structure_Factor;
            
            %uncertainty from bootstrap
            obj.n.doBootstrapErrorAnalysis(10,1000);
            obj.nup.doBootstrapErrorAnalysis(10,1000);
            obj.ndn.doBootstrapErrorAnalysis(10,1000);
            StructFactSzUnc = sqrt(4*obj.nup.Static_Structure_Factor_BootstrapUnc.^2+...
                4*obj.ndn.Static_Structure_Factor_BootstrapUnc.^2+...
                obj.n.Static_Structure_Factor_BootstrapUnc.^2);
            obj.Static_Structure_FactorSz =StructFactSz;
            obj.Static_Structure_FactorSzUnc = StructFactSzUnc;
            
            Qxs = obj.n.Qxs;
            Qys = obj.n.Qys;
            obj.Qxs = Qxs;
            obj.Qys = Qys;
        end
            
        function [XForCuts,Cuts,CutsUnc,FigHandle] = showStructureFactor(obj,Bin)
            if ~exist('Bin','var')
                Bin = 1;
            end
            
            FigHandle = figure('name',obj.Desc_Str);
            subplot(1,2,1);
            imagesc(real(obj.Static_Structure_FactorSz(:,:,Bin(1))));
            obj.n.getColorMap();
            colormap bone;
            %             colormap(obj.ColorMap)
            ax = gca;
            ax.XTickLabel = obj.Qxs(1,:);
            ax.YTickLabel = obj.Qys(:,1);
            xlabel('qx')
            ylabel('qy')
            axis equal; axis image;
            colorbar;
            title('S(q)');
            
            subplot(1,2,2);
            N = (size(obj.Static_Structure_FactorSz,1)-1)/2;
            for ii = Bin
                %cut from (0,0)->(pi,0)
                Cut1 = abs(transpose(obj.Static_Structure_FactorSz(N+1,N+1:end,ii)));
                Cut1_Unc = abs(transpose(obj.Static_Structure_FactorSzUnc(N+1,N+1:end,ii)));
                %(pi,0)->(pi,pi)
                Cut2 = abs(obj.Static_Structure_FactorSz(N+2:end,end,ii));
                Cut2_Unc = abs(obj.Static_Structure_FactorSzUnc(N+2:end,end,ii));
                %(pi,pi)->(0,0)
                Cut3 = flip(abs(diag(obj.Static_Structure_FactorSz(N+1:end-1,N+1:end-1,ii))));
                Cut3_Unc = flip(abs(diag(obj.Static_Structure_FactorSzUnc(N+1:end-1,N+1:end-1,ii))));
                
                
                Cuts = [Cut1;Cut2,;Cut3];
                XForCuts = 1:length(Cuts);
                XForCuts(end-length(Cut3):end) = (XForCuts(end-length(Cut3):end)-XForCuts(end-length(Cut3)))*sqrt(2)+XForCuts(end-length(Cut3));
                CutsUnc = [Cut1_Unc;Cut2_Unc;Cut3_Unc];
                
                errorbar(XForCuts,Cuts,CutsUnc);
                hold on;
            end
            
                
            grid on;
            ax = gca;
            ax.XTick = [1,N+1,2*N+1,XForCuts(end)];
            ax.XTickLabel = {'(0,0)','(0,pi)','(pi,pi)','(0,0)'};
            title('S(q)');
            
            suptitle(obj.Desc_Str);
        end
        
        
        
        function FigHandle = showCorr(obj)
            Title = obj.Desc_Str;
            FigHandle = figure('name',Title);
            NRows = 1;
            NCols = 2;
            
            subplot(NRows,NCols,1)
            errorbar(obj.MeanBinDist,obj.SzSz_NNCorr,obj.SzSz_NNCorrUnc);
            hold on;
            SinglesCorr = squeeze(obj.ns_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            SinglesCorrUnc = squeeze(obj.ns_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            errorbar(obj.MeanBinDist,SinglesCorr,SinglesCorrUnc);
            
            UpsCorr = squeeze(obj.SzSz_Up_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            UpsCorrUnc = squeeze(obj.SzSz_Up_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            errorbar(obj.MeanBinDist,UpsCorr,UpsCorrUnc);
            
            DnsCorr = squeeze(obj.SzSz_Down_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            DnsCorrUnc = squeeze(obj.SzSz_Down_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            errorbar(obj.MeanBinDist,DnsCorr,DnsCorrUnc);
            
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('Spin Correlator')
            legend({'<S^zS^z>_c','<n^s n^s>_c','<n^sup n^sup>_c','<n^sdn n^sdn>_c'})
            
            subplot(NRows,NCols,2)
            errorbar(obj.ns_Dens,obj.SzSz_NNCorr,obj.SzSz_NNCorrUnc);
            hold on;
            errorbar(obj.ns_Dens,SinglesCorr,SinglesCorrUnc);
            errorbar(obj.ns_Dens,UpsCorr,UpsCorrUnc);
            errorbar(obj.ns_Dens,DnsCorr,DnsCorrUnc);
            grid on;
            xlabel('Singles Density, n^s')
            ylabel('Spin Correlator')
            legend({'<S^zS^z>_c','<n^s n^s>_c','<n^sup n^sup>_c','<n^sdn n^sdn>_c'})         
            
            suptitle(Title); 
        end
        
        function FigHandle = showCorrMat(obj,Bin)
            if ~exist('Bin','var')
                Bin = 1;
            end
            
            NRows = 2;
            NCols = 2;
            
            Title = obj.Desc_Str;
            FigHandle = figure('name',Title);
            
            subplot(NRows,NCols,1)
            imagesc(squeeze(obj.SzSz_Corr(:,:,Bin)),[-0.05,0.05]);
            axis equal; axis image;
            colorbar;
            obj.n.getColorMap()
            colormap(obj.n.ColorMap)
            title('Symm Corr');
            
            subplot(NRows,NCols,2)
            imagesc(abs(squeeze(obj.SzSz_Corr(:,:,Bin)))./squeeze(obj.SzSz_CorrUnc(:,:,Bin)),[-5,5]);
            axis equal; axis image;
            colorbar;
            title('Corr/Unc');
            
            subplot(NRows,NCols,3)
            [Xs,Ys] = meshgrid(1:size(obj.SzSz_Corr,2),1:size(obj.SzSz_Corr,1));
            Negs = (-1).^(Xs+Ys);
            imagesc(squeeze(obj.SzSz_Corr(:,:,Bin)).*Negs,[-0.05,0.05]);
            axis equal; axis image;
            colorbar;
            title('(-1)^{(i+j)}*C(i,j)')
            
            subplot(NRows,NCols,4);
            imagesc(squeeze(obj.SzSz_CorrUnSymm(:,:,Bin)),[-0.05,0.05]);
            axis equal; axis image;
            colorbar;
            obj.n.getColorMap()
            colormap(obj.n.ColorMap)
            title('UnSymm Corr');
            
            suptitle(Title);
            
        end
        
        function [CorrLen,Err] = showCorrelationRange(obj,BinIndex,RCutoff,CorrMatStack,CorrMatUncStack)
            %showCorrelationRange(obj,BinIndex)
            %Plot correlators vs. distance.
            %i.e. Plot C(r) for a given bin.
            if ~exist('BinIndex','var')
                BinIndex = 1;
            end
            
            if ~exist('RCutoff','var')
                RCutoff = 0;
            end
            
            if ~exist('CorrMatStack','var')
                CorrMatStack = obj.SzSz_Corr;
                CorrMatUncStack = obj.SzSz_CorrUnc;
            end
            
            
            NNeighbors = (size(CorrMatStack,1)-1)/2;
            CenterIndex = NNeighbors + 1;
            QuadrantCorrMatrix =CorrMatStack(CenterIndex:end,CenterIndex:end,BinIndex);
            QuadrantCorrMatrixUnc = CorrMatUncStack(CenterIndex:end,CenterIndex:end,BinIndex);
            [X,Y] = meshgrid(0:NNeighbors,0:NNeighbors);
            R = sqrt(X.^2+Y.^2);
            
            if ~exist('Axes','var')
                FigName = sprintf('Correlation Range');
                figure('name',FigName);
                Axes = axes;
            end
            
            %fit range
            InitP = [0,0.05,1,0];
            FixedP = [1,0,0,1];
            
%             RsNoZero = R(R~=0);
%             AbsCorrNoZero = abs(QuadrantCorrMatrix(R~=0));
            RsNoZero = R(R>RCutoff);
            AbsCorrNoZero = abs(QuadrantCorrMatrix(R>RCutoff));
            
            [Fp,~,FFH,SE] = fit1D(RsNoZero,AbsCorrNoZero,[],{'exp1D'},InitP,FixedP);
            CorrLen = Fp(3);
            Err = SE(3);
            
            errorbar(Axes,R(:),abs(QuadrantCorrMatrix(:)),QuadrantCorrMatrixUnc(:),'bo');
%             semilogy(Axes,R(:),abs(QuadrantCorrMatrix(:)),'bo');
            hold on;
            InterpRs = linspace(min(RsNoZero),max(RsNoZero),300);
            plot(Axes,InterpRs,FFH(InterpRs),'b')
            grid on;
            xlabel('Correlator Distance (Lattice Sites)')
            ylabel('|C(r)|')
            title(sprintf('Correlation length = %0.2f +/- %0.2f',CorrLen,Err));
        end
        
        function FigHandle = showAzAvg(obj)
            %FigHandle = showAzAvg(obj)
            %Display all data sets azimuthal averages.
            
            Title = obj.Desc_Str;
            FigHandle = figure('name',Title);
            NRows = 2;
            NCols = 2;
            
            subplot(NRows,NCols,1)
            errorbar(obj.MeanBinDist,obj.Sz_Up_Dens,obj.Sz_Up_DensUnc);
            hold on;
            errorbar(obj.MeanBinDist,obj.Sz_Down_Dens,obj.Sz_Down_DensUnc);
            errorbar(obj.MeanBinDist,obj.ns_Dens,obj.ns_DensUnc);
            ylim([0,1])
            grid on;
            legend({'nup_s','ndn_s','ns'});
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
            %nups
            hups = plot(obj.nup.PictureNumberInFolder,obj.nup.AtomNumbers,'bo');
            hold on;
            plot([0,max(obj.nup.PictureNumberInFolder)],[obj.nup.MeanAtomNum,obj.nup.MeanAtomNum],'b');
            plot([0,max(obj.nup.PictureNumberInFolder)],[obj.nup.MeanAtomNum+obj.nup.AtomNumSD,obj.nup.MeanAtomNum+obj.nup.AtomNumSD],'b-.');
            plot([0,max(obj.nup.PictureNumberInFolder)],[obj.nup.MeanAtomNum-obj.nup.AtomNumSD,obj.nup.MeanAtomNum-obj.nup.AtomNumSD],'b-.');
            %ndns
            hdns = plot(obj.ndn.PictureNumberInFolder,obj.ndn.AtomNumbers,'ro');
            plot([0,max(obj.ndn.PictureNumberInFolder)],[obj.ndn.MeanAtomNum,obj.ndn.MeanAtomNum],'r');
            plot([0,max(obj.ndn.PictureNumberInFolder)],[obj.ndn.MeanAtomNum+obj.ndn.AtomNumSD,obj.ndn.MeanAtomNum+obj.ndn.AtomNumSD],'r-.');
            plot([0,max(obj.ndn.PictureNumberInFolder)],[obj.ndn.MeanAtomNum-obj.ndn.AtomNumSD,obj.ndn.MeanAtomNum-obj.ndn.AtomNumSD],'r-.');
            %n singles
            hsingles = plot(obj.n.PictureNumberInFolder,obj.n.AtomNumbers,'go');
            plot([0,max(obj.n.PictureNumberInFolder)],[obj.n.MeanAtomNum,obj.n.MeanAtomNum],'g');
            plot([0,max(obj.n.PictureNumberInFolder)],[obj.n.MeanAtomNum+obj.n.AtomNumSD,obj.n.MeanAtomNum+obj.n.AtomNumSD],'g-.');
            plot([0,max(obj.n.PictureNumberInFolder)],[obj.n.MeanAtomNum-obj.n.AtomNumSD,obj.n.MeanAtomNum-obj.n.AtomNumSD],'g-.');
            ax = gca;
            ax.XLim(1) = 1;
            ax.YLim(1) = 0;
            grid on;
            xlabel('Picture # in Folder')
            ylabel('Atom Number')
            legend([hups,hdns,hsingles],{sprintf('%s=%0.0f +/- %0.0f','n^sup',obj.nup.MeanAtomNum,obj.nup.AtomNumSD),...
                sprintf('%s=%0.0f +/- %0.0f','n^sdn',obj.ndn.MeanAtomNum,obj.ndn.AtomNumSD),...
                sprintf('%s=%0.0f +/- %0.0f','n^sdn',obj.n.MeanAtomNum,obj.n.AtomNumSD)});
            title('Atom Number Stability')


%             title(sprintf('NAtoms = %0.1f +/- %0.1f ',mean(NAtomsList),std(NAtomsList)));
            
            subplot(NRows,NCols,4)
            errorbar(obj.MeanBinDist,obj.ns_Dens,obj.ns_DensUnc);
            hold on;
            errorbar(obj.MeanBinDist,obj.Sz_Dens,obj.Sz_DensUnc);
            grid on;
            xlabel('Lattice Sites')
            ylabel('Density')
            title(sprintf('DetEff = %0.2f',obj.DetectionEff));
            legend({'ns','nup_s+ndn_s'})
            
            suptitle(Title);
        end
        
%         function shrinkToSave(obj)
%             %shrink all data sets to a reasonable size so they can be
%             %saved...
%             obj.n.shrinkToSave();
%             obj.nup.shrinkToSave();
%             obj.ndn.shrinkToSave();
%         end
        
        function reconstituteAllData(obj)
            %for use reconstituting all data sets, bc correlator full
            %matrices are too large to save (~300mb each), but might need
            %these to do more processing. They can easily be obtained
            %again.
            obj.n.processData();
            obj.nup.processData();
            obj.ndn.processData();
%             obj.Up_Sx.processData();
%             obj.Down_Sx.processData();
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
            
            MI_NNCorr = squeeze(obj.ns_Corr(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            MI_NNCorrUnc = squeeze(obj.ns_CorrUnc(obj.n.CenterIndex_CorrMatrix,obj.n.CenterIndex_CorrMatrix+1,:));
            MI_NNNCorr = squeeze(obj.ns_Corr(obj.n.CenterIndex_CorrMatrix+1,obj.n.CenterIndex_CorrMatrix+1,:));
            MI_NNNCorrUnc = squeeze(obj.ns_CorrUnc(obj.n.CenterIndex_CorrMatrix+1,obj.n.CenterIndex_CorrMatrix+1,:));
            
            Mat = cat(2,...
                obj.MeanBinDist,obj.MeanBinDistUnc,...
                obj.Sz_LocalPolarization,obj.Sz_LocalPolarizationUnc,...
                obj.Sz_Up_Dens,obj.Sz_Up_DensUnc,...
                obj.Sz_Down_Dens,obj.Sz_Down_Dens,...
                SzSinglesDensity,SzSinglesDensityUnc,...
                obj.Sx_Up_Dens,obj.Sx_Up_DensUnc,...
                obj.Sx_Down_Dens,obj.Sx_Down_Dens,...
                SxSinglesDensity,SxSinglesDensityUnc,...
                obj.ns_Dens,obj.ns_DensUnc,...
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