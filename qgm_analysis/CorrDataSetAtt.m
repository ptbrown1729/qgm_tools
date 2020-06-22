classdef CorrDataSetAtt < handle
    %
    
    properties
        
        Date
        Folders
             
        %useful experimental data 
        UNominal
        LattServoV
        ScienceFieldV
        LattMonAt10V
        
        %DataFolder classes
        Ones
        Threes
        Singles
        Doubles
        
        %analysis settings
        ExcludePictures
        BinMode = 'spatial' %'spatial', 'density', or 'external'
        
        %densities corrected for efficiencies.
        OnesDensCorr
        OnesDensCorrUnc
        ThreesDensCorr
        ThreesDensCorrUnc
        SinglesDens
        SinglesDensUnc
        DoublesDensCorr
        DoublesDensCorrUnc
        
        %Correlators corrected for efficiencies
        OnesCorr
        OnesCorrUnc
        ThreesCorr
        ThreesCorrUnc
        SinglesCorr
        SinglesCorrUnc
        DoublesCorr
        DoublesCorrUnc
        
        %physical constants
        mLi = 9.9883e-27; %kg
        h = 6.6261e-34; %J.s
        hbar = 1.0546e-34; %J.s
        kb =  1.3806e-23;
        a = 1064e-9/sqrt(2); %m
        
        %positions in terms of distance grid...
        BinEdges
        BinCenters
        MeanBinDist
        MeanBinDistUnc
        
        %mean radial position of bin...
        RadialPos
        RadialPosUnc
        
        %total filling
        Filling_SandD
        Filling_SandD_Unc
        Filling_1and3
        Filling_1and3_Unc
        
        %Efficiencies
        UseEfficiencies = 1
        DetectionEff = 1
        BlowDoublesEff13 = 0.93
        RFAndBlowDoublesEff23 = 0.955
        
        %chemical potentials
        OmegaWeak
        OmegaStrong
        OmegaMean %(2pi)*490-500 Hz??
        CloudAspectRatio
        %
        R_HalfFilling
        ChemPot_0
        Mu1MinusMu3
        ChemPot
        ChemPotUnc
        
        %Other useful quantities
        SinglesFrac
        SinglesFracUnc
        
        LocalPol
        LocalPolUnc
        
        GlobalPol
        GlobalPolUnc
        
        %Non-int FG quantities
        InitFG = 0;
        NonIntFG
        
        %DQMC results      
%         DQMC2_Path = 'DQMC_T=0.5_Mu=-8to8_U=-6to-4.mat'
%         DQMC2_Path = 'DQMC_T=0.5_Mu=-8to8_U=-8to-4.mat'
        DQMC2_Path = 'DQMC_T=0.3-0.7_Mu=-8to8_U=-8to-4.mat'
        DQMC2
        
    end
    
    methods
        
        function obj = CorrDataSetAtt(DateCell,DataSets,ExcludePictures,BinEdges,AzAvgMode,DistGrid)
            %constructor...
            if exist('DateCell','var')
                %if no arguments, instantiate class but do not run this code...
                obj.initialize(DateCell,DataSets,ExcludePictures,BinEdges,AzAvgMode,DistGrid);
            end
        end
        
        function initialize(obj,DateCell,DataSets,ExcludePictures,BinEdges,AzAvgMode,DistGrid)
            try
                obj.loadDQMC();
            catch
                fprintf('failed to load DQMC data \n');
            end
            
            if ~exist('AzAvgMode','var')
                AzAvgMode = 'spatial';
            end
            obj.BinMode = AzAvgMode;
            
            if ~exist('DistGrid','var')
                DistGrid = 0;
            end
            
                obj.loadData(DateCell,DataSets,ExcludePictures,BinEdges,AzAvgMode,DistGrid)
                %do analysis...i.e. get local pol and etc.
                obj.doAnalysis();
        end
        
        function loadData(obj,DateCell,DataSets,ExcludePictures,BinEdges,AzAvgMode,DistGrid)
                            %obj = CorrDataSetAtt(DateCell,DataSets,ExcludePictures,BinEdges)
                %%%Arguments
                %DateCell = {Year Day Month}
                %DataSets = [Doubles,Threes,Ones,Singles]
                %Excluded files = {[],[],[],[]}
                %BinEdges = []
                %
                obj.Date = DateCell;
                obj.Folders = DataSets;

                DataSetList = {DateCell,DateCell,DateCell,DateCell};
                
                if length(DataSets)==4
                    for ii = 1:4
                        DataSetList{ii}{4} = DataSets(ii);
                        DataSetList{ii}{5} = 1;
                        DataSetList{ii}{6} = 1;
                    end
                elseif length(DataSets)==2
                    %probably get rid of this case...too much trouble...
                    DataSetList{1} = {};
                    DataSetList{4} = {};
                    DataSetList{2}{4} = DataSets(1);
                    DataSetList{2}{5} = 1;
                    DataSetList{2}{6} = 1;
                    DataSetList{3}{4} = DataSets(2);
                    DataSetList{3}{5} = 1;
                    DataSetList{3}{6} = 1;
                else
                    exception('Wrong number of datasets in CorrDatasetAtt.m')
                end


                if ~exist('ExcludePictures','var') || isempty(ExcludePictures)
                    ExcludePictures = {[],[],[],[]};
                end

                obj.ExcludePictures = ExcludePictures;
                obj.BinEdges = BinEdges;


                %import data sets
%                 obj.Ones = DataFolder(DataSetList{3},ExcludePictures{3},BinEdges,'density');
                

                
                obj.Ones = DataFolder(DataSetList{3},ExcludePictures{3},BinEdges,obj.BinMode,DistGrid);
                obj.CloudAspectRatio = obj.Ones.CloudAspectRatio;
                DistGrid = obj.Ones.DistGrid;
                
                obj.Doubles = DataFolder(DataSetList{1},ExcludePictures{1},BinEdges,'external',DistGrid);
                obj.Singles = DataFolder(DataSetList{4},ExcludePictures{4},BinEdges,'external',DistGrid);
                obj.Threes = DataFolder(DataSetList{2},ExcludePictures{2},BinEdges,'external',DistGrid);

                obj.BinCenters = 0.5*(BinEdges(1:end-1)+BinEdges(2:end)); %can use the center of the bins
                
%                 obj.RadialPos = transpose(obj.Ones.BinDist);
%                 obj.RadialPosUnc = transpose(obj.Ones.BinDistUnc);
              
                obj.RadialPos = transpose(obj.Ones.BinAvg);
                obj.RadialPosUnc = transpose(obj.Ones.BinUnc);

                obj.MeanBinDist = transpose(obj.Ones.BinAvg); %better to use the average distance of all points in the bin.
                obj.MeanBinDistUnc = transpose(obj.Ones.BinUnc); %and its sdm
                %transposes make same size as BinCenters...this because I
                %replaced this variable for BinCenters in most functions,
                %so they were expecting something that size. In general
                %would prefer to use nx1 vector than 1xn.
                
                 %efficiencies
                if ~obj.UseEfficiencies
                    obj.DetectionEff = 1;
                    obj.BlowDoublesEff13 = 1;
                    obj.RFAndBlowDoublesEff23 = 1;
                end

                %get densities corrected for efficiencies...should only uses
                %these and not directly the Doubles,Ones,Threes,etc.
                obj.SinglesDens = obj.Singles.Occs_AzAvg/obj.DetectionEff;
                obj.SinglesDensUnc = obj.Singles.Occs_AzAvgUnc/obj.DetectionEff;
                obj.DoublesDensCorr = obj.Doubles.Occs_AzAvg/obj.RFAndBlowDoublesEff23/obj.DetectionEff;
                obj.DoublesDensCorrUnc = obj.Doubles.Occs_AzAvgUnc/obj.RFAndBlowDoublesEff23/obj.DetectionEff;

                %1s(corrected for doubles eff) = 1s - Eff*Doubles+Eff*Doubles/Eff
                %i.e. subtract the part we measure for doubles, correct it, and add it back
                % = 1s + (1/Eff-1)*Eff*Doubles = 1s + (1/Eff-1)*(1s+3s-singles)/2
                %As expected this reduces to the measured value if Eff = 1.
                %Reason that I do this, instead of using measured doubles
                %density is that that has a different efficiency due to
                %RF/different doubles blowing.
                obj.OnesDensCorr = obj.Ones.Occs_AzAvg/obj.DetectionEff + (1/obj.BlowDoublesEff13-1)*0.5*(obj.Ones.Occs_AzAvg+obj.Threes.Occs_AzAvg-obj.Singles.Occs_AzAvg)/obj.DetectionEff;
                obj.OnesDensCorrUnc = obj.Ones.Occs_AzAvgUnc/obj.DetectionEff; %can also correct this with the above expression...
                %same thing for threes.
                obj.ThreesDensCorr = obj.Threes.Occs_AzAvg/obj.DetectionEff + (1/obj.BlowDoublesEff13-1)*0.5*(obj.Ones.Occs_AzAvg+obj.Threes.Occs_AzAvg-obj.Singles.Occs_AzAvg)/obj.DetectionEff;
                obj.ThreesDensCorrUnc = obj.Threes.Occs_AzAvgUnc/obj.DetectionEff;
                
                
                %get correlations and correct for efficiencies
                obj.OnesCorr = obj.Ones.Density_Corr_AzAvg/obj.DetectionEff^2;
                obj.OnesCorrUnc = obj.Ones.Density_Corr_AzAvgUnc/obj.DetectionEff^2;
                obj.ThreesCorr = obj.Threes.Density_Corr_AzAvg/obj.DetectionEff^2;
                obj.ThreesCorrUnc = obj.Threes.Density_Corr_AzAvgUnc/obj.DetectionEff^2;
                obj.SinglesCorr = obj.Singles.Density_Corr_AzAvg/obj.DetectionEff^2;
                obj.SinglesCorrUnc = obj.Singles.Density_Corr_AzAvgUnc/obj.DetectionEff^2;
                obj.DoublesCorr = obj.Doubles.Density_Corr_AzAvg/obj.DetectionEff^2/obj.RFAndBlowDoublesEff23^2;
                obj.DoublesCorrUnc = obj.Doubles.Density_Corr_AzAvgUnc/obj.DetectionEff^2/obj.RFAndBlowDoublesEff23^2;
                
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
            obj.ChemPot = obj.RadialPos.^2;
            
            obj.LocalPol = (obj.OnesDensCorr-obj.ThreesDensCorr)./(obj.OnesDensCorr+obj.ThreesDensCorr);
            obj.LocalPolUnc = sqrt((2*obj.OnesDensCorr./(obj.ThreesDensCorr+obj.OnesDensCorr).^2.*obj.ThreesDensCorrUnc).^2+(2*obj.ThreesDensCorr./(obj.ThreesDensCorr+obj.OnesDensCorr).^2.*obj.OnesDensCorrUnc).^2);
            
            obj.SinglesFrac = obj.SinglesDens./(obj.SinglesDens+2*obj.DoublesDensCorr);
            obj.SinglesFracUnc = sqrt((2*obj.DoublesDensCorr./(obj.SinglesDens+2*obj.DoublesDensCorr).^2.*obj.SinglesDensUnc).^2+(2*obj.SinglesDens./(obj.SinglesDens+2*obj.DoublesDensCorr).^2.*obj.DoublesDensCorrUnc).^2);
            
            obj.Filling_SandD = obj.SinglesDens+2*obj.DoublesDensCorr;
            obj.Filling_SandD_Unc = sqrt((obj.SinglesDensUnc).^2+(2*obj.DoublesDensCorrUnc).^2);
            obj.Filling_1and3 = obj.OnesDensCorr + obj.ThreesDensCorr;
            obj.Filling_1and3_Unc = sqrt((obj.OnesDensCorrUnc).^2+(obj.ThreesDensCorrUnc).^2);
            
           %can also correct these for efficiencies...
            N1 = obj.Ones.MeanAtomNum;
            N1_Unc = obj.Ones.AtomNumSD;
            N3 = obj.Threes.MeanAtomNum;
            N3_Unc = obj.Threes.AtomNumSD;
            obj.GlobalPol = (obj.Ones.MeanAtomNum-obj.Threes.MeanAtomNum)/(obj.Ones.MeanAtomNum+obj.Threes.MeanAtomNum);
            obj.GlobalPolUnc = sqrt((2*N1./(N3+N1).^2.*N3_Unc).^2+(2*N3./(N1+N3).^2.*N1_Unc).^2);
            
        end
        
        function Dstring = getDescriptionString(obj)
            Dstring = sprintf('%d_%02d_%d_Folders=%03d-%03d_Pg=%0.2f',obj.Date{1},obj.Date{2},obj.Date{3},min(obj.Folders),max(obj.Folders),obj.GlobalPol);
        end
        
        function FitParams = fitNonIntLattice_Simultaneous(obj,tHz,InitParams,FixedParams)
            %P = [Beta,Omega,Mu1,Mu2]
            %Simultaneously fit two non-interacting lattice fermi gas
            %profiles to the two spin components.
            A = 1; %density correction factor...
            t = obj.h*tHz;
            
            %Generate guesses for parameters from single fits...
            if ~exist('InitParams','var')
                Fp1 = obj.Ones.fitNonIntLattice(tHz);
                Beta1 = Fp1(1); Omega1 = Fp1(2); Mu1 = Fp1(3);
                Fp2 = obj.Threes.fitNonIntLattice(tHz);
                Beta2 = Fp2(1); Omega2 = Fp2(2); Mu2 = Fp2(3);
                Beta = 0.5*(Beta1+Beta2); %units of t
                Omega = 0.5*(Omega1+Omega2);
                InitParams = [Beta,Omega,Mu1,Mu2];
            end

            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            

            %Determine temperature, chemical potentials, and trapping
            %frequency from fit.
            %one trick to let me fix some parameters...
            OnesDensity =  @(P,X)latticeFGRadial1D([P(1),P(2),P(3),A],X);
            ThreesDensity = @(P,X) latticeFGRadial1D([P(1),P(2),P(4),A],X);
            
            %Maybe there is a better way to estimate the uncertainty of
            %these points?
            OnesUnc = transpose(obj.OnesDensCorrUnc);
            OnesUnc(OnesUnc == 0 ) = 1;
            ThreesUnc = transpose(obj.ThreesDensCorrUnc);
            ThreesUnc(ThreesUnc == 0) = 1;
            %
            OnesFit = @(P) (OnesDensity(P,obj.RadialPos)-transpose(obj.OnesDensCorr))./OnesUnc;
            ThreesFit = @(P) (ThreesDensity(P,obj.RadialPos)-transpose(obj.ThreesDensCorr))./ThreesUnc;
            FitFn = @(P) [OnesFit(P.*(1-FixedParams)+InitParams.*FixedParams),ThreesFit(P.*(1-FixedParams)+InitParams.*FixedParams)];
            
            FitParams = lsqnonlin(FitFn,InitParams); 
            
            OmegaMeanHz = FitParams(2)/sqrt(obj.mLi/abs(t))/obj.a;
            OmegaWeak = OmegaMeanHz/sqrt(obj.CloudAspectRatio);
            OmegaStrong = OmegaMeanHz*sqrt(obj.CloudAspectRatio);
            
            T = 1/FitParams(1);
            MuOnes = FitParams(3);
            MuThrees = FitParams(4);
            MuAvg = 0.5*(FitParams(3)+FitParams(4));
            DeltaMu = 0.5*(FitParams(3)-FitParams(4));
            
            %display results
            fprintf('t = %0.2f Hz \n',tHz);
            fprintf('T = %0.2f t = %0.2f nK \n',T,T*t/obj.kb/1e-9);
            fprintf('Chemical potential is zero at half filling \n')
            fprintf('Bandwidth in tight-binding model is 8t \n')
            fprintf('Mu1 = %0.2f t = (2pi) %0.2f Hz \n',MuOnes,MuOnes*t/obj.h);
            fprintf('Mu3 = %0.2f t = (2pi) %0.2f Hz \n',MuThrees,MuThrees*t/obj.h);
            fprintf('MuAvg = %0.2f t = (2pi) %0.2f Hz \n',MuAvg,MuAvg*t/obj.h);
            fprintf('DeltaMu = %0.2f t = (2pi) %0.2f Hz \n',DeltaMu,DeltaMu*t/obj.h);
            fprintf('Omega = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi));
            fprintf('Omega Strong = (2pi) %0.2f Hz \n',OmegaStrong/(2*pi));
            fprintf('Omega Weak = (2pi) %0.2f Hz \n \n',OmegaWeak/(2*pi));
            
            RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            FigName = sprintf('Simultaneous Fit Non-Int Latt FG %03d-%03d, T = %0.2f t, Omega = (2pi) %0.0f Hz',min(obj.Folders),max(obj.Folders),T,OmegaMeanHz/(2*pi));
            fh = figure('name',FigName);
            subplot(2,2,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp),'b');
                title(sprintf('Ones, Mu1 = %0.2f t',MuOnes));
                grid on;
                ylim([0,1])
            subplot(2,2,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesDensity(FitParams,RInterp),'b');
                title(sprintf('Threes, Mu3 = %0.2f t',MuThrees));
                grid on;
                ylim([0,1])
            subplot(2,2,3)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Doubles Profile')
                ylim([0,1])
            subplot(2,2,4)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp)+ThreesDensity(FitParams,RInterp)-2*OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Singles Profile')
                ylim([0,2])
             suptitle(FigName);
        end
        
         function FitParams = fitNonIntLattice_Simultaneous_TrapVariation(obj,tHz,InitParams,FixedParams)
            %P = [Beta,Omega,Mu1,Mu2,Vo,tEr]
            %Simultaneously fit two non-interacting lattice fermi gas
            %profiles to the two spin components.
            A = 1; %density correction factor...
            t = obj.h*tHz;
            
            %Generate guesses for parameters from single fits...
            if ~exist('InitParams','var')
                Fp1 = obj.Ones.fitNonIntLattice(tHz);
                Beta1 = Fp1(1); Omega1 = Fp1(2); Mu1 = Fp1(3);
                Fp2 = obj.Threes.fitNonIntLattice(tHz);
                Beta2 = Fp2(1); Omega2 = Fp2(2); Mu2 = Fp2(3);
                Beta = 0.5*(Beta1+Beta2); %units of t
                Omega = 0.5*(Omega1+Omega2);
                Vo = 120; %in t
                tEr = tHz/14.66e3;
                InitParams = [Beta,Omega,Mu1,Mu2,Vo,tEr];
            end

            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            

            %Determine temperature, chemical potentials, and trapping
            %frequency from fit.
            %one trick to let me fix some parameters...
            OnesDensity =  @(P,X)latticeFGRadial_VaryingLatticeDepth1D([P(1),P(2),P(3),P(5),P(6),A],X);
            ThreesDensity = @(P,X) latticeFGRadial_VaryingLatticeDepth1D([P(1),P(2),P(4),P(5),P(6),A],X);
            
            %Maybe there is a better way to estimate the uncertainty of
            %these points?
            OnesUnc = transpose(obj.OnesDensCorrUnc);
            OnesUnc(OnesUnc == 0 ) = 1;
            ThreesUnc = transpose(obj.ThreesDensCorrUnc);
            ThreesUnc(ThreesUnc == 0) = 1;
            %
            OnesFit = @(P) (OnesDensity(P,obj.RadialPos)-transpose(obj.OnesDensCorr))./OnesUnc;
            ThreesFit = @(P) (ThreesDensity(P,obj.RadialPos)-transpose(obj.ThreesDensCorr))./ThreesUnc;
            FitFn = @(P) [OnesFit(P.*(1-FixedParams)+InitParams.*FixedParams),ThreesFit(P.*(1-FixedParams)+InitParams.*FixedParams)];

            
            FitParams = lsqnonlin(FitFn,InitParams);
            
            
            OmegaMeanHz = FitParams(2)/sqrt(obj.mLi/abs(t))/obj.a;
            OmegaWeak = OmegaMeanHz/sqrt(obj.CloudAspectRatio);
            OmegaStrong = OmegaMeanHz*sqrt(obj.CloudAspectRatio);
            
            T = 1/FitParams(1);
            MuOnes = FitParams(3);
            MuThrees = FitParams(4);
            MuAvg = 0.5*(FitParams(3)+FitParams(4));
            DeltaMu = 0.5*(FitParams(3)-FitParams(4));
            Vo = FitParams(5);
            tEr = FitParams(6);
            
            %display results
            fprintf('t = %0.2f Hz \n',tHz);
            fprintf('T = %0.2f t = %0.2f nK \n',T,T*t/obj.kb/1e-9);
            fprintf('Chemical potential is zero at half filling \n')
            fprintf('Bandwidth in tight-binding model is 8t \n')
            fprintf('Mu1 = %0.2f t = (2pi) %0.2f Hz \n',MuOnes,MuOnes*t/obj.h);
            fprintf('Mu3 = %0.2f t = (2pi) %0.2f Hz \n',MuThrees,MuThrees*t/obj.h);
            fprintf('MuAvg = %0.2f t = (2pi) %0.2f Hz \n',MuAvg,MuAvg*t/obj.h);
            fprintf('DeltaMu = %0.2f t = (2pi) %0.2f Hz \n',DeltaMu,DeltaMu*t/obj.h);
            fprintf('Omega = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi));
            fprintf('Omega Strong = (2pi) %0.2f Hz \n',OmegaStrong/(2*pi));
            fprintf('Omega Weak = (2pi) %0.2f Hz \n',OmegaWeak/(2*pi));
            fprintf('Central Lattice Depth = %0.2f t = %0.2f Er \n',Vo,Vo*tEr);
            fprintf('Central Hopping t = %0.2f Er = %0.2f Hz \n \n',tEr,tEr*14.66e3);
            
            RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            FigName = sprintf('Simultaneous Fit Non-Int Latt FG %03d-%03d, T = %0.2f t, Omega = (2pi) %0.0f Hz',min(obj.Folders),max(obj.Folders),T,OmegaMeanHz/(2*pi));
            fh = figure('name',FigName);
            subplot(2,2,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp),'b');
                title(sprintf('Ones, Mu1 = %0.2f t',MuOnes));
                grid on;
                ylim([0,1])
            subplot(2,2,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesDensity(FitParams,RInterp),'b');
                title(sprintf('Threes, Mu3 = %0.2f t',MuThrees));
                grid on;
                ylim([0,1])
            subplot(2,2,3)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Doubles Profile')
                ylim([0,1])
            subplot(2,2,4)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp)+ThreesDensity(FitParams,RInterp)-2*OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Singles Profile')
                ylim([0,2])
             suptitle(FigName);
        end
        
         function FitParams = fitNonIntLattice_AndCorr_Simultaneous(obj,tHz,InitParams,FixedParams)
            %P = [Beta,Omega,Mu1,Mu2]
            %Simultaneously fit two non-interacting lattice fermi gas
            %profiles to the two spin components.
            A = 1; %density correction factor...
            t = obj.h*tHz;
            
            %Generate guesses for parameters from single fits...
            if ~exist('InitParams','var')
                Fp1 = obj.Ones.fitNonIntLattice(tHz);
                Beta1 = Fp1(1); Omega1 = Fp1(2); Mu1 = Fp1(3);
                Fp2 = obj.Threes.fitNonIntLattice(tHz);
                Beta2 = Fp2(1); Omega2 = Fp2(2); Mu2 = Fp2(3);
                Beta = 0.5*(Beta1+Beta2); %units of t
                Omega = 0.5*(Omega1+Omega2);
                InitParams = [Beta,Omega,Mu1,Mu2];
            end

            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            

            %Determine temperature, chemical potentials, and trapping
            %frequency from fit.
            %one trick to let me fix some parameters...
            OnesDensity =  @(P,X)latticeFGRadial1D([P(1),P(2),P(3),A],X);
            ThreesDensity = @(P,X) latticeFGRadial1D([P(1),P(2),P(4),A],X);
            OnesCorr = @(P,X) latticeFGRadial_NNCorr1D([P(1),P(2),P(3),A^2],X);
            ThreesCorr = @(P,X) latticeFGRadial_NNCorr1D([P(1),P(2),P(4),A^2],X);
            
            %Maybe there is a better way to estimate the uncertainty of
            %these points?
            OnesDensUnc = transpose(obj.OnesDensCorrUnc);
            OnesDensUnc(OnesDensUnc == 0 ) = 1;
            ThreesDensUnc = transpose(obj.ThreesDensCorrUnc);
            ThreesDensUnc(ThreesDensUnc == 0) = 1;
            %correlator uncertainties
            OnesCorrUnc = squeeze(obj.OnesCorrUnc(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            OnesCorrUnc(OnesCorrUnc == 0) = 1;
            ThreesCorrUnc = squeeze(obj.ThreesCorrUnc(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
            ThreesCorrUnc(ThreesCorrUnc == 0) = 1;
            %correlators
            OnesCorrExp = squeeze(obj.OnesCorr(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            ThreesCorrExp = squeeze(obj.ThreesCorr(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
            
            %
            OnesDensFit = @(P) (OnesDensity(P,obj.RadialPos)-transpose(obj.OnesDensCorr))./OnesDensUnc;
            ThreesDensFit = @(P) (ThreesDensity(P,obj.RadialPos)-transpose(obj.ThreesDensCorr))./ThreesDensUnc;
            OnesCorrFit = @(P) (OnesCorr(P,obj.RadialPos) - transpose(OnesCorrExp))./OnesDensUnc;
            ThreesCorrFit = @(P) (ThreesCorr(P,obj.RadialPos) - transpose(ThreesCorrExp))./ThreesDensUnc;
            
            DensWeight = 1;
            CorrWeight = 1; %1/0.04;
            
            FitFn = @(P) [DensWeight*OnesDensFit(P.*(1-FixedParams)+InitParams.*FixedParams),DensWeight*ThreesDensFit(P.*(1-FixedParams)+InitParams.*FixedParams),CorrWeight*OnesCorrFit(P.*(1-FixedParams)+InitParams.*FixedParams),CorrWeight*ThreesCorrFit(P.*(1-FixedParams)+InitParams.*FixedParams)];

            
            FitParams = lsqnonlin(FitFn,InitParams);
            
            
            OmegaMeanHz = FitParams(2)/sqrt(obj.mLi/abs(t))/obj.a;
            OmegaWeak = OmegaMeanHz/sqrt(obj.CloudAspectRatio);
            OmegaStrong = OmegaMeanHz*sqrt(obj.CloudAspectRatio);
            
            T = 1/FitParams(1);
            MuOnes = FitParams(3);
            MuThrees = FitParams(4);
            MuAvg = 0.5*(FitParams(3)+FitParams(4));
            DeltaMu = 0.5*(FitParams(3)-FitParams(4));
            
            %display results
            fprintf('t = %0.2f Hz \n',tHz);
            fprintf('T = %0.2f t = %0.2f nK \n',T,T*t/obj.kb/1e-9);
            fprintf('Chemical potential is zero at half filling \n')
            fprintf('Bandwidth in tight-binding model is 8t \n')
            fprintf('Mu1 = %0.2f t = (2pi) %0.2f Hz \n',MuOnes,MuOnes*t/obj.h);
            fprintf('Mu3 = %0.2f t = (2pi) %0.2f Hz \n',MuThrees,MuThrees*t/obj.h);
            fprintf('MuAvg = %0.2f t = (2pi) %0.2f Hz \n',MuAvg,MuAvg*t/obj.h);
            fprintf('DeltaMu = %0.2f t = (2pi) %0.2f Hz \n',DeltaMu,DeltaMu*t/obj.h);
            fprintf('Omega = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi));
            fprintf('Omega Strong = (2pi) %0.2f Hz \n',OmegaStrong/(2*pi));
            fprintf('Omega Weak = (2pi) %0.2f Hz \n \n',OmegaWeak/(2*pi));
            
            RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            FigName = sprintf('Simultaneous Fit Non-Int Latt FG %03d-%03d, T = %0.2f t, Omega = (2pi) %0.0f Hz',min(obj.Folders),max(obj.Folders),T,OmegaMeanHz/(2*pi));
            fh = figure('name',FigName);
            NRows = 2;
            NCols = 4;
            subplot(NRows,NCols,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp),'b');
                title(sprintf('Ones, Mu1 = %0.2f t',MuOnes));
                grid on;
                ylim([0,1])
            subplot(NRows,NCols,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesDensity(FitParams,RInterp),'b');
                title(sprintf('Threes, Mu3 = %0.2f t',MuThrees));
                grid on;
                ylim([0,1])
            subplot(NRows,NCols,5)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Doubles Profile')
                ylim([0,1])
            subplot(NRows,NCols,6)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp)+ThreesDensity(FitParams,RInterp)-2*OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Singles Profile')
                ylim([0,2])
            %correlators.
            subplot(NRows,NCols,3)
                errorbar(obj.RadialPos,OnesCorrExp,OnesCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesCorr(FitParams,RInterp),'b');
                ylim([-0.05,0])
                title('Ones, NN Corr');
                grid on;
            subplot(NRows,NCols,4)
                errorbar(obj.RadialPos,ThreesCorrExp,ThreesCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesCorr(FitParams,RInterp),'b');
                title('Threes, NN Corr');
                ylim([-0.05,0])
                grid on;
            subplot(NRows,NCols,7)
                DoublesCorr = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                DoublesCorrUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,DoublesCorr,DoublesCorrUnc,'ro')
                hold on;
                %some combinatorial Wick's thrm manipulation leads to...
                %<d_i d_j>_c = n_up^2*<n_up_i n_up_j>_c + n_down^2*<n_down_i n_down_j>_c + n_down^2*n_up^2
                DoublesCorrInferred = OnesCorr(FitParams,RInterp).*ThreesDensity(FitParams,RInterp).^2 + ThreesCorr(FitParams,RInterp).*OnesDensity(FitParams,RInterp).^2 + OnesCorr(FitParams,RInterp).*ThreesCorr(FitParams,RInterp);
                plot(RInterp,DoublesCorrInferred);
                title('Inferred Doubles Correlator')
                grid on;
            subplot(NRows,NCols,8)
                SinglesCorr = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                SinglesCorrUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,SinglesCorr,SinglesCorrUnc,'ro');
                hold on;
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
                SinglesCorrInferred = OnesCorr(FitParams,RInterp)+ThreesCorr(FitParams,RInterp) + 4*DoublesCorrInferred -4*ThreesCorr(FitParams,RInterp).*OnesDensity(FitParams,RInterp) - 4*OnesCorr(FitParams,RInterp).*ThreesDensity(FitParams,RInterp);
                plot(RInterp,SinglesCorrInferred);
                grid on;
                title('Inferred Singles Correlator');

             suptitle(FigName);
         end
         
          function FitParams = fitNonIntLattice_AndCorr_Simultaneous_TrapVariation(obj,tHz,InitParams,FixedParams)
            %P = [Beta,Omega,Mu1,Mu2]
            %Simultaneously fit two non-interacting lattice fermi gas
            %profiles to the two spin components.
            A = 1; %density correction factor...
            t = obj.h*tHz;
            
            %Generate guesses for parameters from single fits...
            if ~exist('InitParams','var')
                Fp1 = obj.Ones.fitNonIntLattice(tHz);
                Beta1 = Fp1(1); Omega1 = Fp1(2); Mu1 = Fp1(3);
                Fp2 = obj.Threes.fitNonIntLattice(tHz);
                Beta2 = Fp2(1); Omega2 = Fp2(2); Mu2 = Fp2(3);
                Beta = 0.5*(Beta1+Beta2); %units of t
                Omega = 0.5*(Omega1+Omega2);
                Vo = 120; %in t
                tEr = tHz/14.66e3;
                InitParams = [Beta,Omega,Mu1,Mu2,Vo,tEr];
            end

            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            

            %Determine temperature, chemical potentials, and trapping
            %frequency from fit.
            %one trick to let me fix some parameters...
            OnesDensity =  @(P,X)latticeFGRadial_VaryingLatticeDepth1D([P(1),P(2),P(3),P(5),P(6),A],X);
            ThreesDensity = @(P,X) latticeFGRadial_VaryingLatticeDepth1D([P(1),P(2),P(4),P(5),P(6),A],X);
            OnesCorr = @(P,X) latticeFGRadial_NNCorr_VaryingLatticeDepth1D([P(1),P(2),P(3),P(5),P(6),A^2],X);
            ThreesCorr = @(P,X) latticeFGRadial_NNCorr_VaryingLatticeDepth1D([P(1),P(2),P(4),P(5),P(6),A^2],X);
            
            %Maybe there is a better way to estimate the uncertainty of
            %these points?
            OnesDensUnc = transpose(obj.OnesDensCorrUnc);
            OnesDensUnc(OnesDensUnc == 0 ) = 1;
            ThreesDensUnc = transpose(obj.ThreesDensCorrUnc);
            ThreesDensUnc(ThreesDensUnc == 0) = 1;
            %correlator uncertainties
            OnesCorrUnc = squeeze(obj.OnesCorrUnc(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            OnesCorrUnc(OnesCorrUnc == 0) = 1;
            ThreesCorrUnc = squeeze(obj.ThreesCorrUnc(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
            ThreesCorrUnc(ThreesCorrUnc == 0) = 1;
            %correlators
            OnesCorrExp = squeeze(obj.OnesCorr(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            ThreesCorrExp = squeeze(obj.ThreesCorr(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
            
            %
            OnesDensFit = @(P) (OnesDensity(P,obj.RadialPos)-transpose(obj.OnesDensCorr))./OnesDensUnc;
            ThreesDensFit = @(P) (ThreesDensity(P,obj.RadialPos)-transpose(obj.ThreesDensCorr))./ThreesDensUnc;
            OnesCorrFit = @(P) (OnesCorr(P,obj.RadialPos) - transpose(OnesCorrExp))./OnesDensUnc;
            ThreesCorrFit = @(P) (ThreesCorr(P,obj.RadialPos) - transpose(ThreesCorrExp))./ThreesDensUnc;
            
            DensWeight = 1;
            CorrWeight = 1; %1/0.04;
            
            FitFn = @(P) [DensWeight*OnesDensFit(P.*(1-FixedParams)+InitParams.*FixedParams),DensWeight*ThreesDensFit(P.*(1-FixedParams)+InitParams.*FixedParams),CorrWeight*OnesCorrFit(P.*(1-FixedParams)+InitParams.*FixedParams),CorrWeight*ThreesCorrFit(P.*(1-FixedParams)+InitParams.*FixedParams)];

            
            FitParams = lsqnonlin(FitFn,InitParams);
            
            
            OmegaMeanHz = FitParams(2)/sqrt(obj.mLi/abs(t))/obj.a;
            OmegaWeak = OmegaMeanHz/sqrt(obj.CloudAspectRatio);
            OmegaStrong = OmegaMeanHz*sqrt(obj.CloudAspectRatio);
            
            T = 1/FitParams(1);
            MuOnes = FitParams(3);
            MuThrees = FitParams(4);
            MuAvg = 0.5*(FitParams(3)+FitParams(4));
            DeltaMu = 0.5*(FitParams(3)-FitParams(4));
            Vo = FitParams(5);
            tEr = FitParams(6);
            
            %display results
            fprintf('t = %0.2f Hz \n',tHz);
            fprintf('T = %0.2f t = %0.2f nK \n',T,T*t/obj.kb/1e-9);
            fprintf('Chemical potential is zero at half filling \n')
            fprintf('Bandwidth in tight-binding model is 8t \n')
            fprintf('Mu1 = %0.2f t = (2pi) %0.2f Hz \n',MuOnes,MuOnes*t/obj.h);
            fprintf('Mu3 = %0.2f t = (2pi) %0.2f Hz \n',MuThrees,MuThrees*t/obj.h);
            fprintf('MuAvg = %0.2f t = (2pi) %0.2f Hz \n',MuAvg,MuAvg*t/obj.h);
            fprintf('DeltaMu = %0.2f t = (2pi) %0.2f Hz \n',DeltaMu,DeltaMu*t/obj.h);
            fprintf('Omega = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi));
            fprintf('Omega Strong = (2pi) %0.2f Hz \n',OmegaStrong/(2*pi));
            fprintf('Omega Weak = (2pi) %0.2f Hz \n',OmegaWeak/(2*pi));
            fprintf('Central Lattice Depth = %0.2f t = %0.2f Er \n',Vo,Vo*tEr);
            fprintf('Central Hopping t = %0.2f Er = %0.2f Hz \n \n',tEr,tEr*14.66e3);
            
            RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            FigName = sprintf('Simultaneous Fit Non-Int Latt FG %03d-%03d, T = %0.2f t, Omega = (2pi) %0.0f Hz',min(obj.Folders),max(obj.Folders),T,OmegaMeanHz/(2*pi));
            fh = figure('name',FigName);
            NRows = 2;
            NCols = 4;
            subplot(NRows,NCols,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp),'b');
                title(sprintf('Ones, Mu1 = %0.2f t',MuOnes));
                grid on;
                ylim([0,1])
            subplot(NRows,NCols,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesDensity(FitParams,RInterp),'b');
                title(sprintf('Threes, Mu3 = %0.2f t',MuThrees));
                grid on;
                ylim([0,1])
            subplot(NRows,NCols,5)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Doubles Profile')
                ylim([0,1])
            subplot(NRows,NCols,6)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp)+ThreesDensity(FitParams,RInterp)-2*OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Singles Profile')
                ylim([0,2])
            %correlators.
            subplot(NRows,NCols,3)
                errorbar(obj.RadialPos,OnesCorrExp,OnesCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesCorr(FitParams,RInterp),'b');
                ylim([-0.05,0])
                title('Ones, NN Corr');
                grid on;
            subplot(NRows,NCols,4)
                errorbar(obj.RadialPos,ThreesCorrExp,ThreesCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesCorr(FitParams,RInterp),'b');
                title('Threes, NN Corr');
                ylim([-0.05,0])
                grid on;
            subplot(NRows,NCols,7)
                DoublesCorr = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                DoublesCorrUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,DoublesCorr,DoublesCorrUnc,'ro')
                hold on;
                %some combinatorial Wick's thrm manipulation leads to...
                %<d_i d_j>_c = n_up^2*<n_up_i n_up_j>_c + n_down^2*<n_down_i n_down_j>_c + n_down^2*n_up^2
                DoublesCorrInferred = OnesCorr(FitParams,RInterp).*ThreesDensity(FitParams,RInterp).^2 + ThreesCorr(FitParams,RInterp).*OnesDensity(FitParams,RInterp).^2 + OnesCorr(FitParams,RInterp).*ThreesCorr(FitParams,RInterp);
                plot(RInterp,DoublesCorrInferred);
                title('Inferred Doubles Correlator')
                grid on;
            subplot(NRows,NCols,8)
                SinglesCorr = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                SinglesCorrUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,SinglesCorr,SinglesCorrUnc,'ro');
                hold on;
                %This correlator is even worse combinatorally...
                SinglesCorrInferred = OnesCorr(FitParams,RInterp)+ThreesCorr(FitParams,RInterp) + 4*DoublesCorrInferred -4*ThreesCorr(FitParams,RInterp).*OnesDensity(FitParams,RInterp) - 4*OnesCorr(FitParams,RInterp).*ThreesDensity(FitParams,RInterp);
                plot(RInterp,SinglesCorrInferred);
                grid on;
                title('Inferred Singles Correlator');

             suptitle(FigName);
         end
        
         function showNonIntLattice_Inferred(obj)
         %show inferred values for singles and doubles density and
         %correlators for a non-interacting gas, using the measured values
         %of ones and threes density and correlators.
          
         %correlator uncertainties
            OnesCorrUnc = squeeze(obj.OnesCorrUnc(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            OnesCorrUnc(OnesCorrUnc == 0) = 1;
            ThreesCorrUnc = squeeze(obj.ThreesCorrUnc(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
            ThreesCorrUnc(ThreesCorrUnc == 0) = 1;
            %correlators
            OnesCorrExp = squeeze(obj.OnesCorr(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            ThreesCorrExp = squeeze(obj.ThreesCorr(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
      
         
         OnesDensity = @(X) interp1(obj.RadialPos,obj.OnesDensCorr,X);
         ThreesDensity = @(X) interp1(obj.RadialPos,obj.ThreesDensCorr,X);
         OnesCorr = @(X) interp1(obj.RadialPos,OnesCorrExp,X);
         ThreesCorr = @(X) interp1(obj.RadialPos,ThreesCorrExp,X);
         
         
         RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            FigName = sprintf('Non-Int Latt Inferred Doubles and Singles, FG %03d-%03d',min(obj.Folders),max(obj.Folders));
            FigHandle = figure('name',FigName);
            NRows = 2;
            NCols = 4;
            subplot(NRows,NCols,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                title(sprintf('Ones'));
                grid on;
                ylim([0,1])
            subplot(NRows,NCols,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                title(sprintf('Threes'));
                grid on;
                ylim([0,1])
            subplot(NRows,NCols,5)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(RInterp).*ThreesDensity(RInterp),'b');
                grid on;
                title('Inferred Doubles Profile')
                ylim([0,1])
            subplot(NRows,NCols,6)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(RInterp)+ThreesDensity(RInterp)-2*OnesDensity(RInterp).*ThreesDensity(RInterp),'b');
                grid on;
                title('Inferred Singles Profile')
                ylim([0,2])
            %correlators.
            subplot(NRows,NCols,3)
                errorbar(obj.RadialPos,OnesCorrExp,OnesCorrUnc,'ro');
                ylim([-0.05,0])
                title('Ones, NN Corr');
                grid on;
            subplot(NRows,NCols,4)
                errorbar(obj.RadialPos,ThreesCorrExp,ThreesCorrUnc,'ro');
                title('Threes, NN Corr');
                ylim([-0.05,0])
                grid on;
            subplot(NRows,NCols,7)
                DoublesCorr = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                DoublesCorrUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,DoublesCorr,DoublesCorrUnc,'ro')
                hold on;
                %some combinatorial Wick's thrm manipulation leads to...
                %<d_i d_j>_c = n_up^2*<n_up_i n_up_j>_c + n_down^2*<n_down_i n_down_j>_c + n_down^2*n_up^2
                DoublesCorrInferred = OnesCorr(RInterp).*ThreesDensity(RInterp).^2 + ThreesCorr(RInterp).*OnesDensity(RInterp).^2 + OnesCorr(RInterp).*ThreesCorr(RInterp);
                plot(RInterp,DoublesCorrInferred);
                title('Inferred Doubles Correlator')
                grid on;
            subplot(NRows,NCols,8)
                SinglesCorr = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                SinglesCorrUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,SinglesCorr,SinglesCorrUnc,'ro');
                hold on;
                %see fitNonIntLattice_AndCorr_Simultaneous fn for
                %explanation of this correlator...
                SinglesCorrInferred = OnesCorr(RInterp)+ThreesCorr(RInterp) + 4*DoublesCorrInferred -4*ThreesCorr(RInterp).*OnesDensity(RInterp) - 4*OnesCorr(RInterp).*ThreesDensity(RInterp);
                plot(RInterp,SinglesCorrInferred);
                grid on;
                title('Inferred Singles Correlator');

             suptitle(FigName);
         
         end
         
         function FitParams = fitNonIntLattice_AndCorr_NoTrap(obj,InitParams,FixedParams)
             %FitParams = fitNonIntLattice_AndCorr_NoTrap(obj,InitParams,FixedParams)
             %FitParams = [DeltaMu,T]
             %use class instance of NonIntFG to fit this.
             
             if ~exist('InitParams','var')
                InitParams = [2,0.3];
            end
            
            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            
            if ~obj.InitFG
                obj.NonIntFG = NonIntFG(1);
            end
            
            
                        %experimental quantities...
            Density = 0.5*(obj.Filling_1and3+obj.Filling_SandD);
            LocalPol = obj.LocalPol;
            SpinDiffExp = obj.OnesDensCorr-obj.ThreesDensCorr;
            
            SpinDiffUnc = sqrt((obj.OnesDensCorrUnc).^2+(obj.ThreesDensCorrUnc).^2);
            SpinDiffUnc(isnan(SpinDiffUnc)) = 1;
            SpinDiffUnc(SpinDiffUnc == 0) = 1;
            
            LocalPol(isnan(LocalPol)) = 1;
            PolUnc = obj.LocalPolUnc;
            PolUnc(isnan(PolUnc)) = 1;
            PolUnc(PolUnc == 0 ) = 1;
            
            OnesCorrUnc = squeeze(obj.OnesCorrUnc(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            OnesCorrUnc(OnesCorrUnc == 0) = 1;
            ThreesCorrUnc = squeeze(obj.ThreesCorrUnc(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
            ThreesCorrUnc(ThreesCorrUnc == 0) = 1;
            DoublesCorrUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
            %correlators
            OnesCorrExp = squeeze(obj.OnesCorr(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
            ThreesCorrExp = squeeze(obj.ThreesCorr(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
            DoublesCorrExp = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExp = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
            
            %initial approach turned out to be quite slow.
%             %non int fermi gas matrices...
%             [ns,delmus,Ts,Pols,Corr1s,Corr2s,SpinDiff] = obj.generateNonIntData();
%             
%             %make functions from these matrices...will fit to these.
%             %note that these functions require all three arguments be the
%             %same size...so you can't have one be vectorized and the other
%             %two be singles.
%             %couldn't figure out a way to get rid of Nans in anonymous
%             %functions...so added some ugly fns to class.
%             pfn = @(n,dmu,T) obj.nonIntFg_PolFn(ns,delmus,Ts,Pols,n,dmu,T);
%             SpinDiffFn = @(n,dmu,T) obj.nonIntFg_SpinDiffFn(ns,delmus,Ts,SpinDiff,n,dmu,T);
%             C1fn = @(n,dmu,T) obj.nonIntFg_CorrFn(ns,delmus,Ts,Corr1s,n,dmu,T);
%             C2fn = @(n,dmu,T) obj.nonIntFg_CorrFn(ns,delmus,Ts,Corr2s,n,dmu,T);
            pfn = @(n,dmu,T) obj.NonIntFG.pfn_n_dmu_T(n,dmu,T);
            SpinDiffFn = @(n,dmu,T) obj.NonIntFG.n1fn_n_dmu_T(n,dmu,T)-obj.NonIntFG.n2fn_n_dmu_T(n,dmu,T);
            C1fn = @(n,dmu,T) obj.NonIntFG.c1fn_n_dmu_T(n,dmu,T);
            C2fn = @(n,dmu,T) obj.NonIntFG.c2fn_n_dmu_T(n,dmu,T);
            

            %correct here...
            Pfit = @(P)(LocalPol - pfn(Density,P(1)*ones(size(Density)),P(2)*ones(size(Density))))./PolUnc;
            Difffit = @(P) (SpinDiffExp - SpinDiffFn(Density,P(1)*ones(size(Density)),P(2)*ones(size(Density))))./SpinDiffUnc;
            C1fit = @(P)(OnesCorrExp-C1fn(Density,P(1)*ones(size(Density)),P(2)*ones(size(Density))))./OnesCorrUnc;
            C3fit = @(P)(ThreesCorrExp-C2fn(Density,P(1)*ones(size(Density)),P(2)*ones(size(Density))))./ThreesCorrUnc;
            
            FitFn = @(P) [Difffit(P.*(1-FixedParams)+InitParams.*FixedParams),C1fit(P.*(1-FixedParams)+InitParams.*FixedParams),C3fit(P.*(1-FixedParams)+InitParams.*FixedParams)];
            FitParams = lsqnonlin(FitFn,InitParams);

            FigHandle = figure('name','Non-Int Gas Fit, No Trap')
            NRows = 2;
            NCols = 3;
            
            subplot(NRows,NCols,1)
            errorbar(Density,SpinDiffExp,SpinDiffUnc,'r-o');
            hold on;
            DensInterp = linspace(min(Density),max(Density),300);
            plot(DensInterp,SpinDiffFn(DensInterp,FitParams(1)*ones(size(DensInterp)),FitParams(2)*ones(size(DensInterp))),'b');
            grid on;
            xlim([0,2]);
            ylim([0,1]);
            xlabel('Density')
            ylabel('n1-n3')
            title('Fit n1-n3');
            
            subplot(NRows,NCols,2)
            errorbar(Density,OnesCorrExp,OnesCorrUnc,'r-o')
            hold on;
            plot(DensInterp,C1fn(DensInterp,FitParams(1)*ones(size(DensInterp)),FitParams(2)*ones(size(DensInterp))),'b');
            grid on;
            ylim([1.5*min(OnesCorrExp),0])
            xlabel('Density')
            ylabel('Fit NN Corr, State 1')
            title('Fit NN Corr, State 1');
            
            subplot(NRows,NCols,3)
            errorbar(Density,ThreesCorrExp,ThreesCorrUnc,'r-o')
            hold on;
            plot(DensInterp,C2fn(DensInterp,FitParams(1)*ones(size(DensInterp)),FitParams(2)*ones(size(DensInterp))),'b');
            ylim([1.5*min(ThreesCorrExp),0])
            grid on;
            xlabel('Density')
            ylabel('NN Corr, State 3')
            title('Fit NN Corr, State 3');
           
            subplot(NRows,NCols,4)
            errorbar(Density,LocalPol,PolUnc,'r-o')
            hold on;
            plot(DensInterp,pfn(DensInterp,FitParams(1),FitParams(2)),'b');
            ylim([0,1])
            grid on;
            xlabel('Density')
            ylabel('Polarization')
            title('Inferred Polarization');
            
            
            subplot(NRows,NCols,5)
            errorbar(Density,DoublesCorrExp,DoublesCorrUnc,'r-o');
            hold on;
            plot(DensInterp,obj.NonIntFG.cdfn_n_dmu_T(DensInterp,FitParams(1),FitParams(2)));
            xlabel('Density')
            ylabel('NN Corr, Doubles');
            title('Inferred NN Corr, Doubles');
            grid on;
            
            subplot(NRows,NCols,6)
            errorbar(Density,SinglesCorrExp,SinglesCorrUnc,'r-o');
            hold on;
            plot(DensInterp,obj.NonIntFG.csfn_n_dmu_T(DensInterp,FitParams(1),FitParams(2)));
            xlabel('Density')
            ylabel('NN Corr, Singles');
            title('Inferred NN Corr, Singles');
            grid on;
            

            suptitle(sprintf('Folders %03d-%03d, Pg = %0.2f, (Mu1-Mu3)/2 = %0.2f t, T = %0.2f t',min(obj.Folders),max(obj.Folders),obj.GlobalPol,FitParams(1),FitParams(2)));

            
         end
         
        function FitParams = fitNonIntLattice_Simultaneous_TwoPancakes(obj,tHz,InitParams,FixedParams)
            %P = [Beta,Omega,Mu1A,Mu2A,Mu1B,Mu2B]
            %Simultaneously fit two non-interacting lattice fermi gas
            %profiles to the two spin components.
            A = 1; %density correction factor...
            t = obj.h*tHz;
            if ~exist('InitParams','var')
                Mu1A = 0;
                Mu2A = 0;
                Mu1B = -2;
                Mu2B = -10;
                Beta = 2; %units of t
                Omega = 2*pi*170*sqrt(obj.mLi/t)*obj.a; %units of sqrt(m/t)*a
                InitParams = [Beta,Omega,Mu1A,Mu2A,Mu1B,Mu2B];
            end
            
            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            
            if isequal(size(FixedParams),size(InitParams))
                FixedParams = zeros(size(InitParams));
            end
            

            %Determine temperature, chemical potentials, and trapping
            %frequency from fit.
            %For pancakes A and B
            %3s = 3sA+3sB - 2*3sA*3sB ... i.e. get all single threes and
            %lose all doubles threes

            OnesDensity =  @(P,X)latticeFGRadial1D([P(1),P(2),P(3),A],X)+latticeFGRadial1D([P(1),P(2),P(5),A],X) - 2*latticeFGRadial1D([P(1),P(2),P(3),A],X).*latticeFGRadial1D([P(1),P(2),P(5),A],X);
            ThreesDensity = @(P,X) latticeFGRadial1D([P(1),P(2),P(4),A],X); %+latticeFGRadial1D([P(1),P(2),P(6),A],X) -2*latticeFGRadial1D([P(1),P(2),P(4),A],X).*latticeFGRadial1D([P(1),P(2),P(6),A],X);
            
                       
            OnesUnc = transpose(obj.OnesDensCorrUnc);
            OnesUnc(OnesUnc == 0 ) = 1;
            ThreesUnc = transpose(obj.ThreesDensCorrUnc);
            ThreesUnc(ThreesUnc == 0) = 1;
            
            OnesFit = @(P) (OnesDensity(P,obj.RadialPos)-transpose(obj.OnesDensCorr))./OnesUnc;
            ThreesFit = @(P) (ThreesDensity(P,obj.RadialPos)-transpose(obj.ThreesDensCorr))./ThreesUnc;
            
           FitFn = @(P) [OnesFit(P.*(1-FixedParams)+InitParams.*FixedParams),ThreesFit(P.*(1-FixedParams)+InitParams.*FixedParams)];

            
            FitParams = lsqnonlin(FitFn,InitParams);
            
            OmegaMeanHz = FitParams(2)/sqrt(obj.mLi/abs(t))/obj.a;

            
            %display results
            fprintf('t = %0.2f Hz \n',tHz);
            fprintf('T = %0.2f t = %0.2f nK \n',1/FitParams(1),t/FitParams(1)/obj.kb/1e-9);
            fprintf('Chemical potential is zero at half filling \n')
            fprintf('Bandwidth in tight-binding model is 8t \n')
            fprintf('Mu1A = %0.2f t = (2pi) %0.2f Hz \n',FitParams(3),FitParams(3)*t/obj.h);
            fprintf('Mu3A = %0.2f t = (2pi) %0.2f Hz \n',FitParams(4),FitParams(4)*t/obj.h);
            fprintf('Mu1B = %0.2f t = (2pi) %0.2f Hz \n',FitParams(5),FitParams(5)*t/obj.h);
            fprintf('Mu3B = %0.2f t = (2pi) %0.2f Hz \n',FitParams(6),FitParams(6)*t/obj.h);
%             fprintf('MuAvg = %0.2f t = (2pi) %0.2f Hz \n',0.5*(FitParams(3)+FitParams(4)),0.5*(FitParams(3)+FitParams(4))*t/obj.h);
%             fprintf('DeltaMu = %0.2f t = (2pi) %0.2f Hz \n',0.5*(FitParams(3)-FitParams(4)),0.5*(FitParams(3)-FitParams(4))*t/obj.h);
            fprintf('Omega = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi));
            fprintf('Omega Strong = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi)*sqrt(obj.CloudAspectRatio));
            fprintf('Omega Weak = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi)/sqrt(obj.CloudAspectRatio));
                     
            RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            fh = figure('name','Simultaneous fit to non-interacting lattice fermi gas. Two pancakes.');
            subplot(2,2,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp),'b');
                plot(RInterp,latticeFGRadial1D([FitParams(1),FitParams(2),FitParams(3),A],RInterp),'g');
                plot(RInterp,latticeFGRadial1D([FitParams(1),FitParams(2),FitParams(5),A],RInterp),'m');
                title('Ones')
                grid on;
                ylim([0,1])
                legend({'Expt','Total','PCakeA','PCakeB'});
            subplot(2,2,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesDensity(FitParams,RInterp),'b');
                plot(RInterp,latticeFGRadial1D([FitParams(1),FitParams(2),FitParams(4),A],RInterp),'g');
                plot(RInterp,latticeFGRadial1D([FitParams(1),FitParams(2),FitParams(6),A],RInterp),'m');
                title('Threes')
                grid on;
                ylim([0,1])
            subplot(2,2,3)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                %Doubles = Doubles_PcakeA + Doubles_PCakeB -Doubles_PcakeA*DoublesPcakeB
                P = FitParams;
                Dbles = latticeFGRadial1D([P(1),P(2),P(3),A],RInterp).*latticeFGRadial1D([P(1),P(2),P(4),A],RInterp) + ...
                    latticeFGRadial1D([P(1),P(2),P(5),A],RInterp).*latticeFGRadial1D([P(1),P(2),P(6),A],RInterp) ...
                    - 2*latticeFGRadial1D([P(1),P(2),P(3),A],RInterp).*latticeFGRadial1D([P(1),P(2),P(4),A],RInterp).*latticeFGRadial1D([P(1),P(2),P(5),A],RInterp).*latticeFGRadial1D([P(1),P(2),P(6),A],RInterp);
%                 OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp) -;
                plot(RInterp,Dbles,'b');
                grid on;
                title('Inferred Doubles Profile')
                ylim([0,1])
            subplot(2,2,4)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                Sngles = OnesDensity(FitParams,RInterp)+ThreesDensity(FitParams,RInterp)-2*Dbles;
                plot(RInterp,Sngles);
                grid on;
                title('Inferred Singles Profile')
                ylim([0,2])
             suptitle('Simultaneous fit to non-interacting lattice fermi gas. Two Pancakes.');
        end
        
        function FitParams = fitNonIntLattice_GaussTrap_Simultaneous(obj,tHz)
            %P = [Beta,Omega,Mu1,Mu2]
            %Simultaneously fit two non-interacting lattice fermi gas
            %profiles to the two spin components.
            A = 1; %density correction factor...
            t = obj.h*tHz;
            Mu1 = 0;
            Mu2 = 0;
            Beta = 0.5; %units of t
            V0 = 40; %this is what it fit too...should estimate if that's reasonable...
            Waist = 80e-6/obj.a;
%             Omega = 2*pi*400*sqrt(obj.mLi/t)*obj.a; %units of sqrt(m/t)*a
            
%             %Fixed Temp
%             OnesDensity =  @(P,X)latticeFGRadial1D([Beta,P(2),P(3),A],X);
%             ThreesDensity = @(P,X) latticeFGRadial1D([Beta,P(2),P(4),A],X);

            %Determine temperature, chemical potentials, and trapping
            %frequency from fit.
            OnesDensity =  @(P,X)latticeFGRadial_GaussPot1D([P(1),P(2),P(3),P(4),A],X);
            ThreesDensity = @(P,X) latticeFGRadial_GaussPot1D([P(1),P(2),P(3),P(5),A],X);
            
            OnesUnc = transpose(obj.OnesDensCorrUnc);
            OnesUnc(OnesUnc == 0 ) = 1;
            ThreesUnc = transpose(obj.ThreesDensCorrUnc);
            ThreesUnc(ThreesUnc == 0) = 1;
            
            OnesFit = @(P) (OnesDensity(P,obj.RadialPos)-transpose(obj.OnesDensCorr))./OnesUnc;
            ThreesFit = @(P) (ThreesDensity(P,obj.RadialPos)-transpose(obj.ThreesDensCorr))./ThreesUnc;
            
%             DensityFn = @(P,X) [OnesDensity(P,X),ThreesDensity(P,X)];
%             FitFn = @(P) (DensityFn(P,obj.RadialPos) - [transpose(obj.OnesDensCorr),transpose(obj.Threes.Occs_AzAvg)])./[transpose(obj.OnesDensCorrUnc),transpose(obj.Threes.Occs_AzAvgUnc)];
            FitFn = @(P) [OnesFit(P),ThreesFit(P)];

            InitP = [Beta,Waist,V0,Mu1,Mu2];
            FitParams = lsqnonlin(FitFn,InitP);
            OmegaMeanHz = sqrt(2*FitParams(3)/FitParams(2)^2)/sqrt(obj.mLi/abs(t))/obj.a;

            
            %display results
            fprintf('t = %0.2f Hz \n',tHz);
            fprintf('T = %0.2f t = %0.2f nK \n',1/FitParams(1),t/FitParams(1)/obj.kb/1e-9);
            fprintf('Chemical potential is zero at half filling \n')
            fprintf('Bandwidth in tight-binding model is 8t \n')
            fprintf('Mu1 = %0.2f t = (2pi) %0.2f Hz \n',FitParams(4),FitParams(4)*t/obj.h);
            fprintf('Mu3 = %0.2f t = (2pi) %0.2f Hz \n',FitParams(5),FitParams(5)*t/obj.h);
            fprintf('MuAvg = %0.2f t = (2pi) %0.2f Hz \n',0.5*(FitParams(4)+FitParams(5)),0.5*(FitParams(4)+FitParams(5))*t/obj.h);
            fprintf('DeltaMu = %0.2f t = (2pi) %0.2f Hz \n',0.5*(FitParams(4)-FitParams(5)),0.5*(FitParams(4)-FitParams(5))*t/obj.h);
            fprintf('Waist = %0.2f a = %0.2f um \n',FitParams(2),FitParams(2)*obj.a/1e-6);
            fprintf('V0 = %0.2f t \n',FitParams(3));
            fprintf('Omega = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi));
            fprintf('Omega Strong = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi)*sqrt(obj.CloudAspectRatio));
            fprintf('Omega Weak = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi)/sqrt(obj.CloudAspectRatio));
            
            RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            fh = figure('name','Simultaneous fit to non-interacting lattice fermi gas');
            subplot(2,2,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp),'b');
                title('Ones')
                grid on;
                ylim([0,1])
            subplot(2,2,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesDensity(FitParams,RInterp),'b');
                title('Threes')
                grid on;
                ylim([0,1])
            subplot(2,2,3)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Doubles Profile')
                ylim([0,1])
            subplot(2,2,4)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp)+ThreesDensity(FitParams,RInterp)-2*OnesDensity(FitParams,RInterp).*ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Inferred Singles Profile')
                ylim([0,2])
             suptitle('Simultaneous fit to non-interacting lattice fermi gas');
        end
        
        function [FitParams,ErrPerPoint] = fitNonIntEdge(obj,Cutoff,FitMode,FitType,InitParams,FixedParams,Display)
            %FitParams = fitNonIntEdge(obj,Cutoff,FitMode,FitType,InitParams,FixedParams,Display)
            %%%
            %Cutoff is either a density or radius value. For radius value,
            %this is the smallest radius included in the fit. For density,
            %this is the largest value of the minority density included in
            %fit.
            %%%FitMode = 'density' or 'radial', to decide what kind of
            %cutoff has been given
            %%%FitType = 'harmonic' or 'gaussian' for two different types
            %of trap profile. Harmonic by default.
            %%%InitParams are the initial parameters to be used for the
            %fit. For the Harmonic trap they are = [Beta,Omega,Mu_Up]
            %for the gaussian trap they are [Beta,Waist,ScaleDepth,Mu_Up]
            %%%FixedParams is a logical array the same size as InitParams.
            %A one fixes the respective initial parameter for the fit. A
            %zero allows that parameter to be fitted.
            %%%Display is a boolean specifying whether or not to plot the
            %results.
            %%%NOTE: to get the chemical potential so that it corresponds
            %%%to the Hubbard model with half-filling taken as zero
            %%%chemical potential, take the chemical potential fit here,
            %%%mu_fit = mu_hubb + U/2
            
            if ~exist('Display','var')
                Display = 1;
            end
            
            if ~exist('FitType','var')
                FitType = 'harmonic'; %or 'gaussian'
            end
            
            if ~strcmp(FitType,'harmonic') && ~strcmp(FitType,'gaussian')
                error('FitType argument to CorrDataSetAtt/fitNonIntEdge must be harmonic or gaussian');
            end
            
            if ~exist('FitMode','var')
                FitMode = 'density'; %or 'radial'
            end
            
            if ~strcmp(FitMode,'density') && ~strcmp(FitMode,'radial')
                error('FitMode argument to CorrDataSetAtt/fitNonIntEdge must be density or radial');
            end
            
            A = 1; %density correction factor...
            if ~exist('InitParams','var')
                if strcmp(FitType,'harmonic')
                    MuOnes = 0;
                    Beta = 2; %units of t
                    Omega = 0.1; %2*pi*200*sqrt(obj.mLi/t)*obj.a; %units of sqrt(m/t)*a
                    %Omega = OmegaHz*sqrt(obj.mLi/t)*obj.a;
                    InitParams = [Beta,Omega,MuOnes];
                elseif strcmp(FitType,'gaussian')
                    MuOnes = 0;
                    Beta = 2; %units of t
                    Vo = 16;
                    Waist = 73;
                    InitParams = [Beta,Waist,Vo,MuOnes];
                end
            end
            
            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            
            if strcmp(FitType,'harmonic')
                OnesDensity =  @(P,X)latticeFGRadial1D([P.*(1-FixedParams)+InitParams.*FixedParams,A],X);
                OnesCorr = @(P,X) latticeFGRadial_NNCorr1D([P.*(1-FixedParams)+InitParams.*FixedParams,A^2],X);
                lb = [0,0,-inf];
                ub = [inf,1,inf];
            elseif strcmp(FitType,'gaussian')
                OnesDensity =  @(P,X)latticeFGRadial_GaussPot1D([P.*(1-FixedParams)+InitParams.*FixedParams,A],X);
                OnesCorr = @(P,X) latticeFGRadial_NNCorr_GaussPot1D([P.*(1-FixedParams)+InitParams.*FixedParams,A^2],X);
                lb = [0,0,0,-inf];
                ub = [inf,inf,inf,inf];
            end
            
            if strcmp(FitMode,'density')
                CutoffR = min(obj.RadialPos(obj.ThreesDensCorr<Cutoff));
            elseif strcmp(FitMode,'radial')
                CutoffR = Cutoff;
            end
            
            Rs = obj.RadialPos(obj.RadialPos>CutoffR);
            OnesOccs = transpose(obj.OnesDensCorr(obj.RadialPos>CutoffR));
            OnesUnc = transpose(obj.OnesDensCorrUnc(obj.RadialPos>CutoffR));
            OnesUnc(OnesUnc == 0 ) = 1;
            
            
            %correlator uncertainties
            OnesCorrUnc = transpose(squeeze(obj.OnesCorrUnc(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:)));
            OnesCorrUnc = OnesCorrUnc(obj.RadialPos>CutoffR);
            OnesCorrUnc(OnesCorrUnc == 0) = 1;
            %correlators
            OnesCorrExp = transpose(squeeze(obj.OnesCorr(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:)));
            OnesCorrExp = OnesCorrExp(obj.RadialPos>CutoffR);
            %
            OnesDensFit = @(P) (OnesDensity(P,Rs)-OnesOccs)./OnesUnc;
            OnesCorrFit = @(P) (OnesCorr(P,Rs) - OnesCorrExp)./OnesCorrUnc;
            OnesFit = @(P) [OnesDensFit(P),OnesCorrFit(P)]; %(OnesDensity(P,Rs)-OnesOccs)./OnesUnc;
            %do fitting.
            if Display
                opts1=  optimset('display','on');
            else
                opts1=  optimset('display','off');
            end
            [FitParams,ResidualNorm] = lsqnonlin(OnesFit,InitParams,lb,ub,opts1);
            ErrPerPoint = ResidualNorm/length(Rs);
            
            if strcmp(FitType,'harmonic')
                Omega_FitUnits = FitParams(2);
                OmegaMean = FitParams(2)/sqrt(obj.mLi/obj.h)/obj.a;
                OmegaStrong = OmegaMean*sqrt(obj.CloudAspectRatio);
                OmegaWeak = OmegaMean/sqrt(obj.CloudAspectRatio);
                T = 1/FitParams(1);
                MuOnes = FitParams(3);
            elseif strcmp(FitType,'gaussian')
                Omega_FitUnits = sqrt(2*FitParams(3)/FitParams(2)^2);
                OmegaMean = sqrt(2*FitParams(3)/FitParams(2)^2)/sqrt(obj.mLi/obj.h)/obj.a;
                OmegaStrong = OmegaMean*sqrt(obj.CloudAspectRatio);
                OmegaWeak = OmegaMean/sqrt(obj.CloudAspectRatio);
                T = 1/FitParams(1);
                MuOnes = FitParams(4);
            end
            
            if Display
                %display results
                fprintf('T = %0.2f t \n',T);
                fprintf('Mu Ones = %0.2f t \n',MuOnes);
                fprintf('Omega = %0.2f = (2pi) %0.2f sqrt(tHz) \n',Omega_FitUnits,OmegaMean/(2*pi));
                fprintf('Omega Strong = (2pi) %0.2f sqrt(tHz) \n',OmegaStrong/(2*pi));
                fprintf('Omega Weak = (2pi) %0.2f sqrt(tHz) \n',OmegaWeak/(2*pi));
                
                RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
                %plot results
                FigHandle = figure('name','Edge fit to non-interacting lattice fermi gas');
                subplot(1,2,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ko');
                plot(RInterp(RInterp>CutoffR),OnesDensity(FitParams,RInterp(RInterp>CutoffR)),'b');
                plot(RInterp(RInterp<=CutoffR),OnesDensity(FitParams,RInterp(RInterp<=CutoffR)),'m');
                grid on;
                legend({'Ones','Threes','Fit to non-int FG','Non-int FG extended'});
                ylim([0,1])
                xlabel('Position (Lattice Sites)')
                ylabel('Density')
                
                subplot(1,2,2)
                UpsCorr = squeeze(obj.OnesCorr(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
                UpsCorrUnc = squeeze(obj.OnesCorrUnc(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,UpsCorr,UpsCorrUnc,'ro');
                hold on;
                plot(RInterp(RInterp>CutoffR),OnesCorr(FitParams,RInterp(RInterp>CutoffR)),'b');
                plot(RInterp(RInterp<=CutoffR),OnesCorr(FitParams,RInterp(RInterp<=CutoffR)),'m');
                grid on;
                legend({'OnesNNCorr','Fit to Non-int FG','Non int FG extended'});
                xlabel('Position (Lattice Sites)')
                ylabel('Correlator')
                
                suptitle(sprintf('%03d - %03d, Edge fit to non-int FG\n T = %0.2ft, Mu1 = %0.2f t, Omega = %0.2f = (2pi) %0.2f sqrt(tHz)',min(obj.Folders),max(obj.Folders),T,MuOnes,Omega_FitUnits,OmegaMean/(2*pi)))
            end
            
        end
        
        function [AllFitParams,CutoffRminList,ErrorsPerPt] = fitNonIntEdge_Incrementally(obj,tHz,UoverT,FitMode,FitType,InitParams,FixedParams,Display)
            %fitNonIntEdge_Incrementally(obj,tHz,UoverT,FitMode,FitType,InitParams,FixedParams)
            %P = [Beta,Omega,Mu1] for harmonic
            %P = [Beta,Waist,DepthScale,Mu1] for gaussian
            %Fit edge of the majority component to non-interacting lattice fermi gas
            %profile.
            
            t = obj.h*tHz;
            if ~exist('FitType','var')
                FitType = 'harmonic'; %or 'gaussian'
            end
            
            if ~exist('FitMode','var')
                FitMode = 'density'; %or 'radial'
            end
            
            if ~exist('UoverT','var')
                UoverT = -6;
                fprintf('WARNING, no UoverT entered to fitNonIntEdge. Using UoverT = -6');
            end
            
            if ~exist('Display','var')
                Display = 0;
            end
            
            if ~exist('InitParams','var')
                if strcmp(FitType,'harmonic')
                    MuOnes = 0;
                    Beta = 2; %units of t
                    Omega = 2*pi*200*sqrt(obj.mLi/t)*obj.a; %units of sqrt(m/t)*a
                    %Omega = OmegaHz*sqrt(obj.mLi/t)*obj.a;
                    InitParams = [Beta,Omega,MuOnes];
                elseif strcmp(FitType,'gaussian')
                    MuOnes = 0;
                    Beta = 2; %units of t
                    Vo = 16;
                    Waist = 73;
                    InitParams = [Beta,Waist,Vo,MuOnes];
                end
            end
            
            if ~exist('FixedParams','var')
                FixedParams = zeros(size(InitParams));
            end
            
            %do edge fits
            Density = 0.5*(obj.Filling_SandD + obj.Filling_1and3);
            CutoffRminList = obj.RadialPos(Density>0.03 & obj.ThreesDensCorr<0.02);
            AllFitParams = [];
            ErrorsPerPt = [];
            for ii = 1:length(CutoffRminList)
                CutoffRmin = CutoffRminList(ii);
                [FitParams,ErrPerPtTemp] = obj.fitNonIntEdge(CutoffRmin,'radial',FitType,InitParams,FixedParams,Display);
                AllFitParams = cat(1,AllFitParams,FitParams);
                ErrorsPerPt = cat(1,ErrorsPerPt,ErrPerPtTemp);
            end
            
            %plot results
            if size(AllFitParams,1)>1
                if strcmp(FitType,'harmonic')
                    FigHandle2 = figure('name','Fit parameters Vs. Cutoff');
                    subplot(1,3,1)
                    plot(CutoffRminList,1./AllFitParams(:,1))
                    grid on;
                    xlabel('R Cutoff (sites)');
                    ylabel('T (t)')
                    title('Temperature')
                    
                    subplot(1,3,2)
                    plot(CutoffRminList,AllFitParams(:,2));
                    grid on;
                    xlabel('R Cutoff (sites)');
                    ylabel('Omega')
                    title('Omega')
                    
                    subplot(1,3,3)
                    plot(CutoffRminList,AllFitParams(:,3));
                    grid on;
                    xlabel('R Cutoff (sites)');
                    ylabel('Mu (t)')
                    title('Mu majority')
                    
                elseif strcmp(FitType,'gaussian')
                    FigHandle2 = figure('name','Fit parameters Vs. Cutoff');
                    subplot(1,4,1)
                    plot(CutoffRminList,1./AllFitParams(:,1))
                    grid on;
                    xlabel('R Cutoff (sites)');
                    ylabel('T (t)')
                    ylim([0,1])
                    title('Temperature')
                    
                    subplot(1,4,2)
                    plot(CutoffRminList,AllFitParams(:,2));
                    grid on;
                    xlabel('R Cutoff (sites)');
                    ylabel('Waist')
                    ylim([0,200])
                    title('Waist')
                    
                    subplot(1,4,3)
                    plot(CutoffRminList,AllFitParams(:,3));
                    grid on;
                    xlabel('R Cutoff (sites)');
                    ylabel('Vo (t)')
                    title('Vo')
                    
                    subplot(1,4,4)
                    plot(CutoffRminList,AllFitParams(:,4));
                    grid on;
                    xlabel('R Cutoff (sites)');
                    ylabel('Mu (t)')
                    ylim([-4,4]);
                    title('Mu majority')
                    
                end
            end      
            
        end
        
        function [MuAvg,MuOneMinusMuThree] = estimateChemPot(obj,tHz,UoverT,OmegaHz)
            %estimate MuAvg and DeltaMu in units of t.
            t = obj.h*tHz;
            Omega = OmegaHz*sqrt(obj.mLi/abs(t))*obj.a;
            Fp = obj.fitNonIntEdge(tHz,UoverT,[1/0.5,Omega,0],[0,1,0]);
            
            MuOne = Fp(3)-UoverT/2; %in units of t
            Omega = Fp(2)/sqrt(obj.mLi/abs(t))/obj.a;
                  
            RHalfFilling = obj.getHalfFilling();
            MuAvg = 0.5*obj.mLi*Omega^2*obj.a^2*RHalfFilling^2/t; %in units of t
            MuOneMinusMuThree = 2*(MuOne - MuAvg); %Mu1-Mu2     
            
            fprintf('Mu1 = %0.2f t \n',MuOne);
            fprintf('Mu3 = %0.2f t \n',2*MuAvg-MuOne);
            fprintf('MuAvg = %0.2f t \n',MuAvg);
            fprintf('Mu1-Mu3 = %0.2f t \n',MuOneMinusMuThree);
            
            obj.Mu1MinusMu3 = MuOneMinusMuThree;
            obj.ChemPot_0 = MuAvg;
            
        end
        
         function [MuAvg,DeltaMu] = estimateChemPot_GaussPot(obj,tHz,Vo_t,Waist_m)
            %estimate MuAvg and DeltaMu in units of t.
            Fp = obj.fitNonIntEdge_GaussPot(tHz,Vo_t,Waist_m);
            
            t = obj.h*tHz;
            MuOne = Fp(4); %in units of t
%             Omega = Fp(2)/sqrt(obj.mLi/abs(t))/obj.a;
            Waist = Fp(2); %units of lattice sites
            Pot0 = Fp(3); %units of t
            RHalfFilling = obj.getHalfFilling();
           
            MuAvg = Pot0*(1-exp(-2*RHalfFilling^2/Waist^2)); %Waist from fit is in units of lattice sites...
%             MuAvg = 0.5*obj.mLi*Omega^2*obj.a^2*RHalfFilling^2/t; %in units of t
            DeltaMu = MuOne - MuAvg;     
            
            fprintf('Mu1 = %0.2f t \n',MuOne);
            fprintf('Mu3 = %0.2f t \n',2*MuAvg-MuOne);
            fprintf('MuAvg = %0.2f t \n',MuAvg);
            fprintf('DeltaMu = %0.2f t \n',DeltaMu);
            
        end
        
        function FitParams = fitAtomicLimit_Simultaneous(obj,UKHz)
            %P = [Beta,Omega,AvgMu,DeltaMu,A,SignOfU];
            %Simultaneously fit to all profiles (ones,threes,singles,and
            %doubles) for atomic limit.
            A = 0.98; %density correction factor...
            U = obj.h*UKHz*1e3;
            SignOfU = sign(UKHz);
            Mu1 = 0;
            Mu2 = 0;
%             DeltaMu = 0.5*(Mu1-Mu2); %fixed this to zero
            AvgMu = 0.5*(Mu1+Mu2);
            Beta = 0.5; %units of t
            Omega = 2*pi*400*sqrt(obj.mLi/abs(U))*obj.a; %units of sqrt(m/t)*a
            
            
            OnesDensity =  @(P,X)atomicLimitRadial_SingleComponentDensity1D([P(1),P(2),P(3),0,A,SignOfU],X);
            ThreesDensity = @(P,X) atomicLimitRadial_SingleComponentDensity1D([P(1),P(2),P(3),0,A,SignOfU],X);
            DoublesDensity = @(P,X) atomicLimitRadial_DoublesDensity1D([P(1),P(2),P(3),0,A,SignOfU],X);
            SinglesDensity = @(P,X) atomicLimitRadial_SinglesDensity1D([P(1),P(2),P(3),0,A,SignOfU],X);
            
            %get rid of zeros in uncertainty that will cause our fit to
            %throw errors
            OnesUnc = transpose(obj.OnesDensCorrUnc);
            OnesUnc(OnesUnc == 0) = 1;
            ThreesUnc = transpose(obj.ThreesDensCorrUnc);
            ThreesUnc(ThreesUnc == 0) = 1;
            SinglesUnc = transpose(obj.SinglesDensUnc);
            SinglesUnc(SinglesUnc == 0 ) = 1;
            DoublesUnc = transpose(obj.DoublesDensCorrUnc);
            DoublesUnc(DoublesUnc == 0 ) = 1;
            
            OnesFit = @(P) (OnesDensity(P,obj.RadialPos)-transpose(obj.OnesDensCorr))./OnesUnc;
            ThreesFit = @(P) (ThreesDensity(P,obj.RadialPos)-transpose(obj.ThreesDensCorr))./ThreesUnc;
            DoublesFit = @(P) (DoublesDensity(P,obj.RadialPos)-transpose(obj.DoublesDensCorr))./DoublesUnc;
            SinglesFit = @(P) (SinglesDensity(P,obj.RadialPos)-transpose(obj.SinglesDens))./SinglesUnc;
            %             DensityFn = @(P,X) [OnesDensity(P,X),ThreesDensity(P,X),DoublesDensity];
            FitFn = @(P) [OnesFit(P),ThreesFit(P),DoublesFit(P),SinglesFit(P)];
            
            InitP = [Beta,Omega,AvgMu];
            FitParams = lsqnonlin(FitFn,InitP);
            
            %print results
            OmegaMeanHz = abs(FitParams(2))/sqrt(obj.mLi/abs(U))/obj.a;
            fprintf('U = %0.2f KHz \n',UKHz);
            fprintf('T = %0.2f U = %0.2f nK \n',1/FitParams(1),abs(U)/FitParams(1)/obj.kb/1e-9);
            fprintf('Chemical potential is zero at half filling \n')
            fprintf('AvgMu = %0.2f U = (2pi) %0.2f Hz \n',FitParams(3),FitParams(3)*abs(U)/obj.h);
%             fprintf('DeltaMu = %0.2f U = (2pi) %0.2f Hz \n',FitParams(4),FitParams(4)*abs(U)/obj.h);
            fprintf('Omega = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi));
            fprintf('Omega Strong = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi)*sqrt(obj.CloudAspectRatio));
            fprintf('Omega Weak = (2pi) %0.2f Hz \n',OmegaMeanHz/(2*pi)/sqrt(obj.CloudAspectRatio));
            
            RInterp = linspace(min(obj.RadialPos),max(obj.RadialPos),100);
            %plot results
            fh = figure('name','Simultaneous Atomic Limit Fit');
            subplot(2,2,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,OnesDensity(FitParams,RInterp));
                grid on;
                title('Ones')
            subplot(2,2,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,ThreesDensity(FitParams,RInterp),'b');
                grid on;
                title('Threes')
            subplot(2,2,3)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc,'ro');
                hold on;
                plot(RInterp,DoublesDensity(FitParams,RInterp),'b');
                grid on;
                title('Doubles')
            subplot(2,2,4)
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc,'ro');
                hold on;
                plot(RInterp,SinglesDensity(FitParams,RInterp),'b');
                grid on;
                title('Singles')
             suptitle('Simultaneous fit to atomic limit');
            
        end
        
                
        function loadDQMC(obj)          
            obj.DQMC2 = load(obj.DQMC2_Path);
        end
        
        function f = fn_no_nans(obj,n,dmu,U,T,fn_handle)
            %helper function for fitDQMC2
                f = fn_handle(n,dmu,U,T);
                f(isnan(f)) = 0;
        end
        
        function [FitParams,FigHandle1,FigHandle2] = fitDQMC_NoTrap(obj,RCutoff,InitP,FixedP)
            %FitParams = [DeltaMu,U,T] in units of the tunneling t
            %simultaneously fit n1-n3,nSingles,nDoubles,Singles Correlator,
            %and Doubles Correlator to determine deltaMu and U. Make this
            %trap independent by fitting everything in the variables
            %(DeltaMu,n,U). At the end, extract expected trap parameters
            %using mu values extract from fit.
            
            %argument checking.
            if ~exist('RCutoff','var')
                RCutoff = 23;
            end
            
            if ~exist('InitP','var')
                %InitP = [0,0,0.1,-6,1e3,14];
                %InitP = [0,0,0.1,-6,1e2,1/16];
                InitP = [0,-6,0.5];
                %InitP = [0,0,0.1,-6,1e2,1];
            end
            
            if ~exist('FixedP','var')
                FixedP = zeros(size(InitP));
            end
            
            if ~isequal(size(FixedP),size(InitP))
                warning('Size of FixedP was not the same as size of InitP. Set FixedP to zeros.')
                FixedP = zeros(size(InitP));
            end
            
            RExp = obj.RadialPos;
            RExpUnc = obj.RadialPosUnc;
            PExp = obj.LocalPol;
            PExpUnc = obj.LocalPolUnc;
            nDiffExp = obj.OnesDensCorr-obj.ThreesDensCorr;
            nDiffExpUnc = sqrt(obj.OnesDensCorrUnc.^2+obj.ThreesDensCorrUnc.^2);
            nExp = 0.5*(obj.Filling_SandD+obj.Filling_1and3);
            nExpUnc = 0.5*sqrt(obj.Filling_SandD_Unc.^2+obj.Filling_1and3_Unc.^2);
            sExp = obj.SinglesDens;
            sExpUnc = obj.SinglesDensUnc;
            dExp = obj.DoublesDensCorr;
            dExpUnc = obj.DoublesDensCorrUnc;
            
            DoublesCorrExp = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            DoublesCorrExpUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExp = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExpUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));

%             RCutoff = 23;
            RExpUnc = RExpUnc(RExp<RCutoff);
            PExp = PExp(RExp<RCutoff);
            PExpUnc = PExpUnc(RExp<RCutoff);
            nExp = nExp(RExp<RCutoff);
            nExpUnc = nExpUnc(RExp<RCutoff);
            nDiffExp = nDiffExp(RExp<RCutoff);
            nDiffExpUnc = nDiffExpUnc(RExp<RCutoff);
            sExp = sExp(RExp<RCutoff);
            sExpUnc = sExpUnc(RExp<RCutoff);
            dExp = dExp(RExp<RCutoff);
            dExpUnc = dExpUnc(RExp<RCutoff);
            
            DoublesCorrExp = DoublesCorrExp(RExp<RCutoff);
            DoublesCorrExpUnc = DoublesCorrExpUnc(RExp<RCutoff);
            SinglesCorrExp = SinglesCorrExp(RExp<RCutoff);
            SinglesCorrExpUnc = SinglesCorrExpUnc(RExp<RCutoff);
            
            RExp = RExp(RExp<RCutoff);

            %get rid of any problem points in uncertainties
            RExpUnc(RExpUnc==0)=1;
            RExpUnc(isnan(RExpUnc)) = 1;
            PExpUnc(PExpUnc==0) = 1;
            PExpUnc(isnan(PExpUnc)) = 1;
            nExpUnc(nExpUnc==0) = 1;
            nExpUnc(isnan(nExpUnc)) = 1;
            nDiffExpUnc(nDiffExpUnc==0) = 1;
            nDiffExpUnc(isnan(nDiffExpUnc)) = 1;
            sExpUnc(sExpUnc==0) = 1;
            sExpUnc(isnan(sExpUnc)) = 1;
            dExpUnc(dExpUnc==0) = 1;
            dExpUnc(isnan(dExpUnc)) = 1;
            DoublesCorrExpUnc(DoublesCorrExpUnc==0) = 1;
            DoublesCorrExpUnc(isnan(DoublesCorrExpUnc)) = 1;
            SinglesCorrExpUnc(SinglesCorrExpUnc==0) = 1;
            SinglesCorrExpUnc(isnan(SinglesCorrExpUnc)) = 1;
            
            %DQMC results from file
            n = obj.DQMC2.nUp+obj.DQMC2.nDn;
            n(isnan(n)) = 0;
            ndiff = obj.DQMC2.nUp-obj.DQMC2.nDn;
            ndiff(isnan(ndiff)) = 0;
            p = obj.DQMC2.Pol;
            p(isnan(p)) = 0;
            s = obj.DQMC2.SingOcc;
            d = 0.5*(n - s);
            corrD = 0.25*obj.DQMC2.ddNN;
            corrD(isnan(corrD)) = 0;
            corrS = 0.25*obj.DQMC2.ssNN;
            corrS(isnan(corrS)) = 0;
            
            
            DelMus = 0.5*(obj.DQMC2.Mu_Ups-obj.DQMC2.Mu_Dns);
            Mus = 0.5*(obj.DQMC2.Mu_Ups+obj.DQMC2.Mu_Dns);
            Us = obj.DQMC2.Us;
            Ts = obj.DQMC2.Ts;
            
            %but we want them at Mus and DeltaMus, instead of Mu1 and Mu2s
            DelMuList = -3.9:0.1:3.9;
            %MuList = -3.9:0.1:3.9;
            nList = 0.01:0.025:1.5;
            UList = squeeze(Us(1,1,:,1));
            TList = squeeze(Ts(1,1,1,:));
%             [DelMuGrid,nGrid,UsGrid] = meshgrid(DelMuList,nList,UList);
            [nGrid,DelMuGrid,UsGrid,TsGrid] = ndgrid(nList,DelMuList,UList,TList);
%   don't want to do too much work. Still a grid along one direction!

            %Generate desired quantites on grid of DeltaMus and ns
%             SimulationPts = [DelMus(:),n(:),Us(:),Ts(:)];
%             SamplePts = [DelMuGrid(:),nGrid(:),UsGrid(:),TsGrid(:)];
%             pOnGrid = griddatan(SimulationPts,p(:),SamplePts);
            pOnGrid = zeros(size(nGrid));
            ndiffOnGrid = zeros(size(nGrid));
            sOnGrid = zeros(size(nGrid));
            dOnGrid = zeros(size(nGrid));
            corrDOnGrid = zeros(size(nGrid));
            corrSOnGrid = zeros(size(nGrid));
            muOnGrid = zeros(size(nGrid));
            
            %may need to change this if start doing lots of temperature
            %points?
            warning('off','all');
            for ii = 1:length(TList)
                %believe duplicate point warnings due to chemical potential
                %regions where n -> 0.
                %simulation points
                DelMusTemp = DelMus(:,:,:,ii);
                nTemp = n(:,:,:,ii);
                UTemp = Us(:,:,:,ii); 
                %sample points
                DelMuGridTemp = DelMuGrid(:,:,:,1);
                nGridTemp = nGrid(:,:,:,1);
                UsGridTemp = UsGrid(:,:,:,1);
                pOnGrid(:,:,:,ii) = griddata(DelMusTemp,nTemp,UTemp,p(:,:,:,ii),DelMuGridTemp,nGridTemp,UsGridTemp);
                ndiffOnGrid(:,:,:,ii) = griddata(DelMusTemp,nTemp,UTemp,ndiff(:,:,:,ii),DelMuGridTemp,nGridTemp,UsGridTemp);
                sOnGrid(:,:,:,ii) = griddata(DelMusTemp,nTemp,UTemp,s(:,:,:,ii),DelMuGridTemp,nGridTemp,UsGridTemp);
                dOnGrid(:,:,:,ii) = griddata(DelMusTemp,nTemp,UTemp,d(:,:,:,ii),DelMuGridTemp,nGridTemp,UsGridTemp);
                corrDOnGrid(:,:,:,ii) = griddata(DelMusTemp,nTemp,UTemp,corrD(:,:,:,ii),DelMuGridTemp,nGridTemp,UsGridTemp);
                corrSOnGrid(:,:,:,ii) = griddata(DelMusTemp,nTemp,UTemp,corrS(:,:,:,ii),DelMuGridTemp,nGridTemp,UsGridTemp);
                muOnGrid(:,:,:,ii) = griddata(DelMusTemp,nTemp,UTemp,Mus(:,:,:,ii),DelMuGridTemp,nGridTemp,UsGridTemp);
            end
            warning('on','all')
            
            pOnGrid(isnan(pOnGrid)) = 1;
            ndiffOnGrid(isnan(ndiffOnGrid)) = 1;
            sOnGrid(isnan(sOnGrid)) = 1;
            dOnGrid(isnan(dOnGrid)) = 1;
            corrDOnGrid(isnan(corrDOnGrid)) = 1;
            corrSOnGrid(isnan(corrSOnGrid)) = 1;
            muOnGrid(isnan(muOnGrid)) = 1;
            
            %then create interpolating function. This is much faster using
            %interp3, which assumes on grid points, than griddata, which
            %does not use gridpoints.
            pfn_interp = @(n,dmu,U,T) interpn(nGrid,DelMuGrid,UsGrid,TsGrid,pOnGrid,n,dmu,U,T);
%             pfn_interp = griddedInterpolant(nGrid,DelMuGrid,UsGrid,TsGrid,pOnGrid);
            pfn = @(n,dmu,U,T) obj.fn_no_nans(n,dmu,U,T,pfn_interp);
            
            ndifffn_interp = @(n,dmu,U,T) interpn(nGrid,DelMuGrid,UsGrid,TsGrid,ndiffOnGrid,n,dmu,U,T);
            ndifffn = @(n,dmu,U,T) obj.fn_no_nans(n,dmu,U,T,ndifffn_interp);
            
            sfn_interp = @(n,dmu,U,T) interpn(nGrid,DelMuGrid,UsGrid,TsGrid,sOnGrid,n,dmu,U,T);
            sfn = @(n,dmu,U,T) obj.fn_no_nans(n,dmu,U,T,sfn_interp);
            
            dfn_interp = @(n,dmu,U,T) interpn(nGrid,DelMuGrid,UsGrid,TsGrid,dOnGrid,n,dmu,U,T);
            dfn = @(n,dmu,U,T) obj.fn_no_nans(n,dmu,U,T,dfn_interp);
            
            corrDfn_interp = @(n,dmu,U,T) interpn(nGrid,DelMuGrid,UsGrid,TsGrid,corrDOnGrid,n,dmu,U,T);
            corrDfn = @(n,dmu,U,T) obj.fn_no_nans(n,dmu,U,T,corrDfn_interp);
            
            corrSfn_interp = @(n,dmu,U,T) interpn(nGrid,DelMuGrid,UsGrid,TsGrid,corrSOnGrid,n,dmu,U,T);
            corrSfn = @(n,dmu,U,T) obj.fn_no_nans(n,dmu,U,T,corrSfn_interp);

            mufn_interp = @(n,dmu,U,T)interpn(nGrid,DelMuGrid,UsGrid,TsGrid,muOnGrid,n,dmu,U,T);
            mufn = @(n,dmu,U,T)obj.fn_no_nans(n,dmu,U,T,mufn_interp);
            
            pfitfn_prelim = @(P,n) pfn(n,P(1)*ones(size(n)),P(2)*ones(size(n)),P(3)*ones(size(n)));
            ndifffitfn_prelim = @(P,n) ndifffn(n,P(1)*ones(size(n)),P(2)*ones(size(n)),P(3)*ones(size(n)));
            sfitfn_prelim = @(P,n) sfn(n,P(1)*ones(size(n)),P(2)*ones(size(n)),P(3)*ones(size(n)));
            dfitfn_prelim = @(P,n) dfn(n,P(1)*ones(size(n)),P(2)*ones(size(n)),P(3)*ones(size(n)));
            corrDfitfn_prelim = @(P,n) corrDfn(n,P(1)*ones(size(n)),P(2)*ones(size(n)),P(3)*ones(size(n)));
            corrSfitfn_prelim = @(P,n) corrSfn(n,P(1)*ones(size(n)),P(2)*ones(size(n)),P(3)*ones(size(n)));

            %construct fit function 
            PrelimFixedP = FixedP;
            PrelimInitP = InitP;
            Prelim_pfit = @(P) (pfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,nExp)-(PExp))./(PExpUnc);
            Prelim_ndifffit = @(P)(ndifffitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,nExp)-(nDiffExp))./(nDiffExpUnc);
            Prelim_sfit = @(P) (sfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,nExp)-(sExp))./(sExpUnc);
            Prelim_dfit = @(P) (dfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,nExp)-(dExp))./(dExpUnc);
            Prelim_corrDfit = @(P) (corrDfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,nExp)-(DoublesCorrExp))./(DoublesCorrExpUnc);
            Prelim_corrSfit = @(P) (corrSfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,nExp)-(SinglesCorrExp))./(SinglesCorrExpUnc);
            
            Prelimfullfit = @(P) [Prelim_ndifffit(P),Prelim_sfit(P),Prelim_dfit(P),Prelim_corrDfit(P),Prelim_corrSfit(P)];

            
            LowerBounds_prelim = [min(DelMuList),min(UList),min(TList)];
            UpperBounds_prelim = [max(DelMuList),max(UList),max(TList)];
            options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1000);
            FitParams = lsqnonlin(Prelimfullfit,PrelimInitP,LowerBounds_prelim,UpperBounds_prelim,options);
            
            nInterp = linspace(0,max(nExp),100);
            Rinterp = linspace(min(RExp),max(RExp),100);
            NRows = 2;
            NCols = 4;
            
            %no trap variation
            FigHandle1 = figure('name','No U/t Variation Across Trap');
            subplot(NRows,NCols,1)
            errorbar(nExp,PExp,PExpUnc,'b.');
            hold on;
            plot(nInterp,pfitfn_prelim(FitParams,nInterp),'r');
            grid on;
            ylim([0,1])
            title('Local Polarization, Inferred')
            
            subplot(NRows,NCols,2)
            errorbar(nExp,nDiffExp,nDiffExpUnc,'b.');
            hold on;
            plot(nInterp,ndifffitfn_prelim(FitParams,nInterp),'r');
            grid on;
            ylim([0,1])
            title('n1-n3, Fit')
            ylabel('n1-n3')
            
            subplot(NRows,NCols,3)
            errorbar(nExp,sExp,sExpUnc,'b.');
            hold on;
            plot(nInterp,sfitfn_prelim(FitParams,nInterp),'r');
            grid on;
            ylim([0,1.5])
            title('Singles Density, Fit')
            
            subplot(NRows,NCols,4)
            errorbar(nExp,dExp,dExpUnc,'b.');
            hold on;
            plot(nInterp,dfitfn_prelim(FitParams,nInterp),'r');
            grid on;
            ylim([0,1.5])
            title('Doubles Density, Fit')
            
            
            subplot(NRows,NCols,8)
            errorbar(nExp,DoublesCorrExp,DoublesCorrExpUnc,'b.');
            hold on;
            plot(nInterp,corrDfitfn_prelim(FitParams,nInterp),'r');
            grid on;
            ylim([-0.05,0.02])
            title('Doubles Correlator, Fit')
            
            subplot(NRows,NCols,7)
            errorbar(nExp,SinglesCorrExp,SinglesCorrExpUnc,'b.');
            hold on;
            plot(nInterp,corrSfitfn_prelim(FitParams,nInterp),'r');
            grid on;
            title('Singles Correlator, Fit')
            
            %if fg initiated, also do that
            if obj.InitFG
                
                %pol
                subplot(NRows,NCols,1)
                plot(nInterp,obj.NonIntFG.pfn_n_dmu_T(nInterp,FitParams(1),FitParams(3)));
                legend({'Expt','DQMC','NonIntFG'});
                %n1-n3
                subplot(NRows,NCols,2)
                 plot(nInterp,obj.NonIntFG.n1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3))-obj.NonIntFG.n2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)));
                %nsingles
                subplot(NRows,NCols,3)
                 plot(nInterp,obj.NonIntFG.n1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3))+obj.NonIntFG.n2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3))-obj.NonIntFG.n1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).*obj.NonIntFG.n2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)));
                %ndoubles
                subplot(NRows,NCols,4)
                plot(nInterp,obj.NonIntFG.n1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).*obj.NonIntFG.n2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)));
                %doubles corr
                subplot(NRows,NCols,8)
                FGDoublesCorr = obj.NonIntFG.c1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).*obj.NonIntFG.n2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).^2 ...
                + obj.NonIntFG.c2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).*obj.NonIntFG.n1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).^2 ...
                + obj.NonIntFG.c1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).*obj.NonIntFG.c2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3));
                plot(nInterp,FGDoublesCorr);
                %singles corr
                subplot(NRows,NCols,7)
                FGSinglesCorr = obj.NonIntFG.c1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)) ...
                    + obj.NonIntFG.c2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)) ...
                    + 4*FGDoublesCorr ...
                    - 4*obj.NonIntFG.c1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).*obj.NonIntFG.n2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)) ...
                    - 4*obj.NonIntFG.c2fn_n_dmu_T(nInterp,FitParams(1),FitParams(3)).*obj.NonIntFG.n1fn_n_dmu_T(nInterp,FitParams(1),FitParams(3));
                plot(nInterp,FGSinglesCorr);
            end
            
            th = suptitle(sprintf('%s \nDMu = %0.2ft, U/t = %0.2f, T = %0.2ft',obj.getDescriptionString(),FitParams(1),FitParams(2),FitParams(3)));
            set(th,'interpreter','none');
             
            FigHandle2 = figure('name','Inferred chemical potential across trap');
            plot(RExp,mufn(nExp,FitParams(1),FitParams(2),FitParams(3)));
            %fit harmonic trap frq
            MuTrapFn = @(P,R) P(1) - 0.5*P(2)^2*R.^2;
            ExcludePts = ones(size(RExp));
            ExcludePts(mufn(nExp,FitParams(1),FitParams(2),FitParams(3)) == 0) = 0; %exclude points where got nans from interpolation. These were set to zero earlier.
            MuTrapFit = @(P) (MuTrapFn(P,RExp)-transpose(mufn(nExp,FitParams(1),FitParams(2),FitParams(3))))./RExpUnc.*ExcludePts;
            FpTrap = lsqnonlin(MuTrapFit,[0.4,0.08]);
            OmegaMeanSqrtT = FpTrap(2)/sqrt(obj.mLi/obj.h)/obj.a;
            %need to multiply this by sqrt(TinHz) to get real trapping frq.
            %Also fit to gaussian potential
            %Omega = sqrt(2*Pot0/Waist^2) in the same units...
            %Params = [mu0,Pot0,Waist)
            MuTrapGaussFn = @(P,R) P(1) -P(2)*(1-exp(-2*R.^2/P(3)^2));
            MuTrapFitGauss = @(P) (MuTrapGaussFn(P,RExp)-transpose(mufn(nExp,FitParams(1),FitParams(2),FitParams(3))))./RExpUnc.*ExcludePts;
            FpTrapGauss = lsqnonlin(MuTrapFitGauss,[0.4,40,90e-6/obj.a]);
            OmegaGauss = sqrt(2*FpTrapGauss(2)/FpTrapGauss(3).^2);
            OmegaGaussSqrtT = OmegaGauss/sqrt(obj.mLi/obj.h)/obj.a;
            
            hold on;
            plot(RExp,MuTrapFn(FpTrap,RExp),'r--');
            plot(RExp,MuTrapGaussFn(FpTrapGauss,RExp),'g.-');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('Mu (t)')
            legend({'Inferred Mu','Harmonic Fit','Gaussian Fit'});
            th = title(sprintf('Harmonic: Mu0 = %0.2f t, Omega = %0.2f = (2*pi) %0.2f sqrt(t) Hz \n Gaussian: Mu0 = %0.2f t, Pot0 = %0.2f, Waist = %0.2f a, Omega = %0.2f = (2pi) %0.2f sqrt(t)',...
                FpTrap(1),FpTrap(2),OmegaMeanSqrtT/(2*pi),FpTrapGauss(1),FpTrapGauss(2),FpTrapGauss(3),OmegaGauss,OmegaGaussSqrtT/(2*pi)));
            set(th,'interpreter','none');
        end
        
        function Fp = fitDQMC2(obj,RCutoff,InitP,FixedP)
            %Fit to DQMC data including variation of U/t across the trap.
            %InitP = [DeltaMu,Mu0,Vo,Omega,T]
            %DeltaMu, Mu0, and T are measured in units of t, the hopping
            %Omega is the trapping frequency for a harmonic trap.
            
            %%%argument checking
            if ~exist('RCutoff','var')
                RCutoff = 30;
            end
            
            if ~exist('InitP','var')
                InitP = [2.4,0,0.1,0.11,0.5];
            end
            
            if ~exist('FixedP','var')
                FixedP = zeros(size(InitP));
            end
            
            %%%get experimental measurements
            RExp = obj.RadialPos;
            PExp = obj.LocalPol;
            PExpUnc = obj.LocalPolUnc;
            nDiffExp = (obj.OnesDensCorr-obj.ThreesDensCorr);
            nDiffExpUnc = sqrt(obj.OnesDensCorrUnc.^2+obj.ThreesDensCorrUnc.^2);
            nExp = 0.5*(obj.Filling_SandD+obj.Filling_1and3);
            nExpUnc = 0.5*sqrt(obj.Filling_SandD_Unc.^2+obj.Filling_1and3_Unc.^2);
            sExp = obj.SinglesDens;
            sExpUnc = obj.SinglesDensUnc;
            dExp = obj.DoublesDensCorr;
            dExpUnc = obj.DoublesDensCorrUnc;
            
            DoublesCorrExp = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            DoublesCorrExpUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExp = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExpUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));

            %remove points larger than cutoff radius
            PExp = PExp(RExp<RCutoff);
            PExpUnc = PExpUnc(RExp<RCutoff);
            nDiffExp = nDiffExp(RExp<RCutoff);
            nDiffExpUnc = nDiffExpUnc(RExp<RCutoff);
            nExp = nExp(RExp<RCutoff);
            nExpUnc = nExpUnc(RExp<RCutoff);
            sExp = sExp(RExp<RCutoff);
            sExpUnc = sExpUnc(RExp<RCutoff);
            dExp = dExp(RExp<RCutoff);
            dExpUnc = dExpUnc(RExp<RCutoff);
            
            DoublesCorrExp = DoublesCorrExp(RExp<RCutoff);
            DoublesCorrExpUnc = DoublesCorrExpUnc(RExp<RCutoff);
            SinglesCorrExp = SinglesCorrExp(RExp<RCutoff);
            SinglesCorrExpUnc = SinglesCorrExpUnc(RExp<RCutoff);
            
            RExp = RExp(RExp<RCutoff);

            %get rid of any problem points in uncertainties
            PExpUnc(PExpUnc==0) = 1;
            PExpUnc(isnan(PExpUnc)) = 1;
            nDiffExpUnc(nDiffExpUnc==0) = 1;
            nDiffExpUnc(isnan(nDiffExpUnc)) = 1;
            nExpUnc(nExpUnc==0) = 1;
            nExpUnc(isnan(nExpUnc)) = 1;
            sExpUnc(sExpUnc==0) = 1;
            sExpUnc(isnan(sExpUnc)) = 1;
            dExpUnc(dExpUnc==0) = 1;
            dExpUnc(isnan(dExpUnc)) = 1;
            DoublesCorrExpUnc(DoublesCorrExpUnc==0) = 1;
            DoublesCorrExpUnc(isnan(DoublesCorrExpUnc)) = 1;
            SinglesCorrExpUnc(SinglesCorrExpUnc==0) = 1;
            SinglesCorrExpUnc(isnan(SinglesCorrExpUnc)) = 1;
            
            %%%DQMC results from file
            n = obj.DQMC2.nUp+obj.DQMC2.nDn;
            n(isnan(n)) = 0;
            ndiff = obj.DQMC2.nUp-obj.DQMC2.nDn;
            ndiff(isnan(ndiff)) = 0;
            p = obj.DQMC2.Pol;
            p(isnan(p)) = 0;
            s = obj.DQMC2.SingOcc;
            d = 0.5*(n - s);
            corrD = 0.25*obj.DQMC2.ddNN;
            corrD(isnan(corrD)) = 0;
            corrS = 0.25*obj.DQMC2.ssNN;
            corrS(isnan(corrS)) = 0;
            
            DelMus = 0.5*(obj.DQMC2.Mu_Ups-obj.DQMC2.Mu_Dns);
            Mus = 0.5*(obj.DQMC2.Mu_Ups+obj.DQMC2.Mu_Dns);
            Us = obj.DQMC2.Us;
            Ts = obj.DQMC2.Ts;
            
            %%%Generate everything as functions Mus,DeltaMus, Us, and Ts
            %instead of Mu_Ups, Mu_Dns, Us, and Ts
            %but we want them at Mus and DeltaMus, instead of Mu1 and Mu2s
            DelMuList = -3.9:0.1:3.9;
            MuList = -3.9:0.1:3.9;
            UList = Us(1,1,:,1);
            TList = Ts(1,1,1,:);
%             [DelMuGrid,MuGrid,UsGrid] = meshgrid(DelMuList,MuList,UList);
            [MuGrid,DelMuGrid,UsGrid,TsGrid] = ndgrid(MuList,DelMuList,UList,TList);

            pOnGrid = zeros(size(MuGrid));
            ndiffOnGrid = zeros(size(MuGrid));
            nOnGrid = zeros(size(MuGrid)); 
            sOnGrid = zeros(size(MuGrid));
            dOnGrid = zeros(size(MuGrid));
            corrDOnGrid = zeros(size(MuGrid));
            corrSOnGrid = zeros(size(MuGrid)); 
            warning('off','all');
            for ii = 1:length(TList)
                %believe duplicate point warnings due to chemical potential
                %regions where n -> 0.
                %simulation points
                DelMusTemp = DelMus(:,:,:,ii);
                MuTemp = Mus(:,:,:,ii);
                UTemp = Us(:,:,:,ii); 
                %sample points
                DelMuGridTemp = DelMuGrid(:,:,:,1);
                MuGridTemp = MuGrid(:,:,:,1);
                UsGridTemp = UsGrid(:,:,:,1);
                pOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,p(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                ndiffOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,ndiff(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                nOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,n(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                sOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,s(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                dOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,d(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                corrDOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,corrD(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                corrSOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,corrS(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
            end
            warning('on','all')
            
            %Generate desired quantites on grid of DeltaMus and Mus
            pOnGrid(isnan(pOnGrid)) = 1;
            ndiffOnGrid(isnan(ndiffOnGrid)) = 1;
            nOnGrid(isnan(nOnGrid)) = 1;
            sOnGrid(isnan(sOnGrid)) = 1;
            dOnGrid(isnan(dOnGrid)) = 1;
            corrDOnGrid(isnan(corrDOnGrid)) = 1;
            corrSOnGrid(isnan(corrSOnGrid)) = 1;
            
            %then create interpolating function. This is much faster using
            %interp3, which assumes on grid points, than griddata, which
            %does not use gridpoints.
            pfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,pOnGrid,mu,dmu,U,T);
            pfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,pfn_interp);
            
            ndifffn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,ndiffOnGrid,mu,dmu,U,T);
            ndifffn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,ndifffn_interp);
            
            nfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,nOnGrid,mu,dmu,U,T);
            nfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,nfn_interp);
            
            sfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,sOnGrid,mu,dmu,U,T);
            sfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,sfn_interp);
            
            dfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,dOnGrid,mu,dmu,U,T);
            dfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,dfn_interp);
            
            corrDfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,corrDOnGrid,mu,dmu,U,T);
            corrDfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,corrDfn_interp);
            
            corrSfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,corrSOnGrid,mu,dmu,U,T);
            corrSfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,corrSfn_interp);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %first fit with no U/t variation across the trap. Can do this
            %in a trap independent way.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            PrelimInitP = [3,-7.5,0.5];
            PrelimFixedP = [0,0,0];
            PrelimFp = obj.fitDQMC_NoTrap(RCutoff,PrelimInitP,PrelimFixedP);
            
            %delmu,mu,U
%             mufn = @(P,R)P(2) - 0.5*P(3)^2*R.^2;
%             pfitfn_prelim = @(P,R) pfn(P(1)*ones(size(R)),mufn(P,R),P(4)*ones(size(R)));
%             nfitfn_prelim = @(P,R) nfn(P(1)*ones(size(R)),mufn(P,R),P(4)*ones(size(R)));
%             sfitfn_prelim = @(P,R) sfn(P(1)*ones(size(R)),mufn(P,R),P(4)*ones(size(R)));
%             dfitfn_prelim = @(P,R) dfn(P(1)*ones(size(R)),mufn(P,R),P(4)*ones(size(R)));
%             corrDfitfn_prelim = @(P,R) corrDfn(P(1)*ones(size(R)),mufn(P,R),P(4)*ones(size(R)));
%             corrSfitfn_prelim = @(P,R) corrSfn(P(1)*ones(size(R)),mufn(P,R),P(4)*ones(size(R)));
% 
%             PrelimFixedP = FixedP(1:4);
%             PrelimInitP = InitP(1:4);
%             
%             %delmu,mu0,omega,U
%             LowerBounds_prelim = [-4,-4,0.02,-8];
%             UpperBounds_prelim = [4,4,1,-4];
%             %construct fit function
%             % Prelimfullfit = @(P) (fitfn_prelim(P,RExp)-transpose(PExp))./transpose(PExpUnc);
%             Prelim_pfit = @(P) (pfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,RExp)-transpose(PExp))./transpose(PExpUnc);
%             Prelim_nfit = @(P) (nfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,RExp)-transpose(nExp))./transpose(nExpUnc);
%             Prelim_sfit = @(P) (sfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,RExp)-transpose(sExp))./transpose(sExpUnc);
%             Prelim_dfit = @(P) (dfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,RExp)-transpose(dExp))./transpose(dExpUnc);
%             Prelim_corrDfit = @(P) (corrDfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,RExp)-transpose(DoublesCorrExp))./transpose(DoublesCorrExpUnc);
%             Prelim_corrSfit = @(P) (corrSfitfn_prelim(P.*(1-PrelimFixedP)+PrelimInitP.*PrelimFixedP,RExp)-transpose(SinglesCorrExp))./transpose(SinglesCorrExpUnc);
%             
%             Prelimfullfit = @(P) [Prelim_pfit(P),Prelim_nfit(P),Prelim_corrDfit(P),Prelim_corrSfit(P)];
%             
%             
%             % options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',4000);
%             PrelimFp = lsqnonlin(Prelimfullfit,PrelimInitP,LowerBounds_prelim,UpperBounds_prelim);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now fit including trap variations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            %this expression comes from assuming a gaussian wavefunctionand
            %computing t and U integrals
%             uovertfn = @(P,R) P(4)*(1-0.5*P(3)^2/P(5)*R.^2).^(-3/4).*exp(2*sqrt((P(5)-0.5*P(3)^2*R.^2)))*exp(-2*sqrt(P(5)));
            %but this expression, need the factors in the exponential to be
            %in units of the lattice recoil, not the hopping.
            %need an extra fit parameter for this.
%             uovertfn = @(P,R) P(4)*(1-0.5*P(3)^2/P(5)*R.^2).^(-3/4).*exp(2*sqrt(P(6))*sqrt((P(5)-0.5*P(3)^2*R.^2)))*exp(-2*sqrt(P(6))*sqrt(P(5)));

            %can instead get U/t from trap depth measurements + band
            %structure calculations.
            LattAtten = 0.48; %determined from fit 04/04/2017, for retro PD = 170mV @10V latt servo
            DepthList = linspace(2,20,100); %in Er
            [ts,txs,tys,tdiags,us,~] = Lattice_2D_Asymmetric_Interpolated(DepthList,LattAtten);
            uFn = @(Depth) interp1(DepthList,us,Depth);
            tFn = @(Depth) interp1(DepthList,ts,Depth);
            uOvert_BandStruct = @(Depth) uFn(Depth)./tFn(Depth);
        
            %previous trap measurements...added zero for convenience
            rads = [0,5.3465,9.7489,12.6254,14.9491,16.9302,18.7301,20.4012,...
                21.9388,23.3554,24.6893,25.9374,27.1337,28.2754,29.3925,30.4709,...
                31.5201,32.5339,33.5149,34.4593,35.3653,36.2474,37.1212,37.9911,...
                38.8321,39.6420];
            vs = [57.7016,57.7016,57.2359,56.7780,56.5207,56.1764,55.9290,...
                55.4548,54.9821,54.6439,54.1269,53.9811,53.5650,53.1263,52.5598,...
                52.4731,52.1904,51.9089,51.5164,51.1449,50.9224,50.6223,50.3907,...
                49.8677,49.8024,49.2935];
            pot = @(R) interp1(rads,vs,R);
            uovertfn = @(P,R) uOvert_BandStruct(P(3)*pot(R));
            
            %to generate initial guess from the other fit, have to do some
            %gymnastics...
            
            GuessForPotentialScale = fzero(@(d) uOvert_BandStruct(d)-PrelimFp(2),10)/vs(1);
            %or from U/t inferred at this depth directly...for this one
            %might need more points between -7 and -6
%             uoverts_meas_infer = [-6.2141,-6.1051,-6.0003,-5.9425,-5.8661,-5.8120,-5.7100,-5.6106,-5.5408,-5.4361,-5.4071,-5.3249,-5.2266,-5.1031,-5.0845];
%             rads = [5.3465,9.7489,12.6254,14.9491,16.9302,18.7301,20.4012,21.9388,23.3554,24.6893,25.9374,27.1337,28.2754,29.3925,30.4709];
%             uovertfn = @(P,R) interp1(rads,uoverts_meas_infer,R);
            
            %harmonic potential
            mufn = @(P,R)P(2) - 0.5*P(4)^2*R.^2;
            %use measured potential
%             mufn = @(P,R) P(2) + P(4)*P(3)*(pot(R)-vs(1));
            %but actually already get t from band structure
%             mufn = @(P,R) P(2) + 1/tFn(P(3)*vs(1))*P(3)*(pot(R)-vs(1));
            
%fitfn not evaluating to size I expect...
            pfitfn = @(P,R) pfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(5)*ones(size(R)));
            ndifffitfn = @(P,R) ndifffn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(5)*ones(size(R)));
            nfitfn = @(P,R) nfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(5)*ones(size(R)));
            sfitfn = @(P,R) sfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(5)*ones(size(R)));
            dfitfn = @(P,R) dfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(5)*ones(size(R)));
            corrDfitfn = @(P,R) corrDfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(5)*ones(size(R)));
            corrSfitfn = @(P,R) corrSfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(5)*ones(size(R)));

            %P = [DelMu,Mu_o,Omega,Uovert_o,Vo,Erec_Per_t]
            %these are all measured in units of t
            if ~FixedP(1)
                InitP(1) = PrelimFp(1);
            end
%             InitP(3) = GuessForPotentialScale;
            if ~FixedP(5)
                InitP(5) = PrelimFp(3);
            end
            
            %construct fit function
            pfullfit = @(P) (pfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(PExp))./transpose(PExpUnc);
            ndifffit = @(P) (ndifffitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(nDiffExp))./transpose(nDiffExpUnc);
            nfullfit = @(P) (nfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(nExp))./transpose(nExpUnc);
            sfullfit = @(P) (sfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(sExp))./transpose(sExpUnc);
            dfullfit = @(P) (dfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(dExp))./transpose(dExpUnc);
            corrDfullfit = @(P) (corrDfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(DoublesCorrExp))./transpose(DoublesCorrExpUnc);
            corrSfullfit = @(P) (corrSfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(SinglesCorrExp))./transpose(SinglesCorrExpUnc);
           
%             fullfit = @(P) [pfullfit(P),nfullfit(P),corrDfullfit(P),corrSfullfit(P)];
            fullfit = @(P) [ndifffit(P),nfullfit(P),sfullfit(P),dfullfit(P),corrDfullfit(P),corrSfullfit(P)];

            LowerBounds = [-4,-4,0,0,0.31];
            UpperBounds = [4,4,1,Inf,0.7];
            Fp = lsqnonlin(fullfit,InitP,LowerBounds,UpperBounds);
%             Fp = lsqnonlin(fullfit,InitP);
            %inferred parameters
%             t_fit_Hz = 14.66e3/Fp(4);
%             OmegaMeanHz = Fp(3)/sqrt(obj.mLi/abs(t_fit_Hz*obj.h))/obj.a;
            
            Rinterp = linspace(min(RExp),max(RExp),300);

            NRows = 2;
            NCols = 4;
            
            %with trap variation
            figure('name','U/t Variation Across Trap')
            subplot(NRows,NCols,1)
            errorbar(RExp,nDiffExp,nDiffExpUnc,'b.');
            hold on;
            plot(Rinterp,ndifffitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1])
            title('n1-n3')
            
            subplot(NRows,NCols,5)
            errorbar(RExp,PExp,PExpUnc,'b.');
            hold on;
            plot(Rinterp,pfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1])
            title('Local Polarization')
            
            subplot(NRows,NCols,2)
            errorbar(RExp,nExp,nExpUnc,'b.');
            hold on;
            plot(Rinterp,nfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Total Density')
            
            subplot(NRows,NCols,3)
            errorbar(RExp,sExp,sExpUnc,'b.');
            hold on;
            plot(Rinterp,sfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Singles Density')
            
            subplot(NRows,NCols,4)
            errorbar(RExp,dExp,dExpUnc,'b.');
            hold on;
            plot(Rinterp,dfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Doubles Density')
           
            
            subplot(NRows,NCols,7)
            errorbar(RExp,SinglesCorrExp,SinglesCorrExpUnc,'b.');
            hold on;
            plot(Rinterp,corrSfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([-0.05,0.01]);
            title('Singles Correlator')
            
            subplot(NRows,NCols,8)
            errorbar(RExp,DoublesCorrExp,DoublesCorrExpUnc,'b.');
            hold on;
            plot(Rinterp,corrDfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([-0.05,0.01]);
            title('Doubles Correlator')
            
%             suptitle(sprintf('U/t Variation across trap \nDMu = %0.2ft, Mu = %0.2ft, Omega = %0.3f, U/t = %0.2f, Vo = %0.2ft, Er/t = %0.2f \n t = %0.1f Hz, Omega = (2pi) %0.1fHz',Fp(1),Fp(2),Fp(3),uovertfn(Fp,0),Fp(5),Fp(6),t_fit_Hz,OmegaMeanHz/(2*pi)));
            th = suptitle(sprintf('%s \nU/t Variation across trap \nDMu = %0.2ft, T = %0.2ft, Mu = %0.2ft, Vo = %0.2f, Omega = %0.2f',obj.getDescriptionString(),Fp(1),Fp(5),Fp(2),Fp(3),Fp(4)));
            set(th,'interpreter','none')
            
            %trap variation
            figure('name','Variation of U/t')
            subplot(2,2,1)
            plot(Rinterp,uovertfn(Fp,Rinterp));
            grid on;
            title('U/t variation across trap')
            xlabel('Position (Lattice Sites)')
            ylabel('U/t')
            
            subplot(2,2,2)
            plot(Rinterp,Fp(3)*pot(Rinterp),'r');
            grid on;
             xlabel('Position (Lattice Sites)')
            ylabel('Potential (Er)')
            
          
            subplot(2,2,3)
            plot(Rinterp,mufn(Fp,Rinterp),'r');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('Mu (t)')
            
            th = suptitle(sprintf('%s \n Trap Variation',obj.getDescriptionString()));
            set(th,'interpreter','none');
        end
              
        function Fp = fitDQMC3(obj,RCutoff,InitP,FixedP)
            %Fit to DQMC data including variation of U/t across the trap.
            %InitP = [DeltaMu,T,U/t @ center,Mu0,Er/t]
            %DeltaMu, Mu0, and T are measured in units of t, the hopping
            %Omega is the trapping frequency for a harmonic trap.
            %cleaned up version of fitDQMC2 in more natural fit parameters.
            
            %%%argument checking
            if ~exist('RCutoff','var')
                RCutoff = 30;
            end
            
            if ~exist('InitP','var')
                InitP = [2.4,0.5,-6,0,0.11];
            end
            
            if ~exist('FixedP','var')
                FixedP = zeros(size(InitP));
            end
            
            if ~isequal(size(FixedP),size(InitP));
                warning('sizes of FixedP and InitP were not the same in fitDQMC3');
                FixedP = zeros(size(InitP));
            end
            
            %%%get experimental measurements
            RExp = obj.RadialPos;
            PExp = obj.LocalPol;
            PExpUnc = obj.LocalPolUnc;
            nDiffExp = (obj.OnesDensCorr-obj.ThreesDensCorr);
            nDiffExpUnc = sqrt(obj.OnesDensCorrUnc.^2+obj.ThreesDensCorrUnc.^2);
            nExp = 0.5*(obj.Filling_SandD+obj.Filling_1and3);
            nExpUnc = 0.5*sqrt(obj.Filling_SandD_Unc.^2+obj.Filling_1and3_Unc.^2);
            sExp = obj.SinglesDens;
            sExpUnc = obj.SinglesDensUnc;
            dExp = obj.DoublesDensCorr;
            dExpUnc = obj.DoublesDensCorrUnc;
            
            DoublesCorrExp = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            DoublesCorrExpUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExp = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExpUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));

            %remove points larger than cutoff radius
            PExp = PExp(RExp<RCutoff);
            PExpUnc = PExpUnc(RExp<RCutoff);
            nDiffExp = nDiffExp(RExp<RCutoff);
            nDiffExpUnc = nDiffExpUnc(RExp<RCutoff);
            nExp = nExp(RExp<RCutoff);
            nExpUnc = nExpUnc(RExp<RCutoff);
            sExp = sExp(RExp<RCutoff);
            sExpUnc = sExpUnc(RExp<RCutoff);
            dExp = dExp(RExp<RCutoff);
            dExpUnc = dExpUnc(RExp<RCutoff);
            
            DoublesCorrExp = DoublesCorrExp(RExp<RCutoff);
            DoublesCorrExpUnc = DoublesCorrExpUnc(RExp<RCutoff);
            SinglesCorrExp = SinglesCorrExp(RExp<RCutoff);
            SinglesCorrExpUnc = SinglesCorrExpUnc(RExp<RCutoff);
            
            RExp = RExp(RExp<RCutoff);

            %get rid of any problem points in uncertainties
            PExpUnc(PExpUnc==0) = 1;
            PExpUnc(isnan(PExpUnc)) = 1;
            nDiffExpUnc(nDiffExpUnc==0) = 1;
            nDiffExpUnc(isnan(nDiffExpUnc)) = 1;
            nExpUnc(nExpUnc==0) = 1;
            nExpUnc(isnan(nExpUnc)) = 1;
            sExpUnc(sExpUnc==0) = 1;
            sExpUnc(isnan(sExpUnc)) = 1;
            dExpUnc(dExpUnc==0) = 1;
            dExpUnc(isnan(dExpUnc)) = 1;
            DoublesCorrExpUnc(DoublesCorrExpUnc==0) = 1;
            DoublesCorrExpUnc(isnan(DoublesCorrExpUnc)) = 1;
            SinglesCorrExpUnc(SinglesCorrExpUnc==0) = 1;
            SinglesCorrExpUnc(isnan(SinglesCorrExpUnc)) = 1;
            
            %%%DQMC results from file
            n = obj.DQMC2.nUp+obj.DQMC2.nDn;
            n(isnan(n)) = 0;
            ndiff = obj.DQMC2.nUp-obj.DQMC2.nDn;
            ndiff(isnan(ndiff)) = 0;
            p = obj.DQMC2.Pol;
            p(isnan(p)) = 0;
            s = obj.DQMC2.SingOcc;
            d = 0.5*(n - s);
            corrD = 0.25*obj.DQMC2.ddNN;
            corrD(isnan(corrD)) = 0;
            corrS = 0.25*obj.DQMC2.ssNN;
            corrS(isnan(corrS)) = 0;
            
            DelMus = 0.5*(obj.DQMC2.Mu_Ups-obj.DQMC2.Mu_Dns);
            Mus = 0.5*(obj.DQMC2.Mu_Ups+obj.DQMC2.Mu_Dns);
            Us = obj.DQMC2.Us;
            Ts = obj.DQMC2.Ts;
            
            %%%Generate everything as functions Mus,DeltaMus, Us, and Ts
            %instead of Mu_Ups, Mu_Dns, Us, and Ts
            %but we want them at Mus and DeltaMus, instead of Mu1 and Mu2s
            DelMuList = -3.9:0.1:3.9;
            MuList = -3.9:0.1:3.9;
            UList = Us(1,1,:,1);
            TList = Ts(1,1,1,:);
%             [DelMuGrid,MuGrid,UsGrid] = meshgrid(DelMuList,MuList,UList);
            [MuGrid,DelMuGrid,UsGrid,TsGrid] = ndgrid(MuList,DelMuList,UList,TList);

            pOnGrid = zeros(size(MuGrid));
            ndiffOnGrid = zeros(size(MuGrid));
            nOnGrid = zeros(size(MuGrid)); 
            sOnGrid = zeros(size(MuGrid));
            dOnGrid = zeros(size(MuGrid));
            corrDOnGrid = zeros(size(MuGrid));
            corrSOnGrid = zeros(size(MuGrid)); 
            warning('off','all');
            for ii = 1:length(TList)
                %believe duplicate point warnings due to chemical potential
                %regions where n -> 0.
                %simulation points
                DelMusTemp = DelMus(:,:,:,ii);
                MuTemp = Mus(:,:,:,ii);
                UTemp = Us(:,:,:,ii); 
                %sample points
                DelMuGridTemp = DelMuGrid(:,:,:,1);
                MuGridTemp = MuGrid(:,:,:,1);
                UsGridTemp = UsGrid(:,:,:,1);
                pOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,p(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                ndiffOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,ndiff(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                nOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,n(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                sOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,s(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                dOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,d(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                corrDOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,corrD(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                corrSOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,corrS(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
            end
            warning('on','all')
            
            %Generate desired quantites on grid of DeltaMus and Mus
            pOnGrid(isnan(pOnGrid)) = 1;
            ndiffOnGrid(isnan(ndiffOnGrid)) = 1;
            nOnGrid(isnan(nOnGrid)) = 1;
            sOnGrid(isnan(sOnGrid)) = 1;
            dOnGrid(isnan(dOnGrid)) = 1;
            corrDOnGrid(isnan(corrDOnGrid)) = 1;
            corrSOnGrid(isnan(corrSOnGrid)) = 1;
            
            %then create interpolating function. This is much faster using
            %interp3, which assumes on grid points, than griddata, which
            %does not use gridpoints.
            pfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,pOnGrid,mu,dmu,U,T);
            pfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,pfn_interp);
            
            ndifffn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,ndiffOnGrid,mu,dmu,U,T);
            ndifffn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,ndifffn_interp);
            
            nfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,nOnGrid,mu,dmu,U,T);
            nfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,nfn_interp);
            
            sfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,sOnGrid,mu,dmu,U,T);
            sfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,sfn_interp);
            
            dfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,dOnGrid,mu,dmu,U,T);
            dfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,dfn_interp);
            
            corrDfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,corrDOnGrid,mu,dmu,U,T);
            corrDfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,corrDfn_interp);
            
            corrSfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,corrSOnGrid,mu,dmu,U,T);
            corrSfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,corrSfn_interp);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %first fit with no U/t variation across the trap. Can do this
            %in a trap independent way.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            PrelimInitP = InitP([1,3,2]);
            PrelimFixedP = FixedP([1,3,2]);
            PrelimFp = obj.fitDQMC_NoTrap(RCutoff,PrelimInitP,PrelimFixedP);
            
            %get guesses from preliminary fit, unless the parameters are fixed.
            if ~FixedP(1)
                InitP(1) = PrelimFp(1);
            end
            if ~FixedP(2)
                InitP(2) = PrelimFp(3);
            end
            if ~FixedP(3)
                InitP(3) = PrelimFp(2);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now fit including trap variations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            %this expression comes from assuming a gaussian wavefunctionand
            %computing t and U integrals
%             uovertfn = @(P,R) P(4)*(1-0.5*P(3)^2/P(5)*R.^2).^(-3/4).*exp(2*sqrt((P(5)-0.5*P(3)^2*R.^2)))*exp(-2*sqrt(P(5)));
            %but this expression, need the factors in the exponential to be
            %in units of the lattice recoil, not the hopping.
            %need an extra fit parameter for this.
%             uovertfn = @(P,R) P(4)*(1-0.5*P(3)^2/P(5)*R.^2).^(-3/4).*exp(2*sqrt(P(6))*sqrt((P(5)-0.5*P(3)^2*R.^2)))*exp(-2*sqrt(P(6))*sqrt(P(5)));

            %can instead get U/t from trap depth measurements + band
            %structure calculations.
            LattAtten = 0.48; %determined from fit 04/04/2017, for retro PD = 170mV @10V latt servo
            DepthList = linspace(2,20,100); %in Er
            [ts,txs,tys,tdiags,us,~] = Lattice_2D_Asymmetric_Interpolated(DepthList,LattAtten);
            uFn = @(Depth) interp1(DepthList,us,Depth);
            tFn = @(Depth) interp1(DepthList,ts,Depth);
            uOvert_BandStruct = @(Depth) uFn(Depth)./tFn(Depth);
            %from initial parameter, guess central lattice depth
            %GuessLatticeDepthAtCenter = fzero(@(d) uOvert_BandStruct(d)-InitP(3),10);
            
            %previous trap measurements...added zero for convenience
            rads = [0,5.3465,9.7489,12.6254,14.9491,16.9302,18.7301,20.4012,...
                21.9388,23.3554,24.6893,25.9374,27.1337,28.2754,29.3925,30.4709,...
                31.5201,32.5339,33.5149,34.4593,35.3653,36.2474,37.1212,37.9911,...
                38.8321,39.6420];
            vs = [57.7016,57.7016,57.2359,56.7780,56.5207,56.1764,55.9290,...
                55.4548,54.9821,54.6439,54.1269,53.9811,53.5650,53.1263,52.5598,...
                52.4731,52.1904,51.9089,51.5164,51.1449,50.9224,50.6223,50.3907,...
                49.8677,49.8024,49.2935];
            vs = vs/max(vs(:));
            pot = @(R) interp1(rads,vs,R);
            %some fun to generate a function which uses the most natural
            %parameter, the band structure at the center of the trap.
            uovertfn = @(P,R) uOvert_BandStruct(pot(R)*fzero(@(d) uOvert_BandStruct(d)-P(3),10));
           
            %harmonic potential
%             mufn = @(P,R)P(4) - 0.5*P(5)^2*R.^2;
            %use measured potential. Now Omega becomes scale factor between
            %t and Er. i.e. P(5) = Er/t
            mufn = @(P,R)P(4) - 0.5*P(5)*(pot(0) - pot(R))*fzero(@(d) uOvert_BandStruct(d)-P(3),10);
            
            %fitfns
            pfitfn = @(P,R) pfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(2)*ones(size(R)));
            ndifffitfn = @(P,R) ndifffn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(2)*ones(size(R)));
            nfitfn = @(P,R) nfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(2)*ones(size(R)));
            sfitfn = @(P,R) sfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(2)*ones(size(R)));
            dfitfn = @(P,R) dfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(2)*ones(size(R)));
            corrDfitfn = @(P,R) corrDfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(2)*ones(size(R)));
            corrSfitfn = @(P,R) corrSfn(mufn(P,R),P(1)*ones(size(R)),real(uovertfn(P,R)),P(2)*ones(size(R)));
            
            %construct fit function
            pfullfit = @(P) (pfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(PExp))./transpose(PExpUnc);
            ndifffit = @(P) (ndifffitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(nDiffExp))./transpose(nDiffExpUnc);
            nfullfit = @(P) (nfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(nExp))./transpose(nExpUnc);
            sfullfit = @(P) (sfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(sExp))./transpose(sExpUnc);
            dfullfit = @(P) (dfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(dExp))./transpose(dExpUnc);
            corrDfullfit = @(P) (corrDfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(DoublesCorrExp))./transpose(DoublesCorrExpUnc);
            corrSfullfit = @(P) (corrSfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(SinglesCorrExp))./transpose(SinglesCorrExpUnc);
           
            fullfit = @(P) [ndifffit(P),nfullfit(P),sfullfit(P),dfullfit(P),corrDfullfit(P),corrSfullfit(P)];

            LowerBounds = [min(DelMuList),min(TList),min(UList),min(MuList),0];
            UpperBounds = [max(DelMuList),max(TList),max(UList),max(MuList),inf];
            Fp = lsqnonlin(fullfit,InitP,LowerBounds,UpperBounds);
            %inferred parameters
%             OmegaSqrtT = Fp(5)/sqrt(obj.mLi/obj.h)/obj.a;
%             OmegaMeanHz = Fp(3)/sqrt(obj.mLi/abs(t_fit_Hz*obj.h))/obj.a;
            
            Rinterp = linspace(min(RExp),max(RExp),300);

            NRows = 2;
            NCols = 4;
            
            %with trap variation
            figure('name','U/t Variation Across Trap')
            subplot(NRows,NCols,1)
            errorbar(RExp,nDiffExp,nDiffExpUnc,'b.');
            hold on;
            plot(Rinterp,ndifffitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1])
            title('n1-n3')
            
            subplot(NRows,NCols,5)
            errorbar(RExp,PExp,PExpUnc,'b.');
            hold on;
            plot(Rinterp,pfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1])
            title('Local Polarization')
            
            subplot(NRows,NCols,2)
            errorbar(RExp,nExp,nExpUnc,'b.');
            hold on;
            plot(Rinterp,nfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Total Density')
            
            subplot(NRows,NCols,3)
            errorbar(RExp,sExp,sExpUnc,'b.');
            hold on;
            plot(Rinterp,sfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Singles Density')
            
            subplot(NRows,NCols,4)
            errorbar(RExp,dExp,dExpUnc,'b.');
            hold on;
            plot(Rinterp,dfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Doubles Density')
           
            subplot(NRows,NCols,7)
            errorbar(RExp,SinglesCorrExp,SinglesCorrExpUnc,'b.');
            hold on;
            plot(Rinterp,corrSfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([-0.05,0.02]);
            title('Singles Correlator')
            
            subplot(NRows,NCols,8)
            errorbar(RExp,DoublesCorrExp,DoublesCorrExpUnc,'b.');
            hold on;
            plot(Rinterp,corrDfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([-0.05,0.02]);
            title('Doubles Correlator')
            
            th = suptitle(sprintf('%s \nU/t Variation across trap \nDMu = %0.2ft, T = %0.2ft, central U/t = %0.2f, Mu = %0.2ft, Er/t = %0.2f -> t = %0.1f Hz',...
                obj.getDescriptionString(),Fp(1),Fp(2),Fp(3),Fp(4),Fp(5),14.66e3/Fp(5)));
            set(th,'interpreter','none')
            
            %trap variation
            figure('name','Variation of U/t')
            subplot(2,2,1)
            plot(Rinterp,uovertfn(Fp,Rinterp));
            grid on;
            title('U/t variation across trap')
            xlabel('Position (Lattice Sites)')
            ylabel('U/t')
            
            subplot(2,2,2)
            plot(Rinterp,pot(Rinterp)*fzero(@(d) uOvert_BandStruct(d)-Fp(3),10),'r');
            grid on;
             xlabel('Position (Lattice Sites)')
            ylabel('Potential (Er)')
            
            subplot(2,2,3)
            plot(Rinterp,mufn(Fp,Rinterp),'r');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('Mu (t)')
            
            th = suptitle(sprintf('%s \n Trap Variation',obj.getDescriptionString()));
            set(th,'interpreter','none');
        end
        
        function Fp = fitDQMC4(obj,RCutoff,InitP,FixedP)
            %Fit to DQMC data including variation of U/t across the trap.
            %InitP = [DeltaMu,T,U/t @ center,Mu0,Er/t]
            %DeltaMu, Mu0, and T are measured in units of t, the hopping
            %Omega is the trapping frequency for a harmonic trap.
            %Realized that I am not handling the DQMC data correctly in
            %fitDQMC3. I am accounting for the trap variation by letting
            %U/t vary spatially. However, I am not correctly accounting for
            %the fact that t varies across the trap in my treatiment of
            %DeltaMu,T, and Mu.
            
            %%%argument checking
            if ~exist('RCutoff','var')
                RCutoff = 30;
            end
            
            if ~exist('InitP','var')
                InitP = [2.4,0.5,-6,0,8.89];
            end
            
            if ~exist('FixedP','var')
                FixedP = zeros(size(InitP));
            end
            
            if ~isequal(size(FixedP),size(InitP))
                warning('sizes of FixedP and InitP were not the same in fitDQMC3');
                FixedP = zeros(size(InitP));
            end
            
            %%%get experimental measurements
            RExp = obj.RadialPos;
            PExp = obj.LocalPol;
            PExpUnc = obj.LocalPolUnc;
            nDiffExp = (obj.OnesDensCorr-obj.ThreesDensCorr);
            nDiffExpUnc = sqrt(obj.OnesDensCorrUnc.^2+obj.ThreesDensCorrUnc.^2);
            nExp = 0.5*(obj.Filling_SandD+obj.Filling_1and3);
            nExpUnc = 0.5*sqrt(obj.Filling_SandD_Unc.^2+obj.Filling_1and3_Unc.^2);
            sExp = obj.SinglesDens;
            sExpUnc = obj.SinglesDensUnc;
            dExp = obj.DoublesDensCorr;
            dExpUnc = obj.DoublesDensCorrUnc;
            
            DoublesCorrExp = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            DoublesCorrExpUnc = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix,obj.Doubles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExp = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));
            SinglesCorrExpUnc = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix,obj.Singles.CenterIndex_CorrMatrix+1,:));

            %remove points larger than cutoff radius
            PExp = PExp(RExp<RCutoff);
            PExpUnc = PExpUnc(RExp<RCutoff);
            nDiffExp = nDiffExp(RExp<RCutoff);
            nDiffExpUnc = nDiffExpUnc(RExp<RCutoff);
            nExp = nExp(RExp<RCutoff);
            nExpUnc = nExpUnc(RExp<RCutoff);
            sExp = sExp(RExp<RCutoff);
            sExpUnc = sExpUnc(RExp<RCutoff);
            dExp = dExp(RExp<RCutoff);
            dExpUnc = dExpUnc(RExp<RCutoff);
            
            DoublesCorrExp = DoublesCorrExp(RExp<RCutoff);
            DoublesCorrExpUnc = DoublesCorrExpUnc(RExp<RCutoff);
            SinglesCorrExp = SinglesCorrExp(RExp<RCutoff);
            SinglesCorrExpUnc = SinglesCorrExpUnc(RExp<RCutoff);
            
            RExp = RExp(RExp<RCutoff);

            %get rid of any problem points in uncertainties
            PExpUnc(PExpUnc==0) = 1;
            PExpUnc(isnan(PExpUnc)) = 1;
            nDiffExpUnc(nDiffExpUnc==0) = 1;
            nDiffExpUnc(isnan(nDiffExpUnc)) = 1;
            nExpUnc(nExpUnc==0) = 1;
            nExpUnc(isnan(nExpUnc)) = 1;
            sExpUnc(sExpUnc==0) = 1;
            sExpUnc(isnan(sExpUnc)) = 1;
            dExpUnc(dExpUnc==0) = 1;
            dExpUnc(isnan(dExpUnc)) = 1;
            DoublesCorrExpUnc(DoublesCorrExpUnc==0) = 1;
            DoublesCorrExpUnc(isnan(DoublesCorrExpUnc)) = 1;
            SinglesCorrExpUnc(SinglesCorrExpUnc==0) = 1;
            SinglesCorrExpUnc(isnan(SinglesCorrExpUnc)) = 1;
            
            %%%DQMC results from file
            n = obj.DQMC2.nUp+obj.DQMC2.nDn;
            n(isnan(n)) = 0;
            ndiff = obj.DQMC2.nUp-obj.DQMC2.nDn;
            ndiff(isnan(ndiff)) = 0;
            p = obj.DQMC2.Pol;
            p(isnan(p)) = 0;
            s = obj.DQMC2.SingOcc;
            d = 0.5*(n - s);
            corrD = 0.25*obj.DQMC2.ddNN;
            corrD(isnan(corrD)) = 0;
            corrS = 0.25*obj.DQMC2.ssNN;
            corrS(isnan(corrS)) = 0;
            
            DelMus = 0.5*(obj.DQMC2.Mu_Ups-obj.DQMC2.Mu_Dns);
            Mus = 0.5*(obj.DQMC2.Mu_Ups+obj.DQMC2.Mu_Dns);
            Us = obj.DQMC2.Us;
            Ts = obj.DQMC2.Ts;
            
            %%%Generate everything as functions Mus,DeltaMus, Us, and Ts
            %instead of Mu_Ups, Mu_Dns, Us, and Ts
            %but we want them at Mus and DeltaMus, instead of Mu1 and Mu2s
            DelMuList = -3.9:0.1:3.9;
            MuList = -3.9:0.1:3.9;
            UList = Us(1,1,:,1);
            TList = Ts(1,1,1,:);
%             [DelMuGrid,MuGrid,UsGrid] = meshgrid(DelMuList,MuList,UList);
            [MuGrid,DelMuGrid,UsGrid,TsGrid] = ndgrid(MuList,DelMuList,UList,TList);

            pOnGrid = zeros(size(MuGrid));
            ndiffOnGrid = zeros(size(MuGrid));
            nOnGrid = zeros(size(MuGrid)); 
            sOnGrid = zeros(size(MuGrid));
            dOnGrid = zeros(size(MuGrid));
            corrDOnGrid = zeros(size(MuGrid));
            corrSOnGrid = zeros(size(MuGrid)); 
            warning('off','all');
            for ii = 1:length(TList)
                %believe duplicate point warnings due to chemical potential
                %regions where n -> 0.
                %simulation points
                DelMusTemp = DelMus(:,:,:,ii);
                MuTemp = Mus(:,:,:,ii);
                UTemp = Us(:,:,:,ii); 
                %sample points
                DelMuGridTemp = DelMuGrid(:,:,:,1);
                MuGridTemp = MuGrid(:,:,:,1);
                UsGridTemp = UsGrid(:,:,:,1);
                pOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,p(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                ndiffOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,ndiff(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                nOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,n(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                sOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,s(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                dOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,d(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                corrDOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,corrD(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
                corrSOnGrid(:,:,:,ii) = griddata(DelMusTemp,MuTemp,UTemp,corrS(:,:,:,ii),DelMuGridTemp,MuGridTemp,UsGridTemp);
            end
            warning('on','all')
            
            %Generate desired quantites on grid of DeltaMus and Mus
            pOnGrid(isnan(pOnGrid)) = 1;
            ndiffOnGrid(isnan(ndiffOnGrid)) = 1;
            nOnGrid(isnan(nOnGrid)) = 1;
            sOnGrid(isnan(sOnGrid)) = 1;
            dOnGrid(isnan(dOnGrid)) = 1;
            corrDOnGrid(isnan(corrDOnGrid)) = 1;
            corrSOnGrid(isnan(corrSOnGrid)) = 1;
            
            %then create interpolating function. This is much faster using
            %interp3, which assumes on grid points, than griddata, which
            %does not use gridpoints.
            pfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,pOnGrid,mu,dmu,U,T);
            pfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,pfn_interp);
            
            ndifffn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,ndiffOnGrid,mu,dmu,U,T);
            ndifffn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,ndifffn_interp);
            
            nfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,nOnGrid,mu,dmu,U,T);
            nfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,nfn_interp);
            
            sfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,sOnGrid,mu,dmu,U,T);
            sfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,sfn_interp);
            
            dfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,dOnGrid,mu,dmu,U,T);
            dfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,dfn_interp);
            
            corrDfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,corrDOnGrid,mu,dmu,U,T);
            corrDfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,corrDfn_interp);
            
            corrSfn_interp = @(mu,dmu,U,T) interpn(MuGrid,DelMuGrid,UsGrid,TsGrid,corrSOnGrid,mu,dmu,U,T);
            corrSfn = @(mu,dmu,U,T) obj.fn_no_nans(mu,dmu,U,T,corrSfn_interp);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %first fit with no U/t variation across the trap. Can do this
            %in a trap independent way.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            PrelimInitP = InitP([1,3,2]);
            PrelimFixedP = FixedP([1,3,2]);
            PrelimFp = obj.fitDQMC_NoTrap(RCutoff,PrelimInitP,PrelimFixedP);
            
            %get guesses from preliminary fit, unless the parameters are fixed.
            if ~FixedP(1)
                InitP(1) = PrelimFp(1);
            end
            if ~FixedP(2)
                InitP(2) = PrelimFp(3);
            end
            if ~FixedP(3)
                InitP(3) = PrelimFp(2);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now fit including trap variations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            %this expression comes from assuming a gaussian wavefunctionand
            %computing t and U integrals
%             uovertfn = @(P,R) P(4)*(1-0.5*P(3)^2/P(5)*R.^2).^(-3/4).*exp(2*sqrt((P(5)-0.5*P(3)^2*R.^2)))*exp(-2*sqrt(P(5)));
            %but this expression, need the factors in the exponential to be
            %in units of the lattice recoil, not the hopping.
            %need an extra fit parameter for this.
%             uovertfn = @(P,R) P(4)*(1-0.5*P(3)^2/P(5)*R.^2).^(-3/4).*exp(2*sqrt(P(6))*sqrt((P(5)-0.5*P(3)^2*R.^2)))*exp(-2*sqrt(P(6))*sqrt(P(5)));

            %can instead get U/t from trap depth measurements + band
            %structure calculations.
            LattAtten = 0.48; %determined from fit 04/04/2017, for retro PD = 170mV @10V latt servo
            DepthList = linspace(2,20,100); %in Er
            [ts,txs,tys,tdiags,us,~] = Lattice_2D_Asymmetric_Interpolated(DepthList,LattAtten);
            uFn = @(Depth) interp1(DepthList,us,Depth); %-6.28/14.66; %spectroscopy value.%
            tFn = @(Depth) interp1(DepthList,ts,Depth);
            uOvert_BandStruct = @(Depth) uFn(Depth)./tFn(Depth);
            %from initial parameter, guess central lattice depth
            %GuessLatticeDepthAtCenter = fzero(@(d) uOvert_BandStruct(d)-InitP(3),10);
            
            %previous trap measurements...added zero for
            %convenience...maybe should extrapolate to get better value?
            rads = [0,5.3465,9.7489,12.6254,14.9491,16.9302,18.7301,20.4012,...
                21.9388,23.3554,24.6893,25.9374,27.1337,28.2754,29.3925,30.4709,...
                31.5201,32.5339,33.5149,34.4593,35.3653,36.2474,37.1212,37.9911,...
                38.8321,39.6420];
            vs = [57.7016,57.7016,57.2359,56.7780,56.5207,56.1764,55.9290,...
                55.4548,54.9821,54.6439,54.1269,53.9811,53.5650,53.1263,52.5598,...
                52.4731,52.1904,51.9089,51.5164,51.1449,50.9224,50.6223,50.3907,...
                49.8677,49.8024,49.2935];
            vs = vs/max(vs(:));
%             pot = @(R) interp1(rads,vs,R);
            %some fun to generate a function which uses the most natural
            %parameter, the band structure at the center of the trap.
%             uoverttrapfn = @(P,R) uOvert_BandStruct(pot(R)*fzero(@(d) uOvert_BandStruct(d)-P(3),10));
%             relative_tunneling_trapfn = @(P,R) tFn(pot(R)*fzero(@(d) uOvert_BandStruct(d)-P(3),10))/tFn(pot(0)*fzero(@(d) uOvert_BandStruct(d)-P(3),10));
            
            %harmonic potential
%             mufn = @(P,R)P(4) - 0.5*P(5)^2*R.^2;
            %use measured potential. Now Omega becomes scale factor between
            %t and Er. i.e. P(5) = Er/t
            %mutrapfn = @(P,R)(P(4) - P(5)*(pot(0) - pot(R))*fzero(@(d) uOvert_BandStruct(d)-P(3),10))./tunneling_scaled(P,R);
            pottrapfn = @(P,R)interp1(rads,vs,R)*fzero(@(d) uOvert_BandStruct(d)-P(3),10);
            uoverttrapfn = @(P,R) uOvert_BandStruct(pottrapfn(P,R));
            relative_tunneling_trapfn = @(P,R) tFn(pottrapfn(P,R))/tFn(pottrapfn(P,0));
            %mutrapfn = @(P,R)(P(4) - P(5)*(pottrapfn(P,0) - pottrapfn(P,R)))./relative_tunneling_trapfn(P,R);
            %no need to fit for the Er/t in the center. We know it already.
            mutrapfn = @(P,R)(P(4) - 1./tFn(pottrapfn(P,0))*(pottrapfn(P,0) - pottrapfn(P,R)))./relative_tunneling_trapfn(P,R);
            Ttrapfn = @(P,R) P(2)./relative_tunneling_trapfn(P,R);
            DelMutrapfn = @(P,R) P(1)./relative_tunneling_trapfn(P,R);
            
            %fitfns
            pfitfn = @(P,R) pfn(mutrapfn(P,R),DelMutrapfn(P,R),real(uoverttrapfn(P,R)),Ttrapfn(P,R));
            ndifffitfn = @(P,R) ndifffn(mutrapfn(P,R),DelMutrapfn(P,R),real(uoverttrapfn(P,R)),Ttrapfn(P,R));
            nfitfn = @(P,R) nfn(mutrapfn(P,R),DelMutrapfn(P,R),real(uoverttrapfn(P,R)),Ttrapfn(P,R));
            sfitfn = @(P,R) sfn(mutrapfn(P,R),DelMutrapfn(P,R),real(uoverttrapfn(P,R)),Ttrapfn(P,R));
            dfitfn = @(P,R) dfn(mutrapfn(P,R),DelMutrapfn(P,R),real(uoverttrapfn(P,R)),Ttrapfn(P,R));
            corrDfitfn = @(P,R) corrDfn(mutrapfn(P,R),DelMutrapfn(P,R),real(uoverttrapfn(P,R)),Ttrapfn(P,R));
            corrSfitfn = @(P,R) corrSfn(mutrapfn(P,R),DelMutrapfn(P,R),real(uoverttrapfn(P,R)),Ttrapfn(P,R));
            
            %construct fit function
            pfullfit = @(P) (pfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(PExp))./transpose(PExpUnc);
            ndifffit = @(P) (ndifffitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(nDiffExp))./transpose(nDiffExpUnc);
            nfullfit = @(P) (nfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(nExp))./transpose(nExpUnc);
            sfullfit = @(P) (sfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(sExp))./transpose(sExpUnc);
            dfullfit = @(P) (dfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(dExp))./transpose(dExpUnc);
            corrDfullfit = @(P) (corrDfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(DoublesCorrExp))./transpose(DoublesCorrExpUnc);
            corrSfullfit = @(P) (corrSfitfn(P.*(1-FixedP)+InitP.*FixedP,RExp)-transpose(SinglesCorrExp))./transpose(SinglesCorrExpUnc);
           
            fullfit = @(P) [ndifffit(P),nfullfit(P),sfullfit(P),dfullfit(P),corrDfullfit(P),corrSfullfit(P)];

            LowerBounds = [min(DelMuList),min(TList),min(UList),min(MuList),0];
            UpperBounds = [max(DelMuList),max(TList),max(UList),max(MuList),inf];
            Fp = lsqnonlin(fullfit,InitP,LowerBounds,UpperBounds);
            %inferred parameters
%             OmegaSqrtT = Fp(5)/sqrt(obj.mLi/obj.h)/obj.a;
%             OmegaMeanHz = Fp(3)/sqrt(obj.mLi/abs(t_fit_Hz*obj.h))/obj.a;
            
            Rinterp = linspace(min(RExp),max(RExp),300);

            NRows = 2;
            NCols = 4;
            
            %with trap variation
            figure('name','U/t Variation Across Trap')
            subplot(NRows,NCols,1)
            errorbar(RExp,nDiffExp,nDiffExpUnc,'b.');
            hold on;
            plot(Rinterp,ndifffitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1])
            title('n1-n3')
            legend({'Expt','Fit'});
            
            subplot(NRows,NCols,5)
            errorbar(RExp,PExp,PExpUnc,'b.');
            hold on;
            plot(Rinterp,pfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1])
            title('Local Polarization')
            
            subplot(NRows,NCols,2)
            errorbar(RExp,nExp,nExpUnc,'b.');
            hold on;
            plot(Rinterp,nfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Total Density')
            
            subplot(NRows,NCols,3)
            errorbar(RExp,sExp,sExpUnc,'b.');
            hold on;
            plot(Rinterp,sfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Singles Density')
            
            subplot(NRows,NCols,4)
            errorbar(RExp,dExp,dExpUnc,'b.');
            hold on;
            plot(Rinterp,dfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([0,1.5])
            title('Doubles Density')
           
            subplot(NRows,NCols,7)
            errorbar(RExp,SinglesCorrExp,SinglesCorrExpUnc,'b.');
            hold on;
            plot(Rinterp,corrSfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([-0.05,0.02]);
            title('Singles Correlator')
            
            subplot(NRows,NCols,8)
            errorbar(RExp,DoublesCorrExp,DoublesCorrExpUnc,'b.');
            hold on;
            plot(Rinterp,corrDfitfn(Fp,Rinterp),'r');
            grid on;
            ylim([-0.05,0.02]);
            title('Doubles Correlator')
            
            th = suptitle(sprintf('%s \nU/t Variation across trap \nDMu = %0.2ft, T = %0.2ft, central U/t = %0.2f, Mu = %0.2ft, Er/t = %0.2f -> t = %0.1f Hz',...
                obj.getDescriptionString(),Fp(1),Fp(2),Fp(3),Fp(4),Fp(5),14.66e3/Fp(5)));
            set(th,'interpreter','none')
            
            %trap variation
            NRows2 = 2;
            NCols2 = 3;
            figure('name','Variation of U/t')
            subplot(NRows2,NCols2,1)
            plot(Rinterp,uoverttrapfn(Fp,Rinterp));
            grid on;
            title('U/t variation across trap')
            xlabel('Position (Lattice Sites)')
            ylabel('U/t')
            
            subplot(NRows2,NCols2,2)
            plot(Rinterp,pottrapfn(Fp,Rinterp),'r');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('Potential (Er)')
            title('Lattice Depth')
            
            
            subplot(NRows2,NCols2,3)
            plot(Rinterp,relative_tunneling_trapfn(Fp,Rinterp),'r');
            hold on;
            plot(Rinterp,uFn(pottrapfn(Fp,Rinterp))/uFn(pottrapfn(Fp,0)),'b');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('t/t(r) or u/u(r)')
            title('t')
            legend({'t/t(r)','U/U(r)'});
            
            subplot(NRows2,NCols2,4)
            plot(Rinterp,mutrapfn(Fp,Rinterp),'r');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('Mu/t(r))')
            title('Mu')
            
            subplot(NRows2,NCols2,5)
            plot(Rinterp,DelMutrapfn(Fp,Rinterp),'r');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('DeltaMu/t(r)')
            title('DeltaMu')
            
            subplot(NRows2,NCols2,6)
            plot(Rinterp,Ttrapfn(Fp,Rinterp),'r');
            grid on;
            xlabel('Position (Lattice Sites)')
            ylabel('T/t(r)')
            title('T')
            
            th = suptitle(sprintf('%s \n Trap Variation',obj.getDescriptionString()));
            set(th,'interpreter','none');
        end
        
        function Radius = findFilling(obj,Filling)
            %Find radius where we reach a certain filling.
            [XPts,I] = sort([obj.RadialPos(:);-obj.RadialPos(:)]);
            YPts = [obj.Filling_SandD(:);obj.Filling_SandD(:)];
            YPts = YPts(I);
            fn = @(x) interp1(XPts,YPts,x)-Filling;
            RFilling_SD = abs(fzero(fn,0));
            
            YPts = [obj.Filling_1and3(:);obj.Filling_1and3(:)];
            YPts = YPts(I);
            fn2 = @(x) interp1(XPts,YPts,x)-Filling;
            RFilling_1and3 = abs(fzero(fn2,0));
            
            if ~isnan(RFilling_SD) && ~isnan(RFilling_1and3)
                RFilling = 0.5*(RFilling_SD+RFilling_1and3);
            elseif ~isnan(RFilling_1and3)
                RFilling = 0.5*(RFilling_1and3+0);
            elseif ~isnan(RFilling_SD)
                RFilling = 0.5*(0+RFilling_SD);
            else
                RFilling = 0;
            end
            
            
            Radius = RFilling; 
           
            fprintf('%0.2f filling point from singles-doubles inferred density = %0.2f sites \n',Filling,RFilling_SD);
            fprintf('%0.2f filling point from 1-3 inferred density = %0.2f sites \n',Filling,RFilling_1and3);
            fprintf('%0.2f filling point estimate = %0.2f sites \n',Filling,RFilling);
        end
        
        function RHalfFilling = getHalfFilling(obj)
            %RHalfFilling = getHalfFilling(obj)
            RHalfFilling = obj.findFilling(1);
        end
        
        function setOmega(obj,OmegaMean,AspectRatio)
            obj.OmegaMean = OmegaMean;
            obj.AspectRatio = max(AspectRatio,1/AspectRatio);
            obj.OmegaStrong = OmegaMean*sqrt(AspectRatio);
            obj.OmegaWeak = OmegaMean/sqrt(AspectRatio);
        end
        
        function ChemPot = setChemPotFromHalfFilling(obj)
            %taking half filling point as having average \bar{\mu} = 0;
            %This amounts to a U/2 offset compared with the typical hubbard
            %hamiltonian.
            obj.getHalfFilling();
%             mLi = 9.9883e-27; %kg
%             a = 750e-9; %m
            obj.ChemPot_0 = 0.5*obj.mLi*obj.OmegaMean^2*(ojb.a*obj.R_HalfFilling).^2;
            obj.ChemPot = obj.ChemPot_0 - 0.5*obj.mLi*obj.OmegaMean^2*(obj.a*obj.RadialPos).^2;
        end
        
        function ChemPot = setChemPotFromCenter(obj,Mu0)
%             mLi = 9.9883e-27; %kg
%             a = 750e-9; %m
            obj.ChemPot_0 = Mu0;
            obj.ChemPot = obj.ChemPot_0 - 0.5*obj.mLi*obj.OmegaMean^2*(obj.a*obj.RadialPos).^2;
        end
        
        function Compressibility = getCompressibility(obj)
            
            %from ones and threes
            Density = 0.5*(obj.Filling_SandD+obj.Filling_1and3);
            DensityUnc = 0.5*sqrt((obj.Filling_SandD_Unc).^2+(obj.Filling_1and3_Unc).^2);
            
            DensityBtw = 0.5*(obj.RadialPos(1:end-1)+obj.RadialPos(2:end));
            
%             DensityBetweenBins = 0.5*(obj.Filling_SandD(1:end-1)+obj.Filling_SandD(2:end));
%             ChemPotBetweenBins = 0.5*(obj.ChemPot(1:end-1)+obj.ChemPot(2:end));
%             RBetweenBins = 0.5*(obj.BinCenters(1:end-1)+obj.BinCenters(2:end));
            %             Compressibility = 1./DensityBetweenBins.^2.*(obj.Occs_AzAvg(2:end)-obj.Occs_AzAvg(1:end-1))./(obj.ChemPot(2:end)-obj.ChemPot(1:end-1));
%             Compressibility = transpose((obj.Filling_SandD(2:end)-obj.Filling_SandD(1:end-1)))./(obj.ChemPot(2:end)-obj.ChemPot(1:end-1));
            dn_dr = transpose(Density(2:end)-Density(1:end-1))./(obj.RadialPos(2:end)-obj.RadialPos(1:end-1));
            Compressibility = -dn_dr./DensityBtw;

            %dn/du = (-omega^2*r)*dn/dr
%             Compressibility = -obj.mLi*obj.OmegaMean^2*obj.a*RBetweenBins.*transpose(obj.Filling_SandD(2:end)-obj.Filling_SandD(1:end-1))./(obj.BinCenters(2:end)-obj.BinCenters(1:end-1));
            
            fh = figure('name','Compressibility');
%             subplot(2,1,1)
%             plot(ChemPotBetweenBins,Compressibility,'ro');
%             grid on;
%             xlabel('\mu');
%             ylabel('\kappa');
%             subplot(2,1,2)
            plot(DensityBtw,Compressibility,'ro');
            grid on;
            xlabel('n');
            ylabel('\kappa');
        end
        
        function FigHandle = showAzAvg(obj)
            %FigHandle = showAzAvg(obj)
            %Display all data sets azimuthal averages.
            
            Title = sprintf('%d/%d/%d, Folders %d, %d, %d, and %d',obj.Date{1},obj.Date{2},obj.Date{3},obj.Folders(1),obj.Folders(2),obj.Folders(3),obj.Folders(4));
            FigHandle = figure('name',Title);
            NRows = 2;
            NCols = 2;
            
            subplot(NRows,NCols,1)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc);
                hold on;
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc);
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc);
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc);
                ylim([0,1])
                grid on;
                legend({'Doubles','Singles','Threes','Ones'});
                xlabel('Lattice Sites')
                ylabel('Density')
                title('Az Avg');
            subplot(NRows,NCols,2);
                errorbar(obj.RadialPos,obj.LocalPol,obj.LocalPolUnc,'k');
                grid on;
                hold on;
                plot([obj.RadialPos(1),obj.RadialPos(end)],[obj.GlobalPol,obj.GlobalPol],'b-')
                plot([obj.RadialPos(1),obj.RadialPos(end)],[obj.GlobalPol+obj.GlobalPolUnc,obj.GlobalPol+obj.GlobalPolUnc],'b--')
                plot([obj.RadialPos(1),obj.RadialPos(end)],[obj.GlobalPol-obj.GlobalPolUnc,obj.GlobalPol-obj.GlobalPolUnc],'b--')
                ylim([0,1])
                xlabel('Lattice Sites')
                ylabel('Polarization')
                title(sprintf('P_g = %0.2f +/- %0.2f',obj.GlobalPol,obj.GlobalPolUnc));
            subplot(NRows,NCols,3);
                plot(obj.Doubles.PictureNumberInFolder,obj.Doubles.AtomNumbers,'o')
                hold on;
                plot(obj.Singles.PictureNumberInFolder,obj.Singles.AtomNumbers,'o');
                plot(obj.Threes.PictureNumberInFolder,obj.Threes.AtomNumbers,'o');
                plot(obj.Ones.PictureNumberInFolder,obj.Ones.AtomNumbers,'o');
                ax = gca;
                ax.YLim(1) = 0;
                grid on;
                xlabel('Picture # in Folder')
                ylabel('Atom Number')
                legend({'Doubles','Singles','Threes','Ones'});
                title('Atom Number Stability')
                
            subplot(NRows,NCols,4)
                errorbar(obj.RadialPos,obj.SinglesFrac,obj.SinglesFracUnc);
                hold on;
                errorbar(obj.RadialPos,obj.Filling_SandD,obj.Filling_SandD_Unc);
                errorbar(obj.RadialPos,obj.Filling_1and3,obj.Filling_1and3_Unc);
                errorbar(obj.RadialPos,obj.OnesDensCorr-obj.ThreesDensCorr,sqrt((obj.OnesDensCorrUnc).^2+(obj.ThreesDensCorrUnc).^2));
                grid on;
                xlabel('Lattice Sites')
                ylabel('Density')
                title(sprintf('DetEff = %0.2f, 13Dbl Eff = %0.2f, Rf+23Dbl Eff = %0.2f',obj.DetectionEff,obj.BlowDoublesEff13,obj.RFAndBlowDoublesEff23))
                legend({'Ns/(Ns+2Nd)','Ns+2Nd','N1+N3','N1-N3'})
                
            suptitle(Title);
        end
        
        function FigHandle = showNNCorr(obj,Indices)
            
            if ~exist('Indices','var')
                Indices = [0,1];
            end
            
            FigName = sprintf('Correlators Vs. Distance, C(%d,%d)',Indices(1),Indices(2));
            FigHandle = figure('name',FigName);
            NRows = 2;
            NCols = 2;
            
            Axes = subplot(NRows,NCols,1);
            obj.Doubles.showNNCorr(Indices,Axes)
            title(sprintf('Doubles, %d shots',length(obj.Doubles.AtomNumbers)))
            Axes = subplot(NRows,NCols,2);
            obj.Singles.showNNCorr(Indices,Axes)
            title(sprintf('Singles, %d shots',length(obj.Singles.AtomNumbers)))
            Axes = subplot(NRows,NCols,3);
            obj.Threes.showNNCorr(Indices,Axes)
            title(sprintf('Threes, %d shots',length(obj.Threes.AtomNumbers)))
            Axes = subplot(NRows,NCols,4);
            obj.Ones.showNNCorr(Indices,Axes)
            title(sprintf('Ones, %d shots',length(obj.Ones.AtomNumbers)))
             
            suptitle(sprintf('Folders %d, %d, %d, and %d',obj.Folders(1),obj.Folders(2),obj.Folders(3),obj.Folders(4)));

        end
        
        function saveDataSet(obj,OverWrite,Path,FName)
            %save useful info as either mat file or text file, depending on
            %extension...TODO finish implementing the .txt file...
            
            if ~exist('Path','var')
%                 Path = '\\128.112.86.75\lithium\Publications\04_Attractive_Hubb_EoS\Data';
                  Path = '';
            end
            
            
            if ~exist('FName','var')
                FName = sprintf('%s_LattServoV=%0.2f_ScienceFieldV=%0.2f_LattMonAt10V=%0.2f.mat',obj.getDescriptionString(),obj.LattServoV,obj.ScienceFieldV,obj.LattMonAt10V);
            end
            
            if ~(strcmp(FName(end-3:end),'.mat')||strcmp(FName(end-3:end),'.txt'))
                FName = strcat(FName,'.mat');
            end
            FPath = fullfile(Path,FName);
            
            RadPos = obj.RadialPos;
            RadPosUnc = obj.RadialPosUnc;
            
            ThreesDensity = obj.ThreesDensCorr;
            ThreesDensityUnc = obj.ThreesDensCorrUnc;
            OnesDensity = obj.OnesDensCorr;
            OnesDensityUnc = obj.OnesDensCorrUnc;
            DoublesDensity = obj.DoublesDensCorr;
            DoublesDensityUnc = obj.DoublesDensCorrUnc;
            SinglesDensity = obj.SinglesDens;
            SinglesDensityUnc = obj.SinglesDensUnc;
            
            TotalDensity_13s = obj.Filling_1and3;
            TotalDensity_13sUnc = obj.Filling_1and3_Unc;
            TotalDensity_SD = obj.Filling_SandD;
            TotalDensity_SDUnc = obj.Filling_SandD_Unc;
            
            GlobalPolarization = obj.GlobalPol;
            GlobalPolarizationUnc = obj.GlobalPolUnc;
            LocalPolarization = obj.LocalPol;
            LocalPolarizationUnc = obj.LocalPolUnc;
            
            UsingEfficiencies = obj.UseEfficiencies;
            DetEff = obj.DetectionEff;
            EffRFAndBlow23 = obj.RFAndBlowDoublesEff23;
            EffBlow13 = obj.BlowDoublesEff13;
            
            CenterIndex = obj.Threes.CenterIndex_CorrMatrix;
            %for convenience, nearest neighbor correlators
            ThreesNNCorr_Symmetrized = squeeze(obj.ThreesCorr(CenterIndex,CenterIndex+1,:));
            ThreesNNCorr_SymmetrizedUnc = squeeze(obj.ThreesCorrUnc(CenterIndex,CenterIndex+1,:));
            OnesNNCorr_Symmetrized = squeeze(obj.OnesCorr(CenterIndex,CenterIndex+1,:));
            OnesNNCorr_SymmetrizedUnc = squeeze(obj.OnesCorrUnc(CenterIndex,CenterIndex+1,:));
            DoublesNNCorr_Symmetrized = squeeze(obj.DoublesCorr(CenterIndex,CenterIndex+1,:));
            DoublesNNCorr_SymmetrizedUnc = squeeze(obj.DoublesCorrUnc(CenterIndex,CenterIndex+1,:));
            SinglesNNCorr_Symmetrized = squeeze(obj.SinglesCorr(CenterIndex,CenterIndex+1,:));
            SinglesNNCorr_SymmetrizedUnc = squeeze(obj.SinglesCorrUnc(CenterIndex,CenterIndex+1,:));
            
            
            %full correlators, vs. bin.
            ThreesCorrMatrix_Symmetrized = obj.ThreesCorr;
            ThreesCorrMatrix_Symmetrized_Unc = obj.ThreesCorrUnc;
            OnesCorrMatrix_Symmetrized = obj.OnesCorr;
            OnesCorrMatrix_Symmetrized_Unc = obj.OnesCorrUnc;
            DoubleCorrMatrix_Symmetrized = obj.DoublesCorr;
            DoubleCorrMatrix_Symmetrized_Unc = obj.DoublesCorrUnc;
            SinglesCorrMatrix_Symmetrized = obj.SinglesCorr;
            SinglesCorrMatrix_Symmetrized_Unc = obj.SinglesCorrUnc;
            
            %ExcludedFiles
            ExcludedPictures = obj.ExcludePictures;
            
            %Other useful information
            RHalfFilling = obj.getHalfFilling();
            R0p9Filling = obj.findFilling(0.9);
            MuAvgEstimate = obj.ChemPot_0;
            Mu1MinusMu3Estimate_EdgeFit = obj.Mu1MinusMu3;
            
            if ~exist(FPath,'file') || OverWrite
                if strcmp(FPath(end-3:end),'.mat')
                    save(FPath,'RadPos','RadPosUnc','ThreesDensity','ThreesDensityUnc',...
                    'OnesDensity','OnesDensityUnc','DoublesDensity','DoublesDensityUnc',...
                    'SinglesDensity','SinglesDensityUnc','LocalPolarization',...
                    'LocalPolarizationUnc','GlobalPolarization','GlobalPolarizationUnc',...
                    'TotalDensity_13s','TotalDensity_13sUnc','TotalDensity_SD','TotalDensity_SDUnc',...
                    'ThreesNNCorr_Symmetrized','ThreesNNCorr_SymmetrizedUnc',...
                    'OnesNNCorr_Symmetrized','OnesNNCorr_SymmetrizedUnc','DoublesNNCorr_Symmetrized',...
                    'DoublesNNCorr_SymmetrizedUnc','SinglesNNCorr_Symmetrized',...
                    'SinglesNNCorr_SymmetrizedUnc','ThreesCorrMatrix_Symmetrized',...
                    'ThreesCorrMatrix_Symmetrized_Unc','OnesCorrMatrix_Symmetrized',...
                    'OnesCorrMatrix_Symmetrized_Unc','DoubleCorrMatrix_Symmetrized',...
                    'DoubleCorrMatrix_Symmetrized_Unc','SinglesCorrMatrix_Symmetrized',...
                    'SinglesCorrMatrix_Symmetrized_Unc','ExcludedPictures','RHalfFilling',...
                    'R0p9Filling','EffRFAndBlow23','EffBlow13','DetEff','UsingEfficiencies',...
                    'MuAvgEstimate','Mu1MinusMu3Estimate_EdgeFit');
                    fprintf('Saved data to %s \n',FPath);
                elseif strcmp(FPath(end-3:end),'.txt')
                    Delimiter = ',';
                    
                    Names = {'RadPos','RadPosUnc','ThreesDensity','ThreesDensityUnc',...
                    'OnesDensity','OnesDensityUnc','DoublesDensity','DoublesDensityUnc',...
                    'SinglesDensity','SinglesDensityUnc','LocalPolarization',...
                    'LocalPolarizationUnc','TotalDensity_13s','TotalDensity_13sUnc',...
                    'TotalDensity_SD','TotalDensity_SDUnc',...
                    'ThreesNNCorr_Symmetrized','ThreesNNCorr_SymmetrizedUnc',...
                    'OnesNNCorr_Symmetrized','OnesNNCorr_SymmetrizedUnc','DoublesNNCorr_Symmetrized',...
                    'DoublesNNCorr_SymmetrizedUnc','SinglesNNCorr_Symmetrized',...
                    'SinglesNNCorr_SymmetrizedUnc'};
                    FID = fopen(FPath,'w');
                    for ii = 1:length(Names)-1
                        fprintf(FID,'%s%s',Names{ii},Delimiter);
                    end
                    fprintf(FID,'%s\n',Names{end});
                    fclose(FID);
                    Mat = [transpose(RadPos),transpose(RadPosUnc),ThreesDensity,ThreesDensityUnc,...
                        OnesDensity,OnesDensityUnc,DoublesDensity,DoublesDensityUnc,SinglesDensity,...
                        SinglesDensityUnc,LocalPolarization,LocalPolarizationUnc,...
                        TotalDensity_13s,TotalDensity_13sUnc,TotalDensity_SD,TotalDensity_SDUnc,...
                        ThreesNNCorr_Symmetrized,ThreesNNCorr_SymmetrizedUnc,...
                        OnesNNCorr_Symmetrized,OnesNNCorr_SymmetrizedUnc,...
                        DoublesNNCorr_Symmetrized,DoublesNNCorr_SymmetrizedUnc,...
                        SinglesNNCorr_Symmetrized,SinglesNNCorr_SymmetrizedUnc];
                    dlmwrite(FName,Mat,'-append','delimiter',Delimiter);
                    fprintf('Saved data to %s \n',FPath);
                else
                    error('Unsupported file type supplied to save in CorrDataSetAtt');
                end
            end
        end
        
        function saveClass(obj,FName)
            %Save this class object as a .mat file.
             if ~exist('FPath','var')
                FName = sprintf('FullClass_Pg=%0.2f_Folders%03d-%03d_Date=%d_%02d_%02d.mat',obj.GlobalPol,min(obj.Folders),max(obj.Folders),obj.Date{1},obj.Date{2},obj.Date{3});
            end
            
            if ~strcmp(FName(end-3:end),'.mat')
                FName = strcat(FName,'.mat');
            end
            
            save(FName,'obj');
        end
        
        function saveSummaryGraphs(obj,OverWrite,Path)
            %saveSummaryGraphs(obj,OverWrite)
            %If you omit OverWrite or specify 0, will not overwrite files
            %that already exist.
            if ~exist('Path','var')
                Path = getOtherDataPath(datenum(obj.Date{1},obj.Date{2},obj.Date{3}));
            end
            
            if ~exist('OverWrite','var')
                OverWrite = 0;
            end
            
            DensFName = sprintf('%s_Density_Summary_LattServoV=%0.2f_ScienceFieldV=%0.2f_LattMonAt10V=%0.2f.fig',obj.getDescriptionString(),obj.LattServoV,obj.ScienceFieldV,obj.LattMonAt10V);
%             DensFName = sprintf('%03d-%03d_Density_Summary.fig',min(obj.Folders),max(obj.Folders));
            DensPath = fullfile(Path,DensFName);
            DensHandle = obj.showAzAvg();
            
            if ~exist(DensPath,'file')||OverWrite
                savefig(DensHandle,DensPath);
                fprintf('Saved density summary at %s \n',DensPath);
            end
            
             CorrFName = sprintf('%s_NN_Correlator_Summary_LattServoV=%0.2f_ScienceFieldV=%0.2f_LattMonAt10V=%0.2f.fig',obj.getDescriptionString(),obj.LattServoV,obj.ScienceFieldV,obj.LattMonAt10V);
%             CorrFName = sprintf('%03d-%03d_NN_Correlator_Summary.fig',min(obj.Folders),max(obj.Folders));
            CorrPath = fullfile(Path,CorrFName);
            CorrHandle = obj.showNNCorr();
            
            if ~exist(CorrPath,'file')||OverWrite
                savefig(CorrHandle,CorrPath);
                fprintf('Saved NN correlator summary at %s \n',CorrPath);
            end
        end
        
        function compareDensities(obj,objStack)
            %compareDensities(obj,objStack)
            %accepts a stack of CorrDataSetAtt class objects.
            %will compare these all.
            
            FigName = sprintf('Compare Sets');
            Title = sprintf('Set1: %d/%d/%d, %03d-%03d, Pg = %0.2f;',obj.Date{2},obj.Date{3},obj.Date{1},min(obj.Folders),max(obj.Folders),obj.GlobalPol);
            Leg = {'Set1'};
            for ii = 1:length(objStack)
                obj2 = objStack(ii);
                TitleAddition = sprintf(' Set%d: %d/%d/%d, %03d-%03d, Pg = %0.2f',ii+1,obj2.Date{2},obj2.Date{3},obj2.Date{1},min(obj2.Folders),max(obj2.Folders),obj2.GlobalPol);
                Title = strcat(Title,TitleAddition);
                Leg = cat(2,Leg,sprintf('Set%d',ii+1));
            end
            
            FigHandle = figure('name',FigName);
            NRows = 2;
            NCols = 2;
            
            
            subplot(NRows,NCols,1)
                errorbar(obj.RadialPos,obj.OnesDensCorr,obj.OnesDensCorrUnc);
                hold on;
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    errorbar(obj2.RadialPos,obj2.OnesDensCorr,obj2.OnesDensCorrUnc);
                end
                ylim([0,1])
                grid on;
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Ones');
            subplot(NRows,NCols,2)
                errorbar(obj.RadialPos,obj.ThreesDensCorr,obj.ThreesDensCorrUnc);
                hold on;
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    errorbar(obj2.RadialPos,obj2.ThreesDensCorr,obj2.ThreesDensCorrUnc);
                end
                ylim([0,1])
                grid on;
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Threes');
            subplot(NRows,NCols,3)                
                errorbar(obj.RadialPos,obj.SinglesDens,obj.SinglesDensUnc);
                hold on;
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    errorbar(obj2.RadialPos,obj2.SinglesDens,obj2.SinglesDensUnc);
                end
                ylim([0,2])
                grid on;
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Singles');
            subplot(NRows,NCols,4)
                errorbar(obj.RadialPos,obj.DoublesDensCorr,obj.DoublesDensCorrUnc);
                hold on;
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    errorbar(obj2.RadialPos,obj2.DoublesDensCorr,obj2.DoublesDensCorrUnc);
                end
                grid on;
                ylim([0,1]);
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Doubles');
            suptitle(Title);  
        end
        
        function compareCorrelators(obj,objStack)
            %compareCorrelators(obj,objStack)
            %accepts a stack of CorrDataSetAtt class objects.
            %will compare these all.
            FigName = sprintf('Compare Sets');
            Title = sprintf('Set1: %d/%d/%d, %03d-%03d, Pg = %0.2f;',obj.Date{2},obj.Date{3},obj.Date{1},min(obj.Folders),max(obj.Folders),obj.GlobalPol);
            Leg = {'Set1'};
            for ii = 1:length(objStack)
                obj2 = objStack(ii);
                TitleAddition = sprintf(' Set%d: %d/%d/%d, %03d-%03d, Pg = %0.2f',ii+1,obj2.Date{2},obj2.Date{3},obj2.Date{1},min(obj2.Folders),max(obj2.Folders),obj2.GlobalPol);
                Title = strcat(Title,TitleAddition);
                Leg = cat(2,Leg,sprintf('Set%d',ii+1));
            end
            
            FigHandle = figure('name',FigName);
            NRows = 2;
            NCols = 2;
            
            
            subplot(NRows,NCols,1)
                OnesCorr1 = squeeze(obj.OnesCorr(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
                OnesCorrUnc1 = squeeze(obj.OnesCorrUnc(obj.Ones.CenterIndex_CorrMatrix+0,obj.Ones.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,OnesCorr1,OnesCorrUnc1);
                hold on;
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    OnesCorr2 = squeeze(obj2.OnesCorr(obj2.Ones.CenterIndex_CorrMatrix+0,obj2.Ones.CenterIndex_CorrMatrix+1,:));
                    OnesCorrUnc2 = squeeze(obj2.OnesCorrUnc(obj2.Ones.CenterIndex_CorrMatrix+0,obj2.Ones.CenterIndex_CorrMatrix+1,:));
                    errorbar(obj2.RadialPos,OnesCorr2,OnesCorrUnc2);
                end
%                 ylim([0,1])
                grid on;
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Ones');
            subplot(NRows,NCols,2)
                ThreesCorr1 = squeeze(obj.ThreesCorr(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
                ThreesCorrUnc1 = squeeze(obj.ThreesCorrUnc(obj.Threes.CenterIndex_CorrMatrix+0,obj.Threes.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,ThreesCorr1,ThreesCorrUnc1);
                hold on;
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    ThreesCorr2 = squeeze(obj2.ThreesCorr(obj2.Threes.CenterIndex_CorrMatrix+0,obj2.Threes.CenterIndex_CorrMatrix+1,:));
                    ThreesCorrUnc2 = squeeze(obj2.ThreesCorrUnc(obj2.Threes.CenterIndex_CorrMatrix+0,obj2.Threes.CenterIndex_CorrMatrix+1,:));
                    errorbar(obj2.RadialPos,ThreesCorr2,ThreesCorrUnc2);
                end
%                 ylim([0,1])
                grid on;
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Threes');
            subplot(NRows,NCols,3)
                SinglesCorr1 = squeeze(obj.SinglesCorr(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                SinglesCorrUnc1 = squeeze(obj.SinglesCorrUnc(obj.Singles.CenterIndex_CorrMatrix+0,obj.Singles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,SinglesCorr1,SinglesCorrUnc1);
                hold on;
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    SinglesCorr2 = squeeze(obj2.Singles.Density_Corr_AzAvg(obj2.Singles.CenterIndex_CorrMatrix+0,obj2.Singles.CenterIndex_CorrMatrix+1,:));
                    SinglesCorrUnc2 = squeeze(obj2.SinglesCorrUnc(obj2.Singles.CenterIndex_CorrMatrix+0,obj2.Singles.CenterIndex_CorrMatrix+1,:));
                    errorbar(obj2.RadialPos,SinglesCorr2,SinglesCorrUnc2);
                end
%                 ylim([0,2])
                grid on;
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Singles');
            subplot(NRows,NCols,4)
                DoublesCorr1 = squeeze(obj.DoublesCorr(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                DoublesCorrUnc1 = squeeze(obj.DoublesCorrUnc(obj.Doubles.CenterIndex_CorrMatrix+0,obj.Doubles.CenterIndex_CorrMatrix+1,:));
                errorbar(obj.RadialPos,DoublesCorr1,DoublesCorrUnc1);
                hold on;
                
                for ii = 1:length(objStack)
                    obj2 = objStack(ii);
                    DoublesCorr2 = squeeze(obj2.DoublesCorr(obj2.Doubles.CenterIndex_CorrMatrix+0,obj2.Doubles.CenterIndex_CorrMatrix+1,:));
                    DoublesCorrUnc2 = squeeze(obj2.DoublesCorrUnc(obj2.Doubles.CenterIndex_CorrMatrix+0,obj2.Doubles.CenterIndex_CorrMatrix+1,:));
                    errorbar(obj2.RadialPos,DoublesCorr2,DoublesCorrUnc2);
                end
                grid on;
%                 ylim([0,1]);
                legend(Leg);
                xlabel('Lattice Sites')
                ylabel('Filling')
                title('Doubles');
            suptitle(Title);  
        end
    end
    
end

