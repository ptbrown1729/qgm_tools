%Generate functions to use for fitting non-interacting two component fermi
%gas properties.

classdef NonIntFG < handle
    
    properties
        
        %
        IsInitialized = 0
        
        %parameters
        GridSize = 50; %number of points in each grid
        NSites = 100;
        MuEndPts = [-30,30];
        TEndPts = [0.1,10];
        NumNeighborCorrs = 1
        MaxSFactKs = 1;
        
        Kxs_SFact
        Kys_SFact
        %to help keep track of correlators
        Dxs
        Dys
        %Limitation: temperature should be large compared to spacing
        %between discrete energy levels. Can decrease this spacing by
        %increasing NSites.
        
        
        %f(mu1,mu2,T)
        Mu1sGrid %y-like grid
        Mu2sGrid %x-like
        TsGrid %z-like
        n1s
        n2s
        c1s
        c2s
        cmat1s
        cmat2s
        
        %f(MuAvg,DelMu,T)
        MuAvgGrid %y-like
        DelMuGrid %x-like, DelMu = 0.5*(Mu1-Mu2)
        %TsGrid %z-like
        n1_muavg_dmu_T
        n2_muavg_dmu_T
        %n_muavg_dmu_T %should work a little more to eliminate redundant variables...
        %p_muavg_dmu_T %redundant
        c1_muavg_dmu_T
        c2_muavg_dmu_T
        cmat1s_muavg_dmu_T
        cmat2s_muavg_dmu_T
        k_muavg_dmu_T
        chi_muavg_dmu_T
        
        %f(n,DelMu,T)
        nsGrid %y-like
        DelMuGrid2 %x-like
        %TsGrid %z-like
        n1_n_dmu_T
        n2_n_dmu_T
        %p_n_dmu_T %redundant
        c1_n_dmu_T
        c2_n_dmu_T
        cmat1s_n_dmu_T
        cmat2s_n_dmu_T
        k_n_dmu_T
        chi_n_dmu_T
        mu_n_dmu_T
        
        %f(Singles,DelMu,T)
        SinglesGrid %y-like
        %DelMuGrid2 %x-like
        %TsGrid %z-like
        n1_singles_dmu_T
        n2_singles_dmu_T
        %p_singles_dmu_T %redundant
        c1_singles_dmu_T
        c2_singles_dmu_T
        mu_singles_dmu_T
        k_singles_dmu_T
        chi_singles_dmu_T
        
    end
    
    methods
        function obj = NonIntFG(Run)
            if exist('Run','var')
                obj.initialize();
            end
        end
        
        function initialize(obj)
            %Generate data for single FG
            Mus = linspace(obj.MuEndPts(1),obj.MuEndPts(2),obj.GridSize);
            Ts = linspace(obj.TEndPts(1),obj.TEndPts(2),obj.GridSize);
            
            %make full grid of quantities
            %probably some of thsi can be replaced by meshgrid in higher
            %dimension???
            NSites = obj.NSites;
            [Kxs,Kys] = meshgrid((2*pi)/(NSites)*(0:NSites-1),(2*pi)/(NSites)*(0:NSites-1));
            [MusGrid,BetasGrid] = meshgrid(Mus,1./Ts);
            MusAndKGrid = repmat(MusGrid,[ones(1,ndims(MusGrid)),NSites,NSites]);
            BetasAndKGrid = repmat(BetasGrid,[ones(1,ndims(BetasGrid)),NSites,NSites]);
            
            Kxxs = permute(repmat(Kxs,[1,1,size(MusGrid)]),[3:ndims(MusGrid)+2,1,2]);
            Kyys = permute(repmat(Kys,[1,1,size(MusGrid)]),[3:ndims(MusGrid)+2,1,2]);
            
            %calculate quantities
            Es = -2*(cos(Kxxs)+cos(Kyys));
            Occ = 1./(exp(BetasAndKGrid.*(Es - MusAndKGrid))+1);
            NNExp = exp(1i*Kxxs*1+1i*Kyys*0);
            ExpFn = @(dx,dy) exp(1i*Kxxs*dx+1i*Kyys*dy);
            %sum over k's to get useful quantities.
            n = sum(sum(Occ,ndims(MusAndKGrid)-1),ndims(MusAndKGrid))/NSites^2;
            C = -abs(sum(sum(NNExp.*Occ,ndims(MusAndKGrid)-1),ndims(MusAndKGrid))/NSites^2).^2;
            
            [Dxs,Dys] = meshgrid(0:obj.NumNeighborCorrs,0:obj.NumNeighborCorrs);
            obj.Dxs = Dxs;
            obj.Dys = Dys;
            
            CMat = zeros([size(n),size(Dxs)]);
            for ii = 1:numel(Dxs)
                CMat(:,:,ii) = -abs(sum(sum(ExpFn(Dxs(ii),Dys(ii)).*Occ,ndims(MusAndKGrid)-1),ndims(MusAndKGrid))/NSites^2).^2;
            end
            
            %S(q) = -sum_k <n_(k-q) n_k>
%             StructFactMat = zeros(size(Kxxs));
            %is circshift correct here???
            if size(Kxs,2)>obj.MaxSFactKs
                XStep = floor(size(Kxs,2)/obj.MaxSFactKs);
            end
            if size(Kys,1)>obj.MaxSFactKs
                YStep = floor(size(Kys,1)/obj.MaxSFactKs);
            end
            
            [KShiftsX,KShiftsY] = meshgrid(0:XStep:size(Kxs,2),0:YStep:size(Kys,2));
%             
%             Kxs_SFact = 0;
%             Kys_SFact = 0;
%             obj.Kxs_SFact = Kxs_SFact;
%             obj.Kys_SFact = Kys_SFact;
            
            StructFact = zeros([size(n),size(KShiftsX)]);
            for ii = 1:numel(KShiftsX)
                StructFact(:,:,ii) = -sum(sum(Occ.*circshift(Occ,[0,0,-KShiftsX(ii),-KShiftsY(ii)]),ndims(MusAndKGrid)-1),ndims(MusAndKGrid))/NSites^2;
            end
                          
            %get rid of memory intensive variables.
            clear Kxs;
            clear Kys;
            clear Kxxs;
            clear Kyys;
            clear MusGrid;
            clear BetasGrid;
            clear Occ;
            clear NNExp;
            clear Es;
            clear BetasAndKGrid;
            clear MusAndKGrid;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %f(Mu1,Mu2,T)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Recast this in shape for two non-interacting fermi gases.
            %Full matrix for interpolation...
            [Mu2sGridd,Mu1sGridd,TsGridd] = meshgrid(Mus,Mus,Ts);
            n1s = permute(repmat(n,[1,1,length(Mus)]),[2,3,1]);
            c1s = permute(repmat(C,[1,1,length(Mus)]),[2,3,1]);
            cmat1s = permute(repmat(CMat,[1,1,1,1,length(Mus)]),[2,5,1,3,4]);
            
            n2s = permute(repmat(n,[1,1,length(Mus)]),[3,2,1]);
            c2s = permute(repmat(C,[1,1,length(Mus)]),[3,2,1]);
            cmat2s = permute(repmat(CMat,[1,1,1,1,length(Mus)]),[5,2,1,3,4]);
            
            ps = (n1s-n2s)./(n1s+n2s);
            ns = n1s+n2s;
            
            obj.Mu1sGrid = Mu1sGridd;
            obj.Mu2sGrid = Mu2sGridd;
            obj.TsGrid = TsGridd;
            obj.n1s = n1s;
            obj.n2s = n2s;
            obj.c1s = c1s;
            obj.c2s = c2s;
            obj.cmat1s = cmat1s;
            obj.cmat2s = cmat2s;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %f(MuAvg,DelMu,T)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %in mu/delmu, no longer have the full range of mu1 and mu2. Have the full half
            %range, and then some extra extension (think of this as rotating
            %coordinates by 45deg, losing the corners).
            %create a larger grid than want initially, then shrink it. This
            %way I can get quantities which are derivatives at the same
            %points as quantities which aren't.
            AvgMusList = linspace(min(Mus+1)/2,max(Mus-1)/2,obj.GridSize+1);
            DelMusList = linspace(min(Mus+1)/2,max(Mus-1)/2,obj.GridSize+1);
            [ExtDelMuGridd,ExtMuGridd,ExtTsGridd] = meshgrid(DelMusList,AvgMusList,Ts);
            DelMuGridd = 0.5*(ExtDelMuGridd(2:end,2:end,:)+ExtDelMuGridd(1:end-1,1:end-1,:));
            MuGridd = 0.5*(ExtMuGridd(2:end,2:end,:)+ExtMuGridd(1:end-1,1:end-1,:));
            TsGridd2 = 0.5*(ExtTsGridd(2:end,2:end,:)+ExtTsGridd(1:end-1,1:end-1,:));
            
            
            obj.DelMuGrid = DelMuGridd;
            obj.MuAvgGrid = MuGridd;
            if ~isequal(obj.TsGrid,TsGridd2)
                error('TGridd and obj.TsGrid not equal');
            end
            
            %interpolate to get points on the regular grid. This sort of
            %interpolation is slow.
            obj.n1_muavg_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,n1s,DelMuGridd,MuGridd,TsGridd);
            obj.n2_muavg_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,n2s,DelMuGridd,MuGridd,TsGridd);
            %obj.n_muavg_dmu_T = obj.n1_muavg_dmu_T + obj.n2_muavg_dmu_T;
            %obj.n_muavg_dmu_T = griddata(0.5*(Mu2sGridd+Mu1sGridd),0.5*(Mu1sGridd-Mu2sGridd),TsGridd,ns,MuGridd,DelMuGridd,TsGridd);
            %obj.p_muavg_dmu_T = (obj.n1_muavg_dmu_T - obj.n2_muavg_dmu_T)./obj.n_muavg_dmu_T;
            %obj.p_muavg_dmu_T = griddata(0.5*(Mu1sGridd+Mu2sGridd),0.5*(Mu1sGridd-Mu2sGridd),TsGridd,ps,MuGridd,DelMuGridd,TsGridd);
            obj.c1_muavg_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,c1s,DelMuGridd,MuGridd,TsGridd);
            obj.c2_muavg_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,c2s,DelMuGridd,MuGridd,TsGridd);
            
            cmat1s_muavg_dmu_T = zeros([size(Mu1sGridd),size(Dxs)]);
            cmat2s_muavg_dmu_T = zeros([size(Mu1sGridd),size(Dxs)]);
            for ii = 1:numel(Dxs)
                cmat1s_muavg_dmu_T(:,:,:,ii) = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,cmat1s(:,:,:,ii),DelMuGridd,MuGridd,TsGridd);
                cmat2s_muavg_dmu_T(:,:,:,ii) = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,cmat2s(:,:,:,ii),DelMuGridd,MuGridd,TsGridd);
            end
            
            
            obj.cmat1s_muavg_dmu_T = cmat1s_muavg_dmu_T;
            obj.cmat2s_muavg_dmu_T = cmat2s_muavg_dmu_T;
            
            %derived functions. First need some expanded quantities to take
            %derivatives.

            %TODO: start with n1 and n2 expanded...use those to generate
            %both this and the desired obj.n1_n_dmu_T's above. Speed up
            %code some.
            nexpanded = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,ns,ExtDelMuGridd,ExtMuGridd,ExtTsGridd);
            n1_minus_n2_expanded = griddata(0.5*(Mu1sGridd-Mu2sGridd),0.5*(Mu1sGridd+Mu2sGridd),TsGridd,n1s-n2s,ExtDelMuGridd,ExtMuGridd,ExtTsGridd);
            
            %dndmu at constant delta mu and T. Effectively evaluated at
            %final points for T and Mu. But still one extra point in DelMu
            %direction.
            dndmuExpanded = (nexpanded(2:end,:,:)-nexpanded(1:end-1,:,:))./(ExtMuGridd(2:end,:,:)-ExtMuGridd(1:end-1,:,:));
            %to get correct DelMu points, average two nearest ones.
            obj.k_muavg_dmu_T = (1./(obj.n1_muavg_dmu_T+obj.n2_muavg_dmu_T)).^2.*(dndmuExpanded(:,2:end,:)+dndmuExpanded(:,1:end-1,:))*0.5;
            %for chi the situation is opposite. First get at final pionts
            %for Ta nd delMu, but still one extra point in Mu directino.
            dmdhExpanded = (n1_minus_n2_expanded(:,2:end,:)-n1_minus_n2_expanded(:,1:end-1,:))./(ExtDelMuGridd(:,2:end,:)-ExtDelMuGridd(:,1:end-1,:));
            %to get correct Mu points, average two nearest ones.
            obj.chi_muavg_dmu_T = 0.5*(dmdhExpanded(2:end,:,:)+dmdhExpanded(1:end-1,:,:));
            %now we are evaluating k and chi at the same points as the other quantities
            %             obj.k_muavg_dmu_T = griddata(0.5*(Mu1sGridd+Mu2sGridd),0.5*(Mu1sGridd-Mu2sGridd),TsGridd,k,MuGridd,DelMuGridd,TsGridd);
            %             obj.chi_muavg_dmu_T = griddata(0.5*(Mu1sGridd+Mu2sGridd),0.5*(Mu1sGridd-Mu2sGridd),TsGridd,chi,MuGridd,DelMuGridd,TsGridd);
            clear nexpanded;
            clear n1_minus_n2_expanded;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %f(n,DelMu,T)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            delmusList = linspace(min(Mus+1)/2,max(Mus-1)/2,obj.GridSize+1);%linspace(min(Mus),max(Mus),obj.GridSize+1);
            nsList = linspace(0.01,1.99,obj.GridSize);
            [ExtDelMuGrid2,ExtnsGrid,ExtTsGrid3] = meshgrid(delmusList,nsList,Ts);
            DelMuGridd2 = 0.5*(ExtDelMuGrid2(:,2:end,:)+ExtDelMuGrid2(:,1:end-1,:));
            nsGridd = 0.5*(ExtnsGrid(:,2:end,:)+ExtnsGrid(:,1:end-1,:));
            TsGrid3 = 0.5*(ExtTsGrid3(:,2:end,:)+ExtTsGrid3(:,1:end-1,:));
            
            obj.DelMuGrid2 = DelMuGridd2;
            obj.nsGrid = nsGridd;
            if ~isequal(TsGridd,TsGrid3)
                error('obj.TsGrid not equal to TsGrid3');
            end
            clear TsGrid3;
            
            %First interpolation. Very slow because points are not equally spaced.
            obj.n1_n_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,n1s,DelMuGridd2,nsGridd,TsGridd);
            obj.n2_n_dmu_T = nsGridd - obj.n1_n_dmu_T;
            %obj.n2_n_dmu_T = griddata(ns,0.5*(Mu1sGridd-Mu2sGridd),TsGridd,n2s,NsGrid,DelMusGrid,TsGrid2);
            %obj.p_n_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,ps,DelMuGridd2,nsGridd,TsGridd);
            obj.c1_n_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,c1s,DelMuGridd2,nsGridd,TsGridd);
            obj.c2_n_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,c2s,DelMuGridd2,nsGridd,TsGridd);
            obj.k_n_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,obj.k_muavg_dmu_T,DelMuGridd2,nsGridd,TsGridd); %Compressfn(NsGrid,DelMusGrid,TsGrid2);
            obj.mu_n_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,0.5*(Mu1sGridd+Mu2sGridd),DelMuGridd2,nsGridd,TsGridd);
            
            cmat1s_n_dmu_T = zeros([size(Mu1sGridd),size(Dxs)]);
            cmat2s_n_dmu_T = zeros([size(Mu1sGridd),size(Dxs)]);
            for ii = 1:numel(Dxs)
                cmat1s_n_dmu_T(:,:,:,ii) = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,cmat1s(:,:,:,ii),DelMuGridd2,nsGridd,TsGridd);
                cmat2s_n_dmu_T(:,:,:,ii) = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,cmat2s(:,:,:,ii),DelMuGridd2,nsGridd,TsGridd);
            end
            obj.cmat1s_n_dmu_T = cmat1s_n_dmu_T;
            obj.cmat2s_n_dmu_T = cmat2s_n_dmu_T;
            
            %TODO: start with n1 and n2 expanded...use those to generate
            %both this and the desired obj.n1_n_dmu_T's above. Speed up
            %code some.
            %chi
            n1_minus_n2_expanded = griddata(0.5*(Mu1sGridd-Mu2sGridd),ns,TsGridd,n1s-n2s,ExtDelMuGrid2,ExtnsGrid,ExtTsGrid3);
            %for spin 1/2...
            obj.chi_n_dmu_T = 0.5*(n1_minus_n2_expanded(:,2:end,:) - n1_minus_n2_expanded(:,1:end-1,:))./(ExtDelMuGrid2(:,2:end,:)-ExtDelMuGrid2(:,1:end-1,:));
            
            clear n1_minus_n2_expanded
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %f(Singles,DelMu,T)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Actually this may be quite problematic, because Singles is not
            %single valued over chemical potential range.
            %Maybe some way to fix this by restricting values???
            
            delmusList = linspace(min(Mus+1)/2,max(Mus-1)/2,obj.GridSize+1);%linspace(min(Mus),max(Mus),obj.GridSize+1);
            %delmusList = linspace(min(Mus),max(Mus),obj.GridSize+1);
            nSinglesList = linspace(0,1,obj.GridSize);
            [ExtDelMuGridd3,ExtSinglesGridd,ExtTsGrid4] = meshgrid(delmusList,nSinglesList,Ts);
            DelMuGridd3 = 0.5*(ExtDelMuGridd3(:,2:end,:)+ExtDelMuGridd3(:,1:end-1,:));
            SinglesGridd = 0.5*(ExtSinglesGridd(:,2:end,:)+ExtSinglesGridd(:,1:end-1,:));
            TsGrid4 = 0.5*(ExtTsGrid4(:,2:end,:)+ExtTsGrid4(:,1:end-1,:));
          
            
            obj.SinglesGrid = SinglesGridd;
            if ~isequal(TsGridd,TsGrid4)
                error('obj.TsGrid not equal to TsGrid4');
            end
            if ~isequal(DelMuGridd2,DelMuGridd3)
                error('obj.DelMuGrid2 and DelMuGridd3 not equal');
            end
            clear TsGrid4;
            clear DelMuGridd3;
            
            obj.n1_singles_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),n1s+n2s-2*n1s.*n2s,TsGridd,n1s,DelMuGridd2,SinglesGridd,TsGridd);
            obj.n2_singles_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),n1s+n2s-2*n1s.*n2s,TsGridd,n2s,DelMuGridd2,SinglesGridd,TsGridd);
            %obj.p_singles_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),n1s+n2s-2*n1s.*n2s,TsGridd,ps,DelMuGridd2,nsGridd,TsGridd);
            obj.c1_singles_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),n1s+n2s-2*n1s.*n2s,TsGridd,c1s,DelMuGridd2,SinglesGridd,TsGridd);
            obj.c2_singles_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),n1s+n2s-2*n1s.*n2s,TsGridd,c2s,DelMuGridd2,SinglesGridd,TsGridd);
            obj.mu_singles_dmu_T = griddata(0.5*(Mu1sGridd-Mu2sGridd),n1s+n2s-2*n1s.*n2s,TsGridd,0.5*(Mu1sGridd+Mu2sGridd),DelMuGridd2,SinglesGridd,TsGridd);
            obj.k_singles_dmu_T = 0;
            obj.chi_singles_dmu_T = 0;
            
            %TODO, finish implement k and chi. Basically need to copy
            %f(n,DelMu,T) code here with a few name changes.
            
            obj.IsInitialized = 1;
            
            
        end
        
        %now, get the real functions we want, which interpolate on a regular grid,
        %and hence are faster.
        
        function save(obj)
            FName = 'NonIntFG.mat';
            data = obj;
            save(FName,'data');
        end
        
        %Functions (mu1,mu2,T)
        function n1 = n1fn_mu1_mu2_T(obj,mu1,mu2,T)
            n1 = interp3(obj.Mu2sGrid,obj.Mu1sGrid,obj.TsGrid,obj.n1s,mu2,mu1,T);
%             n1 = interp3(obj.Mu1sGrid,obj.Mu2sGrid,obj.TsGrid,obj.n1s,mu1,mu2,T);
            n1(isnan(n1)) = 0;
        end
        
        function n2 = n2fn_mu1_mu2_T(obj,mu1,mu2,T)
            n2 = interp3(obj.Mu2sGrid,obj.Mu1sGrid,obj.TsGrid,obj.n2s,mu2,mu1,T);
%             n2 = interp3(obj.Mu1sGrid,obj.Mu2sGrid,obj.TsGrid,obj.n2s,mu1,mu2,T);
            n2(isnan(n2)) = 0;
        end
        
        function n = nfn_mu1_mu2_T(obj,mu1,mu2,T)
            %n = obj.n1(mu1,mu2,T) + obj.n2(mu1,mu2,T);
            %guessing it is faster to interpolate once, rather than twice.
            npts = obj.n1s + obj.n2s;
            n = interp3(obj.Mu2sGrid,obj.Mu1sGrid,obj.TsGrid,npts,mu2,mu1,T);
%             n = interp3(obj.Mu1sGrid,obj.Mu2sGrid,obj.TsGrid,npts,mu1,mu2,T);
            n(isnan(n)) = 0;
        end
        
        function ns = nsfn_mu1_mu2_T(obj,mu1,mu2,T)
            ns = obj.n1fn_mu1_mu2_T(mu1,mu2,T)+obj.n2fn_mu1_mu2_T(mu1,mu2,T) - 2*obj.n1fn_mu1_mu2_T(mu1,mu2,T).*obj.n2fn_mu1_mu2_T(mu1,mu2,T);
            ns(isnan(ns)) = 0;
        end
        
        function p = pfn_mu1_mu2_T(obj,mu1,mu2,T)
            %p = (obj.n1(mu1,mu2,T) + obj.n2(mu1,mu2,T))./(obj.n1(mu1,mu2,T) - obj.n2(mu1,mu2,T));
            ppts = (obj.n1s - obj.n2s)./(obj.n1s + obj.n2s);
            p = interp3(obj.Mu2sGrid,obj.Mu1sGrid,obj.TsGrid,ppts,mu2,mu1,T);
%             p = interp3(obj.Mu1sGrid,obj.Mu2sGrid,obj.TsGrid,ppts,mu1,mu2,T);
            p(isnan(p)) = 0;
        end
        
        function c1 = c1fn_mu1_mu2_T_All(obj,mu1,mu2,T,dx,dy)
            c1 = squeeze(interp3(obj.Mu2sGrid,obj.Mu1sGrid,obj.TsGrid,obj.cmat1s(:,:,:,dx+1,dy+1),mu2,mu1,T));
%             c1 = interp3(obj.Mu1sGrid,obj.Mu2sGrid,obj.TsGrid,obj.cmat1s(:,:,:,dx+1,dy+1),mu1,mu2,T);
            c1(isnan(c1)) = 0;
        end
        
        function c1 = c1fn_mu1_mu2_T(obj,mu1,mu2,T)
            c1 = interp3(obj.Mu2sGrid,obj.Mu1sGrid,obj.TsGrid,obj.c1s,mu2,mu1,T);
%             c1 = interp3(obj.Mu1sGrid,obj.Mu2sGrid,obj.TsGrid,obj.c1s,mu1,mu2,T);
            c1(isnan(c1)) = 0;
        end
        
        function c2 = c2fn_mu1_mu2_T(obj,mu1,mu2,T)
            c2 = interp3(obj.Mu2sGrid,obj.Mu1sGrid,obj.TsGrid,obj.c2s,mu2,mu1,T);
%             c2 = interp3(obj.Mu1sGrid,obj.Mu2sGrid,obj.TsGrid,obj.c2s,mu1,mu2,T);
            c2(isnan(c2)) = 0;
        end
        
        function cd = cdfn_mu1_mu2_T(obj,mu1,mu2,T)
            %some combinatorial Wick's thrm manipulation leads to...
            %<d_i d_j>_c = n_dn^2*<n_up_i n_up_j>_c + n_up^2*<n_down_in_down_j>_c + <n_down_i n_down_j><n_up_i n_up_j>
            cd = obj.n1fn_mu1_mu2_T(mu1,mu2,T).^2.*obj.c2fn_mu1_mu2_T(mu1,mu2,T)...
                +obj.n2fn_mu1_mu2_T(mu1,mu2,T).^2.*obj.c1fn_mu1_mu2_T(mu1,mu2,T)...
                +obj.c1fn_mu1_mu2_T(mu1,mu2,T).*obj.c2fn_mu1_mu2_T(mu1,mu2,T);
            cd(isnan(cd)) = 0;
        end
        
        function cs =csfn_mu1_mu2_T(obj,mu1,mu2,T)
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
             
            cs = obj.c1fn_mu1_mu2_T(mu1,mu2,T)+obj.c2fn_mu1_mu2_T(mu1,mu2,T)...
                + 4*obj.cdfn_mu1_mu2_T(mu1,mu2,T)...
                - 4*obj.n1fn_mu1_mu2_T(mu1,mu2,T).*obj.c2fn_mu1_mu2_T(mu1,mu2,T)...
                -4*obj.n2fn_mu1_mu2_T(mu1,mu2,T).*obj.c1fn_mu1_mu2_T(mu1,mu2,T);
            cs(isnan(cs)) = 0;
        end
        
        %Functions (mu,dmu,T)
        function n1 = n1fn_mu_dmu_T(obj,mu,dmu,T)
            n1 = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,obj.n1_muavg_dmu_T,dmu,mu,T);
%             n1 = interp3(obj.MuAvgGrid,obj.DelMuGrid,obj.TsGrid,obj.n1_muavg_dmu_T,mu,dmu,T);
            n1(isnan(n1)) = 0;
        end
        function n2 = n2fn_mu_dmu_T(obj,mu,dmu,T)
            n2 = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,obj.n2_muavg_dmu_T,dmu,mu,T);
%             n2 = interp3(obj.MuAvgGrid,obj.DelMuGrid,obj.TsGrid,obj.n2_muavg_dmu_T,mu,dmu,T);
            n2(isnan(n2)) = 0;
        end
        
        function n = nfn_mu_dmu_T(obj,mu,dmu,T)
            npts = obj.n1_muavg_dmu_T+obj.n2_muavg_dmu_T;
            n = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,npts,dmu,mu,T);
%             n = interp3(obj.MuAvgGrid,obj.DelMuGrid,obj.TsGrid,npts,mu,dmu,T);
            n(isnan(n)) = 0;
        end
        
        function ns = nsfn_mu_dmu_T(obj,mu,dmu,T)
            ns = obj.nfn_mu_dmu_T(mu,dmu,T) - 2*obj.n1fn_mu_dmu_T(mu,dmu,T).*obj.n2fn_mu_dmu_T(mu,dmu,T);
            ns(isnan(ns)) = 0;
        end
        
        function p = pfn_mu_dmu_T(obj,mu,dmu,T)
            ppts = (obj.n1_muavg_dmu_T-obj.n2_muavg_dmu_T)./(obj.n1_muavg_dmu_T+obj.n2_muavg_dmu_T);
            p = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,ppts,dmu,mu,T);
%             p = interp3(obj.MuAvgGrid,obj.DelMuGrid,obj.TsGrid,ppts,mu,dmu,T);
            p(isnan(p)) = 0;
        end
        
        function c1 = c1fn_mu_dmu_T(obj,mu,dmu,T)
            c1 = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,obj.c1_muavg_dmu_T,dmu,mu,T);
%             c1 = interp3(obj.MuAvgGrid,obj.DelMuGrid,obj.TsGrid,obj.c1_muavg_dmu_T,mu,dmu,T);
            c1(isnan(c1)) = 0;
        end
        
        function c2 = c2fn_mu_dmu_T(obj,mu,dmu,T)
            c2 = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,obj.c2_muavg_dmu_T,dmu,mu,T);
%             c2 = interp3(obj.MuAvgGrid,obj.DelMuGrid,obj.TsGrid,obj.c2_muavg_dmu_T,mu,dmu,T);
            c2(isnan(c2)) = 0;
        end
        
        function cd = cdfn_mu_dmu_T(obj,mu,dmu,T)
            cd = obj.n1fn_mu_dmu_T(mu,dmu,T).^2.*obj.c2fn_mu_dmu_T(mu,dmu,T)...
                +obj.n2fn_mu_dmu_T(mu,dmu,T).^2.*obj.c1fn_mu_dmu_T(mu,dmu,T)...
                +obj.c1fn_mu_dmu_T(mu,dmu,T).*obj.c2fn_mu_dmu_T(mu,dmu,T);
            cd(isnan(cd)) = 0;
        end
        
        function cs =csfn_mu_dmu_T(obj,mu,dmu,T)
            cs = obj.c1fn_mu_dmu_T(mu,dmu,T)+obj.c2fn_mu_dmu_T(mu,dmu,T)...
                + 4*obj.cdfn_mu_dmu_T(mu,dmu,T)...
                - 4*obj.n1fn_mu_dmu_T(mu,dmu,T).*obj.c2fn_mu_dmu_T(mu,dmu,T)...
                -4*obj.n2fn_mu_dmu_T(mu,dmu,T).*obj.c1fn_mu_dmu_T(mu,dmu,T);
            cs(isnan(cs)) = 0;
        end
        
        function k = kfn_mu_dmu_T(obj,mu,dmu,T)
            k = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,obj.k_muavg_dmu_T,dmu,mu,T);
%             k = interp3(obj.MuAvgGrid,obj.DelMuGrid,obj.TsGrid,obj.k_muavg_dmu_T,mu,dmu,T);
            k(isnan(k)) = 0;
        end
        
        function chi = chifn_mu_dmu_T(obj,mu,dmu,T)
            chi = interp3(obj.DelMuGrid,obj.MuAvgGrid,obj.TsGrid,obj.chi_muavg_dmu_T,dmu,mu,T);
            chi(isnan(chi)) = 0;
        end
        
        %Functions (n,dmu,T)
        function n1 = n1fn_n_dmu_T(obj,n,dmu,T)
            n1 = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.n1_n_dmu_T,dmu,n,T);
            n1(isnan(n1)) = 0;
        end
        
        function n2 = n2fn_n_dmu_T(obj,n,dmu,T)
            n2 = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.n2_n_dmu_T,dmu,n,T);
            n2(isnan(n2)) = 0;
        end
        
        function ns = nsfn_n_dmu_T(obj,n,dmu,T)
            ns = obj.n1fn_n_dmu_T(n,dmu,T) + obj.n2fn_n_dmu_T(n,dmu,T) - ...
                2*obj.n1fn_n_dmu_T(n,dmu,T).*obj.n2fn_n_dmu_T(n,dmu,T);
            ns(isnan(ns)) = 0;
        end
        
        function p = pfn_n_dmu_T(obj,n,dmu,T)
            %p = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.p_n_dmu_T,dmu,n,T);
            ppts = (obj.n1_n_dmu_T-obj.n2_n_dmu_T)./(obj.n1_n_dmu_T+obj.n2_n_dmu_T);
            p = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,ppts,dmu,n,T);
            p(isnan(p)) = 0;
        end
        
        function c1 = c1fn_n_dmu_T(obj,n,dmu,T)
            c1 = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.c1_n_dmu_T,dmu,n,T);
            c1(isnan(c1)) = 0;
        end
        
        function c2 = c2fn_n_dmu_T(obj,n,dmu,T)
            c2 = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.c2_n_dmu_T,dmu,n,T);
            c2(isnan(c2)) = 0;
        end
        
        function cd = cdfn_n_dmu_T(obj,n,dmu,T)
            cd = obj.n1fn_n_dmu_T(n,dmu,T).^2.*obj.c2fn_n_dmu_T(n,dmu,T)...
                +obj.n2fn_n_dmu_T(n,dmu,T).^2.*obj.c1fn_n_dmu_T(n,dmu,T)...
                +obj.c1fn_n_dmu_T(n,dmu,T).*obj.c2fn_n_dmu_T(n,dmu,T);
            cd(isnan(cd)) = 0;
        end
        
        function cs =csfn_n_dmu_T(obj,n,dmu,T)
            cs = obj.c1fn_n_dmu_T(n,dmu,T)+obj.c2fn_n_dmu_T(n,dmu,T)...
                + 4*obj.cdfn_n_dmu_T(n,dmu,T)...
                - 4*obj.n1fn_n_dmu_T(n,dmu,T).*obj.c2fn_n_dmu_T(n,dmu,T)...
                -4*obj.n2fn_n_dmu_T(n,dmu,T).*obj.c1fn_n_dmu_T(n,dmu,T);
            cs(isnan(cs)) = 0;
        end
        
        function muavg = mufn_n_dmu_T(obj,n,dmu,T)
            muavg = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.mu_n_dmu_T,dmu,n,T);
            muavg(isnan(muavg)) = 0;
        end
        
        function k = kfn_n_dmu_T(obj,n,dmu,T)
            k = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.k_n_dmu_T,dmu,n,T);
            k(isnan(k)) = 0;
        end
        
        function chi = chifn_n_dmu_T(obj,n,dmu,T)
            chi = interp3(obj.DelMuGrid2,obj.nsGrid,obj.TsGrid,obj.chi_n_dmu_T,dmu,n,T);
            chi(isnan(chi)) = 0;
        end
        
        
        %SINGLES PORTION DOESN'T WORK PRESENTLY!
        %My understanding is that it's related to the fact singles is not
        %single valued fn.
        %Functions (singles,dmu,T)
%         function n1 = n1fn_singles_dmu_T(obj,singles,dmu,T)
%             n1 = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.n1_singles_dmu_T,dmu,singles,T);
%             n1(isnan(n1)) = 0;
%         end
%         
%         function n2 = n2fn_singles_dmu_T(obj,singles,dmu,T)
%             n2 = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.n2_singles_dmu_T,dmu,singles,T);
%             n2(isnan(n2)) = 0;
%         end
%         
%         function n = nfn_singles_dmu_T(obj,singles,dmu,T)
%             npts = obj.n1_singles_dmu_T+obj.n2_singles_dmu_T;
%             n = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,npts,dmu,singles,T);
%             n(isnan(n)) = 0;
%         end
%         
%         function p = pfn_singles_dmu_T(obj,singles,dmu,T)
%             %p = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.p_singles_dmu_T,dmu,singles,T);
%             ppts = (obj.n1_singles_dmu_T-obj.n2_singles_dmu_T)./(obj.n1_singles_dmu_T+obj.n2_singles_dmu_T);
%             p = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,ppts,dmu,singles,T);
%             p(isnan(p)) = 0;
%         end
%         
%         function psingles = psinglesfn_singles_dmu_T(obj,singles,dmu,T)
%             n1_singles = obj.n1_singles_dmu_T - obj.n1_singles_dmu_T.*obj.n2_singles_dmu_T;
%             n2_singles = obj.n2_singles_dmu_T - obj.n1_singles_dmu_T.*obj.n2_singles_dmu_T;
%             ppts = (n1_singles-n2_singles)./(n1_singles+n2_singles);
%             psingles = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,ppts,dmu,singles,T);
%             psingles(isnan(psingles)) = 0;
%         end
%         
%         function c1 = c1fn_singles_dmu_T(obj,singles,dmu,T)
%             c1 = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.c1_singles_dmu_T,dmu,singles,T);
%             c1(isnan(c1)) = 0;
%         end
%         
%         function c2 = c2fn_singles_dmu_T(obj,singles,dmu,T)
%             c2 = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.c2_singles_dmu_T,dmu,singles,T);
%             c2(isnan(c2)) = 0;
%         end
%         
%         function muavg = mufn_singles_dmu_T(obj,singles,dmu,T)
%             muavg = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.mu_singles_dmu_T,dmu,singles,T);
%             muavg(isnan(muavg)) = 0;
%         end
%         
%         function k = kfn_singles_dmu_T(obj,singles,dmu,T)
%             k = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.k_singles_dmu_T,dmu,singles,T);
%             k(isnan(k)) = 0;
%         end
%         
%         function chi = chifn_singles_dmu_T(obj,singles,dmu,T)
%             chi = interp3(obj.DelMuGrid2,obj.SinglesGrid,obj.TsGrid,obj.chi_singles_dmu_T,dmu,singles,T);
%             chi(isnan(chi)) = 0;
%         end
        
        
    end
end





