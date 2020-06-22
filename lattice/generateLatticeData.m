% Solve lattice band structure at a variety of different depths and
% retrorefelection attenuation factors and store data in text file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
attenuation_factors = [0.025,0.1,0.2];
depths = 2:2:70;
z_trp_frq = 21e3; %green trap frq
scattering_length = -888;
beam_angles = 91.6267; %2*atan(ax / ay) where ax, ay are the lattice spacings

fname_out = 'LatticeBandData_Interpolation5.txt';
ii = 1;
while exist(fname_out, 'file')
    [folders, name, ext] = fileparts(fname_out);
    fname_out = fullfile(folders, sprintf('%s_%d%s', name, ii, ext));
    ii = ii + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : length(depths)
    for jj = 1 : length(attenuation_factors)
        [tTb, tx, ty, tdiag, u,...
            BandWidths, GapFromGroundBand, BandTopToGroundBandBottom,...
            MeanSplittingToGroundBand, EBand,...
            WannierFn, BeamAngle] = ...
            Lattice_2D_Asymmetric(depths(ii), scattering_length, z_trp_frq, attenuation_factors(jj), beam_angles, 1, 0);
        
        Row = transpose(cat(1, depths(ii), attenuation_factors(jj), ...
            beam_angles, scattering_length, z_trp_frq, tTb, tx, ty, tdiag, u,...
            MeanSplittingToGroundBand, GapFromGroundBand,...
            BandTopToGroundBandBottom, BandWidths));
        
        if ii == 1 && jj ==1
            Titles = cat(1,{},'Depth(Er)','Attenuation(Efield)',...
                'BeamAngle(deg)', 'as(abohr)', 'zConf(Hz)',...
                'tTB(Er)', 'tx(Er)', 'ty(Er)', 'tdiag(Er)', 'u(Er)');
            NBands = size(EBand,1);
            NameSplit = {}; 
            NameTp2Btm = {}; 
            NameBtm2Top = {}; 
            NameBw = {};
            for kk = 1:NBands
                NameSplit = cat(1, NameSplit, sprintf('MeanSplittingGndTo%d', kk));
                NameTp2Btm = cat(1, NameTp2Btm, sprintf('GndTopToBandBtm%d', kk));
                NameBtm2Top = cat(1, NameBtm2Top, sprintf('GndBtmToBandTop%d', kk));
                NameBw = cat(1, NameBw, sprintf('BandWidth%d', kk));
            end
            
            Titles = cat(1, Titles, NameSplit, NameTp2Btm, NameBtm2Top, NameBw);
            %write to file titles
            Delimiter = ',';
            Fid = fopen(fname_out, 'w');
            for ll = 1:length(Titles)-1
                fprintf(Fid, sprintf('%s%s', Titles{ll}, Delimiter)); 
            end
            fprintf(Fid, sprintf('%s\n', Titles{end}));
            fclose(Fid);
        end
        dlmwrite(fname_out,Row,'-append','delimiter',',');
    end
end

%%
%combine two files

Path1 = 'Latt_LowVals.txt';
Path2 = 'LatticeBandData_Interpolation_add_4_and_5_to_this.txt';
Output = 'combined.txt';

a = dlmread(Path1,',',1,0);
b = dlmread(Path2,',',1,0);

%in this case, wanted only odd depths...
% depths_a = a(:,1);
% Iodd = find(mod(depths_a,2));
% areduced = a(Iodd,:);
% full = cat(1,areduced,b);
% full_sorted = sortrows(full,[1,2]);

%generally want...
full = cat(1,a,b);
full_sorted = sortrows(full,[1,2]);

dlmwrite(Output,full_sorted,'delimiter',',');