classdef spectralfn < handle
    %spectralfn class for storing and manipulating spectral function data
    
    properties
        % parameters of spectral function
        kxs
        kys
        es
        ts
        ns
        mus
        us
        
        % spectral function and derived properties
        a
        aunc
        af
        afunc
        nk
        dos
        N
        norms
        
        linear_index
        %GXMG
        kxs_gxmg
        kys_gxmg
        es_gxmg
        a_gxmg
        aunc_gxmg
        af_gxmg
        afunc_gxmg
        nk_gxmg
        
        % GYMG
        kxs_gymg
        kys_gymg
        es_gymg
        a_gymg
        aunc_gymg
        af_gymg
        afunc_gymg
        nk_gymg
        
        % data representation information
        representation
        energy_mode
        
    end
    
    methods
        function obj = spectralfn(kxs, kys, es, ns, mus, ts, us, a, af, aunc, afunc, mode)
            % kxs: x momentum, length should be same as size(a, 2)
            %
            % kys: y momentum, length should be same as size(a, 1)
            %
            % es: either a 1D vector the same size as size(a, 3), or a
            % matrix the size of the first three dimensions of a.
            %
            % ts: temperatures, size(a, 4)
            % 
            % ns: If data is at constant ns, then a 1d vector. If data is
            % at constants mus, then will be size nus, nts, nus
            %
            % mus: If data is at constant mus, then a 1d vector. If data is
            % at constant ns, then mus will be size nns x nts x nus
            %
            % us: interactions, size(a, 5)
            %
            % a: array of shape nkys x nkxs x nes x ns x nts x nus. Note
            % that at most one of a and af may be empty.
            %
            % af: Note that at most one of a and af may be empty.
            %
            % mode: 'quad' or 'full'. If 'quad', assumes that the
            % southeast quadrant is provided (i.e. kx>0, ky>0)
            if ~exist('kxs', 'var')
                return;
            end
            
            
            
            obj.a = a;
            if ~isempty(aunc)
                obj.aunc = aunc;
            end
            
            obj.af = af; 
            if ~isempty(afunc)
                obj.afunc = afunc;
            end
            
            if isempty(obj.a) && isempty(obj.af)
                error('at most one of a and af may be empty');
            end
            
            if isempty(ts)
                ts = nan;
            end
            
            if isempty(ns)
                ns = nan;
            end
            
            if isempty(mus)
                mus = nan;
            end
            
            if isempty(us)
                us = nan;
            end
            
            obj.kxs = kxs;
            obj.kys = kys;
            obj.es = es;
            obj.ts = ts;
            obj.ns = ns;
            obj.mus = mus;
            obj.us = us;
            
            % check that all sizes are compatible
            if ~isempty(obj.a)
                size_fn = size(obj.a);
            elseif ~isempty(obj.af)
                size_fn = size(obj.af);
            end
            
            % check kx/ky size
            if length(obj.kxs) ~= size_fn(2)
                error('kxs size was not consistent with spectral function');
            end
            
            if length(obj.kys) ~= size_fn(1)
                error('ky size was not consistent with spectral function');
            end
            
            
            % check energy size
            if ndims(obj.es) > 2
                obj.energy_mode = 'unequal';
                
                if ~isequal(size_fn(1:3), size(obj.es))
                    error('energy size was not consistent with spectral function');
                end
                
            else
                obj.energy_mode = 'equal';
                
                if ~isequal(size_fn(3), length(obj.es))
                    specfn_uneq_es
                end
      
            end
             
            % set symmetry mode
            if strcmp(mode, 'quad') || strcmp(mode, 'full')
               obj.representation = mode;
            else
                error();
            end
            
            obj.get_extended_quantities();
        end
        
        function get_extended_quantities(obj)
            % Extract nk, dos, and gxmg quantities from spectral function
           
            
            if ~isempty(obj.a)
                int = obj.a;
                int( isnan(int) ) = 0;
                obj.dos = squeeze( trapz( obj.kys, trapz(obj.kxs, int, 2), 1) );
                
                % deal with cases of equal or unequal energy
                if strcmp( obj.energy_mode, 'equal')
                    obj.norms = squeeze( trapz( obj.es, int, 3) );
                    
                elseif strcmp( obj.energy_mode, 'unequal')
                    sizea = size(obj.a);
                    obj.norms = zeros( [sizea(1:2), sizea(4:end)]);
                    for ii = 1 : length(obj.kxs)
                        for jj = 1 : length(obj.kys)
                            obj.norms(ii, jj, :) = trapz( squeeze(obj.es(ii, jj, :)), int(ii, jj, :, :), 3);
                        end
                    end
                else
                    error();
                end
            end
            
            if ~isempty(obj.af)
                int = obj.af;
                int( isnan(int) ) = 0;
                
                % unequal vs equal energy
                if strcmp( obj.energy_mode, 'equal')
                    obj.nk = squeeze( trapz( obj.es, int, 3) ); 
                    
                elseif strcmp( obj.energy_mode, 'unequal')
                    sizea = size(obj.af);
                    obj.nk = zeros( [sizea(1:2), sizea(4:end)]);
                    for ii = 1 : length(obj.kys)
                        for jj = 1 : length(obj.kxs)
                            obj.nk(ii, jj, :) = trapz( squeeze(obj.es(ii, jj, :)), int(ii, jj, :, :), 3);
                        end
                    end
                else
                    error();
                end
         
                obj.N = 2 / (2*pi)^2 * squeeze( trapz( obj.kys, trapz(obj.kxs, obj.nk, 2), 1 ) );
            end
            
            % deal with cases of kxs specified as 2D vs. 1D
            [kxkx, kyky] = meshgrid(obj.kxs, obj.kys);

            if strcmp(obj.representation, 'quad')
                obj.dos = 4 * obj.dos;
                obj.N = 4 * obj.N;
                a_int = obj.a;
                aunc_int = obj.aunc;
                af_int = obj.af;
                afunc_int = obj.afunc;
                nk_int = obj.nk; 
                e_int = obj.es;
                
            elseif strcmp(obj.representation, 'full')
                if mod( length(obj.kxs), 2) == 1
                    include_center = 1;
                else
                    include_center = 0;
                end
                
                kxkx = get_quadrant(kxkx, 'southeast', 'include_center');
                kyky = get_quadrant(kyky, 'southeast', 'include_center');
                
                if strcmp(obj.energy_mode, 'unequal')
                    e_int = get_quadrant(obj.es, 'southeast', 'include_center');
                end
                
                if ~isempty(obj.a)
                    a_int = get_quadrant(obj.a, 'southeast', 'include_center');
                    if ~isempty(obj.aunc)
                        aunc_int = get_quadrant(obj.aunc, 'southeast', 'include_center');
                    end
                end
               
                if ~isempty(obj.af)
                    af_int = get_quadrant(obj.af, 'southeast', 'include_center');
                    nk_int = get_quadrant(obj.nk, 'southeast', 'include_center');
                    if ~isempty(obj.afunc)
                        afunc_int = get_quadrant(obj.afunc, 'southeast', 'include_center');
                    end
                end
                        
            else
                error();
                
            end
            
            % do symmetry cuts
            [obj.kxs_gxmg, ~, obj.kxs_gymg, ~, ~] = get_high_symm_cut(kxkx, 'upperleft');
            [obj.kys_gxmg, ~, obj.kys_gymg, ~, ~] = get_high_symm_cut(kyky, 'upperleft');
            
            % 
            if strcmp(obj.energy_mode, 'unequal')
                [obj.es_gxmg, ~, obj.es_gymg, ~, ~] = get_high_symm_cut(e_int, 'upperleft');
            else
                obj.es_gxmg = obj.es;
                obj.es_gymg = obj.es;
            end
            
            if ~isempty(obj.a)
                [obj.a_gxmg, ~, obj.a_gymg, ~, obj.linear_index] = get_high_symm_cut(a_int, 'upperleft');
                if ~isempty(obj.aunc)
                    [obj.aunc_gxmg, ~, obj.aunc_gymg, ~, ~] = get_high_symm_cut(aunc_int, 'upperleft');
                end
            end

            if ~isempty(obj.af)
                [obj.af_gxmg, ~, obj.af_gymg, ~, obj.linear_index] = get_high_symm_cut(af_int, 'upperleft');
                [obj.nk_gxmg, ~, obj.nk_gymg, ~, ~] = get_high_symm_cut(nk_int, 'upperleft');
                if ~isempty(obj.afunc)
                    [obj.afunc_gxmg, ~, obj.afunc_gymg, ~, ~] = get_high_symm_cut(afunc_int, 'upperleft');
                end
            end
                
        end       
        
        function quad2full(obj)
            % convert a quadrant to the full brillouin zone
            
            if ~strcmp(obj.representation, 'quad')
                return;
            end
            
            if mod( length(obj.kxs), 2) == 1
                include_center = 1;
            else
                include_center = 0;
            end
            
            if ~isempty(obj.a)
                [obj.a, kxs, kys] = quad2full(obj.a, 'southeast', include_center, obj.kxs, obj.kys);
                obj.norms = quad2full(obj.norms, 'southeast', include_center);
                
                if ~isempty(obj.afunc)
                    obj.afunc = quad2full(obj.afunc, 'southeast', include_center);
                end
            end
            
            if ~isempty(obj.af)
                [obj.af, kxs, kys] = quad2full(obj.af, 'southeast', include_center, obj.kxs, obj.kys);
                obj.nk = quad2full(obj.nk, 'southeast', include_center);
                
                if ~isempty(obj.aunc)
                    obj.aunc = quad2full(obj.aunc, 'southeast', include_center);
                end
            end
            
            obj.kxs = kxs;
            obj.kys = kys;
            obj.representation = 'full';
        end
        
        function full2quad(obj)
            % convert full brillouin zone to a single quadrant
            
            if ~strcmp(obj.representation, 'full')
                return;
            end
            
            if mod( length(obj.kxs), 2) == 1
                include_center = 1;
            else
                include_center = 0;
            end
            
            if ~isempty(obj.a)
                [obj.a, kxs, kys] = get_quadrant(obj.a, 'southeast', 'include_center', obj.kxs, obj.kys);
                obj.norms = get_quadrant(obj.norms, 'southeast', 'include_center');
                
                if ~isempty(obj.aunc)
                    obj.aunc = get_quadrant(obj.aunc, 'southeast', 'include_center');
                end
            end
            
            if ~isempty(obj.af)
                [obj.af, kxs, kys] = get_quadrant(obj.af, 'southeast', 'include_center', obj.kxs, obj.kys);
                obj.nk = get_quadrant(obj.nk, 'southeast', 'include_center');
                
                if ~isempty(obj.afunc)
                    obj.afunc = get_quadrant(obj.afunc, 'southeast', 'include_center');
                end   
            end
            
            kxs = obj.kxs;
            kys = obj.kys;     
            obj.representation = 'quad';
        end
        
        function [as_struct] = save_struct(obj, fname)
            % save class instance as a structure
            
            if ~exist('compress', 'var') || isempty(compress)
                compress = 0;
            end
            
            as_struct = struct();
            
            fields = fieldnames(obj);
            for ii = 1:length(fields)
                as_struct.(fields{ii}) = obj.(fields{ii});
            end

            
            save(fname, '-struct', 'as_struct');
        end
        
        function load_struct(obj, struct_to_load)
            % load structure to class instance
            
            if ischar(struct_to_load)
                struct_to_load = load(struct_to_load);
            end
            
            % now struct_to_load should be our desired structure
            if ~isstruct(struct_to_load)
                error('struct_to_load was not a structure');
            end
            
            obj_fields = fieldnames(obj);
            for ii = 1:length(obj_fields)
                obj.(obj_fields{ii}) = [];
            end
            
            % load struct and get field names
            struct_fields = fieldnames(struct_to_load);
            % loop over struct field names and assign any shared to our class
            for ii = 1 : length(struct_fields)
                try   
                    obj.(struct_fields{ii}) = struct_to_load.(struct_fields{ii});
                catch
                    fprintf('field %s was present in loaded struct, but it is not a field of arpes_data class. Skipped this field. \n',struct_fields{ii});
                end
            end
            
        end
        
    end
    
end

