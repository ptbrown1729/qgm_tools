function data = convert_qmc_file(filename)
    s = dir(filename);
    if s.bytes < 1000
        data = [];
        warning('file %s corrupted',filename);
        return
    end

    f1 = fopen(filename,'r');
    if (f1 == -1)
        error('failed to open %s', filename);
    end
    data = struct();
    % 2D Periodic Lattice; Nx= 10; Ny= 10; total sites=  100   
    line = fgetl(f1);
    if strfind(line,'General Geometry - Free Format')
        fclose(f1);
        data = convert_qmc_hopfile(filename);
%         data = convert_qmc_gradfile(filename);
        return
    end
    
    vals = strsplit(line,';');
    data.latt.type = strtrim(vals{1});
    data.latt.Nx = sscanf(strtrim(vals{2}),'Nx= %d');
    data.latt.Ny = sscanf(strtrim(vals{3}),'Ny= %d');
    data.latt.total_sites = sscanf(strtrim(vals{4}),'total sites= %d');
    
    while true
        line = fgetl(f1);
        if line(1) ==-1 || line(1) == '='
            break;
        end
        vals = strsplit(line,{':','+-'});
        vals{1} = strtrim(vals{1});
        vals{2} = strtrim(vals{2});
        switch vals{1}
            case {'U','t_up','t_dn','mu_up','mu_dn','dtau','beta','gamma','Number of sites','Number of measurement sweep','Frequency of measurement',...
                    'Frequency of recomputing G','Accept count','Reject count','Global move number of sites',...
                    'Global move accept count','Global move reject count','Global move accept rate','Number of warmup sweep'}
                data.(convert_name(strtrim(vals{1})))  = sscanf(vals{2},'%f');
            case 'Time slice - L'
                data.L = sscanf(vals{2},'%f');
            case 'Random seed'
                data.seed = sscanf(vals{2},'%f');
            case'Approximate accept rate'
                data.accept_rate = sscanf(vals{2},'%f');
            otherwise
                data.(convert_name(strtrim(vals{1})))  = vals{2};
        end
    end
    
    data.('num_greensfunction_errors') = NaN;
    data.('running_time_second') = NaN;
    
    while true
    
        line = fgetl(f1);
        if line(1) == -1
            fprintf('end of file reached in expected way\n');
            break;
        end
        section_name = convert_name(strtrim(line));
       % assert(all(strtrim(line) == 'Sign of equal time measurements:'));
        while true
            line = fgetl(f1);
            if line(1) == -1 || line(1) == '='
                break;
            end
            vals = strsplit(line,{':','+-'});
            field_name = convert_name(strtrim(vals{1}));
            try
                data.(section_name).(field_name) = [parseSingleNum(vals{2}) parseSingleNum(vals{3})];
            catch
                error('cannot parse line "%s" in file ""',line,filename);
            end
        end
        if strcmp(section_name,'pair_field_correlation_function___accumulated')
            line = fgetl(f1);
            if (line ~= -1)
                vals = strsplit(line,{':'});
                if strcmp('Number of Green''s function errors',strtrim(vals{1}))
                    data.('num_greensfunction_errors') = sscanf(vals{2},'%d');
                else
                    warning('unexpected line');    
                end
            end
            line = fgetl(f1);
            if (line ~= -1)
                vals = strsplit(line,{':','('});
                if strcmp('Running time',strtrim(vals{1}))
                    data.('running_time_second') = sscanf(vals{2},'%f');
                else
                    warning('unexpected line');    
                end
                
            end
            break
        end
        
    end
    %data.sign_of_equal_time_measurements
    
    
    
%     data
%     data.latt
%     for i=1:9
%         line = fgetl(f1);
%         fprintf(f2,[line ' 0\n']);
%     end
    fclose(f1);
    % Do some stuff
end

function num = parseSingleNum(s)
    num = sscanf(strtrim(s),'%f');
    if numel(num)==1
        return;
    elseif numel(num)==0 || numel(num)>2
        error('num expected');
    elseif numel(num)==2
        num = num(1)*10^num(2);
    end
end

function ret = convert_name(str)
    if length(str)>1
        ret = lower(str);
    else
        ret = str;
    end
    ret = strrep(ret,' ','_');
    ret = strrep(ret,':','');
    ret = strrep(ret,'<','');
    ret = strrep(ret,'>','');
    ret = strrep(ret,'*','_t_');
    ret = strrep(ret,'''','');
    ret = strrep(ret,';','');
    ret = strrep(ret,'=','');
    ret = strrep(ret,'(','');
    ret = strrep(ret,')','');
    ret = strrep(ret,'-','_');
end