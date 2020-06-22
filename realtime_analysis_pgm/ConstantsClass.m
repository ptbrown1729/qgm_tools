classdef ConstantsClass
    % class for storing useful constants
    
    properties
        units = 'si';
        %SI constants. physics.nist.gov/cuu/Constants/index.html
        K = 1.3806488e-23;
        h = 6.62606957e-34; %Js
        hbar = 6.62606957e-34/(2*pi); %Js
        c = 2.99792458e8; %m/s
        e = 1.602176565e-19; %C
        mu0 = 4*pi*1e-7;
        %epsilon0 = 1/(c^2*mu0);
        epsilon0 = 1/((2.99792458e8)^2*(4*pi*1e-7));
        
        me = 9.10938291e-31; %kg
        mp = 1.672621777e-27; %kg
        alpha = 7.2973525698e-3; %dimensionless
        %abhor = hbar/(me*c*alpha);
        abohr = (6.62606957e-34 / (2*pi)) / (9.10938291e-31 * 2.99792458e8 * 7.2973525698e-3); %m
        
        % atomic data
        mass % kg
        
        % D2 line
        lambda_D2 % nm
        omega_D2 % rad/s
        gamma_D2 % Hz
        Isat_D2 % W/m^2
        
        % D1 line
        lambda_D1 % nm
        omega_D1 % rad/s
        gamma_D1 % Hz
        Isat_D1 % W/m^2
    
    end
    
    methods
        function obj = ConstantsClass(atom)
            if ~exist('atom', 'var') || isempty(atom)
                atom = 'lithium-6';
            end
            
            % initialize different atomic species
            if strcmp(atom, 'lithium-6')
                obj.mass = 9.98834e-27; %Kg;
                obj.lambda_D2 = 670.977e-9; % meters
                obj.gamma_D2 = 2 * pi * 5.8724e6; % rad/s
                obj.Isat_D2 = 25.4; % W/m^2
                
                obj.lambda_D1 = 670.992421e-9; % meters
                obj.gamma_D1 = 2 * pi * 5.8724e6; % rad/s
                obj.Isat_D1 = 75.9; % W/m^2
                
            elseif strcmp(atom, 'sodium-23')
                obj.mass = 0.381754035e-25; % kg
                obj.lambda_D2 = 589.1583264e-9 ; % meters
                obj.gamma_D2 = 2 * pi * 9.7946e6; % rad/s
                obj.Isat_D2 = inf; % W/m^2
                
                obj.lambda_D1 = 589.7558147e-9; % meters
                obj.gamma_D1 = 2 * pi * 9.765e6; % rad/s
                obj.Isat_D1 = inf; % W/m^2
                
            elseif strcmp(atom, 'rubidium-87')
                obj.mass = 1.443160648e-25; % kg
                obj.lambda_D2 = 780.241209686e-9 ; % meters
                obj.gamma_D2 = 2 * pi * 6.0666e6; % rad/s
                obj.Isat_D2 = 16.7; % W/m^2
                
                obj.lambda_D1 = 794.978851156e-9; % meters
                obj.gamma_D1 = 2 * pi * 5.7500e6; % rad/s
                obj.Isat_D1 = inf; % W/m^2
                
            else
                error();
            end
            
            obj.omega_D2 = 2 * pi * obj.c / obj.lambda_D2;
            obj.omega_D1 = 2 * pi * obj.c / obj.lambda_D1;
                
        end
        
    end
    
end

