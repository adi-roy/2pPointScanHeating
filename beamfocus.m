classdef beamfocus
    properties
        NumPhotons
        R
        ZCenter
        FocalLength
        Locations
        Trajectories
    end
    
    methods
        function obj = beamfocus(numPhotons, r, zCenter, fl)
            obj.NumPhotons = numPhotons;
            obj.R = r;
            obj.ZCenter = zCenter;
            obj.FocalLength = fl;
            
            obj.initializePhotons();
        end
        
        function obj = initializePhotons(obj)
%             rinit = sqrt(-0.5*log(rand(1,3*obj.NumPhotons)))*obj.R; %radial location for 2D Gaussian with 1/e^2 radius of r. Generate redundant so that there is still enough after clipping.
%             rinit = rinit(rinit <= constants.f*constants.NA); % clipping by the back aperature
            
            d = 2;
            bx = obj.R;
            by = obj.R;
            sig_x = bx;
            sig_y = by;
            m = [0 0]';
            K0 = [sig_x^2 0;0 sig_y^2];
            eps = 0.001;
            I = [1 0;0 1];
            K = K0+eps*I;
            L = chol(K);
            n = 3*obj.NumPhotons;
            u0 = normrnd(0,1,[d,n]);
            u = reshape(u0,d,n);
            x = m + L*u;
            zeta_x = randsample(x(1,:),n);
            zeta_y = randsample(x(2,:),n);

            ridx = abs(zeta_x) <= bx;
            ridy = abs(zeta_y) <= by;
            zeta_x = zeta_x(ridx);
            zeta_y = zeta_y(ridy);

            zeta_x = zeta_x(1:obj.NumPhotons); zeta_y = zeta_y(1:obj.NumPhotons);
            rinit = sqrt(zeta_x.^2+zeta_y.^2); 
            
            rinit = rinit(1:obj.NumPhotons); 
            thinit = asind(rinit/constants.f/constants.nWater); %sine condition for inifinitely conjugated system, used to calculate theta angle for each photon.
            thinitpos = -180+2*rand(1,obj.NumPhotons)*180; %initial position (angular coordinate)

            %directions are encoded in spherical coordinates:
            % thinit = cosd(th_out)+rand(1,num)*(1-cosd(th_out)); %length of the Z-vector
            % phiinit=-180+2*rand(1,num)*180; %inital angle in the z plane
            %initialize x,y, and z direction matrices such that [x y z] is a unit
            %vector
            zdir=cosd(thinit);
            temp = sqrt(1-zdir.^2);
            xdir=-temp.*cosd(thinitpos);
            ydir=-temp.*sind(thinitpos);

            %simulating scanning by translating the fiber in xy @ zcenter, NOT angular
            %scanning. Give each photon a random translation in xy.
            FOV = constants.FOV;
            dx = (FOV*rand(1,obj.NumPhotons)-0.5*FOV);
            dy = (FOV*rand(1,obj.NumPhotons)-0.5*FOV);
            % dx = 0;
            % dy = 0;

            %initializing x, y, and z coordinates at zcenter
            xcor=obj.FocalLength*tand(thinit).*cosd(thinitpos);
            ycor=obj.FocalLength*tand(thinit).*sind(thinitpos);
            zcor=ones(1,obj.NumPhotons)*obj.ZCenter;
            obj.Locations=[xcor+dx;ycor+dy;zcor];
            obj.Trajectories=[xdir;ydir;zdir];
        end
    end
end
