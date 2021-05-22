function displaySetup3D(filter,varargin)
%DISPLAYSETUP Summary of this function goes here
%   Detailed explanation goes here

%% Variable arguments
 varargin = ieParamFormat(varargin);
 p = inputParser;
 p.addParameter('pixel', pixel3D('width',filter.width));
 p.addParameter('wavepacket', NaN)
 p.addParameter('wavelength', 0.4,@isnumeric)


 p.parse(varargin{:});
 pixel = p.Results.pixel;
 wavepacket = p.Results.wavepacket;



%%
hold on;

% Draw filter
filterwidth=filter.width;
filterthickness=filter.width/10;

%rectangle('Position',[leftbottomcorner filterwidth filterthickness],'FaceColor',[0.7 0.7 0.7])
plotcube([filterwidth filterwidth filterthickness],[-filterwidth/2 -filterwidth/2 0],.8,[1 0 0]);
% Draw pixel
pixelrange_x=pixel.range.x;
pixelrange_y=pixel.range.y;
pixelwidth_x=pixelrange_x(2)-pixelrange_x(1);
pixelwidth_y=pixelrange_y(2)-pixelrange_y(1);
pixelthickness=filterthickness;
plotcube([pixelwidth_x pixelwidth_y pixelthickness],[pixelrange_x(1) pixelrange_y(1) -pixelthickness],.8,[0.3 0.3 0.3]);


%% Plot wavepacket 3D
% %% Plot wavepacket 3D

wl=0.8;

%nu_x = @(th,phi) (1./l).*sin(th).*cos(phi);
%nu_y = @(th,phi) (1./l).*sin(th).*sin(phi);

nu_x = linspace(-1/wl,1/wl,100);
nu_y=nu_x';
F = abs(wavepacket(filterwidth,nu_x,nu_y,wl));

F=F./max(F(:)) .*filterwidth;
for x=1:numel(nu_x)
    for y=1:numel(nu_y)
        azimuth(x,y) = (atan2(nu_y(y),nu_x(x)));
        
        sigma(x,y)=sqrt((2*pi*nu_x(x)).^2 + (2*pi*nu_y(y)).^2);
        zenith(x,y) = real(asin(sigma(x,y)./(2*pi/wl)));
        
        fx(x,y)=-F(x,y).*sin(zenith(x,y)).*cos(azimuth(x,y));
        fy(x,y)=-F(x,y).*sin(zenith(x,y)).*sin(azimuth(x,y));
        fz(x,y)=filterthickness+F(x,y).*cos(zenith(x,y));
    end
end




s=surf(fx,fy,fz);
s.EdgeColor = 'none';
%s.FaceAlpha=0.5000;
%%
axis equal
box on
xlabel('X')
ylabel('Y')
zlabel('Z')


%% Draw normal vector

line([0 0],[0 0], [0 filterwidth],'color','k','linestyle',':','linewidth',2);


end

