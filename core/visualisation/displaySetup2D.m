function [outputArg1,outputArg2] = displaySetup2D(filter,varargin)
%DISPLAYSETUP Summary of this function goes here
%   Detailed explanation goes here

%% Variable arguments
 varargin = ieParamFormat(varargin);
 p = inputParser;
 p.addParameter('pixel', pixel2D('width',filter.width));
 p.addParameter('wavepacket', NaN)


 p.parse(varargin{:});
 pixel = p.Results.pixel;
 wavepacket = p.Results.wavepacket;



%%
hold on;

% Draw filter
filterwidth=filter.width;
filterthickness=filter.width/10;

%rectangle('Position',[leftbottomcorner filterwidth filterthickness],'FaceColor',[0.7 0.7 0.7])
plotcube([filterwidth 0 filterthickness],[-filterwidth/2 0 0],.8,[1 0 0]);
% Draw pixel
pixelrange_x=pixel.range.x;
pixelwidth_x=pixelrange_x(2)-pixelrange_x(1);
pixelthickness=filterthickness;
plotcube([pixelwidth_x 0 pixelthickness],[pixelrange_x(1) 0 -pixelthickness],.8,[0.3 0.3 0.3]);




%% Plot wavepacket 2D
wl=0.4;
f2d = @(phi_x) abs(wavepacket(filterwidth,sin(phi_x)/wl,wl)) ;

azimuth=linspace(-pi/2,pi/2,500);

fmax=1./max(abs(f2d(azimuth)))*filterwidth;

z=filterthickness+fmax*f2d(azimuth).*cos(azimuth);
x=-fmax*f2d(azimuth).*sin(azimuth);
plot3(x,zeros(size(x)),z);





%% Draw normal vector

line([0 0],[0 0], [0 filterwidth],'color','k','linestyle',':','linewidth',2);




%% 2D VIEW
view(0,0)
end

