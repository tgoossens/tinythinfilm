function [outputArg1,outputArg2] = displaySetup(filter,varargin)
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
figure; hold on;

% Draw filter
filterwidth=filter.width;
filterthickness=filter.width/10;;
leftbottomcorner=[-filterwidth/2 0]
rectangle('Position',[leftbottomcorner filterwidth filterthickness],'FaceColor',[0.7 0.7 0.7])

% Draw pixel
pixelrange=pixel.range.x;
pixelwidth=pixelrange(2)-pixelrange(1);
leftbottomcorner=[pixelrange(1) -filterwidth/10]
rectangle('Position',[leftbottomcorner pixelwidth filterthickness],'FaceColor',[0.1 0.1 0.1])


%% Plot wavepacket

f = @(phi) abs(wavepacket(filterwidth,sin(phi)/0.4,0.4)) ;
phi=-pi/2:0.01:pi/2;
fmax=max(abs(f(phi)));

y=filterthickness+f(phi).*cos(phi);
x=-f(phi).*sin(phi);
plot(x,y);

end

