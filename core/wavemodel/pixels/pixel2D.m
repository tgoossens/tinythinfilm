function pixel = pixel2D(varargin)
%  PIXEL_FULLWIDTH  Returns the kernel for a pixel size equal to filter size as an anonymous function which can be evaluated at any spatial frequency.
%
%  Inputs:    
%   - width: width of the pixel (use units consistent with other distances)
%  
%  Outputs:
%    kernel: anonymous function which can be evaluated as kernel(nu)
%    
%  Notes:
%  The pixel kernel is not intended for independent usage. It is used to make TINYTRANSMITTANCE correctly integrate the incidence flux on the pixel area.
%
%  From an implementation perspectrive, this convolution approach reduces memory requirements and facilitates speed ups.
%    
% Copyright Thomas Goossens    

% Kernel for X direction
% If the boundaries are th esame, we can use an 

%% Variable arguments
 varargin = ieParamFormat(varargin);
 p = inputParser;
 p.addParameter('width', NaN, @isnumeric);
 p.addParameter('range',[],@isnumeric);
 p.parse(varargin{:});


width = p.Results.width;
range= p.Results.range;

% Pixel size can only be set in one way, no redundancies 
assert((~isnan(width)  + ~isempty(range))==1,"Insufficient or conflicting settings for pixel size.")


% Resuse code form pixel 3D to avoid redundancy of formulas
if(not(isnan(width)))
    pixel_3d=pixel3D('width',width);
elseif(not(isnan(range)))
    pixel_3d=pixel3D('range x',range,'range y',range);
end

% Select only along one dimension for the 2D pixel
pixel.kernel=pixel_3d.kernel.x;


end