function pixel = pixel2D(varargin)
%  PIXEL_FULLWIDTH  Returns the kernel for a pixel size equal to filter size as an anonymous function which can be evaluated at any spatial frequency.
%
%  Inputs:    
%   - width: width of the pixel (use units consistent with other distances)
%  
%  Outputs:
%    kernel: anonymous function which can be evaluated as kernel(nu)
%            This kernel depends how different spatial frequencies are
%            weighted to effectively impose the finite pixel width. 
%  Notes:
%    
% Copyright Thomas Goossens    


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
pixel.kernel.x=pixel_3d.kernel.x;
pixel.range.x=pixel_3d.range.x;


end