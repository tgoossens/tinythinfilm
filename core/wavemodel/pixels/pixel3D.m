%function pixel = pixel3D(range_x,range_y)
function pixel = pixel3D(varargin)
%  PIXEL_FULLWIDTH  Returns the kernel for a pixel size equal to filter size as an anonymous function which can be evaluated at any spatial frequency.
%
%  Inputs: The pixel size is set with variable parameters.
%    One can either set size or ranges separately for x and y dimensions. The ranges are defined with a
%    coordinate system with the origin in the center of the pixel.
%   Examples
%   "width",2: - Make a square pixel of size 2 (in the units of choice)
%   "width x",2,"range y",[-1 0]  - Make a rectangular 
%     
%  Outputs:
%    pixel , 
%        kernel: anonymous function which can be evaluated as kernel(nu)
%                This kernel depends how different spatial frequencies are
%                weighted to effectively impose the finite pixel width. 
%    
% Copyright Thomas Goossens    

%% Variable arguments
 varargin = ieParamFormat(varargin);
 p = inputParser;
 p.addParameter('width', NaN, @isnumeric);
 p.addParameter('widthx', NaN, @isnumeric);
 p.addParameter('widthy', NaN, @isnumeric);
 p.addParameter('rangex', [], @isnumeric);
 p.addParameter('rangey',[],@isnumeric);
 
 p.parse(varargin{:});
width = p.Results.width;
width_x = p.Results.widthx;
width_y = p.Results.widthy;
range_x= p.Results.rangex;
range_y= p.Results.rangey;

% Pixel size can only be set in one way, no redundancies
conditions(1)=~sum(isnan(width));
conditions(end+1)=and(xor(~isnan(width_x),~isempty(range_x)),xor(~isnan(width_y),~isempty(range_y))); % You can define both a range and width for the same axis.


assert(sum(conditions)==1,"Incomplete or conflicting settings for pixel size.")

%% Bring everything back to ranges:

% Set of width: can be a single number or 2-element vector
if(not(isnan(width)))
    % If its a vector assign corresponding widths to range
   if(numel(width)==2)
       range_x=[-0.5 0.5]*width(1);
       range_y=[-0.5 0.5]*width(2); 
   else
    % If its a single element, both dimensions have the same width
       range_x=[-0.5 0.5]*width;
       range_y=range_x;
   end
end

% Set of width_x or width_y (both need to be single numbers)
if(conditions(2))
    if(~isnan(width_x))
        range_x=[-0.5 0.5]*width_x(1);
    else
        range_x=range_x(1:2);
    end
    
    if(~isnan(width_y))
        range_y=[-0.5 0.5]*width_y(1);
    else
        range_y=range_y(1:2);
    end
end



% Kernel for Y direction
pixel.kernel.x=chooseKernel(range_x);
pixel.kernel.y=chooseKernel(range_y);

% Information for later use
pixel.range.x=range_x;
pixel.range.y=range_y;




    function kernel=  chooseKernel(range)
        
        if(numel(range)==1)
            width=range(1);
            kernel = createCenteredKernel(width);
        elseif(range(1)==-range(2))
            width=abs(range(2)-range(1));
            kernel = createCenteredKernel(width);
        else
            kernel=customRangeKernel(range(1),range(2));
        end
    end

    function kernel = customRangeKernel(startpos,endpos)
        kernel = @K;
        function k = K(nu)
            
            k = 1i *(exp(2*pi*1i*nu*endpos)-exp(2*pi*1i*nu*startpos))./(2*pi*nu);
            
            % When zero, return the asymptotic limit
            k(nu==0) =  (endpos-startpos);
        end
    end


    function kernel = createCenteredKernel(width)
        kernel = @(nu) width*sinc_nopi(pi*width*nu);
    end

end