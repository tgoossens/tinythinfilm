function pixel = pixel3D(range_x,range_y)
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
% varargin = ieParamFormat(varargin);
% p = inputParser;
% p.addParameter('width', NaN, @isnumeric);
% p.addParameter('width_x', NaN, @isnumeric);
% p.addParameter('width_y', NaN, @isnumeric);
% p.addParameter('range_x', NaN, @isnumeric);
% p.addParameter('range_y',NaN,@isnumeric);
% p.parse(varargin{:});





% Kernel for Y direction
pixel.kernel.x=chooseKernel(range_x);
pixel.kernel.y=chooseKernel(range_y);


    function kernel = createCustomKernel(startpos,endpos)
        kernel = @K;
        function k = K(nu)
            
            k = 1i *(exp(-2*pi*1i*nu*endpos)-exp(-2*pi*1i*nu*startpos))./(2*pi*nu);
            
            % When zero, return the asymptotic limit
            k(nu==0) =  (endpos-startpos);
        end
    end


    function kernel = createCenteredKernel(width)
        kernel = @(nu) width*sinc_nopi(pi*width*nu);
    end

    function kernel=  chooseKernel(range)
        
        if(numel(range_x)==1)
            width=range_x(1);
            kernel = createCenteredKernel(width);
        elseif(range_x(1)==-range_x(2))
            width=abs(range_x(2)-range_x(1));
            kernel = createCenteredKernel(width);
        else
            kernel=createCustomKernel(range_x(1),range_x(2));
        end
    end
end