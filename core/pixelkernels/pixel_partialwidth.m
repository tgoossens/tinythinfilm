function  kernel = pixel_partialwidth(startpos,endpos)
%  PIXEL_PARTIALWIDTH  Returns the pixel kernel for a pixel smaller than the integrated filter.
%  
%  Inputs:    
%   - startpos: position where the pixel starts
%   - endpos: position where the pixel ends  
%    
%   The origin of the coordinate system is centered in the width of the filter.
%  
%  Outputs:
%    kernel: anonymous function which can be evaluated as kernel(nu)
%
%  Examples: 
%    If the filter is 4 micron width and the pixel is 2 micron wide and centered:
%    pixel_partialwidth(-1,1)
%   
%   When the pixel is aligned to the left of the filter:
%    pixel_partiawidth(-2,-1); 
%    
%    
%  See also: PIXEL_FULLWIDTH
%
%  Copyright Thomas Goossens
    sinca=@(x)sinc(x/pi);
    kernel = @(nu) 1i *(exp(-2*pi*1i*v*endpos)-exp(-2*pi*1i*v*startpos))./(2*pi*v+eps);
end