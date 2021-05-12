function pixel = pixel2D(range)
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



% Resuse code form pixel 3D
pixel_3d=pixel3D(range,range);
pixel.kernel=pixel_3d.kernel.x;


end