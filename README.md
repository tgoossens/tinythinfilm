# Tiny Thin Film
Calculate the transmittance of pixel-integrated thin-film filters.

Thin-film filter deposition technology has advanced to the stage where filters can be made so tiny that they can be integrated onto separate pixels of image sensors.
This trend has mainly been driven by the need for developing compact and lightweight spectral cameras that aim to combine imaging and spectroscopy. 

The pixel-integrated thin-film filters are used to select specific wavelengths. The spatial width, however, strongly affects the angular sensitivity of thin-film filters.
This is problematic when the sensors are used with an imaging lens and needs to be taken into account at the design stage.

This toolbox aims to provide filter and camera designers with quick estimates of the expected filter transmittance.


<div>
<img src="./doc/img/pixelfilters.png" alt="Pixel integrated thin-film filters" width="45%" >
<img src="./doc/img/tinyfabry.png" alt="Tiny Fabry-Pérot transmittance" width="45%" >
</div>

<img src="./doc/img/polarized.png" alt="Differences in polarization" width="60%" >


# Features


- Efficient transmittance calculation compared to numerical solvers like FDFD or Finite Element
- Explore the effect of the width on the filter performance
- Explore the effect of the angle of incidence of collimated light
- Calculate transmittance for filters that are larger than the pixel
- Simulate for s and p polarization
- Simulate for focused light (under construction)

# How to use

The toolbox is intented to provide quick and good transmittance estimates, not exact results.

It is compatible with MATLAB and the free alternative Octave.

# Publications
Goossens, Thomas. "Tiny thin-film filters from a different angle: Correcting angular dependency for spectral imaging.", PhD Thesis (2020).  

# How to cite this work


# References
