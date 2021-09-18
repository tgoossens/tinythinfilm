ffmpeg -i animation-0deg.avi -gifflags +transdiff -y  -loop 0 -filter_complex  "[0]reverse[r];[0][r]concat=n=2:v=1:a=0,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse"  widthsweep-0deg.gif

ffmpeg -i animation-15deg.avi -gifflags +transdiff -y  -loop 0 -filter_complex  "[0]reverse[r];[0][r]concat=n=2:v=1:a=0,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse"  widthsweep-15deg.gif

