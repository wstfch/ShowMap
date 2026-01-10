# fits-show-map

A small WCS-aware FITS image viewer / helper toolkit.

## Install (pip)

```bash
git clone https://github.com/wstfch/ShowMap.git
cd ShowMap
pip install -e .
Check:
pip show fits-show-map
```
## Uninstall
```
pip uninstall fits-show-map
```
## Usage

```bash
Examp:
1. Quickly in your terminal show a fits image
fits-show-map your.fits

2. In jupyter
from fits_show_map import ShowMap
inf = '/Users/wst/galaxies/NGC2442/NGC2442_EMU_SB59742_I.fits'
header, data = ShowMap.load_fits_image(inf)
ShowMap.show_fits(header=header,data=data,colobar=True,fontsize=22,beam=True,cmap='jet',cb_dedi='%0.4f')
```
