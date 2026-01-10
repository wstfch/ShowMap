# fits-show-map

A small WCS-aware FITS image viewer / helper toolkit.

## Install (pip)

```bash
pip install -e .
```

## Usage

```bash
Examp:
from fits_show_map import ShowMap
path = '/Users/wst/galaxies/NGC2442/spx_index/'
filename = 'NGC2442_EMU_SB59742_I'
inf = path + filename + '.fits'
header, data = ShowMap.load_fits_image(inf)
ShowMap.show_fits(header=header,data=data,colobar=True)
```
