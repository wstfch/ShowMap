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
fig, ax = ShowMap.show_fits(header=header,data=data,colobar=True,fontsize=22,beam=True,cmap='jet',cb_dedi='%0.4f')

### For return the "ax", you can link to matplotlib, for example follow above:
ax.set_title('NGC2442')
ax.set_ylabel('DEC (J2000)',labelpad=1)
ax.set_xlabel('RA (J2000)',labelpad=1)

However, all the functions of image detail adjustment can be fully implemented in "show_fits" without the need for additional links. For example:

fig,ax=ShowMap.show_fits(header=header,data=data,colobar=True,fontsize=23,beam=True,cmap='jet',cb_dedi='%0.4f',\
                        RA='RA (J2000)', DEC='DEC (J2000)', title='NGC2442', beam_color='r',\
                        beam_p_pix=[20,20], beam_fluc=[6,6], savefig='./NGC2442.pdf', dpi=100,\
                        savefits='./NGC2442_new.fits')
```
