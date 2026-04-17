# fits-show-map

A small WCS-aware FITS image viewer and helper toolkit for astronomy.

`fits-show-map` is a lightweight Python package for displaying and manipulating FITS images with WCS information. It supports both quick command-line visualization and flexible Python workflows, and is designed to be simple, practical, and convenient for everyday astronomical image analysis and figure production. The package can be used to load FITS images into a clean 2D format, display images with WCS axes, reproject and rotate images, generate cutouts, make single-panel and multi-panel figures, overlay contours, draw beam ellipses, add colorbars, customize ticks, draw scalebars and arrows, and apply a CARTA-like logarithmic stretch for display.

## To install the package:

```bash
git clone https://github.com/wstfch/ShowMap.git
cd ShowMap
pip install -e .
```

## To check whether the installation was successful:

```bash
pip show fits-show-map
```

## To uninstall the package:

```bash
pip uninstall fits-show-map
```

## Functions
### 1. You can quickly show a FITS image in the terminal with:

```bash
fits-show-map your.fits --cmap inferno --savefig
```

The command-line interface takes the FITS file path as the main input, `--cmap` for the matplotlib colormap name, and `--savefig` to save a PNG image next to the FITS file.

### 2. You can also use the package in Python:

```python
from fits_show_map import ShowMap
```
#### 2.1 load_fits_image

`ShowMap.load_fits_image(filename)` reads a FITS file, reduces it to a 2D image, and cleans the header so the WCS remains usable for plotting. It supports FITS arrays with 2 to 5 axes by extracting the leading 2D image plane, removes higher-dimensional WCS keywords, records the x/y coordinate indices in the header, converts AIPS beam keywords to standard `BMAJ`, `BMIN`, and `BPA` when available, and can derive `CDELT` from `PIXSCAL`.

Example:

```python
from fits_show_map import ShowMap

inf = "/Users/wst/galaxies/NGC2442/NGC2442_EMU_SB59742_I.fits"
header, data = ShowMap.load_fits_image(inf)
```
#### 2.2 project
`ShowMap.project(...)` can reproject, recenter, rotate, and regrid an image to a new WCS. This is useful when you want to compare maps with different image sizes, pixel scales, centers, or orientations.

Example:

```python
from fits_show_map import ShowMap

header, data = ShowMap.load_fits_image("your.fits")

new_header, new_data = ShowMap.project(
    header,
    data,
    RA=114.099166,
    DEC=-69.530833,
    NAXIS1=1200,
    NAXIS2=1200,
    CRPIX1=600,
    CRPIX2=600,
    rotate=30,
    fill_nan=True,
    fill_size=25
)
```
#### 2.3 Show_fits
`ShowMap.show_fits(...)` is the main plotting function for a single image. It supports WCS-aware plotting, beam drawing, colorbars, contour overlays, custom labels, titles, saving figures, saving output FITS files, and logarithmic-style display. In many cases, most figure adjustments can be completed directly inside `show_fits(...)` without additional matplotlib commands.

Example:

```python
from fits_show_map import ShowMap

inf = "/Users/wst/galaxies/NGC2442/NGC2442_EMU_SB59742_I.fits"
header, data = ShowMap.load_fits_image(inf)

fig, ax = ShowMap.show_fits(
    header=header,
    data=data,
    colobar=True,
    fontsize=22,
    beam=True,
    cmap="jet",
    cb_dedi="%0.4f"
)
```

Because `show_fits(...)` returns `fig` and `ax`, you can continue using normal matplotlib commands if needed. For example:

```python
ax.set_title("NGC2442")
ax.set_ylabel("DEC (J2000)", labelpad=1)
ax.set_xlabel("RA (J2000)", labelpad=1)
```

However, the same settings can usually be done directly inside `show_fits(...)`:

```python
fig, ax = ShowMap.show_fits(
    header=header,
    data=data,
    colobar=True,
    fontsize=23,
    beam=True,
    cmap="jet",
    cb_dedi="%0.4f",
    xlabel="RA (J2000)",
    ylabel="DEC (J2000)",
    title="NGC2442",
    beam_color="r",
    beam_p_pix=[20, 20],
    beam_fluc=[6, 6],
    savefig="./NGC2442.pdf",
    dpi=100,
    savefits="./NGC2442_new.fits"
)
```

If you want to use a logarithmic-style display stretch, you can set `log=True`:

```python
fig, ax = ShowMap.show_fits(
    header=header,
    data=data,
    log=True,
    alpha=100,
    colobar=True,
    cb_lab="mJy/beam",
    cmap="inferno"
)
```

You can overlay contours either through `show_fits(...)` directly or by calling `ShowMap.show_contour(...)` yourself.

Example using `show_fits(...)`:

```python
fig, ax = ShowMap.show_fits(
    header=header,
    data=data,
    cont=True,
    cont_data=data,
    cont_levels=[0.001, 0.002, 0.004, 0.008],
    cont_color="cyan"
)
```

Example using `ShowMap.show_contour(...)`:

```python
ShowMap.show_contour(
    ax,
    cont_data=data,
    cont_levels=[0.002, 0.004, 0.008, 0.016],
    cont_color="white",
    cont_alpha=0.7
)
```

`ShowMap.show_beam(...)` can be used to draw the beam ellipse from FITS beam keywords such as `BMAJ`, `BMIN`, and `BPA`.

Example:

```python
fig, ax = ShowMap.show_fits(header=header, data=data, beam=False)

ShowMap.show_beam(
    ax,
    header,
    beam_p_pix=[25, 25],
    beam_fluc=[4, 4],
    beam_color="yellow",
    beam_sque=True
)
```

`ShowMap.show_colobar(...)` adds a colorbar to an existing image and supports custom position, label, font size, and tick formatting.

Example:

```python
fig, ax = ShowMap.show_fits(header=header, data=data, colobar=False)
im = ax.images[-1]

ShowMap.show_colobar(
    im,
    ax,
    data,
    cb_loct="top",
    cb_lab="mJy/beam",
    cb_dedi="%0.3f",
    cb_font=14
)
```

`ShowMap.show_tick(...)` customizes the tick style and spacing for WCS axes.

Example:

```python
ShowMap.show_tick(
    ax,
    set_minor=True,
    direction="in",
    line_color="white",
    tick_maj_len=8,
    tick_minor_len=4,
    tick_ra_spac=0.05,
    tick_dec_spac=0.05,
    tick_freq=5
)
```

`ShowMap.draw_scalebar(...)` draws a scalebar in pixel coordinates.

Example:

```python
ShowMap.draw_scalebar(
    ax,
    scalebar_length=80,
    scalebar_text="5 kpc",
    scale_xstar=50,
    scale_ystar=50,
    color="white",
    linewidth=2,
    rotation_angle=0
)
```

You can also draw a rotated scalebar:

```python
ShowMap.draw_scalebar(
    ax,
    scalebar_length=100,
    scalebar_text="10 kpc",
    scale_xstar=60,
    scale_ystar=60,
    rotation_angle=30,
    linestyle="dashed"
)
```

`ShowMap.draw_arrow(...)` adds arrow annotations to the image.

Example:

```python
ShowMap.draw_arrow(
    ax,
    arrow_x=80,
    arrow_y=120,
    arrow_dx=-20,
    arrow_dy=10,
    arrow_color="cyan"
)
```

#### 2.4 show_fits_multi

`ShowMap.show_fits_multi(...)` is used to display multiple FITS images in a multi-panel layout. This is useful for comparisons between different bands, resolutions, or data products.

Example:

```python
from fits_show_map import ShowMap

h1, d1 = ShowMap.load_fits_image("map1.fits")
h2, d2 = ShowMap.load_fits_image("map2.fits")

fig, im = ShowMap.show_fits_multi(
    headers=[h1, h2],
    data_list=[d1, d2],
    nrows=1,
    ncols=2,
    figsize=(14, 6),
    cmap_list=["gray", "magma"],
    title_list=["Map 1", "Map 2"],
    colobar_list=[True, True],
    beam_list=[True, False]
)
```
### 3. cutout_2D
`ShowMap.cutout_2D(...)` creates a fast 2D cutout from a FITS image using sky coordinates. You must provide either `size_pix` or `size_arcsec`.

Example using pixel size:

```python
out_hdr, out_data = ShowMap.cutout_2D(
    infile="your.fits",
    ra_deg=114.099166,
    dec_deg=-69.530833,
    size_pix=(400, 400)
)
```

Example using angular size:

```python
out_hdr, out_data = ShowMap.cutout_2D(
    infile="your.fits",
    ra_deg=114.099166,
    dec_deg=-69.530833,
    size_arcsec=(300, 300)
)
```
### 4. carta_log_stretch
`ShowMap.carta_log_stretch(...)` applies a CARTA-like logarithmic stretch to an image array.

Example:

```python
stretched = ShowMap.carta_log_stretch(data, alpha=1000)
```

`ShowMap.smooth_edge(...)` can be used to generate a smooth alpha mask for masked-image display.

Example:

```python
mask = data > 0.002
masked_data, alpha = ShowMap.smooth_edge(data, sigma_pix=1.0, mask=mask)

fig, ax = ShowMap.show_fits(header=header, data=masked_data, alpha_lim=alpha)
```

In summary, the most commonly used functions are  `ShowMap.load_fits_image`, `ShowMap.project`, `ShowMap.show_fits`, `ShowMap.show_fits_multi`, `ShowMap.show_contour`, `ShowMap.show_beam`, `ShowMap.show_colobar`, `ShowMap.show_tick`, `ShowMap.draw_scalebar`, `ShowMap.draw_arrow`, `ShowMap.cutout_2D`, `ShowMap.carta_log_stretch`, and `ShowMap.smooth_edge`. For more advanced figure production, `ShowMap.show_fits(...)` and `ShowMap.show_fits_multi(...)` provide most of the functionality needed for astronomical plotting.
