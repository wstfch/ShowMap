
#---------------------------------------------------------------------------------------------------#
# To show the fits image                                                                            #
# Creat data : 2023.07.29                                                                           #
# Author: Shengtao Wang                                                                             #
#---------------------------------------------------------------------------------------------------#

#fits show map
import numpy as np
from astropy.io import fits as pf
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from matplotlib import pylab as plt
from matplotlib.colors import LogNorm
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord
from matplotlib.patches import Ellipse, Rectangle, FancyArrowPatch
from regions import EllipseSkyRegion, EllipsePixelRegion
from matplotlib.lines import Line2D
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from reproject import reproject_interp, reproject_exact
from scipy.ndimage import generic_filter
from copy import deepcopy
from astropy import units as u
from astropy.convolution import convolve, Gaussian2DKernel
import os,sys,re
plt.rcParams.update({'font.size': 18})

def get_args():
    ps = ArgumentParser()
    ps.add_argument("image", type=str)
    ps.add_argument("-c","--cmap", type=str,default='viridis')
    ps.add_argument("-s","--savefig", type=str,default='False')
    args = ps.parse_args()
    return args

#-----------------------------------------------------------------------------#
# Load a fits image and return the data and dimensions                        #
#-----------------------------------------------------------------------------#
def load_fits_image(filename):

    # Read the header and image data from the file
    header=pf.getheader(filename)
    data = pf.getdata(filename)
    naxis = len(data.shape)

    # Strip unused dimensions from the data array
    if naxis == 2:
        xydata = data.copy()
        del data
    elif naxis == 3:
        xydata = data[0].copy()
        del data
    elif naxis == 4:
        xydata = data[0][0].copy()
        del data
    elif naxis == 5:
        xydata = data[0][0][0].copy()
        del data
    else:
        print("Data array contains %s axes" % naxis)
        print("This script supports up to 5 axes only.")
        sys.exit(1)

    # Strip unused dimensions from the header
    stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                     'CD1_','CD2_','CUNIT','PC1_','PC2_','PC3_','PC4_']
    
    stripHeadKeys1 = ['PC1_','PC2_','PC3_','PC4_','PC01_','PC02_','PC03_','PC04_']
    
    for i in range(3,6):
        for key in stripHeadKeys:
            if (key+str(i)) in header: del header[key+str(i)]
    for j in range(1,6):            
        for key1 in stripHeadKeys1:
            if (key1+str(j)) in header: del header[key1+str(j)]
            if (key1+ '0' + str(j)) in header: del header[key1+ '0' + str(j)]
    header['NAXIS'] = 2
    if 'WCSAXES' in header:
        header['WCSAXES']=2

    # Determine the coordinate type and the corresponding keyword index
    # Make a note in the header
    ra_regx = re.compile('^RA')
    dec_regx = re.compile('^DEC')
    glon_regx = re.compile('^GLON')
    glat_regx = re.compile('^GLAT')
    if 'CTYPE1' in header:
        if ra_regx.match(header['CTYPE1']) or glon_regx.match(header['CTYPE1']):
            for i in range(int(header['NAXIS'])):
                keyword = "CTYPE"+str(i+1)
                if ra_regx.match(header[keyword]): coord_type="EQU"; x_indx=i+1
                if dec_regx.match(header[keyword]): y_indx=i+1
                if glon_regx.match(header[keyword]): coord_type="GAL"; x_indx=i+1
                if glat_regx.match(header[keyword]): y_indx=i+1
            if not x_indx or not y_indx:
                print("Failed to find Equatorial or Galactic axis coordinate types.")
                del data; del header
                sys.exit(1)
            else:
                header['XINDEX'] = x_indx
                header['YINDEX'] = y_indx
        else:
            print('Does not have "RA-DEC" and "GLON-GLAT" coordinates!!!')
    else:
        print('key "CTYPE1" does not, please check the header!!!')
    # Convert AIPS clean-beam types to standard BMAJ, BMIN, BPA
    try:
        bmaj = header['CLEANBMJ']
        bmin = header['CLEANBMN']
        bpa = header['CLEANBPA']
        header.update('BMAJ',bmaj)
        header.update('BMIN',bmin)
        header.update('BPA',bpa)
    except Exception:
        print("No AIPS style beam keywords found.")

    # Check for PIXSCAL keywords and write to CDELT standard
    try:
        xdelt=(-1)*(header['PIXSCAL'+str(x_indx)])/3600.0
        ydelt=(header['PIXSCAL'+str(y_indx)])/3600.0
        header['CDELT'+str(x_indx)] = xdelt
        header['CDELT'+str(y_indx)] = ydelt
    except Exception:
        pass

    return [header,xydata]  


def show_map(image,cmap='viridis',savefig='False'):
    filename = re.findall('(.+)(?=.fits)', image)[0] 
    header,data = load_fits_image(image)
    print(header)
    wcs = WCS(header)
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(1,1,1, projection=wcs,slices=('x','y'))
    fmax = np.nanpercentile(data,99)
    fmin = np.nanpercentile(data,1)
    im = ax.imshow(data, origin='lower', cmap = cmap, vmin=fmin, vmax=fmax)
    ax.imshow(data, vmin=fmin, vmax=fmax, cmap = cmap, origin='lower', interpolation='spline36')
   
    plt.rcParams['xtick.direction'] = 'in'  # x-axis scale point to inner
    plt.rcParams['ytick.direction'] = 'in'
    ax.set_xlabel('RIGHT ASCENSION', labelpad = 1)
    ax.set_ylabel('DECLINATION(J2000)',labelpad = 2)
    ax.tick_params(direction='in', width='2',length=8,color='k')
    #Set colorbar
    norm = ImageNormalize(data, interval=MinMaxInterval(), stretch=SqrtStretch())
    cbar = plt.colorbar(im, ax=ax,pad=0.007)
    cbar.ax.set_yticklabels(cbar.ax.get_yticks(), fontsize=20)
    #cbar.set_label(r'mJy/beam',fontsize=22,labelpad=15)
    fmt = ticker.FormatStrFormatter('%0.2f')  #Keep one decimal place
    cbar.ax.yaxis.set_major_formatter(fmt)
    if savefig == 'True':
        plt.savefig('./' + filename + '.png', dpi=300) 
        print('The {}.png has saved in local directory!!!'.format(filename))
    plt.show()
    
def cli(args):
    show_map(image=args.image,
             cmap=args.cmap,
             savefig=args.savefig
    )

class ShowMap:
    """
    Date:2023.7.31
    Author: Shengtao Wang 
    This script is used to show one or more fits images and some convenient functions,
    (including map rotate, filling the Nan, cut large image to smaller, regrid, loading fits to standard 2D images, etc.)
    You can also use it in combination with matplotlib for more convenience
    """
    # def __init__(self, fits=False, header=False, data=False, lim_image=False,fontsize=20, colobar=False, beam=False, cont=False):
    #     self.fits = fits
    #     self.lim_image = lim_image
    #     self.fontsize = fontsize
    #     self.colobar = colobar
    #     self.beam = beam
    #     self.cont = cont
    #     if self.fits:
    #         self.header, self.data = ShowMap.load_fits_image(self.fits)
    #     else:
    #         self.header, self.data = header,data
    
    @staticmethod 
    def load_fits_image(filename):
        # Read the header and image data from the file
        with pf.open(filename) as hdul:
            header = hdul[0].header
            data = hdul[0].data
        #header=pf.getheader(filename)

        #data = pf.getdata(filename)
        naxis = len(data.shape)
        # Strip unused dimensions from the data array
        if naxis == 2:
            xydata = data.copy()
            del data
        elif naxis == 3:
            xydata = data[0].copy()
            del data
        elif naxis == 4:
            xydata = data[0][0].copy()
            del data
        elif naxis == 5:
            xydata = data[0][0][0].copy()
            del data
        else:
            print("Data array contains %s axes" % naxis)
            print("This script supports up to 5 axes only.")
            sys.exit(1)

        # Strip unused dimensions from the header
        stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                        'CD1_','CD2_','CUNIT','PC1_','PC2_','PC01_','PC02_']
        
        stripHeadKeys1 = ['PC1_','PC2_','PC3_','PC4_','PC01_','PC02_','PC03_','PC04_']
        #stripHeadKeys1 = ['PC3_','PC4_','PC03_','PC04_']
        
        for i in range(3,6):
            for key in stripHeadKeys:
                if (key+str(i)) in header: del header[key+str(i)]

        for j in range(1,6):            
            for key1 in stripHeadKeys1:
                if (key1+str(j)) in header: del header[key1+str(j)]
                if (key1+ '0' + str(j)) in header: del header[key1+ '0' + str(j)]
        header['NAXIS'] = 2
        if 'WCSAXES' in header:
            header['WCSAXES']=2
        # Determine the coordinate type and the corresponding keyword index
        # Make a note in the header
        ra_regx = re.compile('^RA')
        dec_regx = re.compile('^DEC')
        glon_regx = re.compile('^GLON')
        glat_regx = re.compile('^GLAT')
        if 'CTYPE1' in header:
            if ra_regx.match(header['CTYPE1']) or glon_regx.match(header['CTYPE1']):
                for i in range(int(header['NAXIS'])):
                    keyword = "CTYPE"+str(i+1)
                    if ra_regx.match(header[keyword]): coord_type="EQU"; x_indx=i+1
                    if dec_regx.match(header[keyword]): y_indx=i+1
                    if glon_regx.match(header[keyword]): coord_type="GAL"; x_indx=i+1
                    if glat_regx.match(header[keyword]): y_indx=i+1
                if not x_indx or not y_indx:
                    print("Failed to find Equatorial or Galactic axis coordinate types.")
                    del data; del header
                    sys.exit(1)
                else:
                    header['XINDEX'] = x_indx
                    header['YINDEX'] = y_indx
            else:
                print('Does not have "RA-DEC" and "GLON-GLAT" coordinates!!!')
        else:
            print('key "CTYPE1" does not, please check the header!!!')
        # Convert AIPS clean-beam types to standard BMAJ, BMIN, BPA
        try:
            bmaj = header['CLEANBMJ']
            bmin = header['CLEANBMN']
            bpa = header['CLEANBPA']
            header.update('BMAJ',bmaj)
            header.update('BMIN',bmin)
            header.update('BPA',bpa)
        except Exception:
            print("No AIPS style beam keywords found.")

        # Check for PIXSCAL keywords and write to CDELT standard
        try:
            xdelt=(-1)*(header['PIXSCAL'+str(x_indx)])/3600.0
            ydelt=(header['PIXSCAL'+str(y_indx)])/3600.0
            header['CDELT'+str(x_indx)] = xdelt
            header['CDELT'+str(y_indx)] = ydelt
        except Exception:
            pass
        return [header,xydata]
        
    @staticmethod
    def project(header0, data0, RA, DEC, NAXIS1, NAXIS2, CRPIX1, CRPIX2, CDELT1=None, CDELT2=None, rotate=None, fill_nan = False, fill_size = 30):
        if isinstance(RA, str) and isinstance(DEC, str):
            coord = SkyCoord(RA, DEC)
            RA, DEC = coord.ra.degree, coord.dec.degree
        if CDELT1 is None:CDELT1 = header0['CDELT1']
        if CDELT2 is None:CDELT2 = header0['CDELT2']
        if rotate is not None:
            # Convert rotate to radians
            theta_rad = np.deg2rad(rotate)
            cos_theta = np.cos(theta_rad)
            sin_theta = np.sin(theta_rad)
            # Create the rotation matrix elements
            PC1_1 = cos_theta
            PC1_2 = -sin_theta
            PC2_1 = sin_theta
            PC2_2 = cos_theta
        # Define the new header
        hdu0 = deepcopy(header0)
        hdu0['NAXIS1'] = NAXIS1
        hdu0['NAXIS2'] = NAXIS2
        hdu0['CRVAL1'] = RA
        hdu0['CRVAL2'] = DEC
        hdu0['CRPIX1'] = CRPIX1
        hdu0['CRPIX2'] = CRPIX2
        hdu0['CDELT1'] = CDELT1
        hdu0['CDELT2'] = CDELT2
        if rotate is not None:
            hdu0['PC1_1'] = PC1_1
            hdu0['PC1_2'] = PC1_2
            hdu0['PC2_1'] = PC2_1
            hdu0['PC2_2'] = PC2_2
        # Reproject the data
        hdu = pf.PrimaryHDU(data0, header=header0)
        data, footprint = reproject_exact(hdu, hdu0, shape_out=(NAXIS2, NAXIS1))
        # To filling the Nan for the image
        if fill_nan:
            def nan_mean_filter(values):
                valid_values = values[np.isfinite(values)]
                if len(valid_values) > 0:
                    return np.mean(valid_values)
                else:
                    return np.nan
            # Create a mask include Nan
            nan_mask = np.isnan(data)
            # only filling the include Nan region
            filled_data = data.copy()
            filled_data[nan_mask] = generic_filter(data, nan_mean_filter, size=fill_size)[nan_mask]    
            data = filled_data   
        return hdu0, data

    @staticmethod
    def show_fits(header, data, lim_image=False,colobar=False,beam=False,cont=False,log=False, fontsize = 20, cmap='viridis', figsize=(12, 9), auto_scaler=True, max=99, min=0,cb_dedi='%0.2f',xpad=1,ypad=2,line_width=2,alpha_lim=None,\
                cb_pad=0.007,cb_loct='right',cb_font=20, cb_aspect=None,cb_shrink=None,cb_percent=False,cb_show_ticks=True, decimals=1,direction='in',tick_maj_len=8,tick_minor_len=4,tick_ra_spac=0.1,tick_dec_spac=0.1,tick_freq=5,set_minor=False,\
                title=None, title_pad=10, cb_lab=None, beam_p_pix=[15, 15], beam_fluc=[2, 2], beam_color='red', line_color='w', ticks=None, beam_sque=True,\
                beam_sque_linw=1.2,RA='RIGHT ASCENSION (J2000)',DEC='DECLINATION (J2000)',NAXIS1=None,NAXIS2=None,savefits=None,cont_data=None,cont_levels=None,cont_color='red',\
                alpha = 100, ticks_n=7, cb_size=0.08, color_pad=0.01, span=(0.1, 0.9),fmt="{:.2f}",cbmin=None,cbmax=None,\
                cont_alpha=0.8,CRPIX1=None,CRPIX2=None,CDELT1=None,CDELT2=None,rotate=None,savefig=None,dpi=300):
        
        plt.rcParams.update({'font.size': 20})
        #header, data = self.header,self.data
        #fontsize=self.fontsize
        if lim_image:
            if RA==None:RA=header['CRVAL1']
            if DEC==None:DEC=header['CRVAL2']  
            if NAXIS1==None:NAXIS1=header['NAXIS1']
            if NAXIS2==None:NAXIS2=header['NAXIS2']
            if CRPIX1==None:CRPIX1=header['CRPIX1']
            if CRPIX2==None:CRPIX2=header['CRPIX2']
            header, data=ShowMap.project(header,data,RA,DEC,NAXIS1,NAXIS2,CRPIX1,CRPIX2,CDELT1,CDELT2,rotate)
        else:
            print('Images are not restricted or clipped!!!')

        if savefits is not None:
            pf.writeto(savefits,data,header,overwrite=True)
            print('The fits has saved in:',savefits)
        wcs = WCS(header)
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=wcs, slices=('x', 'y'))
        if auto_scaler:
            fmax = np.nanpercentile(data, max)
            fmin = np.nanpercentile(data, min)
            if log:
                data_log = ShowMap.carta_log_stretch(data, alpha=alpha, vmin=fmin, vmax=fmax)
        else:
            fmax = max
            fmin = min
        if log:
            ax.imshow(data_log, vmin=0, vmax=1, cmap=cmap, alpha=alpha_lim, origin='lower', interpolation='spline36')
        else:
            im = ax.imshow(data, origin='lower',alpha=alpha_lim, cmap=cmap, vmin=fmin, vmax=fmax)
            ax.imshow(data, vmin=fmin, vmax=fmax, alpha=alpha_lim, cmap=cmap, origin='lower', interpolation='spline36')

        if colobar:
            """
            Show the colorbar with two different type ("top" and "right"), if "top" that these two parameters efficient "aspect, shrink".
            """
            if log:
                ShowMap.show_colobar_liner(fig, ax, alpha, ticks_n,cb_loct, cb_size, color_pad, cb_pad, cb_lab, cmap,span, fmt, cb_font, cbmin=fmin, cbmax=fmax)
            else:
                ShowMap.show_colobar(im, ax, data, cb_loct, cb_pad, cb_dedi, cb_lab, cb_font, cb_aspect, cb_shrink,cb_percent,decimals,cb_show_ticks)

        if beam:
            ShowMap.show_beam(ax, header, beam_p_pix, beam_fluc, beam_color,beam_sque,beam_sque_linw)

        if cont:
            if cont_data is None: cont_data=data  
            ShowMap.show_contour(ax, cont_data=cont_data, cont_levels=cont_levels, cont_color=cont_color, cont_alpha=cont_alpha) 

        plt.rcParams['xtick.direction'] = direction
        plt.rcParams['ytick.direction'] = direction
        # Set the coord interval of RA and DEC
        #ax.coords[0].set_ticks(spacing=0.08 * u.deg, exclude_overlapping=True)  
        #ax.coords[1].set_ticks(spacing=0.06 * u.deg, exclude_overlapping=True)
        ax.set_xlabel(RA, labelpad=xpad)
        ax.set_ylabel(DEC, labelpad=ypad)
        ax.set_title(title,fontsize=fontsize,pad=title_pad)
        #ax.tick_params(direction='in', width=line_width, length=8, color=line_color)
        ShowMap.show_tick(ax,set_minor=set_minor,direction=direction,line_color=line_color,line_width=line_width,tick_maj_len=tick_maj_len,\
            tick_minor_len=tick_minor_len,tick_ra_spac=tick_ra_spac,tick_dec_spac=tick_dec_spac,tick_freq=tick_freq)

        if savefig is not None:
            plt.savefig(savefig, bbox_inches = 'tight', dpi=dpi) 
            #print('The {}.png has saved!!!'.format(self.filename))
            print('The picture has saved in:',savefig)
        #plt.show()
        return [fig, ax]


    @staticmethod
    def show_fits_multi(headers, data_list, \
                        nrows=1, ncols=2, figsize=(12, 9),colobar_list=None, cont_list=None,\
                        fontsize=20, cmap_list=None, auto_scaler=True, x_label_list=None, y_label_list=None,y_label_hide_list=None,wspace=0., hspace=0., \
                        max_list=None, min_list=None, cb_dedi_list=None, \
                        xpad_list=None, ypad_list=None, line_width_list=None, \
                        cb_pad_list=None, cb_loct_list=None, cb_font_list=None, \
                        cb_aspect_list=None, cb_shrink_list=None, cb_percent_list=None, \
                        decimals_list=None, direction_list=None, tick_maj_len_list=None, \
                        tick_minor_len_list=None, tick_ra_spac_list=None, \
                        tick_dec_spac_list=None, tick_freq_list=None, set_minor_list=None,\
                        title_list=None, title_pad_list=None, cb_lab_list=None,\
                        savefig=None, dpi=300, beam_list=None, beam_p_pix_list=None, beam_fluc_list=None, beam_color_list=None, beam_sque_list=None, cont_data_list=None, cont_levels_list=None, beam_sque_linw_list=None,\
                        cont_color_list=None, set_minor=False, direction='in',line_color='k',line_width=1.2,tick_maj_len=6,tick_minor_len=3,tick_ra_spac=0.1,tick_dec_spac=0.1,tick_freq=5):
        
        # # Ensure headers and data are lists
        # assert isinstance(headers, list), "Headers must be a list of FITS headers"
        # assert isinstance(data_list, list), "Data must be a list of FITS data arrays"
        # assert len(headers) == len(data_list), "Headers and data arrays must have the same length"
        
        num_plots = len(headers)  # Number of subplots

        # Set default values if parameter lists are not provided
        cmap_list = cmap_list or ['viridis'] * num_plots
        max_list = max_list or [99] * num_plots
        min_list = min_list or [0] * num_plots
        x_label_list = x_label_list or ['RIGHT ASCENSION (J2000)'] * num_plots
        y_label_list = y_label_list or ['DECLINATION (J2000)'] * num_plots
        xpad_list = xpad_list or [1] * num_plots
        ypad_list = ypad_list or [2] * num_plots
        line_width_list = line_width_list or [2] * num_plots
        cont_list = cont_list or [False] * num_plots
        colobar_list = colobar_list or [False] * num_plots
        cb_pad_list = cb_pad_list or [0.007] * num_plots
        cb_loct_list = cb_loct_list or ['right'] * num_plots
        cb_font_list = cb_font_list or [20] * num_plots
        cb_aspect_list = cb_aspect_list or [None] * num_plots
        cb_shrink_list = cb_shrink_list or [None] * num_plots
        cb_percent_list = cb_percent_list or [False] * num_plots
        cb_dedi_list = cb_dedi_list or ['%0.2f'] * num_plots
        decimals_list = decimals_list or [1] * num_plots      # There used to limt decimal digits if use percent
        direction_list = direction_list or ['in'] * num_plots
        tick_maj_len_list = tick_maj_len_list or [8] * num_plots
        tick_minor_len_list = tick_minor_len_list or [4] * num_plots
        tick_ra_spac_list = tick_ra_spac_list or [0.1] * num_plots
        tick_dec_spac_list = tick_dec_spac_list or [0.1] * num_plots
        tick_freq_list = tick_freq_list or [5] * num_plots
        set_minor_list = set_minor_list or [False] * num_plots
        title_list = title_list or [None] * num_plots
        title_pad_list = title_pad_list or [10] * num_plots
        cb_lab_list = cb_lab_list or [None] * num_plots
        beam_list = beam_list or [False] * num_plots
        beam_color_list = beam_color_list or ['red'] * num_plots
        beam_sque_list = beam_sque_list or [True] * num_plots
        beam_p_pix_list = beam_p_pix_list or [[15,15]] * num_plots
        beam_fluc_list = beam_fluc_list or [[2,2]] * num_plots
        beam_sque_linw_list = beam_sque_linw_list or [1.2] * num_plots
        cont_data_list = cont_data_list or [None] * num_plots
        cont_levels_list = cont_levels_list or [None] * num_plots
        cont_color_list = cont_color_list or ['red'] * num_plots
        y_label_hide_list = y_label_hide_list or [False] * num_plots
        
        # Create figure and subplots
        #fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False)
        plt.rcParams.update({'font.size': fontsize})
        fig = plt.figure(figsize=figsize)
        fig.subplots_adjust(wspace=wspace, hspace=hspace)
        for i in range(len(headers)):
            wcs = WCS(headers[i])
            ax = fig.add_subplot(nrows, ncols, i+1, projection=wcs, slices=('x', 'y'))
            if auto_scaler:
                fmax = np.nanpercentile(data_list[i], max_list[i])
                fmin = np.nanpercentile(data_list[i], min_list[i])
            else:
                fmax = max_list[i]
                fmin = min_list[i]
            im = ax.imshow(data_list[i], origin='lower', cmap=cmap_list[i], vmin=fmin, vmax=fmax)

            # Add Colorbar for Each Subplot
            if colobar_list[i]:
                #fig.colorbar(im, ax=ax, orientation=cb_loct_list[i], pad=cb_pad_list[i])
                ShowMap.show_colobar(im, ax, data_list[i], cb_loct_list[i], cb_pad_list[i], cb_dedi_list[i], cb_lab_list[i], cb_font_list[i], cb_aspect_list[i], cb_shrink_list[i], cb_percent_list[i], decimals_list[i])
            # Apply Beam Settings
            if beam_list[i]:
                 ShowMap.show_beam(ax, headers[i], beam_p_pix_list[i], beam_fluc_list[i], beam_color_list[i],beam_sque_list[i],beam_sque_linw_list[i])

            # Apply Contours
            if cont_list[i]:
                ShowMap.show_contour(ax, cont_data=cont_data_list[i], cont_levels=cont_levels_list[i], cont_color=cont_color_list[i])

            ShowMap.show_tick(ax,set_minor=set_minor,direction=direction,line_color=line_color,line_width=line_width,tick_maj_len=tick_maj_len,\
            tick_minor_len=tick_minor_len,tick_ra_spac=tick_ra_spac,tick_dec_spac=tick_dec_spac,tick_freq=tick_freq)

            # Set Labels and Title
            ax.set_xlabel(x_label_list[i], labelpad=xpad_list[i])
            ax.set_ylabel(y_label_list[i], labelpad=ypad_list[i])
            if y_label_hide_list[i]:
                ax.coords[1].set_ticklabel_visible(False)  # Hide the y-axis label
            if title_list[i] is not None:
                ax.set_title(title_list[i], fontsize=fontsize, pad=title_pad_list[i])
        # Adjust Layout
        #plt.tight_layout()
        # Save Figure
        if savefig is not None:
            plt.savefig(savefig, bbox_inches='tight', dpi=dpi)
            print(f'The figure has been saved at: {savefig}')
        return [fig, im]

    # To draw the colorbar
    @staticmethod
    def show_colobar(im, ax, data, cb_loct='right', cb_pad=0.007, cb_dedi='%0.2f', cb_lab=None, cb_font=20, cb_aspect=None, cb_shrink=None, cb_percent=False, decimals=0, cb_show_ticks=True):
        norm = ImageNormalize(data, interval=MinMaxInterval(), stretch=SqrtStretch())
        colorbar_kwargs = {'pad': cb_pad}
        if cb_aspect is not None:
            colorbar_kwargs['aspect'] = cb_aspect
        if cb_shrink is not None:
            colorbar_kwargs['shrink'] = cb_shrink

        if cb_loct == 'right':
            cbar = plt.colorbar(im, ax=ax, **colorbar_kwargs)
        else:
            colorbar_kwargs.update({'orientation': 'horizontal', 'location': 'top'})
            cbar = plt.colorbar(im, ax=ax, **colorbar_kwargs)

        # Custom tick formatter
        if cb_percent:
            formatter = ticker.PercentFormatter(xmax=1, decimals=decimals)  # Scale values to 100%
        else:
            formatter = ticker.FormatStrFormatter(cb_dedi)
        
        if cb_loct == 'right':
            cbar.ax.yaxis.set_major_formatter(formatter)
        elif cb_loct == 'top':
            cbar.ax.xaxis.set_major_formatter(formatter)

        if not cb_show_ticks:
            if cb_loct == 'right':
                cbar.ax.set_yticklabels([])
            elif cb_loct == 'top':
                cbar.ax.set_xticklabels([])

        #formatter = ticker.FormatStrFormatter(cb_dedi)
        #cbar.ax.xaxis.set_major_formatter(formatter)  # Note that it is changed to xaxis here

        cbar.ax.tick_params(labelsize=cb_font)  # Set the font size of the scale label
        if cb_lab is not None:
            cbar.set_label(cb_lab, fontsize=cb_font, labelpad=15)  # Set the color bar label

    @staticmethod
    def show_colobar1(pos=[0.98, 0.206, 0.03, 0.67],ticks=None,cb_lab=None,cb_pad=20):
        # Set the position and scale of colorbar
        position=fig.add_axes(pos) # [left,down,right,up]
        cbar = fig.colorbar(im, cax = position)
        cbar.set_label(cb_lab,labelpad=cb_pad)
        cbar.set_ticks(ticks)

    # To draw the beam 
    @staticmethod
    def show_beam(ax,header, beam_p_pix=[15,15], beam_fluc=[2,2], beam_color='red', beam_sque=True,beam_sque_linw = 1.2):
        #header = self.header
        pix_size = np.abs(header['CDELT1'])
        p_pix = beam_p_pix
        fluc = beam_fluc
        color = beam_color
        beam_size_x = header['BMAJ'] / pix_size
        beam_size_y = header['BMIN'] / pix_size
        beam = Ellipse((p_pix[0], p_pix[1]), beam_size_x, beam_size_y, angle=header['BPA'],
                       facecolor=color, edgecolor=color, alpha=1, zorder=200)
        ax.add_patch(beam)
        if beam_sque:
            square_center_x, square_center_y = p_pix[0] - fluc[0] / 2, p_pix[1] - fluc[1] / 2
            half_width, half_height = beam_size_x / 2, beam_size_y / 2
            square_x = square_center_x - (beam_size_x / 2)
            square_y = square_center_y - (beam_size_y / 2)
            square = Rectangle(xy=(square_x, square_y), width=beam_size_x + fluc[0],
                            height=beam_size_y + fluc[1], facecolor='none', edgecolor=color, linewidth=beam_sque_linw, alpha=1, zorder=200)
            ax.add_patch(square)

    # To draw the contour
    @staticmethod        
    def show_contour(ax, cont_data=None,cont_levels=None,cont_color='red',cont_alpha=0.8):
        ax.contour(cont_data, levels=cont_levels, colors=cont_color, alpha=cont_alpha) 
    
    @staticmethod
    def show_tick(ax,set_minor=False,direction='in',line_color='w',line_width=1.5,tick_maj_len=8,tick_minor_len=4,tick_ra_spac=0.1,tick_dec_spac=0.1,tick_freq=5):
        """
        set_minor: If Ture, set the minor scale
        direction: The direction of scale ('in' or 'out')
        tick_spac: The spacing of the main scale.(unit:deg)
        tick_freq: The number of minor scale.
        """
        if set_minor:
            # Set major ticks on the WCS coordinates
            lon = ax.coords['ra']
            lat = ax.coords['dec']
            # Set major ticks every 0.1 degree
            lon.set_ticks(spacing=tick_ra_spac * u.deg)
            lat.set_ticks(spacing=tick_dec_spac * u.deg)

            # Customize the appearance of the major ticks
            lon.tick_params(which='major', direction=direction, color=line_color, width=line_width, length=tick_maj_len)
            lat.tick_params(which='major', direction=direction, color=line_color, width=line_width, length=tick_maj_len)

            # Set minor ticks
            lon.set_minor_frequency(tick_freq)  # Increase the number of minor ticks
            lat.set_minor_frequency(tick_freq)

            # Display minor ticks and set their properties
            lon.display_minor_ticks(True)
            lat.display_minor_ticks(True)

            # Set minor tick width, length, and color, note: If set which = 'minor', only support other parameters 'length'
            lon.tick_params(which='minor',length=tick_minor_len)
            lat.tick_params(which='minor',length=tick_minor_len)

            lon.set_ticks_position('bt')   # bottom+top
            lat.set_ticks_position('lr')   # left+right

        else:
            ax.tick_params(direction=direction, width=line_width, length=tick_maj_len, color=line_color)
        
    @staticmethod
    def draw_scalebar(ax, scalebar_length, scalebar_text, scale_xstar, scale_ystar, text_font=20, fluc=4, color='white',
                  linewidth=2, mutation_scale=15, rotation_angle=0, text_yposition=None, linestyle='solid'):
        """
        Draw a scalebar with rotation support.
        Parameters:
        scalebar_length : float
            The length of the scalebar in pixels.
        scalebar_text : str
            Text to label the scalebar (e.g., "10 kpc").
        scale_xstar, scale_ystar : float
            The coordinates of the starting point of the scalebar.
        text_font : int
            Font size for the scalebar text.
        fluc : float
            Offset for the arrow ends (left or right).
        color : str
            Color of the scalebar and arrows.
        linewidth : float
            Line width of the scalebar.
        mutation_scale : float
            Scale of the arrow heads.
        rotation_angle : float
            Rotation angle of the scalebar in degrees.
        text_yposition : float or None
            The y position for the text; if None, a default value is used.
        linestyle : str
            Linestyle of the scalebar ('solid' or 'dashed').
        """
        # Function to apply rotation
        def rotate_point(x, y, cx, cy, angle_deg):
            angle_rad = np.radians(angle_deg)
            cos_theta = np.cos(angle_rad)
            sin_theta = np.sin(angle_rad)
            x_rot = cos_theta * (x - cx) - sin_theta * (y - cy) + cx
            y_rot = sin_theta * (x - cx) + cos_theta * (y - cy) + cy
            return x_rot, y_rot

        # Rotate the endpoints of the scalebar
        x_end = scale_xstar + scalebar_length
        y_end = scale_ystar
        x_start_rot, y_start_rot = rotate_point(scale_xstar, scale_ystar, scale_xstar, scale_ystar, rotation_angle)
        x_end_rot, y_end_rot = rotate_point(x_end, y_end, scale_xstar, scale_ystar, rotation_angle)

        # Draw the main line of the scalebar
        line = Line2D([x_start_rot, x_end_rot], [y_start_rot, y_end_rot], color=color, linewidth=linewidth, linestyle=linestyle)
        ax.add_line(line)

        # Rotate the positions of the arrows
        left_arrow_x, left_arrow_y = rotate_point(scale_xstar - fluc, scale_ystar, scale_xstar, scale_ystar, rotation_angle)
        right_arrow_x, right_arrow_y = rotate_point(x_end + fluc, y_end, scale_xstar, scale_ystar, rotation_angle)

        # Draw the arrows
        ax.add_patch(FancyArrowPatch((left_arrow_x, left_arrow_y), (x_start_rot, y_start_rot),
                                    arrowstyle='<|-', color=color, mutation_scale=mutation_scale))
        ax.add_patch(FancyArrowPatch((right_arrow_x, right_arrow_y), (x_end_rot, y_end_rot),
                                    arrowstyle='<|-', color=color, mutation_scale=mutation_scale))

        # Rotate the text position and add it
        if text_yposition is not None:
            ytex_positon = text_yposition
        else:
            ytex_positon = scale_ystar + 5
        text_x, text_y = rotate_point(scale_xstar + scalebar_length / 2, ytex_positon, scale_xstar, scale_ystar, rotation_angle)
        ax.text(text_x, text_y, scalebar_text, color=color, ha='center', fontsize=text_font, rotation=rotation_angle)
    
    @staticmethod
    def draw_arrow(ax, arrow_x = 20,arrow_y = 50,arrow_dx = -10,arrow_dy = -3,arrow_color='white',arrow_width=1, arrow_headwidth=5, arrow_headlength=10):
        # draw arrow 
        """
        arrow_x: start point (x axis)
        arrow_y: start point (y axis)
        arrow_dx: x add 
        arrow_dy: y add    
        """
        ax.annotate('', xy=(arrow_x + arrow_dx, arrow_y + arrow_dy), xytext=(arrow_x, arrow_y),\
                    arrowprops=dict(facecolor=arrow_color, edgecolor=arrow_color, shrink=0.05, width=arrow_width, headwidth=arrow_headwidth, headlength=arrow_headlength))

    @staticmethod
    def carta_log_stretch(data, alpha=1000, vmin=None, vmax=None):
        if vmin is None: vmin = np.nanmin(data)
        if vmax is None: vmax = np.nanmax(data)
        x = (data - vmin) / (vmax - vmin)
        x = np.clip(x, 0, 1)
        y = np.log1p(alpha * x) / np.log1p(alpha)
        return y

    @staticmethod
    def inv_stretch_scalar(y, alpha, vmin, vmax):
        """
        y: Show (0..1) position → return to liner space value
        """
        y = np.asarray(y, dtype=float)
        # x is the ratio of normalization 
        x = np.expm1(y * np.log1p(alpha)) / alpha
        x = np.clip(x, 0, 1)
        return vmin + x * (vmax - vmin)

    @staticmethod
    def show_colobar_liner(fig, ax, alpha=100, ticks_n=7,cb_loct='right', cb_size=0.08, color_pad=0.01, cb_pad=15, cb_lab=None, cmap=None,span=(0.1, 0.9), fmt="{:.2f}", cb_font=20, cbmin=0.0, cbmax=1.0):
        """
        Under the log display: The colorbar scale positions are at equal intervals in the display space, but the labels show linear values.
        - cbmin/cbmax: linear space vmin/vmax（similar to carta_log_stretch）
        - span: Avoid the scale sticking to the edge (default: 0.1 to 0.9)
        - fmt: Linear label format
        """
        if not ax.images:
            raise RuntimeError("No images found on this axes to build a colorbar.")
        im = ax.images[-1]                 # Take the last imshow layer
        #im.set_interpolation('nearest')    # Optional: Avoid visual deviation
        if cmap is not None:
            im.set_cmap(cmap)              
        cb = fig.colorbar(im, ax=ax, location=cb_loct, fraction=cb_size, pad=color_pad)

        # Display space (0..1) Equal-spaced scale positions
        a, b = span
        if not (0.0 <= a < b <= 1.0):
            raise ValueError("span must be within [0,1] and a < b.")
        ticks_disp = np.linspace(a, b, max(2, int(ticks_n)))
        # Reverse transformation to linear values for labels
        ticks_lin = ShowMap.inv_stretch_scalar(ticks_disp, alpha, cbmin, cbmax)
        cb.set_ticks(ticks_disp)
        cb.set_ticklabels([fmt.format(t) for t in ticks_lin])
        if cb_lab is not None:
            cb.set_label(cb_lab, fontsize=cb_font, labelpad=cb_pad)
        return cb

    @staticmethod
    def smooth_edge(data, sigma_pix, mask):
        """
        This function used to smooth the image edge if you want to limt its data.
        data: You want to show data
        threshold: set the threshould you want to masked, it is masked to Nan if smaller the value.
        sigma_pix: the smooth kernel in pixel (0.6~1.2) 
        data_mask: The shape same as data that you used to limit data.
        """
        data_new_spec = np.ones(data.shape)
        data_new_spec[:,:] = np.nan
        mask_data = np.where(mask, data, np.nan)
        sigma_pix = sigma_pix                  # Between 0.6~1.2
        soft = convolve(mask.astype(float), Gaussian2DKernel(sigma_pix),boundary='extend', nan_treatment='interpolate')
        alpha = np.clip(soft, 0.0, 1.0) # Normalize soft to [0,1] as transparency (alpha)
        return mask_data, alpha

    @staticmethod
    def cutout_2D(infile, ra_deg, dec_deg, size_pix=None, size_arcsec=None):
        """
        Fast cutout from a large FITS image WITHOUT reproject/resampling.
        Parameters
        ----------
        ra_deg, dec_deg : float
            Cutout center in degrees (ICRS).
        size_pix : (nx, ny) or int
            Cutout size in pixels. If int, use (int, int).
        size_arcsec : (sx, sy) or float
            Cutout size in arcsec. Requires valid WCS pixel scale. If float, use (float, float).
        ext : int
            FITS extension index.
        """
        if (size_pix is None) == (size_arcsec is None):
            raise ValueError("Provide exactly one of size_pix or size_arcsec.")

        hdr, data = ShowMap.load_fits_image(infile)

        # Handle possible extra dimensions (e.g., (1,1,Ny,Nx))
        # We cut the last two axes as (y,x)
        w = WCS(hdr)

        center = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame="icrs")

        if size_pix is not None:
            if isinstance(size_pix, (int, np.integer)):
                size = (int(size_pix), int(size_pix))
            else:
                size = (int(size_pix[0]), int(size_pix[1]))
        else:
            if isinstance(size_arcsec, (int, float, np.floating)):
                size_arcsec = (float(size_arcsec), float(size_arcsec))
            size = (size_arcsec[0]*u.arcsec, size_arcsec[1]*u.arcsec)

        # Build Cutout2D on 2D plane
        # If data has more dims, we cut only the last two dims for every leading plane.
        if data.ndim == 2:
            cut = Cutout2D(data, position=center, size=size, wcs=w, mode="trim")
            out_data = cut.data
            out_wcs = cut.wcs
        else:
            print('Only support 2D fits image!')

        out_hdr = out_wcs.to_header()
        #fits.PrimaryHDU(data=out_data, header=out_hdr).writeto(outfile, overwrite=overwrite)
        return out_hdr, out_data
                    




        
