#!/usr/bin/env python
# Routine to check quality of LOFAR images
def main(msin,config_path, python_path, fits_path, tgss_server, nvss_server, first_server, lots_server):
    ''' Main entry function called by the genericpipeline framework.
    Args:
        msin (str): input measurement set being processed.
        config_path (str): path to the pipeline cfg file.
        python_path (str): additional paths to be added to the Python path at runtime.
        tgss_server (str): server to download the TGSS catalogue from.
    Returns:
        0 (int)
    '''
    import os
    import os.path
    import sys


    try:
        import bdsf as bdsm
    except ImportError:
        import lofar.bdsm as bdsm
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    sys.path.append(os.path.abspath(python_path))
    #from quality_parset import option_list  ## This file is not needed anymore
    from options import options,print_options,download_cat
    from astropy.io import fits
    from astropy.table import Table
    from auxcodes import report,run,get_rms,warn,die,sepn
    from crossmatch_utils import match_catalogues,filter_catalogue,select_isolated_sources,bootstrap
    from quality_make_plots import plot_flux_ratios,plot_flux_errors,plot_position_offset
    from getcpus import getcpus

    # Define various angle conversion factors.
    arcsec2deg=1.0/3600
    arcmin2deg=1.0/60
    deg2rad=np.pi/180
    deg2arcsec = 1.0/arcsec2deg
    rad2deg=180.0/np.pi
    arcmin2rad=arcmin2deg*deg2rad
    arcsec2rad=arcsec2deg*deg2rad
    rad2arcmin=1.0/arcmin2rad
    rad2arcsec=1.0/arcsec2rad
    steradians2degsquared = (180.0/np.pi)**2.0
    degsquared2steradians = 1.0/steradians2degsquared
    cat_path= '/'.join(python_path.split('/')[:-2]) + '/catalogues/' # Path to folder with catalogues


    #get all the catalogues
    catlist=[tgss_server,nvss_server,lots_server,first_server]
    for i,cat in enumerate(("TGSS","NVSS","LOTS","FIRST")):
        if catlist[i].upper() != 'NONE':
            catlist[i] = download_cat(cat_path, catlist[i])
        else:
            catlist[i] = None

    #     try:
    #         o[cat]=o['filenames'][i]
    #     except:
    #         pass


    option_list = ( ( 'machine', 'NCPU', int, getcpus() ),
                    ( 'image', 'pbimage', str, fits_path, 'PB-corrected image to use for source finding' ),
                    ( 'image', 'catprefix', str, 'image_full_ampphase1m', 'Prefix to use for output catalogues' ),
                    ( 'control', 'sfind', bool, True, 'Do source finding?' ),
                    ( 'control', 'sfind_pixel_fraction', float, 0.5, 'Source find over what fraction of the image?' ),
                    ( 'control', 'quiet', bool, False ),
                    ( 'control', 'logging', str, 'logs' ),
                    ( 'control', 'dryrun', bool, False ),
                    ( 'control', 'restart', bool, True ),
                    ( 'pybdsm', 'atrous', bool, True ),
                    ( 'comparison_cats', 'list', list, None ),
                    ( 'comparison_cats', 'filenames', list, None),
                    ( 'comparison_cats', 'radii', list, None),
                    ( 'comparison_cats', 'fluxfactor', list, None),
                    ( 'comparison_cats', 'TGSS', str, catlist[0] ),
                    ( 'comparison_cats', 'TGSS_matchrad', float, 10.0 ),
                    ( 'comparison_cats', 'TGSS_match_majkey1', float, 'Maj_1' ),
                    ( 'comparison_cats', 'TGSS_match_majkey2', float, 'Maj_2' ),
                    ( 'comparison_cats', 'TGSS_filtersize', float, 40.0 ),
                    ( 'comparison_cats', 'TGSS_fluxfactor', float, 1000.0 ),
                    ( 'comparison_cats', 'NVSS', str, catlist[1] ),
                    ( 'comparison_cats', 'NVSS_matchrad', float, 10.0 ),
                    ( 'comparison_cats', 'NVSS_match_majkey1', float, 'Maj_1' ),
                    ( 'comparison_cats', 'NVSS_match_majkey2', float, 'Maj_2' ),
                    ( 'comparison_cats', 'NVSS_filtersize', float, 40.0 ),
                    ( 'comparison_cats', 'NVSS_fluxfactor', float, 1000.0 ),
                    ( 'comparison_cats', 'LOTS', str, catlist[2] ),
                    ( 'comparison_cats', 'LOTS_matchrad', float, 10.0 ),
                    ( 'comparison_cats', 'LOTS_match_majkey1', float, 'Maj_1' ),
                    ( 'comparison_cats', 'LOTS_match_majkey2', float, 'Maj_2' ),
                    ( 'comparison_cats', 'LOTS_filtersize', float, 40.0 ),
                    ( 'comparison_cats', 'LOTS_fluxfactor', float, 1000.0 ),
                    ( 'comparison_cats', 'FIRST', str, catlist[3] ),
                    ( 'comparison_cats', 'FIRST_matchrad', float, 10.0 ),
                    ( 'comparison_cats', 'FIRST_match_majkey1', float, 'Maj' ),
                    ( 'comparison_cats', 'FIRST_match_majkey2', float, 'MAJOR' ),
                    ( 'comparison_cats', 'FIRST_filtersize', float, 10.0 ),
                    ( 'comparison_cats', 'FIRST_fluxfactor', float, 1.0 ) )







    def logfilename(s,options=None):
        ''' Returns full path to the log file, using the 'logging' parameter defined in the options.
        Args:
            s (str): log file name.
            options (dict): dictionary containig the options. This is read from quality_parset.py.
        Returns:
            log path (str): full path to the logging file.
            OR
            None: when options['logging'] is not defined.
        '''
        if options is None:
            options=o
        if options['logging'] is not None:
            return options['logging']+'/'+s
        else:
            return None

    def filter_catalog(singlecat,matchedcat,fitsimage,outname,auxcatname,options=None):
        ''' Filter the catalogues by radius, size and isolation.

        Sources are filtered out based on the following conditions:
            - Located further away than 3 degrees from the phase center.
            - Over 10 arcsec in size in the LOFAR catalogue.
            - Isolated.

        Args:
            singlecat (str): (long baseline) catalogue extracted from the fits file.
            matchedcat (str): catalogue containing the matched sources between the LOFAR image and catalogue that was matched against.
            fitsimage (str): the original images.
            outname (str): filename for the output diagnostic plot.
            auxcatname (str): the name of the catalogue that was matched against.
            options (dict): the options dictionary defined in quality_parset.py.
        Returns:
            None
        '''
        if options is None:
            options = o

        if options['restart'] and os.path.isfile(outname):
            warn('File ' + outname +' already exists, skipping source filtering step')
        else:

            matchedcat = Table.read(matchedcat)
            singlecat = Table.read(singlecat)

            fitsimage = fits.open(fitsimage)

            fieldra = fitsimage[0].header['CRVAL1']
            fielddec = fitsimage[0].header['CRVAL2']
            fitsimage.close()

            print 'Originally',len(matchedcat),'sources'
            matchedcat=filter_catalogue(matchedcat,fieldra,fielddec,3.0)

            print '%i sources after filtering for 3.0 deg from centre' % len(matchedcat)

            matchedcat=matchedcat[matchedcat['DC_Maj']<10.0]

            print '%i sources after filtering for sources over 10arcsec in LOFAR' % len(matchedcat)

            # not implemented yet!
            #tooextendedsources_aux = np.array(np.where(matchedcat[1].data[options['%s_match_majkey2'%auxcatname]] > options['%s_filtersize'%auxcatname])).flatten()
            #print '%s out of %s sources filtered out as over %sarcsec in %s'%(np.size(tooextendedsources_aux),len(allsources),options['%s_filtersize'%auxcatname],auxcatname)

            matchedcat=select_isolated_sources(matchedcat,30.0)
            print '%i sources after filtering for isolated sources in LOFAR' % len(matchedcat)

            matchedcat.write(outname)

    def sfind_image(catprefix,pbimage,sfind_pixel_fraction,options=None):
        ''' Find sources in the image using PyBDSF.
        Args:
            catprefix (str): prefix to use for output catalogues.
            pbimage (str): the image in which to find the sources.
            sfind_pixel_fraction (float): a value between 0.0 and 1.0 denoting the fraction of the image over which to look for sources.
            options (dict): the options dictionary defined in quality_parset.py.
        Returns:
            None
        '''
        if options is None:
            options = o
        f = fits.open(pbimage)
        imsizex = f[0].header['NAXIS1']
        imsizey = f[0].header['NAXIS2']
        f.close()
	old_dir = os.getcwd()
	filename = pbimage.split('/')[-1]
	dirname = '/'.join(pbimage.split('/')[:-1]) + '/'
	os.chdir(dirname)
	print('Filename: ', filename)
	print('Dirname: ', dirname)
        kwargs={}
        if o['sfind_pixel_fraction']<1.0:
            lowerx,upperx = int(((1.0-sfind_pixel_fraction)/2.0)*imsizex),int(((1.0-sfind_pixel_fraction)/2.0)*imsizex + sfind_pixel_fraction*imsizex)
            lowery,uppery = int(((1.0-sfind_pixel_fraction)/2.0)*imsizey),int(((1.0-sfind_pixel_fraction)/2.0)*imsizey + sfind_pixel_fraction*imsizey)
            kwargs['trim_box']=(lowerx,upperx,lowery,uppery)

        if options['restart'] and os.path.isfile(catprefix +'.cat.fits'):
            warn('File ' + catprefix +'.cat.fits already exists, skipping source finding step')
        else:
            img = bdsm.process_image(filename, detection_image='', thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=False,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None,**kwargs)
            img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True')
            img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
            img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
            img.export_image(outfile=catprefix +'.pybdsmmask.fits',img_type='island_mask',img_format='fits',clobber=True)
            img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True')
	    #os.chdir(old_dir)

    def crossmatch_image(lofarcat,auxcatname,options=None):
        ''' Cross match the LOFAR catalogue with an auxilliary catalogue.

        Args:
            lofarcat (str): the LOFAR catalogue.
            auxcatname (str): the auxilliary catalogue to match against.
            options (dict): the options dictionary defined in quality_parset.py.

        Returns:
            None
        '''
        if options is None:
            options = o
        auxcat = options[auxcatname]
        if options['restart'] and os.path.isfile(lofarcat + '_' + auxcatname + '_match.fits'):
            warn('File ' + lofarcat + '_' + auxcatname + '_match.fits already exists, skipping source matching step')
        else:
            t=Table.read(lofarcat)
            tab=Table.read(auxcat)
            match_catalogues(t,tab,o[auxcatname+'_matchrad'],auxcatname)
            t=t[~np.isnan(t[auxcatname+'_separation'])]
            t.write(lofarcat+'_'+auxcatname+'_match.fits')
    #----------------------------------------------------------------------------------------------
    #Actual Steps start from here:



    global o
    o=options(config_path,option_list)
    #ingest config data:
    o['']

    print(o)
    if o['pbimage'] is None:
        die('pbimage must be specified')

    # fix up the new list-type options - NO DONT!
    # for i,cat in enumerate(o['list']):
    #     try:
    #         o[cat]=o['filenames'][i]
    #     except:
    #         pass
    #     try:
    #         o[cat+'_matchrad']=o['radii'][i]
    #     except:
    #         pass
    #     try:
    #         o[cat+'_fluxfactor']=o['fluxfactor'][i]
    #     except:
    #         pass

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # pybdsm source finding
    sfind_image(o['catprefix'],o['pbimage'],o['sfind_pixel_fraction'])

    # matching with catalogs
    for cat in o['list']:
        print 'Doing catalogue',cat
        crossmatch_image(o['catprefix'] + '.cat.fits',cat)
        filter_catalog(o['catprefix'] + '.cat.fits',o['catprefix']+'.cat.fits_'+cat+'_match.fits',o['pbimage'],o['catprefix']+'.cat.fits_'+cat+'_match_filtered.fits',cat,options=o)

    # Filter catalogs (only keep isolated compact sources within 3deg of pointing centre)

    # Astrometric plots
    if 'FIRST' in o['list']:
        report('Plotting position offsets')
        plot_position_offset('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_positions.png'%o['catprefix'],'FIRST',options=o)

        t=Table.read(o['catprefix']+'.cat.fits_FIRST_match_filtered.fits')
        bsra=np.percentile(bootstrap(t['FIRST_dRA'],np.mean,10000),(16,84))
        bsdec=np.percentile(bootstrap(t['FIRST_dDEC'],np.mean,10000),(16,84))
        mdra=np.mean(t['FIRST_dRA'])
        mddec=np.mean(t['FIRST_dDEC'])
        print 'Mean delta RA is %.3f arcsec (1-sigma %.3f -- %.3f arcsec)' % (mdra,bsra[0],bsra[1])
        print 'Mean delta DEC is %.3f arcsec (1-sigma %.3f -- %.3f arcsec)' % (mddec,bsdec[0],bsdec[1])

        report('Plotting flux ratios')
        # Flux ratio plots (only compact sources)
        plot_flux_ratios('%s.cat.fits_FIRST_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_FIRST_match_filtered_fluxerrors.png'%o['catprefix'],options=o)

    report('Plotting flux scale comparison')
    # Flux scale comparison plots
    if 'TGSS' in o['list']:
        plot_flux_errors('%s.cat.fits_TGSS_match_filtered.fits'%o['catprefix'],o['pbimage'],'%s.cat.fits_TGSS_match_filtered_fluxratio.png'%o['catprefix'],'TGSS',options=o)
        t=Table.read(o['catprefix']+'.cat.fits_TGSS_match_filtered.fits')
        ratios=t['Total_flux']/(t['TGSS_Total_flux']/o['TGSS_fluxfactor'])
        bsratio=np.percentile(bootstrap(ratios,np.median,10000),(16,84))
        print 'Median LOFAR/TGSS ratio is %.3f (1-sigma %.3f -- %.3f)' % (np.median(ratios),bsratio[0],bsratio[1])
    if 'NVSS' in o['list']:
        t=Table.read(o['catprefix']+'.cat.fits_NVSS_match_filtered.fits')
        t=t[t['Total_flux']>10e-3]
        ratios=t['Total_flux']/t['NVSS_Total_flux']
        bsratio=np.percentile(bootstrap(ratios,np.median,10000),(16,84))
        print 'Median LOFAR/NVSS ratio is %.3f (1-sigma %.3f -- %.3f)' % (np.median(ratios),bsratio[0],bsratio[1])
    # Noise estimate
    hdu=fits.open(o['pbimage'])

    imagenoise = get_rms(hdu)
    print 'An estimate of the image noise is %.3f muJy/beam' % (imagenoise*1E6)

    # my-directory = '/data020/scratch/sean/measurement_sets/'
    #
    # try: # see if the inspection plot directory exists
    #     os.stat(my-directory + 'inspection_plots')
    # except: # create it if it does not exist
    #     os.mkdir(my-directory + 'inspection_plots')
    #
    # # move PNGs into a single directory called inspection_plots and delete any useless files at the end
    # os.rename(my-directory + '*.png', my-directory + 'inspection_plots')
    #
    # try: # delete the intermediate files
    #     os.remove(my-directory + '???')
    # except OSError: # but skip it if they do not exist
    #     pass

    return 0
