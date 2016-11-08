# -----------------------------------------------------------------------------
# Get mean FWHM, SKY MEAN and SKY STANDARD DEVIATION value from several stars
# -----------------------------------------------------------------------------


procedure getdata (image, datamax, ellipticity)

file image {prompt = "Input image name (include extension) or \'*.fits\'"}
real datamax {prompt = "Maximum good data value?"}
real ellipticity {prompt = "Maximum accepted ellipticity?"}
struct gain_key {mode="q", prompt = "GAIN keyword?"}
struct rdnoise_key {mode="q", prompt = "RDNOISE keyword?"}
struct *file_var {mode="h", prompt = "Internal file name variable"}
struct *file_var2 {mode="h", prompt = "Internal file name variable"}
struct *file_var3 {mode="h", prompt = "Internal file name variable"}


begin
  string imname, section, old_imname, g_key, rd_key
  bool check, data_daofind, datapar_search, section_bool, trimming, fitrad_error
  real dmin, dmax, ellip_max, fitrad, sigma, gain, rdnoise, sstd, smean, thresh
  real fannulus, fdannulu, sig, smean2, sstd2
  real diff, diff2, diff3, oldfitrad, oldsmean, oldsstd, oldthresh
  real var1[500], var2[500], var3[500], var4[500], var5[500], var6[500]
  real fitradvar[11], smeanvar[11], sstdvar[11]
  real xcen1, ycen1, xcen2, ycen2, fwh2, ellip, dx, dy, fwhm_sum
  real msky, stdev, msky_aver, stdev_aver, psfmeasure_fwhm
  int iters, max_iters, i, j, k, m, e_r
  struct line, line2

  imname = image
  dmax = datamax
  g_key = gain_key
  rd_key = rdnoise_key
  ellip_max = ellipticity

# Check package
if (! defpac ("daophot")) {
    print ('\n This script must be loaded inside noao/digiphot/daophot')
    bye()
}
else { # Do nothing
}      

# Create 'getdatalist' file
delete ('getdatalist',verify=no,>>&"/dev/null")
if (imname == '*.fits') {
    files ('*.fits', >> 'getdatalist')
}
else {
    print (imname, >> 'getdatalist')
}

# Ask if file(s) should be trimmed
trimming = no
section_bool = no
check = no
print ('\nSearch only a SECTION of the frame(s)? (y/n)')
scan (check)
if (check) {
    section_bool = yes
    print ('\n Input section where stars will be searched for')
    print (' with the format: [x1:x2,y1:y2]')
    scan (section)
}

# This bool variable will indicate later on if the FWHM value was set manually
# to 3 because no FWHM could be calculated
fitrad_error = no
file_var3 = ('getdatalist')

# Getdata file while
while (fscan (file_var3,line2) != EOF) {
    k = strlen(line2)
    if (substr (line2, k-4, k) == ".fits") {
        imname = substr (line2, 1, k-5)
    }
    else {
        print ('\nDid you input the name of the frame WITHOUT the')
        print (' \'fits\' extension? (y/n)')
        scan (check)
        if (check) {
            imname = substr (line2, 1, k)
        }
        else {
            print (' FILENAME ERROR')
            bye()
        }
    }
    
    print ('\nImage: '//imname)  

    hselect.mode = "hl"
    hselect.images = imname
    hselect.fields = g_key
    hselect.expr = yes
    hselect > "tempget"
    file_var = 'tempget'
    while (fscan (file_var,gain) != EOF)
    del ("tempget")
    print ('GAIN = '//gain)

    hselect.mode = "hl"
    hselect.images = imname
    hselect.fields = rd_key
    hselect.expr = yes
    hselect > "tempget"
    file_var = 'tempget'
    while (fscan (file_var,rdnoise) != EOF)
    del ("tempget")
    print ('RDNOISE = '//rdnoise)

    fitrad = 3 # Initial FWHM value
    smean = 1. # Initial SKY MEAN value
    sstd = 1.  # Initial SKY STANDARD DEVIATION value

    oldfitrad = 3.
    oldsmean = 1.
    oldsstd = 1.
    datapar_search = yes
    
    old_imname = imname
    if (section_bool) {
        print (' Trimming...')
        imcopy.input = imname//section
        imcopy.output = imname//'_trim.fits'
        imcopy.verbose = yes
        imcopy.mode = "hl"
        imcopy
        imname = imname//'_trim'
        trimming = yes
    }
    
    max_iters = 5
    iters=1
    sig=1
    # Datapar search while
    while (datapar_search==yes) {
        print ('\n\n***************************************************')
        print ('          Number of iterations left: '//(max_iters-iters))
        print ('***************************************************')
        if ((smean*gain + rdnoise*rdnoise) <= 0. || (smean <= 0.)) {
            print ('\n POSSIBLY A \'LAS CAMPANAS\' OBSERVATORY FRAME')
            print ('       (if not, check results carefully)')
            if (sig==1) {
                sigma = sqrt(1*gain + rdnoise*rdnoise)/gain
                # This structure is meant for 'Las Campanas' observatory frames
                sigma = sigma*5
                sig = 2
            }
        }
        else {
            sigma = sqrt(smean*gain + rdnoise*rdnoise)/gain
        }

        print ("\nDatapars values (fwhmpsf, sigma, datamin, datamax): ")
        dmin = smean-3*sigma
        print(fitrad, sstd, dmin, dmax)
        datapars.fwhmpsf = fitrad
        datapars.sigma = sstd
        datapars.datamin = smean-3*sigma
        datapars.datamax = dmax

        i = 5
        data_daofind=no

        # Daofind while
        while (data_daofind==no && i >= 1) {

            # Use a high threshold so only the brighter stars will be found
            findpars.threshold = i*3.5*sigma
            print ('\n---------------------------------------------------')
            print (' Daofind task (i): '// i)
            print (" Threshold (i*3.5*sigma): "// i*3.5*sigma)
            print ('---------------------------------------------------')

            daofind.verif = no
            daofind.verb = yes
            daofind.interactive = no
            daofind.verbose = no
            daofind.mode = 'hl'
            daofind ((imname), (imname//'.coo.psf.1'))
            # display ((imname), 1)
            # tvmark.interactive = no
            # tvmark.outimage = ""
            # tvmark.mark = 'circle'
            # tvmark.font = "raster"
            # tvmark.txsize = 2
            # tvmark.radii = 10
            # tvmark.color = 204
            # tvmark.number = yes
            # tvmark.label = no
            # tvmark (1, (imname//'.coo.psf.1'))                      

            file_var = (imname//'.coo.psf.1')
            m=0
            while (fscan (file_var,line) != EOF) {
                m = m + 1
            }
            # This is the number of stars found by 'Daofind', the first 41 lines are format text.
            m = m - 41
            print ("\nNumber of stars found by \'Daofind\': "//m)

            # If number of stars found is <= 5 run 'Daofind' again
            # with a lower threshold
            i = i-1
            if (m <= 5) {
                if (i <=0) {
                    print ('\nNot enough stars found using minimum threshold value (<5).')
                    print ('Halting')
                    print ('Not enough stars found using minimum threshold value (<5)', >> imname//'_iter')
                    delete ((imname//'.coo.psf.1'))
                    bye()
                }
                print ("Not enough stars found (<5).")
                print ("Reducing threshold and peforming new \'daofind\'.")
                delete ((imname//'.coo.psf.1'))
            }
            else {
                data_daofind = yes
            }
            
            if (data_daofind == yes) {
                print ('\n---------------------------------------------------')
                print (' Phot task')
                print ('---------------------------------------------------')
                unlearn centerpars
                unlearn fitskypars
                unlearn photpars
                unlearn psf
                # From Massey-Davis guide to stellar CCD photometry
                fitskypars.salgorithm = "mode"
                centerpars.calgorithm = "none" 
                
                # According to 'A Reference Guide to the IRAF-DAOPHOT Package'
                # by L. Davis (page 31): cbox = 2xFWHM (or 5, wichever is greater)
                #                        annulus = 4xFWHM
                #                        dannulu = 2.5-4.0xFWHM
               
                # According to IRAF help: a reasonable value for 'cbox' is:
                # 2.5-4.0 * FWHM
                
                # According to "A User's Guide to Stellar CCD Photometry with IRAF"
                # by Massey-Davis (page 47): cbox = 5 (approx 2.0-3.0xFWHM)
                #                            annulus = 10 (approx 3.0-4.0xFWHM)
                #                            dannulu = 10 (approx 3.0-4.0xFWHM)
                 
                centerpars.cbox = 2.5*fitrad
                fannulus = 4*fitrad
                fdannulu = 3.25*fitrad      
                fitskypars.annulus = fannulus
                fitskypars.annulus = fdannulu              
                
                phot.interactive = no
                phot.radplots = no
                phot.update = yes
                phot.verbose = yes
                phot.verify = no
                phot.verbose = no
                phot.mode = 'hl'
                photpars.apertures = fitrad
                phot ((imname), (imname//'.coo.psf.1'), (imname//'.mag.psf.1'))
                
                txdump.mode = 'hl' 
                txdump.textfile = (imname//'.mag.psf.1')
                txdump.headers = no
                txdump.fields = 'MSKY,STDEV'
                txdump.expr = 'MAG[1]!=INDEF'
                txdump > auxiliar

                file_var = ('auxiliar')
                m=0
                while (fscan (file_var,line) != EOF) {
                    m = m + 1
                }
                print ("\nNumber of stars found (with MAG != INDEF) = "//m)
                print ("Number of stars found (with MAG != INDEF) = "//m, >> imname//'_iter')

                if (m<=5) {
                    if (i <=0) {
                        print ('\nNot enough stars found using minimum threshold value (<5).')
                        print ('Halting')
                        print ('Not enough stars found using minimum threshold value (<5)', >> imname//'_iter')    
                        delete ((imname//'.coo.psf.1'))
                        delete ((imname//'.mag.psf.1'))
                        delete ('auxiliar')
                        delete ('getdatalist')
                        bye()
                    }
                    else {
                        print ('\nNot enough stars found with MAG != INDEF (<5)')
                        print ('The threshold value must be too big (check code)')
                        print ('or the image too saturated (or there\'s something')
                        print ('wrong with the frame)')
                        print ('\nReducing threshold and peforming new \'daofind\'.')
                        data_daofind = no
                        delete ('auxiliar')
                        delete ((imname//'.coo.psf.1'))
                        delete ((imname//'.mag.psf.1'))
                    }
                }
            } # This bracket closes the 'data_daofind' if
        } # This bracket closes the 'data_daofind' while

        # Obtaining SKY MEAN value and SKY STANDARD DEVIATION value
        file_var = "auxiliar"
        msky_aver = 0
        stdev_aver = 0
        m=0
        while (fscan (file_var,msky,stdev) != EOF) {
            msky_aver = msky_aver + msky
            stdev_aver = stdev_aver + stdev
            m = m + 1
        }
        
        smean = msky_aver/m # Final SMEAN value
        sstd = stdev_aver/m # Final SSTD value
        smeanvar[iters] = smean
        sstdvar[iters] = sstd
        
        print ("SMEAN = "//smean)
        print ("STDEV = "//sstd)
        diff2 = smean -oldsmean
        diff3 = sstd - oldsstd
        delete ('auxiliar')
        
        # Obtaining FWHM value
        print ('\n---------------------------------------------------')
        print (' Obtaining average FWHM value')
        print ('---------------------------------------------------')
        
        txdump.mode = 'hl' 
        txdump.textfile = (imname//'.mag.psf.1')
        txdump.headers = yes
        txdump.fields = 'xcenter, ycenter, mag'
        # Cleans the stars with INDEF MAG values
        txdump.expr = 'MAG[1]!=INDEF'
        txdump > auxiliar
        # Sorts stars in descending MAG order
        txsort.ascend = yes
        txsort ('auxiliar', 'MAG')

        txdump.mode = 'hl'
        txdump.textfile = ('auxiliar')
        txdump.headers = no
        txdump.fields = 'xcenter, ycenter'
        # Removes the 'MAG' column
        txdump.expr = 'yes'
        txdump > ('auxiliar2')

        file_var = "auxiliar2"          
        m=0
        while (fscan (file_var,line) != EOF) {
            # Counts number of stars in last file
            m = m + 1
        }

        # If number of stars is >= 100 then keep the first 100
        if (m >= 100.) {
            fields.files = ('auxiliar2')
            fields.fields = "1,2"
            fields.lines = "1-100"
            fields.mode = "h"
            fields > ('auxiliar3')
            print ("\nUse the brightest 100 stars.")
        }
        else {
            print ('')
            cp ('auxiliar2', 'auxiliar3')
        }
             
        print ('q', >> 'cursor.txt')
        print ('Running \'psfmeasure\' task...')


        # This task performs the calculation of FWHM values for multiple stars
        noao
        obsutil
        psfmeasure(coords="markall", wcs="logical", display=no, frame=1,
        level=0.5, size="FWHM", beta=INDEF, scale=1., radius=15, sbuffer=5,
        swidth=5, saturation=62000, ignore_sat=yes, iterations=5, xcenter=INDEF,
        ycenter=INDEF, logfile="", graphcur="cursor.txt",
        images=(imname//".fits"), imagecur="auxiliar3", > "outputpsf")
        # print ("check")
        # scan (check)
        file_var = "outputpsf"
        m=0
        while (fscan (file_var,line) != EOF) {
            m = m + 1
        }

        # fields.mode = "h"
        # fields.files = "auxiliar3"
        # fields.fields = "1,2"
        # fields.lines = "2-"
        # fields >> "output3"

        fields.mode = "h"
        fields.files = "outputpsf"
        fields.fields = "1,2,4,5"
        fields.lines = "5-"//(m-2)
        fields >> "output2"

        fields.mode = "h"
        fields.files = "outputpsf"
        fields.fields = "9"
        fields.lines = m
        fields >> "temp_psfmeasure"
        file_var = "temp_psfmeasure"
        while (fscan (file_var,psfmeasure_fwhm) != EOF)
        del ("temp_psfmeasure")

        # file_var = "output3"
        # m=1
        # while (fscan (file_var, xcen1, ycen1) != EOF) {
        #     var1[m] = xcen1  # XCENTER found by 'phot'
        #     var2[m] = ycen1 # YCENTER found by 'phot'
        #     m = m + 1
        # }     

        file_var = "output2"
        m=1
        while (fscan (file_var, xcen2, ycen2, fwh2, ellip) != EOF) {
            # var3[m] = xcen2 # XCENTER found by 'psfmeasure'
            # var4[m] = ycen2 # YCENTER found by 'psfmeasure'
            var5[m] = fwh2   # FWHM
            var6[m] = ellip  # Ellipticity
            # print (' X,Y,FWHM = '//xcen2, ycen2, fwh2)
            m = m + 1
        }  

        # Rejecting stars with large ellipticity.
        # dx = 0.
        # dy = 0.
        fwhm_sum = 0.
        j = 0
        e_r = 0
        for (k=1; k<=m; k=k+1) {
            # OLD v
            # dx = abs(var1[k] - var3[k])
            # dy = abs(var2[k] - var4[k])
            # Condition to keep stars FWHM value (found by psfmeasure):
            # must be centered within 3 pixels of the center found by 'phot'
            # if ((dx < 3.) && (dy < 3.) && (var6[k] < 0.2)) {
            # OLD ^
            if (var6[k] < ellip_max) {
                fwhm_sum = fwhm_sum + var5[k]
                j = j+1
            }
            else {
                e_r = e_r+1
            }
        }
        print ("Accepted stars: "// j)
        print ("Rejected stars with large ellipticities: "// e_r)

        if (j == 0) {
            fitrad = 3
            fitradvar[iters] = fitrad
            print ('\nNo FWHM could be calculated during this iteration')
            print ('FWHM value set to 3 to avoid division by zero.')
            fitrad_error = yes
        }
        else {
            fitrad = fwhm_sum/j # Final FWHM value
            fitradvar[iters] = fitrad
            fitrad_error = no
        }
                  
        print ('\nFWHM (averaged) = '//fitrad)
        print ('FWHM (psfmeasure) = '//psfmeasure_fwhm)
        diff = fitrad - oldfitrad

        # End of iteration condition
    
        if (smean<0) {
            smean2=-smean
        }
        else {
            smean2=smean
        }
        # This structure accounts for the fact that these values may be negative
        # and so the 'End of iteration condition' below will fail unless they
        # are transformed into positive
        if (sstd<0) {
            sstd2=-ssstd
        }
        else {
            sstd2=sstd
        }

        # End of iteration condition: FWHM, SKY MEAN and STDEV values must
        # ALL have a difference of less than abs(10%) with the previous
        # calculated value; OR after max_iters or more iterations have been executed.
        if ((diff >= -fitrad/10.) && (diff <= fitrad/10.) && (diff2 >= -smean2/10.) && (diff2 <= smean2/10.) && (diff3 >= -sstd2/10.) && (diff3 <= sstd2/10.) || (iters >= max_iters)) {
            datapar_search = no
            # If the script reached the maximum number of iterations,
            # average ALL the values and present this average as the
            # final value.
            if (iters>=max_iters) {
                fitrad = 0.
                smean = 0.
                sstd = 0.
                m = 1
                while (m <=max_iters) {
                    fitrad = fitradvar[m] + fitrad
                    smean = smeanvar[m] + smean
                    sstd = sstdvar[m] + sstd
                    m = m + 1
                }
                fitrad = fitrad/max_iters
                smean = smean/max_iters
                sstd = sstd/max_iters
                print ('\n Maximum number of iterations achieved, using average values')
                print ('Maximum number of iterations achieved, using average values', >> imname//'_iter')
                print ('') 
                if (fitrad_error==yes) {
                    print ('FWHM set to 3 due to error.')
                }    
                else {
                    print ('FWHM = '//fitrad)
                }
                print ('Sky Mean = '//smean)
                print ('STDDEV = '//sstd)
                print ('')                      
            }
            else {
                if (fitrad_error==yes) {
                    print ('\nFWHM set to 3 due to error.')
                }    
                else {
                    print ('\nFWHM = '//fitrad)
                }
                print ('Sky Mean = '//smean)
                print ('STDDEV = '//sstd)
                print ('')
            }
        }
              
        delete ('auxiliar')
        delete ('auxiliar2')
        delete ('auxiliar3')
        delete ((imname//'.coo.psf.1'))
        delete ((imname//'.mag.psf.1'))
        delete ('cursor.txt')
        delete ('outputpsf')
        delete ('output2') 
        # delete ('output3')

        if (fitrad_error==yes) {
            print ('FWHM (manually set to 3 due to error) = '//fitrad, >> imname//'_iter')
        }    
        else {
            print ('FWHM = '//fitrad, >> imname//'_iter')
        }
        print ('STDDEV = '//sstd, >> imname//'_iter')  
        print ('Sky Mean = '//smean, >> imname//'_iter')          
              
        oldfitrad = fitrad
        oldsmean = smean
        oldsstd = sstd
        iters = iters+1
              
        fitrad_error=no # This is to reset this error warning

    } # This bracket closes the 'datapar_search' 'while'
          
    print (fitrad//'  FWHM', >> old_imname//'_data')
    print (sstd//'  STDDEV' , >> old_imname//'_data')  
    print (smean//'  Sky Mean', >> old_imname//'_data')
          
    if (trimming == yes) {
        delete (imname//'.fits')
    }              

} # This bracket closes the 'while' that goes through the 'getdatalist' file
      
delete ('getdatalist')
print ('\n SCRIPT FINISHED SUCCESFULLY') 

end
