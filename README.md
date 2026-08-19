# Ellipsometry
Analyzing ellipsometry from LiF on Si deposits.

Ellipsometry8_3_26.ipynb contains the preliminary analysis of the Woollam data that arrived on August 2nd. It deals with a fused silica blank about 600 nm of LiF that has been scanned three times with a spacing of 2 mm. Scans2 and 4 are with the wafer flat at the 12 o'clock position and Scan3 with it at the 3 o'clock position. The data are in files *Run_6_26_Scan<n>.png* where <n> is the scan number. Scan number 1 was just a measurement of the refractive index of the glass substrate so I did not include it.

The conclusions are:
 - The ellipsometry is very reproducible.
 - The odd far off points do not reproduce.
 - This deposit has the form of a rather sharp ridge (more quartic than quadratic) with a total height of about 1.3%.
 - The difference between two supposedly identical scans is normally distributed with a FWHM of about 0.5 nm.
 - The data from the rotated position match those from the unrotated positions to almost the same precision as the two unrotated positions.

Ellipsometry1.ipynb starts with a summary of our visual inspection of the data in Run_11_15
and some notes on the 7/1/26 Zoom call with Matthew Barber of Blue Ridge and Ron Synowicki of Woollam.

It then contains my ramblings as I developed some analysis code. About half way down the file starting
at the heading SN002 is my first round analysis of each of the four deposits including fits to simple
angled planes, residual plots, and histograms of height noise. The conclusions are:
 - The spot noise on all the deposits accounts for about half of the original height spread
 - All four deposits show a clear linear gradient.
 - With the gradient removed all four show signs of drooping round the edges, though the patterns are
   different in each case.
 - With the gradient removed the central regions are flat with $\pm2\,nm$ of quite normally distributed
   noise plus the outlying spots mentioned above.

