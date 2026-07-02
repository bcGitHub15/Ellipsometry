# Ellipsometry
Analyzing ellipsometry from LiF on Si deposits.

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

