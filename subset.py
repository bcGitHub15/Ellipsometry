#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 10:52:00 2026

subset.py

Down sample a DepositData4 source to match resolution of a pattern.

So long as two data sets share the same dx=dy and the same origin, 
which all so far do, it should be possible to extract the low resolution
subset from a high resolution map. Indeed it should not even need fancy 
searching. If the ratio of the sizes is an integer then we should be 
able to simply index through the arrays at different rates, emitting all 
points that were valid in the larger source array. This does not need 
to be a class method but will need to operate on the (public) instance variables.

We will need two source data sets, the larger one will be the data source and the smaller will serve as the pattern.

Because the arrays are stored linearly and not in a square I am going to work by
translating the linear index in the pattern into an x,y pair and then translating
that to a linear index in the source.

It turns out that there are points on the edge of the image that are in
the low res scan but not in the high res. They have z values of zero and
so are easily detected and removed from the subset map.

@author: bcollett
"""
import numpy as np
from DepositData4 import DepositData4

def subset(src, pattern, debug=False):
    if src.zdata is None:
        raise ValueError('No source zdata.')
    src_dx = src.xdata[1] - src.xdata[0]
    p_dx = pattern.xdata[1] - pattern.xdata[0]
    dratio = p_dx / src_dx
    iratio = int(dratio + 0.5)
    print(dratio, iratio)
    print(dratio - iratio, np.abs(dratio - iratio))
    if np.abs(dratio - iratio) > 1e-3:
        raise ValueError(f'Pattern step must be an integral multiple of src step, not {dratio}')
    #
    # Make sure that the sizes match
    #
    print(src.xdata.shape)
    src_npoint = src.xdata.shape[0]
    src_nrow = int(np.sqrt(src_npoint))
    print(pattern.xdata.shape)
    p_npoint = pattern.xdata.shape[0]
    p_nrow = int(np.sqrt(p_npoint))
    if (src_nrow - 1) != iratio * (p_nrow - 1):
        raise ValueError(f'Src size ({src_nrow}) must be an integral multiple of pattern size ({p_nrow})')
    # Figure out if there is an offset
    p_x0 = pattern.xdata[0]
    s_x0 = src.xdata[0]
    xoffset = p_x0 - s_x0
    print(xoffset)
    yoffset = pattern.ydata[0] - src.ydata[0]
    p_y0 = pattern.ydata[0]
    s_y0 = src.ydata[0]
    yoffset = p_y0 - s_y0
    print(yoffset)
    # build result same size as pattern
    result = DepositData4(pattern.xdata, pattern.ydata, mask=pattern.mask, t1=None, sort=False)
    # Work through arrays copying data, for moment x, y, z
    eps = 1.0e-4
    print(f'Pattern starts at ({pattern.xdata[0], pattern.ydata[0]}')
    print(f'Source starts at ({src.xdata[0], src.ydata[0]}')
    nskipped = 0
    for i in range(p_npoint):
        # Find x,y
#        p_row_no = int(i / p_nrow)
#        p_col_no = i - row_no * p_nrow
#        x = p_x0 + p_col_no * p_dx
#        y = p_y0 + p_row_no * p_dx   # Assumes square!!
        x = pattern.xdata[i]
        y = pattern.ydata[i]
        # Convert to index in src
#        s_row_no = src_nrow - 1 - int((s_y0 - y) / src_dx + eps)
        s_row_no = int((y - s_y0) / src_dx + eps)
#        s_col_no = src_nrow - 1 - int((s_x0 - x) / src_dx + eps) # Assumes square
        s_col_no = int((x - s_x0) / src_dx + eps) # Assumes square
        s_index = s_row_no * src_nrow + s_col_no
        if debug:
            print(f'At index {i} loc {x},{y} srow={s_row_no}, scol={s_col_no}')
        if s_index >= src_npoint:
            print(f'At pattern index {i} s_index={s_index} exceeds {src_npoint-1}')
            break
        if np.abs(x - src.xdata[s_index]) > eps:
            print(f'x error at i = {i}, x={x}, index={s_index}, xdata={src.xdata[s_index]}')
            continue
        if np.abs(y - src.ydata[s_index]) > eps:
            print(f'y error at i = {i}, y={y}, index={s_index}, ydata={src.ydata[s_index]}')
            continue
        # Point is in right place, if it is valid then copy it.
        result.mask[i] = pattern.mask[i]
        if pattern.mask[i]:
            zdat = src.zdata[s_index]
            if zdat > 0:
                result.zdata[i] = zdat
            else:
                result.mask[i] = False
        else:
            nskipped += 1
#            print(f'Skipping invalid point at result index {i}.')
    print(f'Subset complete. Skipped {nskipped} points')
    new_name = f'Subset of {src.name}'
    result.set_name(new_name)
    return result

if __name__ == '__main__':
    g1mm = DepositData4.from_file('Run9_17_WA0904_1mm.txt', 'Glass 1mm')
    g1mm.Ron_plot()
    S4 = DepositData4.from_file('Run6_26_Scan4.txt', 'Scan4')
    g2mm = subset(g1mm, S4)
    g2mm.Ron_plot()