clear; close all;

mex -setup c++

% debug = '-g';
debug = '';

cppfile = 'scan_corr_mex.cpp';
mex(cppfile,debug,'-largeArrayDims','-DNOMINMAX')


cppfile = 'scan_corr_im_mex.cpp';
mex(cppfile,debug,'-largeArrayDims','-DNOMINMAX')

