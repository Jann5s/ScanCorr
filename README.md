# ScanCorr

## About
ScanCorr is a tool to correct image scanning artifacts (typically from SEM) by fusing two images scanned at orthogonal directions into one corrected image. The method uses a specialized DIC method to correlate the scan-lines between the images, allowing each scan-line to move freely but not to deform or rotate. More details about the method are inside the <scan_corr_GUI_help.txt>.


## Organisation
The main function that does all the work is scan_corr.m, which is designed to be called programmatically, for instance to handle large batches of code. However, most users will interact with scan_corr_GUI.m which as the name suggests is a Graphical User Interface. This GUI includes a help, which is exactly <scan_corr_GUI_help.txt>.

This type of algorithm requires a lot of iterations of many pixels, so to keep computation times reasonable, two mex algorithms are used; <scan_corr_mex.cpp> and <scan_corr_im_mex.cpp>. The first handles the DIC part of scan_corr, the second handles the image interpolation, which is based on the excellent B-Spline interpolation library made by Thevenaz [1].

[1] P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited," IEEE Transactions on Medical Imaging, vol. 19, no. 7, pp. 739-758, July 2000.


## Installation
For installation, not much is required, except that the two mex files need to be compiled for your system. To do so, <scan_corr_mex_compiles.m> is provided which does exactly this. Naturally, your system needs to be configured for mex file compilation, for which I refer the user to the Matlab documentation on the matter.

