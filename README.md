# __SPA_ICIPdemosMatlabPackageV1.0__


This MATLAB package contains the implementations of __Spectral Pre-Adaptation (SPA)__ in the following ICIP publications:

>1. J. Portilla, *["Maximum likelihood extension for non-circulant deconvolution"](https://ieeexplore.ieee.org/document/7025868/)*, 2014 IEEE International Conference on  Image Processing (ICIP) , pp. 4276-4279, Oct 2014.
2. C. Dong, J. Portilla, *["Maximum likelihood interpolation for aliasing-aware image restoration"](https://ieeexplore.ieee.org/document/7532420/)*, 2016 IEEE International Conference on Image Processing (ICIP), pp. 564-568, Sept 2016.
3. C. Dong, J. Portilla, *["Spectral pre-adaptation for two-step arbitrary-shape-support image restoration"](https://ieeexplore.ieee.org/document/8296936/)*, 2017 IEEE International Conference on Image Processing (ICIP),pp. 3515-3519, Sept 2017.

__Version 1.0, August 2018.__
- - -
<br/>
## **Description**
**SPA** stands for **S**pectral **P**re-**A**daptation. Given  an  incomplete  (i.e.,  missing  pixels)  uniformly  blurred  image, it gives a likely guess of the unknown pixels by performing a maximum likelihood optimization according to a simple spectral model and based on the consistency of the pixels’ values with the blurring kernel and Gaussian noise. The obtained result is spectrally compatible with a noisy version of an image that has been filtered with the blurring kernel. Then, one can apply any standard deconvolution method to compensate for the noise and blur on the so-completed image. SPA assumes the blurring kernel and the  noise standard deviation are either both known or they have  been  reliably estimated.

It can be used in different problems when doing deconvolution, like:
- Boundary handling (by extending the boundaries and imposing circular boundaries condition) (see ICIP 2014<sup>[1]</sup>).
- Deblurring aliased observations, in which case the known subset is a uniformly sampled version of the complete set of pixels (see ICIP 2016<sup>[2]</sup>).
- Deblurring blurred objects having supports with arbitrary shapes (e.g., out-of-focus background in presence of in-focus foreground) (see ICIP 2017<sup>[3]</sup>).

<br/>

## **Files and Folders**
|Files and Folders            |Description |
|:------------:|:------------:|
|COPYING|GNU AFFERO GENERAL PUBLIC LICENSE Version3|
|core/				     | folder with core function script, SPA.m|
|data/				     | folder with images or data used in all scripts|
|demos/ICIP2014/ICIP2014.m	 |	script for all demos in ICIP 2014<sup>[1]</sup>|
|demos/ICIP2016/ICIP2016.m	 | script for all demos in ICIP 2016<sup>[2]</sup>|
|demos/ICIP2017/	 | folder with scripts for all demos in ICIP 2017<sup>[3]</sup>|
|functions/				 | folder with auxiliary MATLAB functions|
|L2rL0deblurMatlabToolbox2p1/	| folder with the [L2-r-L0deblur toolbox version 2.1](https://www.researchgate.net/publication/325903515_L2rL0deblurMatlabToolbox2p1)|
|README.pdf		   |	this file|


<br/>
## **Important Notes**
- To try the algorithm on
  - all ICIP2014 demos, run *demos/ICIP2014/ICIP2014.m* with proper changes of set up parameters. For details of parameter setting please see the script's comments.
  - all ICIP2016 demos, run *demos/ICIP2016/ICIP2016.m* with proper changes of set up parameters. For details of parameter setting please see the script's comments.
  - demo in Fig.1. ICIP2017, block in-painting with Lena image, run *demos/ICIP2017/ICIP2017_BlockInpainting.m*.
  - demo in Fig.2. ICIP2017, a fore-ground in focus background blurred realistic simulation, run *demos/ICIP2017/ICIP2017_Pirate.m*.
  - demo in Fig.3. ICIP2017, a fore-ground in focus background blurred real high quality photography, run *demos/ICIP2017/ICIP2017_Comb.m*.
  - demo in Fig.4. ICIP2017, a fore-ground in focus background blurred real JPEG-compressed medium quality photography, run *demos/ICIP2017/ICIP2017_Pen.m*.


- Requires the L2-r-L0 Deblur Matlab Toolbox 2.1:
 https://www.researchgate.net/publication/325903515_L2rL0deblurMatlabToolbox2p1
It is included in this ZIP file.  It  is  an  implementation  of  the  method  described  in:
 > J. Portilla, A. Tristan-Vega, I.W. Selesnick, [*"Efficient and robust
   image restoration using multiple-feature L2-relaxed sparse analysis
   priors"*](https://ieeexplore.ieee.org/document/7265041/), IEEE Transactions on Image Processing, vol. 24, no. 12,
   pp. 5046-5059, Dec 2015.

- There is a very good general agreement of the results provided by this code package and the reported results in the corresponding ICIP publications. <br/>Small differences with reported ISNR values are mainly due to here we are using a newer version of the core SPA.m function (*previously referred to
as "**MLE**" or "**MLI**", for Maximum-Likelihood Extension or Interpolation*), which, for uniformity, we enforced here to be common for all reproduced experiments in all three ICIP publications.

  Additionally, a certain variability can be caused due to HW&SW differences:
  - We have observed some differences running in different computers with diverse operating systems and/or MATLAB versions.
  - Some SPA experiments may be sensitive to these differences (like when dealing with heavy blur and low noise).

<br/>
## **Copyright and License**

**Copyright (c) 2018**
**Javier Portilla** <javier.portilla@csic.es\>
**Chaoqun Dong** <cdongae@connect.ust.hk\>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

<br/>
##### **Thank you for your interest in this code. Comments (on bugs, compatibility issues, experimental extensions of the original code, etc.) are all very welcome. Enjoy it!**
