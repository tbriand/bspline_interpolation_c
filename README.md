# Theory and Practice of Image B-Spline Interpolation #

## Summary ##
This program implements the homographic transformation of an image using
B-spline interpolation. It is part of an [IPOL publication](
https://doi.org/10.5201/ipol.2018.221)

## Authors ##

* Thibaud Briand <thibaud.briand@enpc.fr>
* Pascal Monasse <monasse@imagine.enpc.fr>

Laboratoire d'Informatique Gaspard Monge (LIGM)/
Ecole des Ponts ParisTech

## Version ##
Version 1.0, released on 06/27/2018

## License ##
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

You should have received a copy of the GNU General Pulic License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2017-2018, Thibaud Briand <thibaud.briand@enpc.fr>,
                         Pascal Monasse <monasse@imagine.enpc.fr>

All rights reserved.

## Build ##
Required environment: Any unix-like system with a standard compilation
environment (make and C compiler) and [Cmake](https://cmake.org/).
Optional libraries:
[libpng](http://libpng.org/pub/png/libpng.html),
[lipjpeg](http://ijg.org/),
[libtiff](http://simplesystems.org/libtiff/),
[GSL](https://www.gnu.org/software/gsl/)

Build instructions:

    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_BUILD_TYPE=Release ../src
    $ make

It produces library "Splinter" and programs "bspline" and "compute_order".

## Usage ##
The program reads an  homography, an input image, takes some parameters and
produces an homographic transformation of the image using B-spline
interpolation. The meaning of the parameters is thoroughly discussed on the
accompanying IPOL article.

      Usage: ./bspline "homography" in out [order boundary eps larger]
    Homographic transformation of an image using B-spline interpolation

    homography: 9 matrix coefficients ("h11 h12 h13; h21 h22 h23; h31 h32 h33")
    in       : filename of the input image
    out      : filename of the output image
      Options:
      order    : order of interpolation (integer between 0 and 16, default 11)
      boundary : extension (constant, periodic, hsymmetric*, wsymmetric)
      eps      : relative precision (float, default 6) (eps>=1 means 10^-eps)
      larger   : compute on exact (0*) or larger domain (1)
        *default parameters

Execution examples:

  1. Default parameters for an horizontal shift of 0.5 pixel:

      $ bspline "1 0 0.5; 0 1 0; 0 0 1" input.png output.tiff

  2. Using order 3, constant boundary extension and relative precision of 10^-5:

      $ bspline "1 0 0.5; 0 1 0; 0 0 1" input.png output.tiff 3 constant 5

*Remark*: the semicolon row separators in the homography matrix are optional.

*Remark*: the boundary parameter can be abridged (e.g. c=const=constant)

### Test ###
    $ ./bspline "1 0 0.5; 0 1 0; 0 0 1" ../data/lenna.pgm out.pfm
    $ diff -s ../data/lenna_x+0.5.pfm out.pfm

The pfm output format is not standard and not readable by most software.
Displayable output formats include pgm, ppm, *png*, *jpeg*, *tiff*
(*require optional library support*)

### Library usage ###
The core part of the bspline interpolation is included in library "Splinter".
If you want to use bspline interpolation in your own programs, this is the only
part you need. The usage is fairly simple (see function `splinter_homography`):

    #include "splinter.h"
    splinter_plan_t plan = splinter_plan(...); // Prefiltering
    double pixOut[3];                          // Output values (#channels)
    splinter(pixOut, 1.3, 2.4, plan);         // Interpolate at coords (1.3,2.4)
    splinter_destroy_plan(plan);               // Free reserved memory

### Generating HTML documentation ###
    $ cd src
    $ doxygen Doxyfile

Open html/index.html to browse the documentation.

## List of files in the directory src##

* bspline_main.c         : Main program for input/output
* homography_tools.[hc]  : Functions related to homographies
* splinter_transform.[hc]: Compute homographic transformation of image
* bspline.[hc]           : Compute B-spline parameters and kernel (library)
* splinter.[hc]          : Prefilter and indirect B-spline transform (library)

Additional files:

* iio/                   : C library for opening images in any format
* xmtime.h               : Clock with millisecond precision
* compute_bspline.c      : Ccompute the B-spline interpolator parameters
* hom4p.c                : Compute homography from 4 points (for on-line demo)
