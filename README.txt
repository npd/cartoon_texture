% Cartoon+texture decomposition by the TV-L1 model

# ABOUT

* Author    : Vincent Le Guen <vincent.leguen@telecom-paristech.org>
* Copyright : (C) 2013 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides an implementation of the TV-L1
cartoon+texture decomposition, as described in IPOL.

This program reads and writes PNG images, but can be easily
adapted to any other file format.

Only 8bit RGB PNG images are handled. Other TIFF files are implicitly
converted to 8bit color RGB.

# REQUIREMENTS

The code is written in ANSI C and C++, and should compile on any
system with an ANSI C/C++ compiler.

The libpng header and libraries are required on the system for
compilation and execution. The program is parallelized with OPENMP.

# COMPILATION

Simply use the provided makefile, with the command `make`.


# USAGE

`cartoonTexture` takes 4 parameters: `cartoonTexture in.png lambda cartoon.png texture.png`
* `lambda`     : the regularization parameter
* `in.png`   : input image
* `cartoon.png`  : output image without the textures
* `texture.png` : textures (out - in, bounded to [-20, 20])

# ABOUT THIS FILE

Copyright 2013 IPOL Image Processing On Line http://www.ipol.im/
Author: Vincent Le Guen <vincent.leguen@telecom-paristech.org>

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.

