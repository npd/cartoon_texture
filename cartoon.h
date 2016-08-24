/*
 * Copyright (c) 2013, Vincent Le Guen <vincent.leguen@telecom-paristech.org>
 * All rights reserved.
 *  
 *
 * 
 * License:
 *
 * This program is provided for scientific and educational only:
 * you can use and/or modify it for these purposes, but you are
 * not allowed to redistribute this work or derivative works in
 * source or executable form. A license must be obtained from the
 * patent right holders for any other use.
 *
 *
 */


/**
 * @file   cartoon.cpp
 * @brief  Cartoon + texture functions
 *
 * 
 *
 * @author Vincent Le Guen <vincent.leguen@telecom-paristech.org>
 */

#ifndef CARTOON_H
#define CARTOON_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <iostream>
#include <ctime>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

float * Chambolle_TV_L1(float * input_image, int nb_iter, float tau, float sigma, float lambda, float theta, int Width, int Height);
void ComputeImageGradient(float *&image, float*& GradX, float*& GradY, int Width, int Height);
void ComputeDivergence(float*& out_image, float*& p1, float*& p2, int Width, int Height);
void ProximalG(float*& image_out, float*& u, float* g, float lambda, float tau, int Width, int Height);
void ProximalF_star(float *&p1, float *&p2, float*& out1, float*& out2, int Width, int Height);
void imageAddition(float* &out_image, float*& a, float*& b, float lambda1, float lambda2, int Width, int Height);

#endif
