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
 * @mainpage  Cartoon + Texture decomposition
 *
 * README.txt:
 * @verbinclude README.txt
 */


/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 *
 * @author Vincent Le Guen <vincent.leguen@telecom-paristech.org>
 */

#include "cartoon.h"
#include "io_png.h"

extern "C" {
#include "iio.h"
}

using namespace std;

int main(int argc, char **argv) {

    if (argc < 5) {
        printf("usage: cartoonTexture image lambda cartoon texture \n");
        return EXIT_FAILURE;
    }

    int Width, Height, nb_channels;
    float *input_image = iio_read_image_float_split(argv[1], &Width, &Height, &nb_channels);

    if (nb_channels == 2) { nb_channels = 1; }  // we do not use the alpha channel
    if (nb_channels > 3) { nb_channels = 3; }   // we do not use the alpha channel

    int size_image = Width * Height;

    // parameters: input
    float lambda = atof(argv[2]);

    // build cartoon
    float *cartoon;

    // We fix the parameters
    int nb_iter_max = 1000;
    float tau = 0.35;
    float sigma = 0.35;
    float theta = 1;

    float initial_time = clock();

    // Perform the TV-L1 algorithm on the 3 channels for a color image
    float *cartoon1;
    float *cartoon2;
    float *cartoon3;

    if (nb_channels == 1) {
        cartoon = Chambolle_TV_L1(input_image, nb_iter_max, tau, sigma,lambda, theta, Width, Height);
    }
    else {
        cartoon = new float[nb_channels * size_image];
        cartoon1 = Chambolle_TV_L1(input_image, nb_iter_max, tau, sigma,lambda, theta, Width, Height);
        for(int j=0; j < size_image; j++) cartoon[j] = cartoon1[j];

        cartoon2 = Chambolle_TV_L1(&input_image[size_image], nb_iter_max, tau, sigma,lambda, theta, Width, Height);
        for(int j=0; j < size_image; j++) cartoon[size_image+j] = cartoon2[j];

        cartoon3 = Chambolle_TV_L1(&input_image[2*size_image], nb_iter_max, tau, sigma,lambda, theta, Width, Height);
        for(int j=0; j < size_image; j++) cartoon[2*size_image+j] = cartoon3[j];
    }

    float final_time = clock();
    float cpu_time = (final_time - initial_time) / CLOCKS_PER_SEC;
    cout << "time: " << cpu_time << "s   lambda = " <<  lambda  << endl;

    // BV norm
    float * GradX =  new float[Width*Height];
    float * GradY =  new float[Width*Height];

    float BV_norm = 0;
    if(nb_channels == 1) {  // graylevel image
        ComputeImageGradient(cartoon, GradX, GradY, Width, Height);
        for(int i=0;i<size_image;i++) BV_norm += sqrtf(GradX[i]*GradX[i] + GradY[i]*GradY[i]);
        //cout << "graylevel - lambda  = " << lambda  <<  "  BV norm = " << BV_norm/10 << endl;
    }

    if(nb_channels == 3) {  // color image
        ComputeImageGradient(cartoon1, GradX, GradY, Width, Height);
        for(int i=0;i<size_image;i++) BV_norm += sqrtf(GradX[i]*GradX[i] + GradY[i]*GradY[i]);

        ComputeImageGradient(cartoon2, GradX, GradY, Width, Height);
        for(int i=0;i<size_image;i++) BV_norm += sqrtf(GradX[i]*GradX[i] + GradY[i]*GradY[i]);

        ComputeImageGradient(cartoon3, GradX, GradY, Width, Height);
        for(int i=0;i<size_image;i++) BV_norm += sqrtf(GradX[i]*GradX[i] + GradY[i]*GradY[i]);
        //cout << "color - lambda  = " << lambda  <<  "  BV norm = " << BV_norm/10 << endl;
    }

    // compute texture as difference as transform from [-20,20] to [0,255]
    float *texture = new float[nb_channels * size_image];

    for (int i=0; i < nb_channels*size_image; i++) {
        texture[i] = input_image[i]- cartoon[i];
    }

    // save cartoon and texture images
    iio_save_image_float_split(argv[3], cartoon, Width, Height, nb_channels);
    iio_save_image_float_split(argv[4], texture, Width, Height, nb_channels);

    delete[] cartoon1;
    delete[] cartoon2;
    delete[] cartoon3;
    delete[] cartoon;
    delete[] texture;
    return EXIT_SUCCESS;
}
