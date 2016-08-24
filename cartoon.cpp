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

#include "cartoon.h"
using namespace std;

/**
 * @file   cartoon.cpp
 * @brief  Cartoon + texture functions
 *
 *
 *
 * @author Vincent Le Guen <vincent.leguen@telecom-paristeh.org>
 */


/**
 * \brief   Main function to perform the TV-L1 cartoon texture decomposition
 *
 *
 * @param[in]   image : input image
 * @param[in]   nb_iter_max : number of iterations of the algorithm
 * @param[in]   tau, sigma, lambda, theta : parameters of the algorithm
 * @param[in]   Width, Height : size of the image
 * @param[out]  cartoon image
 */
float *Chambolle_TV_L1(float* input_image, int nb_iter_max, float tau, float sigma, float lambda, float theta, int Width, int Height) {
    // Initialization
    float *u = new float[Width*Height];
    float *u_old = new float[Width*Height];
    for(int k=0; k < Width*Height; k++) u[k] = input_image[k];

    float *p1 = new float[Width*Height];
    float *p2 = new float[Width*Height];

    float *div = new float[Width*Height];
    float *im1 = new float[Width*Height];
    float *im2 = new float[Width*Height];

    float *GradX =  new float[Width*Height];
    float *GradY =  new float[Width*Height];

    float *temp1 =  new float[Width*Height];
    float *temp2 =  new float[Width*Height];
    float E = numeric_limits<float>::infinity(), E_old;
    int it;

    for(it=1 ; it <= nb_iter_max ; it++) {

    for(int k=0; k < Width*Height ; k++)  u_old[k] = u[k];
        // First step : apply ProximalF
        ComputeImageGradient(u, GradX, GradY, Width, Height);
        imageAddition(im1,p1,GradX,1,sigma,Width,Height);  // im1 <- p1+sigma*GradX
        imageAddition(im2,p2,GradY,1,sigma,Width,Height);  // im2 <- p2+sigma*GradY
        ProximalF_star(im1, im2, p1, p2, Width, Height);

        // Second step : apply ProximalG
        ComputeDivergence(div, p1, p2, Width,Height);  // div <- divergence(p1,p2)
        imageAddition(im1,u,div,1,tau,Width,Height);  // im1 <- u + tau*div
        ProximalG(u, im1, input_image, lambda, tau, Width,Height);

        // Third step : u <- u + theta*(u-u_old)
        imageAddition(u,u,u_old,1+theta,-theta,Width,Height);

        // Energy calculation and stopping criterion
        if(it%10==0) {
            E_old = E;
            E = 0;
            ComputeImageGradient(u, GradX, GradY, Width, Height);
            for(int k=0; k < Width*Height ; k++) E += lambda * fabs(u[k]-input_image[k]) + sqrtf(GradX[k]*GradX[k] + GradY[k]*GradY[k]);
            E = E/(Width*Height); // normalized energy
            if(fabs(E-E_old) < 0.001) break;
            //cout << "iteration " << it << "  E = " << E << endl;
        }
    }

    //cout << "iterations : " << it << endl;
    delete[] div;
    delete[] im1;
    delete[] im2;
    delete[] p1;
    delete[] p2;
    delete[] temp1;
    delete[] temp2;
    delete[] u_old;
    delete[] GradX;
    delete[] GradY;
    return u;
}

/**
 * \brief   Compute the gradient of an image
 *
 *
 * @param[in]   image : input image
 * @param[in]   Width, Height : size of the image
 * @param[out]  GradX, GradY : coordinates X and Y of the gradient
 *
 * This function computes the discrete gradient of an image. The output gradient coordinates
 * gradX and gradY are given by reference.
 */
void ComputeImageGradient(float*& image, float*& GradX, float*& GradY,int Width, int Height) {
    if (!image) {
        printf("Null input image (ComputeImageGradient)");
        exit(-1);
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Boundary conditions
    for(int j=1 ; j<= Width ; j++)  GradX[(Height-1)*Width+j-1] = 0;
    for(int i=1 ; i<= Height ; i++) GradY[(i-1)*Width+ Height-1] = 0;

    for (int i = 1; i < Height ; i++) {
        for (int j = 1; j < Width ; j++) {
            GradX[(i-1)*Width + j-1] = image[i * Width + j-1] - image[(i-1) * Width + j-1];
            GradY[(i-1)*Width + j-1] = image[(i-1) * Width + j] - image[(i-1) * Width + j-1];
        }
    }
}

/**
 * \brief   Compute the discrete version of the divergence of a vector field
 *
 *
 * @param[in]   p1, p2 : the 2 input coordinates of the vector field
 * @param[in]   Width, Height : size of the image
 * @param[out]  out_image : the output image given by reference
 *
 */
void ComputeDivergence(float*& out_image, float*& p1, float*& p2, int Width, int Height) {
    if ((!p1)||(!p2)) {
        printf("Null input (ComputeDivergence)");
        exit(-1);
    }
    float term_p1 = 0, term_p2 = 0;

#ifdef _OPENMP
#pragma omp parallel for private(term_p1) private(term_p2)
#endif
    for (int i = 1; i <= Height ; i++) {
        for (int j = 1; j <= Width ; j++) {

            if(i==1) term_p1 = p1[(i-1)*Width+j-1];
            else if(i==Height) term_p1 = - p1[(i-2)*Width+j-1];
            else term_p1 = p1[(i-1)*Width+j-1] - p1[(i-2)*Width+j-1];

            if(j==1) term_p2 = p2[(i-1)*Width+j-1];
            else if(j==Width) term_p2 = - p1[(i-1)*Width+j-2];
            else term_p2 = p2[(i-1)*Width+j-1] - p2[(i-1)*Width+j-2];

            out_image[(i-1)*Width+j-1] = term_p1 + term_p2;
        }
    }
}

/**
 * \brief   Compute the ProximalG operator
 *
 *
 * @param[in]   u : input image for the proximal operator
 * @param[in]   g : the original input image of the cartoon texture algorithm
 * @param[in]   lambda, tau : parameters
 * @param[in]   Width, Height : size of the image
 * @param[out]  out_image : the output image given by reference
 *
 */
void ProximalG(float*& image_out, float*& u, float* g, float lambda, float tau, int Width, int Height) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 1; i <= Height ; i++) {
        for (int j = 1; j <= Width ; j++) {
            if(u[(i-1)*Width+j-1] - g[(i-1)*Width+j-1] > tau*lambda)          image_out[(i-1)*Width+j-1] = u[(i-1)*Width+j-1] - tau*lambda;
            if(u[(i-1)*Width+j-1] - g[(i-1)*Width+j-1] < - tau*lambda)        image_out[(i-1)*Width+j-1] = u[(i-1)*Width+j-1] + tau*lambda;
            if(abs(u[(i-1)*Width+j-1] - g[(i-1)*Width+j-1]) <=  tau*lambda)   image_out[(i-1)*Width+j-1] = g[(i-1)*Width+j-1];
        }
    }
}

/**
 * \brief   Compute the ProximalF operator
 *
 *
 * @param[in]   p1,p2 : input vector field for the proximal operator
 * @param[in]   Width, Height : size of the image
 * @param[out]  out1,out2 : output vector field
 *
 */
void ProximalF_star(float*& p1, float*& p2, float*& out1, float*& out2, int Width, int Height) {
    float norm_p = 0;

#ifdef _OPENMP
#pragma omp parallel for private(norm_p)
#endif
    for (int i = 1; i <= Height ; i++) {
        for (int j = 1; j <= Width ; j++) {
            norm_p = sqrtf ( p1[(i-1)*Width+j-1]*p1[(i-1)*Width+j-1] + p2[(i-1)*Width+j-1]*p2[(i-1)*Width+j-1] );
            out1[(i-1)*Width+j-1] = p1[(i-1)*Width+j-1] / fmax(1,norm_p);
            out2[(i-1)*Width+j-1] = p2[(i-1)*Width+j-1] / fmax(1,norm_p);
        }
    }
}

/**
 *
 *\brief   Compute the image addition :  image <- lambda1*a + lambda2*b
 * @param[in]   a,b : input images
 * @param[in]   lambda1, lambda2: parameters
 * @param[in]   Width, Height : size of the image
 * @param[out]  out_image : the output image given by reference
 *
 */
void imageAddition(float* &out_image, float*& a, float*& b, float lambda1, float lambda2, int Width, int Height) {
float value;
#ifdef _OPENMP
#pragma omp parallel for private(value)
#endif
    for (int i = 1; i <= Height ; i++) {
        for (int j = 1; j <= Width ; j++) {
            value = lambda1*a[(i-1)*Width+j-1] + lambda2*b[(i-1)*Width+j-1];
            out_image[(i-1)*Width+j-1] = value;

        }
    }
}
