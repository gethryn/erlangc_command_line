/*
 *  gammadist.h
 *  ErlangC_CommandLine
 *
 *  Created by Gethryn Ghavalas on 19/04/10.
 *  Copyright 2010 Geth Home. All rights reserved.
 *
 */

float gammln(float xx);
long double factrl(int n);
float factln(int n);
float bico(int n, int k);
float beta(float z, float w);
float gammp(float a, float x);
float gammq(float a, float x);
void gcf(float *gammacf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
void nrerror(char error_text[]);