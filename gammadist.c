/*
 *  gammadist.c
 *  ErlangC_CommandLine
 *
 *  Created by Gethryn Ghavalas on 19/04/10.
 *  Copyright 2010 Geth Home. All rights reserved.
 *
 */

#include "gammadist.h"
#include <math.h>
#include <stdio.h>





// returns the value of ln[gamma(xx)] for xx > 0
float gammln(float xx) {
	double x, y, tmp, ser;
	static double cof[6] = {76.18009172947146,-86.50532032941677,
		24.01409824083091, -11.231739572450155,
		0.1208650973866179e-2, -0.5395239384953e-5};
	int j;
	
	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j=0; j<=5; j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


// returns n! as a floating point number
long double factrl(int n) {
	// float gammln(float xx);
	static int  ntop = 4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0}; // fill in only as required.
	int j;
	
	if (n<0) nrerror("Negative factorial in routine factrl");
	if (n>32) return expl(gammln(n+1.0));
	while (ntop<n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}


// returns ln(n!)
float factln(int n) {
	static float a[101];
	
	if (n<0) nrerror("Negative factorial in routine factln");
	if (n<=1) return 0.0;
	if (n<=100) return a[n] ? a[n] : (a[n]=gammln(n+1.0)); // in range of table
	else return gammln(n+1.0);
}


// returns the binomial coefficient (n|k) as a floating point number.
float bico(int n, int k) {
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}


// returns the value of the beta function.
float beta(float z, float w) {
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}


// returns the incomplete gamma function P(a,x).
float gammp(float a, float x) {
	//void gcf(float *gammacf, float a, float x, float *gln);
	//void gser(float *gamser, float a, float x, float *gln);
	//voif nrerror(char error_text[]);
	float gamser, gammcf, gln;
	
	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments for routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}


// returns the incomplete gamma function Q(a,x) = 1 - P(a,x).
float gammq(float a, float x) {
	//void gcf(float *gammacf, float a, float x, float *gln);
	//void gser(float *gamser, float a, float x, float *gln);
	//voif nrerror(char error_text[]);
	float gamser,gammcf,gln;
	
	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments for routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf, a, x, &gln);
		return gammcf;
	}
}

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

// returns the incomplete gamma function P(a,x) evaluated by its series representation as gamser.
// also returns ln[gamma(a)] as gln.
void gser(float *gamser, float a, float x, float *gln) {
	int n;
	float sum, del, ap;
	
	*gln = gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser = 0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1; n<=ITMAX; n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
	}
}

// returns the incomplete gamma function Q(a,x) evaluated by its continued fraction representation
// as gammcf.  Also returns ln[gamma(a)] as gln.
void gcf (float *gammcf, float a, float x, float *gln) {
	int i;
	float an,b,c,d,del,h;
	
	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1; i<=ITMAX; i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void nrerror(char error_text[]) {
	printf("%s",error_text);;
}