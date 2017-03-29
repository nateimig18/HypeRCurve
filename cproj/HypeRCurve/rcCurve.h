#ifndef RCCURVE_H
#define RCCURVE_H

#include <iostream>
#include <math.h>

#include "rcMath.h"

#define CONSTRAIN(X, XL, XH) 	(((X) < (XL))? (XL) : (((X) > (XH))? (XH) : (X)))

struct rcParam1{
	float m0, mw, rmx;
	float w1, w2, g;
};
typedef struct rcParam1 rcParam1_t;

struct rcParam2{
	float m0, mw, np, rmx;
	float w1, w2, g, invNp;
};
typedef struct rcParam2 rcParam2_t;

void initRcParam1(rcParam1_t *pRcParams1, float m0, float mw, float rmx);
float calcRcParam1(float x, rcParam1_t *pRcParams1);

void initRcParam2(rcParam2_t *pRcParams2, float m0, float mw, float np, float rmx);
float calcRcParam2(float x, rcParam2_t *pRcParams2);

/*
	Custom RC Curve 1 (w/ Intuitive params)
	m0  = normalized center stick slope (0.05 < m0 < 0.95)
	mw  = normalized center stick linear width
	rmx = max rate

	tol = 0.2		// static parameter tolerance

	x = normalized stick input  (i.e. -1 < x < 1)
	y = normalized stick output (i.e. -1 < y < 1)
	w = un-normalized rate output (i.e. -rmx < w < rmx)

	g = lambertw(-1, m0 * tol/(1 - m0) * mw .* log(mw)) ./ log(mw);	// define gamma from width
	y = sgn(x)*[(1-m0)*|x|^g + m0*|x|]
	w = rmx * y = sgn(x) * rmx * [(1-m0)*|x|^g + m0*|x|]
*/
void initRcParam1(rcParam1_t *pRcParams1, float m0, float mw, float rmx){
	static const float tol = 0.2;
	float npmin, invNp, Ap, t, logmw, g;

	m0 = CONSTRAIN(m0, 0.05, 0.95);
	mw = CONSTRAIN(mw,  0.0,  1.0);

	logmw = logf(mw);
	t = m0 * tol/(1 - m0) * mw * logmw;
	g = myLambertW(t) / logmw;

	pRcParams1->m0 = m0;
	pRcParams1->mw = mw;
	pRcParams1->rmx = rmx;

	pRcParams1->w1 = (1 - m0);
	pRcParams1->w2 = m0;
	pRcParams1->g = g;

    printf("m0 = %4.2f, mw = %4.2f, rmx = %5.1f, g = %5.3f, w1 = %5.3f, w2 = %5.3f\n", pRcParams1->m0, pRcParams1->mw, pRcParams1->rmx, pRcParams1->g, pRcParams1->w1, pRcParams1->w2);
}

float calcRcParam1(float x, rcParam1_t *pRcParams1){
	float y, w;
	float sgnx = (x < 0) ? -1.0f : 1.0f;

	x *= sgnx;
	// y = myPowf(x, pRcParams2->g);
	y = sgnx * (pRcParams1->w1 * myPowf(x, pRcParams1->g) + pRcParams1->w2 * x);
	w = pRcParams1->rmx * y;
	return(w);
}

/*
	Custom RC Curve 2 (w/ Intuitive params)
	m0  = normalize center stick slope (0.05 < m0 < 0.95)
	mw  = center stick linear width
	p	 = mid-stick linearity control (allows for a much more gradual transition into the end-stick hyper-slope region)
	rmx = max rate

	tol = 0.2		// static parameter tolerance

	x = normalized stick input  (i.e. -1 < x < 1)
	y = normalized stick output (i.e. -1 < y < 1)
	w = un-normalized rate output (i.e. -rmx < w < rmx)


	pmin = powf((((1/m0) - 1.002)), -0.97) + 0.17;	// Estimate minimum p can be.
	p = (p > pmin) ? p : pmin;


	Ap = p.*(exp(1./p)-1);
	g = lambertw(-1, (t0*m0*exp(mw / p))/(1 - m0*Ap) * mw * log(mw)) / log(mw);
	y = sgn(x) * ((1 - m0*Ap) * |x|^g + m0 * p * (exp(|x|/p) - 1));
	w = rmx * y = sgn(x) * rmx * ((1 - m0*Ap) * |x|^g + m0 * p * (exp(|x|/p) - 1));
*/

void initRcParam2(rcParam2_t *pRcParams2, float m0, float mw, float np, float rmx){
	static const float tol = 0.2;
	float npmin, invNp, Ap, t, logmw, g;

	m0 = CONSTRAIN(m0, 0.05, 0.95);
	mw = CONSTRAIN(mw,  0.0,  1.0);

	npmin = myMinRho(m0);
	np = (np < npmin) ? npmin : np;

	invNp = 1.0f / np;
	logmw = logf(mw);
	Ap = np * ( myExp(invNp) - 1);                                          // Ap = np * ( myExp(invNp) - 1);
	t =  ( (tol * m0 * myExp( mw * invNp)) / (1 - m0*Ap) ) * mw * logmw;    // t =  ( (tol * m0 * myExp( mw * invNp)) / (1 - m0*Ap) ) * mw * logmw;
	g = myLambertW(t) / logmw;

	pRcParams2->m0 = m0;
	pRcParams2->mw = mw;
	pRcParams2->np = np;
	pRcParams2->invNp = invNp;
	pRcParams2->rmx = rmx;
	pRcParams2->w1 = (1 - m0 * Ap);
	pRcParams2->w2 = m0*np;
	pRcParams2->g = g;

	printf("m0 = %4.2f, mw = %4.2f, np = %6.4f, rmx = %5.1f, g = %5.3f, w1 = %5.3f, w2 = %5.3f; npmin = %5.3f\n", pRcParams2->m0, pRcParams2->mw, pRcParams2->np, pRcParams2->rmx, pRcParams2->g, pRcParams2->w1, pRcParams2->w2, npmin);
}

float calcRcParam2(float x, rcParam2_t *pRcParams2){
	float y, w;
	float sgnx = (x < 0) ? -1.0f : 1.0f;

	x *= sgnx;
	// y = myPowf(x, pRcParams2->g);
	y = sgnx * (pRcParams2->w1 * myPowf(x, pRcParams2->g) + pRcParams2->w2 * (myExp(pRcParams2->invNp * x) - 1));
	w = pRcParams2->rmx * y;
	return(w);
}

#endif
