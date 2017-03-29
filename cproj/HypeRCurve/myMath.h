#ifndef MYMATH_H
#define MYMATH_H

#define PI      (0x3.243f6cp0)
// Reduce [0, 2*pi] into 8 octants, equates to range reduction x = X + k*C, where C = 2*pi / 2^3 = pi/4
// and an approximation of cos(X) & sin(X) on the interval -pi/8 <= X <= pi/8. Using trigonemtric angle
// sum identities. i.e.
//   cos(X + k*C) = cos(X)*cos(k*C) - sin(X)*sin(k*C)
//   sin(X + k*C) = cos(X)*sin(k*C) + sin(X)*cos(k*C)
// For an 8-point octant this requires 8 known cosine & sine values at each of the 8 angles, however
// since sine is a shifted version of cosine the only one table of 8 values of cosine are really needed.
// Before applying a minimax bounded error polynomial over an interval, its important to look at coeffs
// on the taylor polynomial of each. By restricting the nonzero monomials of the minimax to those
// nonzero monomials of the taylor series a more stable approximation can be formed without .

// Cody & Waite's Range Reduction Constants
#define PI_4    (0x1.921FB6p-01)    //     pi/4 = 7.853981852531433105e-01 = HI(pi/4) + LO(pi/4)
#define PI_4HI  (0x1.920000p-01)    // HI(pi/4) = 7.8515625e-01
#define PI_4LO  (0x1.FB5444p-13)    // LO(pi/4) = 2.4191339616663753986358642578125e-04
#define PI_4INV (0x1.45F306p+00)    //     4/pi = 1.27323949337005615234375e+00

// cos(k*pi/4), sin(k*pi/4)
float cosOff4LUT[] = { 0x1.000000p+00,  0x1.6A09E6p-01,  0x0.000000p+00, -0x1.6A09E6p-01, -0x1.000000p+00, -0x1.6A09E6p-01,  0x0.000000p+00,  0x1.6A09E6p-01 };
float myCoSi(float x, float *pC, float *pS);

float myCoSi(float x, float *pC, float *pS){
    int     m, ms, mc;
    float   xI, xR, xR2;
    float   c, s, cy, sy;

    // Cody & Waite's range reduction Algorithm, [-pi/4, pi/4]
    xI  = floorf(x*PI_4INV + 0.5);
    xR  = (x - xI*PI_4HI) - xI*PI_4LO;
    m   = (int) xI;
    xR2 = xR*xR;

    // Find cosine & sine index for angle offsets indices
    mc = (  m  ) & 0x7;     // two's complement permits upper modulus for negative numbers =P
    ms = (m + 6) & 0x7;     // two's complement permits upper modulus for negative numbers =P

    // Find cosine & sine
    cy = cosOff4LUT[mc];     // Load pointer of matching order horner polynomial coeffcients
    sy = cosOff4LUT[ms];     // Load pointer of matching order horner polynomial coeffcients

    /* Cosine Minimax Approximations */
    // c = cosf(xR);
    // c = 0xf.ff79fp-4 + x2 * (-0x7.e58e9p-4);   // TOL = 1.2786e-4
    c = 0xf.ffffdp-4 + xR2 * (-0x7.ffebep-4 + xR2 * 0xa.956a9p-8);  // TOL = 1.7882e-7
    // c = 0x1.p0 + x2 * (-0x7.fffff8p-4 + x2 * (0xa.aa9bap-8 + x2 * (-0x5.a7948p-12)))   // TOL =  3.06321e-10

    /* Sine Remez Approximation */
    // s = sinf(xR);
    // s = xR * (0xf.ffbf7p-4 + x2 * (-0x2.a41d0cp-4));    // TOL = 4.835251e-6
    s = xR * (0xf.fffffp-4 + xR2 * (-0x2.aaa65cp-4 + xR2 * 0x2.1ea25p-8));  // TOL = 1.1841e-8
    // s = xR * (0x1.p0 + x2 * (-0x2.aaaaa8p-4 + x2 * (0x2.222048p-8 + x2 * (-0xc.f0ce4p-16))));   // TOL = 3.33542e-10
    // s = sinf(xR);

    // Note: c & s are local approximations over [-pi/8, +pi/8], cy & sy are are cosine & sine
    // offset values, thus the offset angle creates a (2x2) rotation matrix to generate over any
    // interval.
    *pC = c*cy - s*sy;
    *pS = c*sy + s*cy;

    return(m);
}

#endif
