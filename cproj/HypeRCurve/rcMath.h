#ifndef RCMATH_H
#define RCMATH_H

#include <stdint.h>     // Needed for standart int types.
#include <math.h>

/*----------------------------------- Type Defines -----------------------------------*/
/* Single Precision Floating Point type defines*/
struct sSFPbitFld {
	uint32_t m : 23;        // mant = 1 + m[22:0]/2^23
	uint32_t e : 8;         // exp  = e[7:0] - (2^7 - 1)
	uint32_t s : 1;         // sign = (-1)^s
};                          // vx = sign * mant * 2^exp = (-1)^s * (1 + m[22:0]/2^23) * 2^{e[7:0] - (2^7 - 1)}
typedef struct sSFPbitFld SFPbitFld_t;

union uSFPtype{
	float ff;
	uint32_t fi;
	SFPbitFld_t SFPbits;
};
typedef union uSFPtype SFPtype_t;


/*------------------------------------ Prototypes -------------------------------------*/
float myExp2(float x);
float myLog2(float x);
float myPowf(float x, float p);
float myExp(float x);
float myLambertW(float x);
float myMinRho(float m0);

/*----------------------------------- Definitions -------------------------------------*/
float myExp2(float x){
	SFPtype_t ret;
	float xf = floorf(x);
	float xr = x - xf;

    if(xf < -127) return(0.0f);

	//--------------- APPROXIMATE 2^x, FOR 0 <= x < 1 -------------------
	// Note: Compounding error still exists but is bounded from [-inf, 0], which is the case for this restricted powf application.
	// Improving would require either a pade approximation which unfortunately includes the un-holy division (14 CLK's) or adding
	// an iterative Halley's or Newton's method to increase the accuracy which is not worth further digression. 8.(

	// Choose the appropriate approximation. TODO: implement using fused multiply accumulate operation (Note: M7 this is a single CLK cycle).
	// ret.ff = 0x1.00a246p0 + xr * (0xa.6abp-4 + xr * 0x5.81075p-4);                                                      // <--- Working with powf. :)
    ret.ff = 0x1.p0 + xr * (0xb.224a4p-4 + xr * (0x3.995f08p-4 + xr * 0x1.4456b8p-4));                               // <--- Working with powf. :)
	// ret.ff = 0x1.00003ep0 + xr * (0xb.1663ap-4 + xr * (0x3.ddc048p-4 + xr * (0xd.3b945p-8 + xr * 0x3.81b0acp-8)));   // <--- Working with powf. :)
	// ret.ff = 0xf.ffffep-4 + xr * (0xb.17297p-4 + xr * (0x3.d79cbcp-4 + xr * (0xe.4d401p-8 + xr * (0x2.4a13c8p-8 + xr * 0x7.c49568p-12))));

	ret.SFPbits.e = ((int8_t) xf) + 127;

	return(ret.ff);
}

float myLog2(float x){
	union { float f; int i; } v;
	float ret, E2, L2;

	v.f = x;
	ret = (int8_t)(v.i >> 23) - 127;					// Take exp. from arguement to initialize output
	v.i = (v.i & ~(0xFF << 23)) | (0x7F << 23);     // Force exponent of float to zero, & value between [1, 2)

	//---------------- APPROXIMATE log2(x), FOR 1 <= x < 2 ----------------
	// Choose the appropriate approximation. TODO: implement using fused multiply accumulate operation (Note: M7 has a single CLK cycle fused multiply and accumulate).
    // L2 = v.f - 1;        // Tol = +-86.1e-3
    // L2 = -0x1.b20438p0 + v.f * (0x2.0b0654p0 + v.f * (-0x5.9021cp-4));      // Tol = +-7.784e-3
    // L2 = -0x2.2a24cp0 + v.f * (0x3.10d81cp0 + v.f * (-0x1.0f83dcp0 + v.f * 0x2.8d07fp-4));  // Tol = +-0.898e-3

    L2 = -0x2.851078p0 + v.f * (0x4.1618bp0 + v.f * (-0x2.2257a4p0 + v.f * (0xa.65721p-4 + v.f * (-0x1.507b62p-4))));   // Tol = +-0.116e-3

    // L2 = -0x2.ce7f28p0 + v.f * (0x5.1c8b68p0 + v.f * (-0x3.933bfp0 + v.f * (0x1.a54af2p0 + v.f * (-0x6.bb68d8p-4 + v.f * 0xb.9b4afp-8))));      // Tol = +-1.6e-5
    // L2 = -0x3.0b16f8p0 + v.f * (0x6.200b78p0 + v.f * (-0x5.5c449p0 + v.f * (0x3.4d0798p0 + v.f * (-0x1.45e1e2p0 + v.f * (0x4.6c513p-4 + v.f * (-0x6.9abb2p-8))))));      // Tol = +-2.17e-6

	return(ret + L2);
}

float myPowf(float x, float p){
	float Exp2PLog2, ret;

	if(x != 0){
		Exp2PLog2 = p * myLog2(x);

		ret = myExp2(Exp2PLog2);
	}else{
		ret = 0;
	}

	return(ret);
}

#define LOG2E   (0x1.715476p0)
float myExp(float x){
	float y;
	// This is an approximation for the Transcendental Exponential Function exp(x) = e^x on a bounded interval
	// of 0 <= x <= 2.71828, e^-1 <= rho < 1
	// of 0 <= x <= 5.6

	//---------------- APPROXIMATE exp(x), FOR 0 <= x <= exp(1); m0 >= 0.1; ----------------
	// y = 0xf.fef0fp-4 + x * (0x1.02712ap0 + x * (0x7.2bb18p-4 + x * (0x4.2f8bd8p-4 + x * (-0x7.cfbe3p-8 + x * 0x7.c5e34p-8))));	// ERR = 3.917756e-3
	// y = 0x1.0001a8p0 + x * (0xf.fafc9p-4 + x * (0x8.24dc2p-4 + x * (0x2.4b0f9cp-4 + x * (0x1.1a257ep-4 + x * (-0x1.c297fp-8 + x * 0x1.4e6aep-8)))));		// ERR = 3.826171e-4
	// y = 0xf.fffdcp-4 + x * (0x1.0008acp0 + x * (0x7.fad638p-4 + x * (0x2.bc83acp-4 + x * (0x8.da923p-8 + x * (0x3.a9b234p-8 + x * (-0x5.164f38p-12 + x * 0x3.01441p-12))))));	// ERR = 3.268547e-5
	// y = 0x1.000002p0 + x * (0xf.fff41p-4 + x * (0x8.00946p-4 + x * (0x2.a81004p-4 + x * (0xb.048e3p-8 + x * (0x1.b923f2p-8 + x * (0xa.04745p-12 + x * (-0xc.3150dp-16 + x * 0x6.07148p-16)))))));	// ERR = 2.517432e-6

	//---------------- APPROXIMATE exp(x), FOR 0 <= x <= 4.4; m0 >= 0.05 or rho >= 0.2275 ----------------
	// sollya command: "display=hexadecimal; E =[||]; B=[0;4.4]; F=exp(x); for i from 1 to 10 do { P=fpminimax(F, i, [|SG...|], B, relative, floating); print(horner(P)); R=dirtyinfnorm(P-F,B); E=E:.round(R, SG, RN);}; display=decimal; print(E);"
	// Choose the appropriate approximation. TODO: implement using fused multiply accumulate operation (Note: M7 this is a single CLK cycle).
	// y = 0x1.002c6p0 + x * (0xf.9dd42p-4 + x * (0xa.0013ap-4 + x * (-0xe.9113p-8 + x * (0x3.75ac9cp-4 + x * (-0xe.04e01p-8 + x * 0x2.c17d94p-8)))));		// 5.514643e-2
	// y = 0xf.ff9d2p-4 + x * (0x1.01115ap0 + x * (0x7.8dd268p-4 + x * (0x3.b962cp-4 + x * (-0x7.c87fe8p-8 + x * (0xc.564fdp-8 + x * (-0x2.80f564p-8 + x * 0x6.6849e8p-12))))));		// 7.6815072e-3
	// y = 0x1.0000c4p0 + x * (0xf.fd6c2p-4 + x * (0x8.1563dp-4 + x * (0x2.6a3098p-4 + x * (0x1.070c8ap-4 + x * (-0x2.484c6p-8 + x * (0x2.3ab838p-8 + x * (-0x6.0747ap-12 + x * 0xc.fdafap-16)))))));		// 9.4887166e-4
	// y = 0xf.fffeap-4 + x * (0x1.00058p0 + x * (0x7.fc8b68p-4 + x * (0x2.b78278p-4 + x * (0x9.373d9p-8 + x * (0x3.94583p-8 + x * (-0x7.c68438p-12 + x * (0x5.6d2fp-12 + x * (-0xc.8f36dp-16 + x * 0x1.75b528p-16))))))));		// 1.05118786e-4
	// y = 0x1.000002p0 + x * (0xf.fff5ep-4 + x * (0x8.007afp-4 + x * (0x2.a8804p-4 + x * (0xa.f7b4fp-8 + x * (0x1.c1ab92p-8 + x * (0xa.3ce89p-12 + x * (-0x1.5028a4p-12 + x * (0xb.62812p-16 + x * (-0x1.705cap-16 + x * 0x2.5aeecp-20)))))))));	// 1.0709154e-5

    y = myExp2(LOG2E * x);
	return(y);
}

#define INVE    	(0x1.78B564p-02)
#define INVSQRT2	(0x1.6A09E6p-01)

float myLambertW(float x){
    // Note x must be between -1/e <= x < 0
	// Absolute Error <= 0.002 "Good Enough"
	// Based on the "Analytical approximations for real values of the W-function" which comes from
	// "Real values of the W-Fucntion" by D.A. Barry, P.J Culligan-Hensley, S.J. Barry

	if((x <= -INVE) || (x > 0)) return(NAN);

	// Based on "Analytical approximations for real values of the W-function"
	static const float M1 = 0.3361, M2 = -0.0042, M3 = -0.0201;
	float sigma, sqrtSigma, expSqrtSigma, w, t;

	sigma = -(1 + logf(-x));
	sqrtSigma = sqrtf(sigma);
	t = M1 * INVSQRT2 * sqrtSigma / (1 + M2 * sigma * expf(M3 * sqrtSigma));

	w = -(1+sigma);
	w += -2*(1 - 1/(1 + t)) / M1;

	return(w);
}

float myMinRho(float m0){
	float m0inv, ret;

	if((0.01 <= m0) && (m0 <= 0.99)){
		m0inv = 1.0f/m0;
		ret = powf(m0inv - 1.002, -0.97) + 0.181;	// Custom rough fit above implicit curve p(x) defined by x = p(x)*(exp(1/p(x)) - 1)
	}else{
		ret = 15;
	}
	return(ret);
}

#endif
