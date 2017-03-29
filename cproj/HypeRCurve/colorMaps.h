#ifndef COLORMAPS_H
#define COLORMAPS_H

void LumToParulaRGB(float L, float *pR, float *pG, float *pB);
void LumToViridisRGB(float L, float *pR, float *pG, float *pB);
void LumToSurreyRGB(float L, float *pR, float *pG, float *pB);

const float ParulaPts[] = {0.0f, 0x2.323234p-4, 0x5.c5c5c8p-4, 0xd.adadbp-4};		// {0, 35/255, 92/255, 218/255}

const float ParulaR1[] = { 0x3.5460acp-4, -0xf.4daf1p-8,  0x3.494ddcp4,  -0xc.0bf93p8, 0x1.25fdc4p16,  -0xe.180e2p16, 0x5.0d37b8p20, -0xa.bfe3cp20};
const float ParulaR2[] = {   0xa.a6533p0,  -0x1.457d4p8, 0x1.0691fep12, -0x7.48c94p12, 0x1.ed5c28p16,  -0x4.da58cp16, 0x6.b4bb18p16, -0x3.e8f7cp16};
const float ParulaR3[] = {   0x3.0e82ep0, -0x3.e9d114p4,  0x2.085b38p8,  -0x8.a7913p8, 0x1.45e914p12, -0x1.aac5eap12, 0x1.23c938p12, -0x5.154a48p8};
const float ParulaR4[] = {  0x2.13092cp8,  -0xc.67d55p8, 0x1.d46676p12, -0x2.2158ap12, 0x1.39c218p12,  -0x4.7738c8p8,             0,             0};

const float ParulaG1[] = { 0x2.a92a3p-4, 0x1.a87502p0, -0x2.456afcp4,    0x8.2257p8, -0xc.4849dp12,  0x9.11157p16, -0x3.3008bcp20, 0x6.d2c548p20};
const float ParulaG2[] = {-0x2.21b69cp0, 0x4.dc8ed8p4,  -0x4.1c408p8,  0x1.f15d2p12, -0x8.b3c6dp12, 0x1.6f543ep16, -0x2.0e5c84p16, 0x1.3b5c3cp16};
const float ParulaG3[] = {-0x2.445964p4, 0x1.f2e4e2p8,  -0xb.4f5fap8, 0x2.4a40acp12, -0x4.83207p12,  0x5.8514bp12,  -0x4.136d6p12, 0x1.a5d6f2p12, -0x4.66c9p8};
const float ParulaG4[] = {-0x1.3ff742p4,  0x4.da23fp4,  -0x6.b8fcap4,  0x4.016838p4,  -0xd.2dd84p0,             0,              0,             0};

const float ParulaB1[] = { 0x8.779a7p-4,    0x3.14637p0, -0x7.a94508p0,   0x2.697ae8p8, -0x4.9da3c8p12, 0x4.45d518p16, -0x1.d69d4cp20, 0x4.9c5d38p20};
const float ParulaB2[] = { -0x4.4535ep0,    0xa.6cc14p4,  -0x8.e6eeep8,  0x4.2a3aa8p12, -0x1.26dadap16, 0x3.006e14p16,  -0x4.4154cp16, 0x2.8714a8p16};
const float ParulaB3[] = { -0x4.3a7a3p0,    0x4.ba65bp4, -0x1.d8594ep8,    0x6.482f3p8,    -0xc.93cep8,   0xe.a7921p8,   -0x9.28b8ep8,   0x2.5e057p8};
const float ParulaB4[] = {  0x8.bc308p8, -0x2.edcf34p12, 0x6.493ab8p12, -0x6.beafa8p12,  0x3.9e395cp12,  -0xc.6b764p8,              0,             0};

void LumToParulaRGB(float L, float *pR, float *pG, float *pB){
	float maccR = 0, maccG = 0, maccB = 0;
	float const *ptrR, *ptrG, *ptrB;

	if(L <= ParulaPts[1]){
		ptrR = &ParulaR1[7];
		ptrG = &ParulaG1[7];
		ptrB = &ParulaB1[7];
	}else if(L <= ParulaPts[2]){
		ptrR = &ParulaR2[7];
		ptrG = &ParulaG2[7];
		ptrB = &ParulaB2[7];
	}else if(L <= ParulaPts[3]){
		ptrR = &ParulaR3[7];
		ptrG = &ParulaG3[8];
		ptrB = &ParulaB3[7];

		maccG = *ptrG--;
	}else{
		ptrR = &ParulaR4[7];
		ptrG = &ParulaG4[7];
		ptrB = &ParulaB4[7];
	}

	for(int N = 8; N > 0; N--){
		maccR = maccR * L + *ptrR--;
		maccG = maccG * L + *ptrG--;
		maccB = maccB * L + *ptrB--;
	}

	*pR = maccR;
	*pG = maccG;
	*pB = maccB;
}

const float ViridisPts[] = {0.0f, 0x5.454548p-4, 0x9.59596p-4, 1.0f};	// {0, 84/255, 149/255, 1.0}

const float ViridisR1[] = {0x4.45a6p-4, 0x6.950e2p-4, -0x2.aef2c8p0, 0x5.fb52f8p-4, 0x5.bedcc8p0};
const float ViridisR2[] = {0x1.9af78p0, -0xc.d3285p0, 0x2.cef05p4, -0x4.7f533p4, 0x2.b23c8p4};
const float ViridisR3[] = {0x3.74e0fcp0, -0xd.e314p0, 0x1.12e098p4, -0x4.f30d7p0, -0xc.e838dp-4};

const float ViridisG1[] = {0x1.3f71bp-8, 0x1.628f54p0, 0xb.5e07bp-4, -0x8.5046cp0, 0xd.c5b4bp0};
const float ViridisG2[] = {-0x1.78ed2ap-4, 0x2.7e88ccp0, -0x4.d1abp0, 0x6.9da7f8p0, -0x3.654cd4p0};
const float ViridisG3[] = {0x1.ef352ep0, -0x9.42bc2p0, 0x1.4bd108p4, -0x1.21c724p4, 0x5.9ae28p0};

const float ViridisB1[] = {0x5.4548ep-4, 0x1.8fe7bep0, -0x2.32e53cp0, -0x5.6170b8p0, 0xc.1fe4p0};
const float ViridisB2[] = {0x7.6cd308p-4, 0xb.f758fp-4, -0x2.a28a88p0, 0x5.073f78p0, -0x4.1c8a4p0};
const float ViridisB3[] = {-0x3.75921p4, 0x1.8605ccp8, -0x4.3a1c88p8, 0x5.d5d328p8, -0x4.02b938p8, 0x1.187c34p8};

void LumToViridisRGB(float L, float *pR, float *pG, float *pB){
	float maccR = 0, maccG = 0, maccB = 0;
	float const *ptrR, *ptrG, *ptrB;

	if(L <= ViridisPts[1]){
		ptrR = &ViridisR1[4];
		ptrG = &ViridisG1[4];
		ptrB = &ViridisB1[4];
	}else if(L <= ViridisPts[2]){
		ptrR = &ViridisR2[4];
		ptrG = &ViridisG2[4];
		ptrB = &ViridisB2[4];
	}else{
		ptrR = &ViridisR3[4];
		ptrG = &ViridisG3[4];
		ptrB = &ViridisB3[5];
		maccB = *ptrB--;
	}

	for(int N = 5; N > 0; N--){
		maccR = maccR * L + *ptrR--;
		maccG = maccG * L + *ptrG--;
		maccB = maccB * L + *ptrB--;
	}

	*pR = maccR;
	*pG = maccG;
	*pB = maccB;
}

const float SurreyR[] = {0x1.a3d07ap-12, -0x7.034638p-8,  0x7.8db23p-4, -0x5.578cd8p-4,  0x7.9d5afp0,  -0x9.f48ebp0,  0x3.3ab888p0};
const float SurreyG[] = {             0,   0x1.013166p0, -0x1.b0a784p0,   0xb.ac0b60p0,  -0x1.42f9p4,   0xc.48893p0, -0x2.158e7cp0};
const float SurreyB[] = {             0,   0x3.b33810p0,   0x7.1b951p0,  -0x3.9aa450p4, 0x5.0bfad8p4, -0x2.202e5cp4,  0x1.1ea4dcp0};

void LumToSurreyRGB(float L, float *pR, float *pG, float *pB){
    float maccR = 0, maccG = 0, maccB = 0;
	float const *ptrR, *ptrG, *ptrB;        // float pointer to constant READONLY data :)

    ptrR = &SurreyR[6];
    ptrG = &SurreyG[6];
    ptrB = &SurreyB[6];

	for(int N = 7; N > 0; N--){
		maccR = maccR * L + *ptrR--;
		maccG = maccG * L + *ptrG--;
		maccB = maccB * L + *ptrB--;
	}

	*pR = maccR;
	*pG = maccG;
	*pB = maccB;
}


#endif
