#include <stdio.h>
#include <stdlib.h>     // Needed for hex float representation...
#include <float.h>      // Needed for hex float representation...
#include <stdint.h>     // Needed for standart int types.
#include <math.h>       // Needed for bounded error testing...

#include "koolplot.h"
#include "rcMath.h"
#include "rcCurve.h"
#include "colorMaps.h"
#include "myMath.h"

#define PI  (0x3.243f6cp0)

Plotdata x, y;

void plotRcCurve1(float m0, float WB0, float WB1, float rmx, int Nc);
void plotRcCurve2(float m0, float WB0, float WB1, float np, float rmx, int Nc);
void plotPowfError(void);
void plotCoSiError(float xStart, float xEnd, float xInc);

int main(){
    int Nc;
    float np, rmx = 1.0f;

    clear(x); clear(y);

    np = 0.900;

    // plotRcCurve1(0.475, 0.15, 0.90, rmx, 30);
    plotRcCurve2(0.15, 0.05, 0.90, np, rmx, 30);
    // plotPowfError();
    // plotCoSiError(-10*PI, 10*PI, 0x1.0p-8);

    return(0);
}

void plotRcCurve1(float m0, float WB0, float WB1, float rmx, int Nc){
    rcParam1_t rcParams;
    float w0, dWB;
    float r, g, b;

    dWB = (WB1 - WB0) / (Nc - 1);

    for(int i = 0; i < Nc; i++){
        // LumToSurreyRGB(((float) i)/ Ni, &r, &g, &b);
        // LumToViridisRGB(((float) i)/ Ni, &r, &g, &b);
        LumToParulaRGB(((float) i)/ Nc, &r, &g, &b);
        setColor(x, y, COLOR(255*r, 255*g, 255*b));

        w0 = dWB * i + WB0;
        printf(" > w0 = %4.3f\n", w0);
        initRcParam1(&rcParams, m0, w0, rmx);

        for(float xx = 0.0, yy; xx <= 1.0; xx += 0x1.0p-10){
            yy = calcRcParam1(xx, &rcParams);
            // yy = myPowf(xx, w0+2);
            // printf(" > x1 = %5.3f, w = %6.4f\n", x1, gg);
            point(x, y, xx, yy);
        }
    }

    setColor(x, y, COLOR(0, 0, 0));
    for(float xx = 0.0; xx <= 1.0; xx += 0x1.0p-10){
        point(x, y, xx, m0*xx);
    }

    axesBotLeft(x, y, 0.0f, 0.0f);
    axesTopRight(x, y, 1.0f, rmx);
    plot(x, y);
}

void plotRcCurve2(float m0, float WB0, float WB1, float np, float rmx, int Nc){
    rcParam2_t rcParams;
    float w0, dWB;
    float r, g, b;

    dWB = (WB1 - WB0) / (Nc - 1);

    for(int i = 0; i < Nc; i++){
        // LumToSurreyRGB(((float) i)/ Ni, &r, &g, &b);
        // LumToViridisRGB(((float) i)/ Ni, &r, &g, &b);
        LumToParulaRGB(((float) i)/ Nc, &r, &g, &b);
        setColor(x, y, COLOR(255*r, 255*g, 255*b));

        w0 = dWB * i + WB0;
        initRcParam2(&rcParams, m0, w0, np, rmx);

        for(float xx = 0.0, yy; xx <= 1.0; xx += 0x1.0p-10){
            yy = calcRcParam2(xx, &rcParams);
            // yy = myPowf(xx, w0+2);
            // printf(" > x1 = %5.3f, w = %6.4f\n", x1, gg);
            point(x, y, xx, yy);
        }
    }

    setColor(x, y, COLOR(0, 0, 0));
    for(float xx = 0.0, yy; xx <= 1.0; xx += 0x1.0p-10){
        yy = rcParams.w2 * (myExp(rcParams.invNp * xx) - 1);

        point(x, y, xx, yy);
    }

    setColor(x, y, COLOR(255, 0, 0));
    for(float xx = 0.0, yy; xx <= 1.0; xx += 0x1.0p-10){
        point(x, y, xx, m0*xx);
    }

    axesBotLeft(x, y, 0.0f, 0.0f);
    axesTopRight(x, y, 1.0f, rmx);
    plot(x, y);
}

void plotPowfError(void){
    rcParam2_t rcParams;
    float w0, dWB, WB1, WB0;
    float r, g, b;
    float v, w;

    int Nc = 10;
    WB0 = 0.05; WB1 = 0.95;
    dWB = (WB1 - WB0) / (Nc - 1);

    for(int i = 0; i < Nc; i++){
        // LumToSurreyRGB(((float) i)/ Ni, &r, &g, &b);
        // LumToViridisRGB(((float) i)/ Ni, &r, &g, &b);
        LumToParulaRGB(((float) i)/ Nc, &r, &g, &b);
        setColor(x, y, COLOR(255*r, 255*g, 255*b));

        w0 = dWB * i + WB0;
        initRcParam2(&rcParams, 0.01, w0, rcParams.np, 1.0);

        for(v = 0.0; v <= 1.0f; v += 0x1.0p-10){
            // w = calcRcParam2(v, &rcParams);
            w = myPowf(v, rcParams.g) - powf(v, rcParams.g);
            point(x, y, v, w);
        }
    }

    plot(x, y);
}

void plotCoSiError(float xStart, float xEnd, float xInc){
    setColor(x, y, COLOR(255, 0, 0));
    for(float xx = xStart, yR, yC, yS; xx <= xEnd; xx += xInc){
        yR = myCoSi(xx, &yC, &yS);

        yC -= cosf(xx);
        yS -= sinf(xx);

        point(x, y, xx, yS);
    }

    setColor(x, y, COLOR(0, 255, 0));
    for(float xx = xStart, yR, yC, yS; xx <= xEnd; xx += xInc){
        yR = myCoSi(xx, &yC, &yS);

        yC -= cosf(xx);
        yS -= sinf(xx);

        point(x, y, xx, yC);
    }

    plot(x, y);
}
