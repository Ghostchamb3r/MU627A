#include "csdl.h"
#include <math.h>

/*
To Compile:

This is the correct command lol
gcc -O2 -dynamiclib -o hyperhermiteosc.dylib hyperhermiteosc.c -DUSE_DOUBLE -I/Library/Frameworks/CsoundLib64.framework/Versions/6.0/Headers -arch x86_64
gcc -O2 -mavx -mavx2 -shared -o hyperhermite.dll hyperhermiteosc.c -DUSE_DOUBLE -I"C:\Program Files\Csound6_x64\include\csound"

*/
/*
Might need to manually define pi
const float pi = 3.141592653589793;
*/
// implement tanh as a table lookup instead of calling tanh every sample
// Compute tanh

/*
Hermite interpolation reconfigured from:
https://people.sc.fsu.edu/~jburkardt/c_src/hermite_interpolant/hermite_interpolant.c
*/

#define TANH_TABLE_SIZE 1024  // Adjust for desired resolution
#define TANH_MIN -1.0  // Range mapping
#define TANH_MAX 1.0  

static MYFLT tanhTable[TANH_TABLE_SIZE];
static MYFLT tanhDerivTable[TANH_TABLE_SIZE]; // Store tanh' (x)

void generateTanhTable() { // we can possibly use AVX here to precompute the tables faster
    for (int i = 0; i < TANH_TABLE_SIZE; i++) {
        MYFLT x = TANH_MIN + (i * (TANH_MAX - TANH_MIN) / (TANH_TABLE_SIZE - 1));
        tanhTable[i] = tanh(x);
        tanhDerivTable[i] = 1 - tanhTable[i] * tanhTable[i]; // tanh'(x) = 1 - tanh^2(x)
    }
}


static MYFLT computeTanhHermite(MYFLT x) {
    if (x <= TANH_MIN) return tanhTable[0];
    if (x >= TANH_MAX) return tanhTable[TANH_TABLE_SIZE - 1];

    MYFLT normIndex = ((x - TANH_MIN) / (TANH_MAX - TANH_MIN)) * (TANH_TABLE_SIZE - 1); //2.98485588666341e-4
    int index = (int)normIndex; // 2
    MYFLT frac = normIndex - index; // .98485588666341e-4
// function values f0, f1
    MYFLT f0 = tanhTable[index];
    MYFLT f1 = tanhTable[index + 1];
// derivative values d0, d1
    MYFLT d0 = tanhDerivTable[index];
    MYFLT d1 = tanhDerivTable[index + 1];

    // Cubic Hermite interpolation formula
    // basis functions
    MYFLT h00 = (1 + 2 * frac) * (1 - frac) * (1 - frac);  // 1.96971177332682e-4 * 0.9999015144113337 * 0.9999015144113337
    MYFLT h10 = frac * (1 - frac) * (1 - frac);
    MYFLT h01 = frac * frac * (3 - 2 * frac);
    MYFLT h11 = frac * frac * (frac - 1);

    return h00 * f0 + h10 * d0 + h01 * f1 + h11 * d1;
}




// DATA
typedef struct _hyperhermite{
OPDS h;
MYFLT *out;
MYFLT *in1, *in2, *in3, *in4;
MYFLT freq;
MYFLT amp;
MYFLT phase;
MYFLT m;
} hyperhermite;

// INITIALIZE 

int hyperhermite_init(CSOUND *csound, hyperhermite *p){
    p->phase = 0.0f;
    p->m = *p->in4;
    if (p->m < 0.0) p->m = 0.0;
    if (p->m > 1.0) p->m = 1.0;
    generateTanhTable();  // Precompute tanh values

return OK;
}

// PROCESS

int hyperhermite_process_audio(CSOUND *csound, hyperhermite *p){

int inNumSamples = CS_KSMPS;
MYFLT *aout = p->out;
MYFLT freq = *p->in1;
MYFLT amp = *p->in2;
MYFLT phase = *p->in3;
MYFLT m = *p->in4;
MYFLT sr = csound->GetSr(csound);

float phaseIncrement = (2 * M_PI * freq) / sr;

for(int i = 0; i < inNumSamples; i++){
    float sinewave = amp * sin(p->phase + phase);
    float k = 12000 / (freq * log10(freq));
    float shapedwave = computeTanhHermite(k * sinewave);
    float ringmodder = shapedwave * cos(p->phase);
    float sawtoothwave = shapedwave * (1 - p->m/2) + ringmodder * p->m;

    aout[i] = sawtoothwave;
    p->phase += phaseIncrement;
    if (p->phase > 2 * M_PI) p->phase -= 2 * M_PI;
}

return OK;
}

// REGISTER

static OENTRY localops[] = {
{ "hyperhermite", sizeof(hyperhermite), 0, 7, "a", "iiii",(SUBR) hyperhermite_init,
(SUBR) hyperhermite_process_audio }
};


LINKAGE
