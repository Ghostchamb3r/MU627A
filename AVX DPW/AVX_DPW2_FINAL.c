#include "csdl.h"
#include <immintrin.h>
#include <math.h>

#define VECTOR_SIZE 4

//Vesa Valimaki, Juhan Nam, Julius O. Smith and Jonathan S. Abel
//Alias-Suppressed Oscillators Based on Differentiated Polynomial Waveforms
//IEEE Transactions on Audio, Speech, and Language Processing 18(4) May 2010, pp 786--798
// Original SuperCollider implementation at:
// https://github.com/supercollider/sc3-plugins/blob/main/source/AntiAliasingOscillators/AntiAliasingOscillators.cpp
// Under DPW4Saw functions
/*
To Compile:

This is the correct command lol
gcc -O2 -mavx -mavx2 -shared -o AVXDPW.dll AVX_DPW2_FINAL.c -DUSE_DOUBLE -I"C:\Program Files\Csound6_x64\include\csound"

*/

/*
DATA STRUCTURE 
************************************/
typedef struct _AVXDPW {
OPDS h;
MYFLT *out;
MYFLT *in1, *in2;
MYFLT freq;
MYFLT amp;
double differentiations[3];
double differentiations2[2];
double phase;
} AVXDPW;


/*
INITIALIZATION
************************************/
int AVXDPW_init(CSOUND *csound, AVXDPW *p){
        
    p->phase = 0.0;
    for(int i =0; i < 3;++i){
        p->differentiations[i]=0.0;
    }
     for(int i =0; i < 2;++i){
        p->differentiations2[i]=0.0;
    }
    return OK;
}

/*
PROCESSING
************************************/

int AVXDPW_process_audio(CSOUND *csound, AVXDPW *p) {
    MYFLT *aout = p->out;
    MYFLT sar = csound->GetSr(csound);
    MYFLT freq = *p->in1;
    MYFLT amp = *p->in2;

    double phasestep = freq / sar;
    double period = 1.0 / phasestep;
    double phase = p->phase;
    double ampcompensation = 0.0052083333333333 * period * period * period;
    double ampcompensation2 = 0.25 * period;
    int inNumSamples = CS_KSMPS;
    int d;

    double g, oneminusg;
    if (freq > 600.0) {
        g = 1.0;
        oneminusg = 0.0;
    } else if (freq > 400.0) {
        g = (freq - 400.0) * 0.005;
        g *= g;
        oneminusg = 1.0 - g;
    } else {
        g = 0.0;
        oneminusg = 1.0;
    }

    __m256d two = _mm256_set1_pd(2.0);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d ampcompensation_vector = _mm256_set1_pd(ampcompensation);
    __m256d ampcompensation_vector2 = _mm256_set1_pd(ampcompensation2);
    __m256d g_vector = _mm256_set1_pd(g);
    __m256d oneminusg_vector = _mm256_set1_pd(oneminusg);

    __m256d diff1_vector = _mm256_set1_pd(p->differentiations[0]);
    __m256d diff2_vector = _mm256_set1_pd(p->differentiations[1]);
    __m256d diff3_vector = _mm256_set1_pd(p->differentiations[2]);
    __m256d diff4_vector = _mm256_set1_pd(p->differentiations2[0]);
    __m256d diff5_vector = _mm256_set1_pd(p->differentiations2[1]);

    for (d = 0; d <= inNumSamples - VECTOR_SIZE; d += VECTOR_SIZE) {
        // Phase vector
        __m256d phase_vector = _mm256_set_pd(
            phase + 3 * phasestep,
            phase + 2 * phasestep,
            phase + 1 * phasestep,
            phase
        );
        phase_vector = _mm256_sub_pd(phase_vector, _mm256_floor_pd(phase_vector));

        // temp = 2 * phase - 1
        __m256d temp = _mm256_sub_pd(_mm256_mul_pd(two, phase_vector), one);
        // temp * temp
        __m256d temp_squared = _mm256_mul_pd(temp, temp);

        // DPW2 path
        __m256d temp2 = temp_squared;
        __m256d dpw2 = _mm256_sub_pd(temp2, diff4_vector);
        diff4_vector = temp2;
        dpw2 = _mm256_mul_pd(dpw2, ampcompensation_vector2); // Apply amp compensation

        // DPW4 path
        __m256d temp_polymorph = _mm256_sub_pd(_mm256_mul_pd(temp_squared, temp_squared), _mm256_mul_pd(two, temp_squared));

        // 3 differentiations
        temp2 = temp_polymorph;
        temp_polymorph = _mm256_sub_pd(temp2, diff1_vector);
        diff1_vector = temp2;

        temp2 = temp_polymorph;
        temp_polymorph = _mm256_sub_pd(temp2, diff2_vector);
        diff2_vector = temp2;

        temp2 = temp_polymorph;
        temp_polymorph = _mm256_sub_pd(temp2, diff3_vector);
        diff3_vector = temp2;

        temp_polymorph = _mm256_mul_pd(temp_polymorph, ampcompensation_vector);

        // Mix
        __m256d out_mix = _mm256_add_pd(
            _mm256_mul_pd(g_vector, temp_polymorph),
            _mm256_mul_pd(oneminusg_vector, dpw2)
        );

        // Send to output buff
        _mm256_storeu_pd(&aout[d], out_mix);

        // Advance phase
        phase += VECTOR_SIZE * phasestep;
        if (phase > 1.0) phase -= floor(phase);
    }

    double tempArraynge[4];
    _mm256_storeu_pd(tempArraynge, diff1_vector); p->differentiations[0] = tempArraynge[3];
    _mm256_storeu_pd(tempArraynge, diff2_vector); p->differentiations[1] = tempArraynge[3];
    _mm256_storeu_pd(tempArraynge, diff3_vector); p->differentiations[2] = tempArraynge[3];
    _mm256_storeu_pd(tempArraynge, diff4_vector); p->differentiations2[0] = tempArraynge[3];
    _mm256_storeu_pd(tempArraynge, diff5_vector); p->differentiations2[1] = tempArraynge[3];

    p->phase = phase;
    // csound->Message(csound, "frequency: %f Hz\n", freq);
    return OK;

}

/*
REGISTRATION
*************************************/

static OENTRY localops[] = {
    {"AVXDPW", sizeof(AVXDPW), 0, 7, "a", "ii",
     (SUBR) AVXDPW_init, (SUBR) AVXDPW_process_audio,}
};

LINKAGE