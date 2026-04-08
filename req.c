/*
 * ReQ — Parametric EQ for PocketDAW
 * 4-band stereo biquad EQ: Low Shelf | Peak | Peak | High Shelf
 *
 * Parameters (0.0–1.0 normalized, 3 per band → matches host band*3+sub mapping):
 *   0  LS Freq     20–500 Hz       (low shelf)
 *   1  LS Gain     -12 to +12 dB
 *   2  LS Slope    0.3–1.0
 *   3  P1 Freq     80–2000 Hz      (peak)
 *   4  P1 Gain     -12 to +12 dB
 *   5  P1 Q        0.3–4.0
 *   6  P2 Freq     800–12000 Hz    (peak)
 *   7  P2 Gain     -12 to +12 dB
 *   8  P2 Q        0.3–4.0
 *   9  HS Freq     2000–20000 Hz   (high shelf)
 *  10  HS Gain     -12 to +12 dB
 *  11  HS Slope    0.3–1.0
 *
 * Build:
 *   gcc -shared -fPIC -O2 -o req.so req.c -lm
 *   aarch64-none-linux-gnu-gcc -shared -fPIC -O2 -o req.so req.c -lm
 *   x86_64-w64-mingw32-gcc -shared -o req.dll req.c -lm
 */

#include <SDL2/SDL.h>      /* for SDL_Renderer, SDL_Rect, SDL_* draw calls */
#include "pocketdaw.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Biquad filter state (Direct Form I) ── */
typedef struct {
    double b0, b1, b2, a1, a2; /* coefficients */
    double x1L, x2L, y1L, y2L; /* left  channel state */
    double x1R, x2R, y1R, y2R; /* right channel state */
} Biquad;

static inline int bad_double(double v) {
    /* Catch NaN, Inf, and extreme values that would corrupt the filter */
    return !(v > -1e15 && v < 1e15);  /* NaN fails all comparisons */
}

static void biquad_process(Biquad* f, const float* inL, const float* inR,
                            float* outL, float* outR, int n) {
    double b0=f->b0, b1=f->b1, b2=f->b2, a1=f->a1, a2=f->a2;

    /* If coefficients are bad, pass through unchanged */
    if (bad_double(b0) || bad_double(a1) || bad_double(a2)) {
        if (inL != outL) for (int i = 0; i < n; i++) outL[i] = inL[i];
        if (inR != outR) for (int i = 0; i < n; i++) outR[i] = inR[i];
        f->x1L=f->x2L=f->y1L=f->y2L=0.0;
        f->x1R=f->x2R=f->y1R=f->y2R=0.0;
        return;
    }

    double x1L=f->x1L, x2L=f->x2L, y1L=f->y1L, y2L=f->y2L;
    double x1R=f->x1R, x2R=f->x2R, y1R=f->y1R, y2R=f->y2R;
    for (int i = 0; i < n; i++) {
        double xL = (double)inL[i];
        double yL = b0*xL + b1*x1L + b2*x2L - a1*y1L - a2*y2L;
        x2L=x1L; x1L=xL; y2L=y1L; y1L=yL;
        outL[i] = (float)yL;

        double xR = (double)inR[i];
        double yR = b0*xR + b1*x1R + b2*x2R - a1*y1R - a2*y2R;
        x2R=x1R; x1R=xR; y2R=y1R; y1R=yR;
        outR[i] = (float)yR;
    }
    /* Sanitize state — if filter went unstable, reset to prevent ongoing clicks */
    if (bad_double(y1L) || bad_double(y2L)) { x1L=x2L=y1L=y2L=0.0; }
    if (bad_double(y1R) || bad_double(y2R)) { x1R=x2R=y1R=y2R=0.0; }
    f->x1L=x1L; f->x2L=x2L; f->y1L=y1L; f->y2L=y2L;
    f->x1R=x1R; f->x2R=x2R; f->y1R=y1R; f->y2R=y2R;
}

/* ── Coefficient computation ── */

/* Low shelf — RBJ cookbook (S = shelf slope, 0.3–1.0; 1.0 = steepest) */
static void coeff_low_shelf(Biquad* f, double sr, double freq, double gainDb, double S) {
    double A  = pow(10.0, gainDb / 40.0);
    double w0 = 2.0 * M_PI * freq / sr;
    double cw = cos(w0), sw = sin(w0);
    if (S < 0.01) S = 0.01;
    double a  = sw / 2.0 * sqrt((A + 1.0/A) * (1.0/S - 1.0) + 2.0);
    double b0 =  A*((A+1.0) - (A-1.0)*cw + 2.0*sqrt(A)*a);
    double b1 =  2.0*A*((A-1.0) - (A+1.0)*cw);
    double b2 =  A*((A+1.0) - (A-1.0)*cw - 2.0*sqrt(A)*a);
    double a0 =  (A+1.0) + (A-1.0)*cw + 2.0*sqrt(A)*a;
    double a1 = -2.0*((A-1.0) + (A+1.0)*cw);
    double a2 =  (A+1.0) + (A-1.0)*cw - 2.0*sqrt(A)*a;
    f->b0=b0/a0; f->b1=b1/a0; f->b2=b2/a0;
    f->a1=a1/a0; f->a2=a2/a0;
}

/* Peaking EQ — RBJ cookbook */
static void coeff_peak(Biquad* f, double sr, double freq, double gainDb, double Q) {
    double A  = pow(10.0, gainDb / 40.0);
    double w0 = 2.0 * M_PI * freq / sr;
    double alpha = sin(w0) / (2.0 * Q);
    double b0 =  1.0 + alpha * A;
    double b1 = -2.0 * cos(w0);
    double b2 =  1.0 - alpha * A;
    double a0 =  1.0 + alpha / A;
    double a1 = -2.0 * cos(w0);
    double a2 =  1.0 - alpha / A;
    f->b0=b0/a0; f->b1=b1/a0; f->b2=b2/a0;
    f->a1=a1/a0; f->a2=a2/a0;
}

/* High shelf — RBJ cookbook (S = shelf slope, 0.3–1.0) */
static void coeff_high_shelf(Biquad* f, double sr, double freq, double gainDb, double S) {
    double A  = pow(10.0, gainDb / 40.0);
    double w0 = 2.0 * M_PI * freq / sr;
    double cw = cos(w0), sw = sin(w0);
    if (S < 0.01) S = 0.01;
    double a  = sw / 2.0 * sqrt((A + 1.0/A) * (1.0/S - 1.0) + 2.0);
    double b0 =  A*((A+1.0) + (A-1.0)*cw + 2.0*sqrt(A)*a);
    double b1 = -2.0*A*((A-1.0) + (A+1.0)*cw);
    double b2 =  A*((A+1.0) + (A-1.0)*cw - 2.0*sqrt(A)*a);
    double a0 =  (A+1.0) - (A-1.0)*cw + 2.0*sqrt(A)*a;
    double a1 =  2.0*((A-1.0) - (A+1.0)*cw);
    double a2 =  (A+1.0) - (A-1.0)*cw - 2.0*sqrt(A)*a;
    f->b0=b0/a0; f->b1=b1/a0; f->b2=b2/a0;
    f->a1=a1/a0; f->a2=a2/a0;
}

/* ── Plugin state ── */
#define NUM_BANDS    4
#define NUM_PARAMS   12
#define REQ_WAVE_LEN 512  /* ring buffer size for waveform display */

typedef struct {
    float sr;
    float p[NUM_PARAMS];    /* normalized params */
    Biquad bands[NUM_BANDS];
    int   dirty;            /* recompute coeffs next process */
    int   selectedNode;     /* for custom draw highlight */
    int   dragging;         /* 1 = node locked during drag */
    int   prevMouseDown;    /* previous frame mouse state — detect press start */
    float dragOffsetX;      /* pixel offset from cursor to node center at grab time */
    float dragOffsetY;

    /* ── Waveform oscilloscope ring buffers (written by pdfx_process) ── */
    float wavePreL [REQ_WAVE_LEN];  /* pre-FX  mono mix (L+R)/2 */
    float wavePostL[REQ_WAVE_LEN];  /* post-FX mono mix */
    volatile int wavePos;           /* next write position (audio thread writes, UI reads) */
    int   showPreFx;                /* 0 = post-FX waveform, 1 = pre-FX */
    char  fmt[64];          /* format_param buffer */

    /* Computed real values (for draw) */
    float bandFreq[NUM_BANDS];
    float bandGain[NUM_BANDS];
    float bandQ[NUM_BANDS];
} ReQ;

/* Map normalized params to real values.
 * Layout: 3 params per band (Freq, Gain, Q/Slope) → matches host band*3+sub mapping.
 *   Band 0 (Low Shelf):  p[0]=Freq  p[1]=Gain  p[2]=Slope
 *   Band 1 (Peak):       p[3]=Freq  p[4]=Gain  p[5]=Q
 *   Band 2 (Peak):       p[6]=Freq  p[7]=Gain  p[8]=Q
 *   Band 3 (High Shelf): p[9]=Freq  p[10]=Gain p[11]=Slope
 */
static void _update_values(ReQ* eq) {
    eq->bandFreq[0] = 20.0f   + eq->p[0] * 480.0f;   /* 20–500 Hz  */
    eq->bandGain[0] = -12.0f  + eq->p[1] * 24.0f;    /* ±12 dB     */
    eq->bandQ[0]    = 0.3f    + eq->p[2] * 0.7f;     /* shelf slope 0.3–1.0 */

    eq->bandFreq[1] = 80.0f   + eq->p[3] * 1920.0f;  /* 80–2000 Hz */
    eq->bandGain[1] = -12.0f  + eq->p[4] * 24.0f;
    eq->bandQ[1]    = 0.3f    + eq->p[5] * 3.7f;     /* 0.3–4.0    */

    eq->bandFreq[2] = 800.0f  + eq->p[6] * 11200.0f; /* 800–12k Hz */
    eq->bandGain[2] = -12.0f  + eq->p[7] * 24.0f;
    eq->bandQ[2]    = 0.3f    + eq->p[8] * 3.7f;

    eq->bandFreq[3] = 2000.0f + eq->p[9] * 18000.0f; /* 2k–20k Hz  */
    eq->bandGain[3] = -12.0f  + eq->p[10] * 24.0f;
    eq->bandQ[3]    = 0.3f    + eq->p[11] * 0.7f;    /* shelf slope 0.3–1.0 */
}

static void _recompute(ReQ* eq) {
    _update_values(eq);
    double sr = (double)eq->sr;
    coeff_low_shelf (&eq->bands[0], sr, eq->bandFreq[0], eq->bandGain[0], eq->bandQ[0]);
    coeff_peak      (&eq->bands[1], sr, eq->bandFreq[1], eq->bandGain[1], eq->bandQ[1]);
    coeff_peak      (&eq->bands[2], sr, eq->bandFreq[2], eq->bandGain[2], eq->bandQ[2]);
    coeff_high_shelf(&eq->bands[3], sr, eq->bandFreq[3], eq->bandGain[3], eq->bandQ[3]);
    eq->dirty = 0;
}

/* ── Required exports ── */

int pdfx_api_version(void) { return 2; }  /* must match PDFX_API_VERSION in host */
const char* pdfx_name(void) { return "ReQ"; }

PdFxInstance pdfx_create(float sr) {
    ReQ* eq = calloc(1, sizeof(ReQ));
    eq->sr = sr;
    /* Defaults: bands centered, all at 0 dB, Q=1 */
    eq->p[0]  = 0.2f;  /* B1 freq  ~116 Hz  */
    eq->p[1]  = 0.5f;  /* B1 gain   0 dB    */
    eq->p[2]  = 1.0f;  /* B1 slope  1.0     */
    eq->p[3]  = 0.2f;  /* B2 freq  ~464 Hz  */
    eq->p[4]  = 0.5f;  /* B2 gain   0 dB    */
    eq->p[5]  = 0.2f;  /* B2 Q     ~1.0     */
    eq->p[6]  = 0.2f;  /* B3 freq  ~3.0 kHz */
    eq->p[7]  = 0.5f;  /* B3 gain   0 dB    */
    eq->p[8]  = 0.2f;  /* B3 Q     ~1.0     */
    eq->p[9]  = 0.5f;  /* B4 freq  ~11 kHz  */
    eq->p[10] = 0.5f;  /* B4 gain   0 dB    */
    eq->p[11] = 1.0f;  /* B4 slope  1.0     */
    eq->selectedNode = -1;
    eq->dragging     = 0;
    eq->prevMouseDown= 0;
    memset(eq->wavePreL,  0, sizeof(eq->wavePreL));
    memset(eq->wavePostL, 0, sizeof(eq->wavePostL));
    eq->wavePos    = 0;
    eq->showPreFx  = 0;
    eq->dirty = 1;
    _recompute(eq);
    return eq;
}

void pdfx_destroy(PdFxInstance inst) { free(inst); }

void pdfx_process(PdFxInstance inst, PdFxAudio* a) {
    ReQ* eq = (ReQ*)inst;
    int i, n;
    int wpos;
    float master;
    float tL[8192], tR[8192];  /* must handle max host buffer (MIX_BUF_SIZE/2) */

    /* Detect sample rate change (host reconfigured audio device) */
    if (a->sampleRate > 0.0f && a->sampleRate != eq->sr) {
        eq->sr = a->sampleRate;
        eq->dirty = 1;
    }
    if (eq->dirty) _recompute(eq);

    master = 1.0f;  /* no dedicated master gain param — band gains control level */
    n = a->bufferSize;

    /* Capture pre-FX audio into ring buffer (mono mix of input) */
    wpos = eq->wavePos;
    for (i = 0; i < n; i++) {
        eq->wavePreL[wpos] = (a->inputL[i] + a->inputR[i]) * 0.5f;
        wpos = (wpos + 1) % REQ_WAVE_LEN;
    }

    /* Chain bands: inputL/R → tL/R → outputL/R */
    biquad_process(&eq->bands[0], a->inputL, a->inputR, tL, tR, n);
    biquad_process(&eq->bands[1], tL, tR, a->outputL, a->outputR, n);
    biquad_process(&eq->bands[2], a->outputL, a->outputR, tL, tR, n);
    biquad_process(&eq->bands[3], tL, tR, a->outputL, a->outputR, n);

    /* Apply master gain + soft-clip to prevent overdriving downstream FX */
    for (i = 0; i < n; i++) {
        float oL = a->outputL[i] * master;
        float oR = a->outputR[i] * master;
        if (oL > 1.0f) oL = 1.0f; else if (oL < -1.0f) oL = -1.0f;
        if (oR > 1.0f) oR = 1.0f; else if (oR < -1.0f) oR = -1.0f;
        a->outputL[i] = oL;
        a->outputR[i] = oR;
    }

    /* Capture post-FX audio and advance write pointer */
    wpos = eq->wavePos;  /* reset to same start position */
    for (i = 0; i < n; i++) {
        eq->wavePostL[wpos] = (a->outputL[i] + a->outputR[i]) * 0.5f;
        wpos = (wpos + 1) % REQ_WAVE_LEN;
    }
    eq->wavePos = wpos;  /* advance once after both buffers written */
}

/* ── Parameters ── */

int pdfx_param_count(void) { return NUM_PARAMS; }

const char* pdfx_param_name(int i) {
    const char* names[] = {
        "LS Freq", "LS Gain", "LS Slope",
        "P1 Freq", "P1 Gain", "P1 Q",
        "P2 Freq", "P2 Gain", "P2 Q",
        "HS Freq", "HS Gain", "HS Slope"
    };
    return (i >= 0 && i < NUM_PARAMS) ? names[i] : NULL;
}

const char* pdfx_param_group(int i) {
    if (i <= 2) return "Band 1 (Low Shelf)";
    if (i <= 5) return "Band 2 (Peak)";
    if (i <= 8) return "Band 3 (Peak)";
    return "Band 4 (High Shelf)";
}

float pdfx_get_param(PdFxInstance inst, int i) {
    return (i >= 0 && i < NUM_PARAMS) ? ((ReQ*)inst)->p[i] : 0.0f;
}

void pdfx_set_param(PdFxInstance inst, int i, float v) {
    ReQ* eq = (ReQ*)inst;
    if (i >= 0 && i < NUM_PARAMS) {
        eq->p[i] = v;
        eq->dirty = 1;
    }
}

const char* pdfx_format_param(PdFxInstance inst, int i) {
    ReQ* eq = (ReQ*)inst;
    int band = i / 3, sub = i % 3;
    _update_values(eq);
    if (sub == 0) { /* Freq */
        float f = eq->bandFreq[band];
        if (f >= 1000.0f) snprintf(eq->fmt, 64, "%.1f kHz", f/1000.0f);
        else snprintf(eq->fmt, 64, "%.0f Hz", f);
        return eq->fmt;
    }
    if (sub == 1) { /* Gain */
        snprintf(eq->fmt, 64, "%+.1f dB", eq->bandGain[band]);
        return eq->fmt;
    }
    if (sub == 2) { /* Q or Slope */
        if (band == 0 || band == 3)
            snprintf(eq->fmt, 64, "S %.0f%%", eq->bandQ[band] * 100.0f);
        else
            snprintf(eq->fmt, 64, "%.2f Q", eq->bandQ[band]);
        return eq->fmt;
    }
    return NULL;
}

void pdfx_get_accent_color(PdFxInstance inst, uint8_t* r, uint8_t* g, uint8_t* b) {
    *r = 0; *g = 255; *b = 65; /* PocketDAW green */
}

/* ── Presets ── */

int pdfx_preset_count(void) { return 5; }

const char* pdfx_preset_name(int i) {
    const char* n[] = {"Flat","Low Cut","Vocal Presence","Bass Boost","Air"};
    return (i >= 0 && i < 5) ? n[i] : NULL;
}

void pdfx_load_preset(PdFxInstance inst, int i) {
    ReQ* eq = (ReQ*)inst;
    /*          B1F  B1G  B1S   B2F  B2G  B2Q   B3F  B3G  B3Q   B4F  B4G  B4S  */
    float pre[5][NUM_PARAMS] = {
        /* Flat */
        {0.2f,0.5f,1.0f, 0.2f,0.5f,0.2f, 0.2f,0.5f,0.2f, 0.5f,0.5f,1.0f},
        /* Low Cut — cut lows, boost presence */
        {0.05f,0.2f,0.5f, 0.15f,0.5f,0.3f, 0.35f,0.65f,0.25f, 0.6f,0.5f,0.5f},
        /* Vocal Presence — dip 300Hz, boost 3kHz */
        {0.2f,0.45f,0.5f, 0.15f,0.35f,0.3f, 0.22f,0.7f,0.4f, 0.6f,0.55f,0.5f},
        /* Bass Boost — low shelf up, high shelf neutral */
        {0.15f,0.7f,0.5f, 0.1f,0.6f,0.2f, 0.2f,0.5f,0.2f, 0.5f,0.5f,0.5f},
        /* Air — high shelf boost */
        {0.2f,0.5f,1.0f, 0.2f,0.5f,0.2f, 0.3f,0.5f,0.2f, 0.65f,0.7f,0.5f},
    };
    if (i >= 0 && i < 5) {
        for (int j = 0; j < NUM_PARAMS; j++) eq->p[j] = pre[i][j];
        eq->dirty = 1;
    }
}

/* ── Custom draw — EQ curve + draggable nodes ── */

/* Convert frequency to X pixel (log scale) */
static int freqToX(float freq, int x, int w) {
    float logMin = log10f(20.0f);
    float logMax = log10f(20000.0f);
    float t = (log10f(freq < 20.0f ? 20.0f : freq) - logMin) / (logMax - logMin);
    return x + (int)(t * w);
}

/* Convert dB gain to Y pixel */
static int gainToY(float dB, int y, int h) {
    float t = (dB + 12.0f) / 24.0f; /* 0=bottom(-12dB), 1=top(+12dB) */
    return y + h - (int)(t * h);
}

/* Compute magnitude response of all 4 bands at a given frequency (dB) */
/* Compute magnitude response for drawing — uses display values (bandFreq/bandGain/bandQ)
 * instead of live biquad coefficients to avoid racing with the audio thread.
 * Recomputes temporary coefficients on the UI thread stack. */
static float responseDb(ReQ* eq, float freq) {
    double sr = (double)eq->sr;
    if (sr < 1.0) sr = 44100.0;
    double w = 2.0 * M_PI * freq / sr;
    double totalDb = 0.0;

    for (int b = 0; b < NUM_BANDS; b++) {
        /* Compute temporary biquad coefficients from display values */
        Biquad tmp;
        memset(&tmp, 0, sizeof(tmp));
        if (b == 0)      coeff_low_shelf (&tmp, sr, eq->bandFreq[b], eq->bandGain[b], eq->bandQ[b]);
        else if (b == 3)  coeff_high_shelf(&tmp, sr, eq->bandFreq[b], eq->bandGain[b], eq->bandQ[b]);
        else              coeff_peak      (&tmp, sr, eq->bandFreq[b], eq->bandGain[b], eq->bandQ[b]);

        double re_n = tmp.b0 + tmp.b1*cos(w) + tmp.b2*cos(2*w);
        double im_n =         -tmp.b1*sin(w) - tmp.b2*sin(2*w);
        double re_d = 1.0   + tmp.a1*cos(w) + tmp.a2*cos(2*w);
        double im_d =         -tmp.a1*sin(w) - tmp.a2*sin(2*w);
        double mag2 = (re_n*re_n + im_n*im_n) / (re_d*re_d + im_d*im_d + 1e-30);
        totalDb += 10.0 * log10(mag2 < 1e-30 ? 1e-30 : mag2);
    }
    return (float)totalDb;
}

static const uint8_t NODE_COLORS[NUM_BANDS][3] = {
    {255, 160,  40},  /* Band 1 — amber      */
    {  0, 200, 255},  /* Band 2 — cyan       */
    {200,  80, 255},  /* Band 3 — purple     */
    {255,  60,  60},  /* Band 4 — red        */
};


/* ─────────────────────────────────────────────────────────
 * Pixel font (3×5 bitmap) — digits, +, -, ., k, H, z, d, B, Q
 * ───────────────────────────────────────────────────────── */
static const uint8_t PIXEL_FONT[19][5] = {
    {0x7,0x5,0x5,0x5,0x7}, /* 0 */ {0x2,0x6,0x2,0x2,0x7}, /* 1 */
    {0x7,0x1,0x7,0x4,0x7}, /* 2 */ {0x7,0x1,0x7,0x1,0x7}, /* 3 */
    {0x5,0x5,0x7,0x1,0x1}, /* 4 */ {0x7,0x4,0x7,0x1,0x7}, /* 5 */
    {0x7,0x4,0x7,0x5,0x7}, /* 6 */ {0x7,0x1,0x1,0x1,0x1}, /* 7 */
    {0x7,0x5,0x7,0x5,0x7}, /* 8 */ {0x7,0x5,0x7,0x1,0x7}, /* 9 */
    {0x0,0x2,0x7,0x2,0x0}, /* + */ {0x0,0x0,0x7,0x0,0x0}, /* - */
    {0x0,0x0,0x0,0x0,0x2}, /* . */ {0x5,0x5,0x6,0x5,0x5}, /* k */
    {0x5,0x5,0x7,0x5,0x5}, /* H */ {0x7,0x1,0x2,0x4,0x7}, /* z */
    {0x7,0x5,0x6,0x5,0x6}, /* d */ {0x7,0x5,0x6,0x5,0x7}, /* B */
    {0x7,0x5,0x7,0x1,0x7}, /* Q = 9 visual */
};

static int req_charToIdx(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c == '+') return 10; if (c == '-') return 11; if (c == '.') return 12;
    if (c == 'k') return 13; if (c == 'H') return 14; if (c == 'z') return 15;
    if (c == 'd') return 16; if (c == 'B') return 17; if (c == 'Q') return 18;
    return -1;
}

static void req_drawStr(SDL_Renderer* r, int px, int py, const char* s,
                        uint8_t cr, uint8_t cg, uint8_t cb) {
    int cx = px;
    for (; *s; s++, cx += 4) {
        int idx = req_charToIdx(*s);
        if (idx < 0) { cx -= 2; continue; }
        for (int row = 0; row < 5; row++) {
            uint8_t bits = PIXEL_FONT[idx][row];
            for (int col = 0; col < 3; col++) {
                if (bits & (4 >> col)) {
                    SDL_SetRenderDrawColor(r, 0, 0, 0, 255);
                    SDL_RenderDrawPoint(r, cx+col+1, py+row+1);
                    SDL_SetRenderDrawColor(r, cr, cg, cb, 255);
                    SDL_RenderDrawPoint(r, cx+col, py+row);
                }
            }
        }
    }
}

static void req_fmtFreq(float p, float fmin, float fmax, char* out, int n) {
    float hz = fmin + p * (fmax - fmin);
    if (hz >= 1000.0f) {
        int w = (int)(hz / 1000.0f), f = (int)((hz / 1000.0f - w) * 10.0f);
        snprintf(out, n, "%d.%dkHz", w, f);
    } else {
        snprintf(out, n, "%.0fHz", hz);
    }
}

static void req_fmtGain(float p, char* out, int n) {
    float db = (p - 0.5f) * 24.0f;
    int w = (int)fabsf(db), f = (int)(fabsf(db) * 10.0f) % 10;
    if (db >= 0) snprintf(out, n, "+%d.%ddB", w, f);
    else         snprintf(out, n, "-%d.%ddB", w, f);
}

static void req_fmtQ(float p, char* out, int n) {
    float q = 0.3f + p * 3.7f;
    snprintf(out, n, "Q%.1f", q);
}

/* Shelf slope label — shared by bands 0 and 3 */
static void req_fmtSlope(float p, char* out, int n) {
    float s = 0.3f + p * 0.7f;
    /* Display as 0–100% steepness */
    int pct = (int)(s * 100.0f + 0.5f);
    snprintf(out, n, "SL%d", pct);
}

/* Per-band param indices and freq ranges — matches 3-per-band layout */
static const int  REQ_FREQ_P[4] = {0, 3, 6, 9};
static const int  REQ_GAIN_P[4] = {1, 4, 7, 10};
static const int  REQ_Q_P[4]    = {2, 5, 8, 11};
static const float REQ_FMIN[4]  = {20.f,  80.f,  800.f, 2000.f};
static const float REQ_FMAX[4]  = {500.f, 2000.f,12000.f,20000.f};

/* Bottom panel: band tabs + Freq/Gain/Q sliders */
#define REQ_PANEL_H 46
static void req_drawPanel(SDL_Renderer* r, const ReQ* eq, int pw, int ph,
                           int selectedParam, int editMode) {
    int panelY = ph - REQ_PANEL_H;
    /* Background */
    SDL_SetRenderDrawColor(r, 8, 10, 18, 255);
    SDL_Rect bg2 = {0, panelY, pw, REQ_PANEL_H};
    SDL_RenderFillRect(r, &bg2);
    SDL_SetRenderDrawColor(r, 35, 35, 55, 255);
    SDL_RenderDrawLine(r, 0, panelY, pw-1, panelY);

    /* Band tabs */
    int tabH = 14, tabY = panelY + 2, tabW = pw / NUM_BANDS;
    int b;
    for (b = 0; b < NUM_BANDS; b++) {
        int tx = b * tabW;
        int sel = (b == eq->selectedNode);
        uint8_t cr = NODE_COLORS[b][0], cg = NODE_COLORS[b][1], cb = NODE_COLORS[b][2];
        SDL_SetRenderDrawColor(r, sel ? cr/4 : 16, sel ? cg/4 : 16, sel ? cb/4 : 24, 255);
        SDL_Rect tr2 = {tx+1, tabY, tabW-2, tabH};
        SDL_RenderFillRect(r, &tr2);
        SDL_SetRenderDrawColor(r, cr, cg, cb, sel ? 255 : 100);
        SDL_Rect strip2 = {tx+1, tabY, tabW-2, 3};
        SDL_RenderFillRect(r, &strip2);
        /* Band label */
        const char* lbl4[] = {"L","B1","B2","H"};
        req_drawStr(r, tx + tabW/2 - 4, tabY + 5, lbl4[b],
                    sel ? cr : cr/2, sel ? cg : cg/2, sel ? cb : cb/2);
        if (sel) {
            SDL_SetRenderDrawColor(r, cr, cg, cb, 180);
            SDL_RenderDrawRect(r, &tr2);
        }
    }

    /* Freq/Gain/Q sliders for selected band */
    int sb = eq->selectedNode;
    if (sb < 0 || sb >= NUM_BANDS) return;
    int barY = tabY + tabH + 3;
    int barH = panelY + REQ_PANEL_H - barY - 3;
    if (barH < 4) return;

    int mx2 = 3;
    int bw2 = (pw - mx2*2 - 4) / 3;
    uint8_t bcr = NODE_COLORS[sb][0], bcg = NODE_COLORS[sb][1], bcb = NODE_COLORS[sb][2];

    int i;
    for (i = 0; i < 3; i++) {
        int bx2 = mx2 + i * (bw2 + 2);
        int qIdx = REQ_Q_P[sb];
        int pIdx = (i == 0) ? REQ_FREQ_P[sb] : (i == 1) ? REQ_GAIN_P[sb] : (qIdx >= 0 ? qIdx : -1);
        float val = (pIdx >= 0) ? eq->p[pIdx] : 0.5f;

        /* Track */
        SDL_SetRenderDrawColor(r, 22, 22, 32, 255);
        SDL_Rect trk = {bx2, barY, bw2, barH};
        SDL_RenderFillRect(r, &trk);

        if (1) { /* all slots active now (shelf bands use slope slider) */
            /* Fill */
            SDL_SetRenderDrawColor(r, bcr*2/3, bcg*2/3, bcb*2/3, 200);
            if (i == 1) {
                int mid = bx2 + bw2/2;
                int fill = (int)((val - 0.5f) * bw2);
                SDL_Rect fr = fill >= 0 ? (SDL_Rect){mid, barY+1, fill, barH-2}
                                        : (SDL_Rect){mid+fill, barY+1, -fill, barH-2};
                SDL_RenderFillRect(r, &fr);
                SDL_SetRenderDrawColor(r, 40, 40, 60, 255);
                SDL_RenderDrawLine(r, mid, barY, mid, barY+barH-1);
            } else {
                int fw = (int)(val * (bw2-2)); if (fw < 1) fw = 1;
                SDL_Rect fr = {bx2+1, barY+1, fw, barH-2};
                SDL_RenderFillRect(r, &fr);
            }
            /* Edge accent */
            SDL_SetRenderDrawColor(r, bcr, bcg, bcb, 255);
            int ex = (i == 1) ? bx2 + bw2/2 + (int)((val-0.5f)*bw2)
                               : bx2 + 1 + (int)(val*(bw2-2));
            ex = ex < bx2+1 ? bx2+1 : (ex > bx2+bw2-2 ? bx2+bw2-2 : ex);
            SDL_RenderDrawLine(r, ex, barY+1, ex, barY+barH-2);
            /* Label */
            char vlbl[14];
            if (i == 0)      req_fmtFreq(val, REQ_FMIN[sb], REQ_FMAX[sb], vlbl, sizeof(vlbl));
            else if (i == 1) req_fmtGain(val, vlbl, sizeof(vlbl));
            /* Q slot: use slope label for shelf bands (0 and 3), Q label for peak bands */
            else if (sb == 0 || sb == 3) req_fmtSlope(val, vlbl, sizeof(vlbl));
            else             req_fmtQ(val, vlbl, sizeof(vlbl));
            int vw = (int)strlen(vlbl)*4 - 1;
            int vx2 = bx2 + (bw2 - vw)/2, vy2 = barY + (barH-5)/2;
            if (vy2 < barY) vy2 = barY;
            req_drawStr(r, vx2, vy2, vlbl, bcr, bcg, bcb);
        }
        /* Border */
        SDL_SetRenderDrawColor(r, 45, 45, 65, 255);
        SDL_RenderDrawRect(r, &trk);

        /* Gamepad highlight: bright border on the selected sub-parameter */
        if (selectedParam >= 0 && sb == selectedParam / 3 && i == selectedParam % 3) {
            SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_BLEND);
            SDL_SetRenderDrawColor(r, 255, 255, 255, editMode == 1 ? 200 : 100);
            SDL_Rect hlr = {bx2-1, barY-1, bw2+2, barH+2};
            SDL_RenderDrawRect(r, &hlr);
            SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_NONE);
        }
    }
}

int pdfx_draw(PdFxInstance inst, void* renderer, PdDrawContext* ctx) {
    SDL_Renderer* r = (SDL_Renderer*)renderer;
    ReQ* eq = (ReQ*)inst;

    int x = ctx->x, y = ctx->y, w = ctx->w, h = ctx->h;
    /* Reserve bottom for panel */
    int graphH = h - REQ_PANEL_H;
    int cy = y + graphH / 2;

    /* ── Sync gamepad selection: host provides selectedParam via D-pad ── */
    if (ctx->selectedParam >= 0 && ctx->selectedParam < NUM_BANDS * 3) {
        int gpBand = ctx->selectedParam / 3;
        if (gpBand != eq->selectedNode && !ctx->mouseDown) {
            eq->selectedNode = gpBand;
        }
    }

    /* ── INPUT FIRST: process mouse before drawing so visuals reflect current state ── */
    {
        const int* nodeFreqParamI = REQ_FREQ_P; /* {0, 3, 6, 9} */
        const int* nodeGainParamI = REQ_GAIN_P; /* {1, 4, 7, 10} */
        int freshPress = ctx->mouseDown && !eq->prevMouseDown;

        if (ctx->mouseDown) {
            float mx = (float)ctx->mouseX;
            float my = (float)ctx->mouseY;

            int inPanel = (my >= (float)(h - REQ_PANEL_H));

            if (freshPress) {
                if (inPanel) {
                    /* Tab click: select band */
                    int panelY2 = h - REQ_PANEL_H;
                    int tabH2 = 14, tabY2 = panelY2 + 2, tabW2 = w / NUM_BANDS;
                    int b2;
                    for (b2 = 0; b2 < NUM_BANDS; b2++) {
                        int tx2 = b2 * tabW2;
                        if (mx >= (float)tx2 && mx < (float)(tx2 + tabW2) &&
                            my >= (float)tabY2 && my < (float)(tabY2 + tabH2)) {
                            eq->selectedNode = b2;
                            break;
                        }
                    }
                    /* Slider drag — mark dragging with encoded slider index */
                    int barY2 = tabY2 + tabH2 + 3;
                    if (my >= (float)barY2 && eq->selectedNode >= 0) {
                        int bw2 = (w - 3*2 - 4) / 3;
                        int si;
                        for (si = 0; si < 3; si++) {
                            int bx2 = 3 + si * (bw2 + 2);
                            if (mx >= (float)bx2 && mx < (float)(bx2 + bw2)) {
                                eq->dragging = 10 + si; /* 10/11/12 = freq/gain/Q slider */
                                break;
                            }
                        }
                    }
                } else if (mx >= (float)x && mx < (float)(x + 30) &&
                           my >= (float)(y + 2) && my < (float)(y + 13)) {
                    /* PRE/PST toggle button click */
                    eq->showPreFx = !eq->showPreFx;
                } else {
                    /* Graph: find closest node, record grab offset */
                    float bestDist = 900.0f;
                    int   bestNode = -1;
                    int b3;
                    if (eq->dirty) _recompute(eq); /* ensure positions are current */
                    for (b3 = 0; b3 < NUM_BANDS; b3++) {
                        int nx2 = freqToX(eq->bandFreq[b3], x, w);
                        int ny2 = gainToY(eq->bandGain[b3], y, graphH);
                        float dx = mx - (float)nx2;
                        float dy = my - (float)ny2;
                        float dist = dx*dx + dy*dy;
                        if (dist < bestDist) { bestDist = dist; bestNode = b3; }
                    }
                    if (bestNode >= 0) {
                        eq->selectedNode = bestNode;
                        eq->dragging = 1; /* graph drag */
                        /* Record offset from cursor to node center so drag doesn't jump */
                        int nx3 = freqToX(eq->bandFreq[bestNode], x, w);
                        int ny3 = gainToY(eq->bandGain[bestNode], y, graphH);
                        eq->dragOffsetX = mx - (float)nx3;
                        eq->dragOffsetY = my - (float)ny3;
                    }
                }
            }

            /* Continue graph drag — compensate for grab offset */
            if (eq->dragging == 1 && eq->selectedNode >= 0) {
                int b4 = eq->selectedNode;
                float ex = mx - eq->dragOffsetX; /* effective x = cursor minus grab offset */
                float ey = my - eq->dragOffsetY;
                float t = (ex - (float)x) / (float)w;
                t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
                eq->p[nodeFreqParamI[b4]] = t;
                float gt = 1.0f - (ey - (float)y) / (float)graphH;
                gt = gt < 0.0f ? 0.0f : (gt > 1.0f ? 1.0f : gt);
                eq->p[nodeGainParamI[b4]] = gt;
                eq->dirty = 1;
            }
            /* Panel slider drag */
            if (eq->dragging >= 10 && eq->selectedNode >= 0) {
                int si2 = eq->dragging - 10;
                int sb2 = eq->selectedNode;
                int bw3 = (w - 3*2 - 4) / 3;
                int bx3 = 3 + si2 * (bw3 + 2);
                float sv = (mx - (float)bx3) / (float)bw3;
                sv = sv < 0.0f ? 0.0f : (sv > 1.0f ? 1.0f : sv);
                if (si2 == 1) sv = (mx - ((float)bx3 + bw3*0.5f)) / (float)bw3 + 0.5f;
                sv = sv < 0.0f ? 0.0f : (sv > 1.0f ? 1.0f : sv);
                int pIdx2 = (si2 == 0) ? REQ_FREQ_P[sb2]
                          : (si2 == 1) ? REQ_GAIN_P[sb2]
                          : (REQ_Q_P[sb2] >= 0 ? REQ_Q_P[sb2] : -1);
                if (pIdx2 >= 0) { eq->p[pIdx2] = sv; eq->dirty = 1; }
            }
        } else {
            /* Mouse released — unlock */
            eq->dragging = 0;
        }

        eq->prevMouseDown = ctx->mouseDown;
    }

    /* ── Update display values so curve/nodes reflect this frame's changes ── */
    /* Only update the mapped freq/gain/Q values for drawing — do NOT recompute
     * biquad coefficients here.  The audio thread's pdfx_process handles _recompute
     * safely.  Calling _recompute from both threads races on the biquad state and
     * can produce NaN coefficients that crash responseDb → gainToY → SDL draw. */
    _update_values(eq);

    /* ── DRAW with fully up-to-date state ── */

    /* Background */
    SDL_SetRenderDrawColor(r, 8, 8, 8, 255);
    SDL_Rect bg = {x, y, w, graphH};
    SDL_RenderFillRect(r, &bg);

    /* ── Waveform oscilloscope (drawn over background, under grid) ── */
    {
        const float* waveBuf = eq->showPreFx ? eq->wavePreL : eq->wavePostL;
        int wstart = eq->wavePos;       /* oldest sample = current write pos */
        int prevWy = -1;
        int px;
        /* Pre=amber, Post=teal */
        uint8_t wcr = eq->showPreFx ? 200 : 0;
        uint8_t wcg = eq->showPreFx ?  90 : 180;
        uint8_t wcb = eq->showPreFx ?   0 :  90;
        SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(r, wcr, wcg, wcb, 55);
        for (px = 0; px < w; px++) {
            int idx = (wstart + px * REQ_WAVE_LEN / w) % REQ_WAVE_LEN;
            float s = waveBuf[idx];
            /* Clamp and scale: max swing = 35% of graphH */
            if (s >  1.0f) s =  1.0f;
            if (s < -1.0f) s = -1.0f;
            int wy = cy - (int)(s * (float)graphH * 0.35f);
            if (wy < y)          wy = y;
            if (wy > y+graphH-1) wy = y + graphH - 1;
            if (prevWy >= 0 && abs(wy - prevWy) <= graphH/2)
                SDL_RenderDrawLine(r, x + px - 1, prevWy, x + px, wy);
            prevWy = wy;
        }
        SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_NONE);
    }

    /* ── PRE / PST toggle button (top-left of graph) ── */
    {
        int btnX = x + 2, btnY = y + 2, btnW = 28, btnH = 11;
        /* Button background: amber when PRE active, teal when POST */
        uint8_t bfr = eq->showPreFx ? 180 :  20;
        uint8_t bfg = eq->showPreFx ?  70 : 130;
        uint8_t bfb = eq->showPreFx ?   0 :  80;
        SDL_SetRenderDrawColor(r, bfr/3, bfg/3, bfb/3, 255);
        SDL_Rect btnBg = {btnX, btnY, btnW, btnH};
        SDL_RenderFillRect(r, &btnBg);
        SDL_SetRenderDrawColor(r, bfr, bfg, bfb, 255);
        SDL_RenderDrawRect(r, &btnBg);
        req_drawStr(r, btnX + 4, btnY + 3,
                    eq->showPreFx ? "PRE" : "PST",
                    bfr, bfg, bfb);
    }

    /* Grid lines — dB markers */
    {
        int db;
        for (db = -9; db <= 9; db += 3) {
            int gy = gainToY((float)db, y, graphH);
            SDL_SetRenderDrawColor(r, db == 0 ? 50 : 22, db == 0 ? 50 : 22, db == 0 ? 50 : 22, 255);
            SDL_RenderDrawLine(r, x, gy, x+w-1, gy);
        }
    }
    /* Grid lines — frequency markers */
    {
        float freqMarks[] = {50,100,200,500,1000,2000,5000,10000};
        int f2;
        for (f2 = 0; f2 < 8; f2++) {
            int gx = freqToX(freqMarks[f2], x, w);
            SDL_SetRenderDrawColor(r, 22, 22, 22, 255);
            SDL_RenderDrawLine(r, gx, y, gx, y+graphH);
        }
    }

    /* EQ response curve */
    {
        int prevPy = -1;
        int px;
        SDL_SetRenderDrawColor(r, 0, 255, 65, 200);
        for (px = 0; px < w; px++) {
            float t    = (float)px / (float)(w - 1);
            float freq = powf(10.0f, 1.301f + t * 2.699f);
            float db   = responseDb(eq, freq);
            db = db < -20.0f ? -20.0f : (db > 20.0f ? 20.0f : db);
            int py = gainToY(db, y, graphH);
            if (prevPy >= 0) SDL_RenderDrawLine(r, x+px-1, prevPy, x+px, py);
            prevPy = py;
        }
        /* Fill under curve */
        for (px = 0; px < w; px++) {
            float t    = (float)px / (float)(w - 1);
            float freq = powf(10.0f, 1.301f + t * 2.699f);
            float db   = responseDb(eq, freq);
            db = db < -20.0f ? -20.0f : (db > 20.0f ? 20.0f : db);
            int py = gainToY(db, y, graphH);
            SDL_SetRenderDrawColor(r, 0, 255, 65, 18);
            if (py < cy) SDL_RenderDrawLine(r, x+px, py, x+px, cy);
            else         SDL_RenderDrawLine(r, x+px, cy, x+px, py);
        }
    }

    /* Band nodes */
    {
        int b;
        for (b = 0; b < NUM_BANDS; b++) {
            int nx = freqToX(eq->bandFreq[b], x, w);
            int ny = gainToY(eq->bandGain[b], y, graphH);
            uint8_t cr = NODE_COLORS[b][0];
            uint8_t cg = NODE_COLORS[b][1];
            uint8_t cb = NODE_COLORS[b][2];
            int selected = (eq->selectedNode == b);
            int rs = selected ? 12 : 10;
            int a;
            SDL_SetRenderDrawColor(r, cr, cg, cb, selected ? 200 : 80);
            for (a = 0; a < 360; a += (selected ? 10 : 15)) {
                double rad = a * M_PI / 180.0;
                SDL_RenderDrawPoint(r, nx + (int)(cos(rad)*rs), ny + (int)(sin(rad)*rs));
            }
            SDL_SetRenderDrawColor(r, cr, cg, cb, 255);
            SDL_Rect dot = {nx-4, ny-4, 8, 8};
            SDL_RenderFillRect(r, &dot);
        }
    }

    /* Panel — drawn last so it reflects current selectedNode */
    req_drawPanel(r, eq, w, h, ctx->selectedParam, ctx->editMode);

    return 1;
}
