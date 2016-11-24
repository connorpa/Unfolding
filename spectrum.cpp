#include <cmath>
#include "spectrum.h"
using namespace std;

double Spectrum (double *x, double *p)
{
    const double X = *x;
    if (X <= 0) return 0;
    double sum;
    for (short i = 0 ; i < NTERMS ; i++)
        sum += p[i] / pow(X, i);
    return sum;
}

#ifdef __CLING__
#include <TF1.h>
#include <TCanvas.h>
TF1 * spectrum ()
{
    TCanvas * c = new TCanvas ("c", "");
    c->SetGridx();
    c->SetGridy();
    c->SetLogy();
    TF1 * spectrum = new TF1 ("spectrum", Spectrum, 0, 1, 100);
    spectrum->SetParameters(0,0,0,0,1,0,0,0,0,0);
    spectrum->Draw();
    return spectrum;
}
#endif
