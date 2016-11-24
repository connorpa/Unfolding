#include "efficiency.h"
#include <cmath>
using namespace std;

double Efficiency (double *x, double *p)
{
    const double pt = *x;
    const double a = p[0], mu = p[1], sigma = p[2];
    return a + 0.5*(1-a)*(1+erf((pt-mu)/sigma));
}

#ifdef __CLING__
TF1 * efficiency ()
{
    TCanvas * c = new TCanvas ("c", "");
    c->SetGridx();
    c->SetGridy();
    c->SetLogy();
    TF1 * efficiency = new TF1 ("efficiency", Efficiency, 0, 1, 100);
    efficiency->SetParameters(0.5, 0, 1);
    efficiency->Draw();
    return spectrum;
}
#endif
