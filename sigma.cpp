#include <cmath>
using namespace std;

double Sigma (double *x, double *p)
{
    const double pt = *x;
    if (p[1] == 0)
        return p[0] + p[4]*pt;
    else 
        return p[0] + p[1]/(pow(pt,p[2]) + p[3]*pt) + p[4]*pt;
}

#ifdef __CLING__
TF1 * spectrum ()
{
    TCanvas * c = new TCanvas ("c", "");
    c->SetGridx();
    c->SetGridy();
    c->SetLogy();
    TF1 * sigma = new TF1("sigma", Sigma, 0.1, fit_ptmax);
    sigma->SetParameters(0.1,0.1,0.1,0.1,0.1);
    sigma->Draw();
    return spectrum;
}
#endif
