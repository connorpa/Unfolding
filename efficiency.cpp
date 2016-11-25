#include "efficiency.h"
#include <cmath>
using namespace std;

double Efficiency (double *x, double *p)
{
    const double pt = *x;
    const double a = p[0], mu = p[1], sigma = p[2];
    if (sigma == 0) return 0;
    return a + 0.5*(1-a)*(1+erf((pt-mu)/sigma));
}

#ifdef __CLING__
#include <TCanvas.h>
#include <TF1.h>
TCanvas * efficiency ()
{
    TCanvas * c = new TCanvas ("c");
    c->SetGridx();
    c->SetGridy();

    double xmin = 30, xmax = 3000, rel_sigma = 0.2;

    c->DrawFrame(xmin, 0, xmax, 1.05, "trigger efficiency");

    TF1 * f = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    f->SetParameters(0.2, xmin, rel_sigma*xmin);
    f->SetLineColor(kRed);
    f->Draw("same");

    TF1 * g = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    g->SetParameters(0.5, xmin, rel_sigma*xmin);
    g->SetLineColor(kBlue);
    g->Draw("same");

    TF1 * h = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    h->SetParameters(0.8, xmin, rel_sigma*xmin);
    h->SetLineColor(kViolet);
    h->Draw("same");


    TF1 * f2 = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    f2->SetParameters(0.2, (xmax-xmin)/2, rel_sigma*(xmax-xmin)/2);
    f2->SetLineColor(kRed);
    f2->SetLineStyle(2);
    f2->Draw("same");

    TF1 * g2 = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    g2->SetParameters(0.5, (xmax-xmin)/2, rel_sigma*(xmax-xmin)/2);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    g2->Draw("same");

    TF1 * h2 = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    h2->SetParameters(0.8, (xmax-xmin)/2, rel_sigma*(xmax-xmin)/2);
    h2->SetLineColor(kViolet);
    h2->SetLineStyle(2);
    h2->Draw("same");


    TF1 * f3 = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    f3->SetParameters(0.2, xmax, rel_sigma*xmax);
    f3->SetLineColor(kRed);
    f3->SetLineStyle(3);
    f3->Draw("same");

    TF1 * g3 = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    g3->SetParameters(0.5, xmax, rel_sigma*xmax);
    g3->SetLineColor(kBlue);
    g3->SetLineStyle(3);
    g3->Draw("same");

    TF1 * h3 = new TF1 ("efficiency", Efficiency, xmin, xmax, 100);
    h3->SetParameters(0.8, xmax, rel_sigma*xmax);
    h3->SetLineColor(kViolet);
    h3->SetLineStyle(3);
    h3->Draw("same");


    return c;
}
#endif
