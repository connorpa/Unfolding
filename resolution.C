#include <cmath>
using namespace std;

double Resolution (double *x, double *p)
{
    static const double sqrtLn4 = 1.17741002251547466;
    double X = x[0];
    double N = p[0],
           mu = p[1],
           sigma = p[2],
           tau = p[3], // Novosibirsk argument
           kL = p[4],
           kR = p[5],
           aL = p[6],
           nL = p[7],
           aR = p[8],
           nR = p[9];
    if (sigma <= 0) return 0;
    bool kNovosibirsk = (tau >= 1e-7);

    auto z = [&kNovosibirsk,&mu,&sigma,&tau](double x)
    {    
        if (kNovosibirsk)
        {
            const double Lambda = asinh( tau * sqrtLn4 ) / ( sigma*tau*sqrtLn4 );
            return sqrt(pow( log(1+Lambda*tau*(x-mu))/tau ,2) + pow(tau,2));
        }
        else
            return (x-mu)/sigma;
    };   

    if (X <= aL)
        return N * exp(  pow(z(kL),2)/2 - z(kL)*z(aL) ) * pow(-nL/z(kL),nL) * pow(z(aL) - nL/z(kL) - z(X), -nL) ;
    else if (X <= kL)
        return N * exp(  pow(z(kL),2)/2 - z(kL)*z( X) );
    else if (X <= kR)
        return N * exp( -pow(z( X),2)/2 );
    else if (X <= aR)
        return N * exp(  pow(z(kR),2)/2 - z(kR)*z( X) );
    else
        return N * exp(  pow(z(kR),2)/2 - z(kR)*z(aR) ) * pow(-nR/z(kR),nR) * pow(z(aR) - nR/z(kR) - z(X), -nR) ;
}

#ifdef __CLING__
TF1 * resolution ()
{
    TCanvas * c = new TCanvas ("c", "");
    c->Divide(2,1);
    c->cd(1);
    c->GetPad(1)->DrawFrame(-1, 0, 1, 1.1, "linear scale");
    c->GetPad(1)->SetGridx();
    c->GetPad(1)->SetGridy();
    c->cd(2);
    c->GetPad(2)->DrawFrame(-1, 1e-4, 1, 1e1, "log scale");
    c->GetPad(2)->SetGridx();
    c->GetPad(2)->SetGridy();
    c->GetPad(2)->SetLogy();
    TF1 * core       = new TF1 ("core"     , Resolution, -1, 1, 10),
        * transition = new TF1 ("deviation", Resolution, -1, 1, 10),
        * CB_tail    = new TF1 ("tail"     , Resolution, -1, 1, 10);
    //                         N  mu sigma tau   kL     kR    aL  nL   aR nR
    core      ->SetParameters( 1,  0, 0.1, 0.05,   -1,    1,   -1,  0,   1, 0);
    transition->SetParameters( 1,  0, 0.1, 0.05,-0.01, 0.05,   -1,  0,   1, 0);
    CB_tail   ->SetParameters( 1,  0, 0.1, 0.05,-0.01, 0.05,-0.20,  1,0.30, 5);
    
    core      ->SetLineColor(kRed   );
    transition->SetLineColor(kBlue  );
    CB_tail   ->SetLineColor(kViolet);

    c->cd(1);
    core      ->DrawCopy("same");
    transition->DrawCopy("same");
    CB_tail   ->DrawCopy("same");

    c->cd(2);
    core      ->DrawCopy("same");
    transition->DrawCopy("same");
    CB_tail   ->DrawCopy("same");

    return CB_tail;
}
#endif
