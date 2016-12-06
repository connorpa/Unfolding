#ifndef __HEADER__
#define __HEADER__

#include <map>
#include <string>
#include <vector>

#include <TUnfold.h>

class TH1;
class TH2;
class TCanvas;
class TString;

struct UnfoldingParameters
{
    const enum Type { kBinByBin, kInversion, kSVD, kBayes, kTUnfold } type;
    UnfoldingParameters(Type t) : type(t) {}
};
struct UnfoldingParametersSVD : public UnfoldingParameters
{
    int kreg;
    UnfoldingParametersSVD (Type t, int k) : UnfoldingParameters(t), kreg(k) {}
};
struct UnfoldingParametersBayes : public UnfoldingParameters
{
    int niterations;
    const enum RMtype { kSquare, kFine } rm_type;
    UnfoldingParametersBayes (Type t, int ni, RMtype rt) : UnfoldingParameters(t), niterations(ni), rm_type(rt) {}
};
struct UnfoldingParametersTUnfold : public UnfoldingParameters
{
    TUnfold::ERegMode regularisation;
    //double tau; // TODO?
    UnfoldingParametersTUnfold (Type t, TUnfold::ERegMode r) : UnfoldingParameters(t), regularisation(r) {}
};


double Spectrum   (double *, double *);
double Resolution (double *, double *);
double Efficiency (double *, double *);
double Sigma      (double *, double *);

void make_RM (TH2 *, TH2 *, TF1 *, double, double, TF1 *, TF1 *, unsigned long, TString);
void make_measurement (TH1 *, TH1 *, TF1 *, /*TF1 *,*/ double, double, TF1 *, TF1 *, unsigned long); 
std::vector<double> find_binning (TH2 *, float, float);
std::vector<TH1 *> make_unfolding (TH2 *, TH1 *, TH2 *, TH1 *);
std::pair<TH1 *, TH1 *> make_miss_fake (TH2 *, TH1 *, TH1 *);
std::vector<TH1 *> make_ABPS (TH2 *);
TCanvas * make_canvas (TH1 *,  TH1 *, TH2 *, TH2 *, std::vector<double>, TString, TString, std::vector<UnfoldingParameters *>, std::vector<TString>, std::vector<TString>, std::vector<TString>, double turnon);

#endif
