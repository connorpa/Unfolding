/********************************************************************************
 * UNFOLDING                                                                    *
 *                                                                              *
 * The aim of the program is to play with RooUnfold and                         *
 * the different unfolding methods                                              *
 * under different conditions, including binning, resolution, etc.              *
 * but NOT to unfold your data!                                                 *
 *                                                                              *
 * How to use it?                                                               *
 *  - make sure to use Root6                                                    *
 *  - download the latest RooUnfold and compile it                              *
 *  - adapt the makefile                                                        *
 *  - browse through the TODOs and adapt the code according to your needs       *
 *  - make it and run it!                                                       *
 *                                                                              *
 * To download RooUnfold:                                                       *
 * > svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk RooUnfold   *
 *                                                                              *
 * Contact: Patrick L.S. Connor (patrick.connor@desy.de)                        *
 ********************************************************************************/

#include <map>
#include <vector>
#include <cstdlib>
// Root
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TString.h>
#include <TLatex.h>
#include <TColor.h>
#include <TStyle.h>
// RooUnfold
#include <RooUnfold.h>
#include <RooUnfoldBinByBin.h>
#include <RooUnfoldErrors.h>
#include <RooUnfoldParms.h>
#include <RooUnfoldSvd.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldDagostini.h>
#include <RooUnfoldInvert.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldTUnfold.h>

using namespace std;

TCanvas * divider (const vector<TString> lines, Color_t fillcolor = kGray)
{
    TCanvas * c = new TCanvas("divider", "Unfolding");
    c->Draw();
    c->cd();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * pave  = new TPaveText(0.2, 0.2 , 0.8, 0.8, "NBNDC");
    pave->SetFillColor(fillcolor);
    pave->SetTextFont(42);
    for (const TString& line: lines)
        pave->AddText(line);
    pave->Draw();
    return c;
}
TCanvas * divider (const TString line, Color_t fillcolor = kGray)
{
    vector<TString> lines = {line};
    return divider(lines, fillcolor);
}

double flat_spectrum(double pt)
{
    return 1;
}
double falling_spectrum(double pt)
{
    if (pt == 0) return 0;
    return pow(pt,-4);
}
double even_more_falling_spectrum(double pt)
{
    if (pt == 0) return 0;
    return pow(pt,-6);
}

//double resolution (double *x, double *p)
//{
//    // TODO: add tails to and deviations from gaussian core
//    double mu = p[0],
//           sigma = p[1];
//    if (sigma > 0)
//        return exp(-pow((*x-mu)/sigma,2)/2);
//    else
//        return (*x == mu) ? 1. : 0.;
//}
double resolution (double *x, double *p)
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
    if (sigma <= 0) return N;
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

void make_RM (TH2 * h_RM,
              TH2 * h_resolution,
              double (* xsec)(double),
              double minpt,
              double maxpt,
              TF1 * f_resolution,
              unsigned long nevents,
              TString sampling)
{
    bool kDiagonalOnly = (f_resolution->GetParameter(1) == 0),
         kPerfectSampling = (sampling == "perfect"),
         kUniformSampling = (sampling == "uniform"),
         kCoreSampling = (sampling == "core");

    double minresolution = f_resolution->GetXmin(),
           maxresolution = f_resolution->GetXmax();

    TRandom3 r;
    cout << ">>> RM and resolution" << endl;
    if (kPerfectSampling)
    {
        for (int ybin = 1 ; ybin <= h_RM->GetNbinsY() ; ybin++)
        {
            double ybinwidth = h_RM->GetYaxis()->GetBinWidth(ybin);
            for (double pt_gen = h_RM->GetYaxis()->GetBinLowEdge(ybin)+1 ; pt_gen <= h_RM->GetYaxis()->GetBinUpEdge(ybin) ; pt_gen++) // TODO: change step
            {
                double current_xsec = xsec(pt_gen);
                for (int xbin = 1 ; xbin <= h_RM->GetNbinsX() ; xbin++)
                {
                    double xbinwidth = h_RM->GetXaxis()->GetBinWidth(xbin);
                    for (double pt_rec = h_RM->GetXaxis()->GetBinLowEdge(xbin)+1 ; pt_rec <= h_RM->GetXaxis()->GetBinUpEdge(xbin) ; pt_rec++) // TODO: change step
                    {
                        double x = (pt_gen-pt_rec)/pt_gen,
                               resolution = f_resolution->Eval(x);
                        h_RM->SetBinContent(xbin, ybin, h_RM->GetBinContent(xbin, ybin) + current_xsec*resolution/(xbinwidth*ybinwidth));
                        h_resolution->SetBinContent(h_resolution->GetXaxis()->FindBin(x), xbin, h_resolution->GetBinContent(h_resolution->GetXaxis()->FindBin(x), xbin) + current_xsec*resolution/(xbinwidth*ybinwidth));
                    }
                }
            }
        }
        // attempt to add an uncertainty TODO
        for (int ybin = 1 ; ybin <= h_RM->GetNbinsY() ; ybin++)
            for (int xbin = 1 ; xbin <= h_RM->GetNbinsX() ; xbin++)
                h_RM->SetBinError(xbin, ybin, 1./sqrt(h_RM->GetBinContent(xbin, ybin)));
    }
    else
    {
        for (unsigned long i = 0 ; i < nevents ; i++)
        {
            // pick up a generated value and fill it with the weight equal to the cross section
            double pt_gen = minpt + r.Rndm()*(maxpt-minpt),
                   weight = xsec(pt_gen);
            // determine pt_rec and x according to sampling method
            double pt_rec, x;
            if (kDiagonalOnly)
            {
                pt_rec = pt_gen*(1-x);
                x = f_resolution->GetParameter(0);
            }
            else if (kUniformSampling)
            {
                pt_rec = minpt + r.Rndm()*(maxpt-minpt);
                x = (pt_gen-pt_rec)/pt_gen;
                if (x < minresolution)
                    continue; // to mimick core sampling
                else if (x > maxresolution) // this should normally never happen if maxresolution is correctly set to 1
                    throw TString::Format("pt_gen=%f && pt_rec=%f => x=%f", pt_gen, pt_rec, x);
                weight *= f_resolution->Eval(x);
            }
            else  if (kCoreSampling)
            {
                x = f_resolution->GetRandom(minresolution, maxresolution);
                pt_rec = pt_gen*(1-x);
            }
            else
                throw TString("Undefined sampling");

            h_RM        ->Fill(pt_rec, pt_gen, weight);
            h_resolution->Fill(x     , pt_gen, weight);
        }
    }
    h_RM->Scale(1./h_RM->Integral());

    // get nonnegative minimum
    double minimum = h_RM->GetMaximum(), maximum = h_RM->GetMinimum();
    for (unsigned short ibin = 1 ; ibin <= h_RM->GetNbinsX() ; ibin++) for (unsigned short jbin = 1 ; jbin <= h_RM->GetNbinsY() ; jbin++)
    {
        const double current_content = h_RM->GetBinContent(ibin,jbin);
        if (current_content > 0 && current_content < minimum) minimum = current_content;
        if (current_content > 0 && current_content > maximum) maximum = current_content;
    }
    h_RM->SetMinimum(minimum);
    h_RM->SetMaximum(maximum);
    h_RM->GetXaxis()->SetNoExponent();
    h_RM->GetXaxis()->SetMoreLogLabels();
    h_RM->GetYaxis()->SetNoExponent();
    h_RM->GetYaxis()->SetMoreLogLabels();
    h_RM->SetStats(0);

    // get nonnegative minimum
    minimum = h_resolution->GetMaximum(); maximum = h_resolution->GetMaximum(); // initial value
    for (unsigned short ibin = 0 ; ibin <= h_resolution->GetNbinsX() ; ibin++) for (unsigned short jbin = 1 ; jbin <= h_resolution->GetNbinsY() ; jbin++)
    {
        const double current_content = h_resolution->GetBinContent(ibin,jbin);
        if (current_content > 0 && current_content < minimum) minimum = current_content;
        if (current_content > 0 && current_content > maximum) maximum = current_content;
        if (current_content < 0) throw TString("Negative bin content!!");
    }
    h_resolution->SetMinimum(minimum);
    h_resolution->SetMaximum(maximum);
    h_resolution->GetYaxis()->SetNoExponent();
    h_resolution->GetYaxis()->SetMoreLogLabels();
    h_resolution->SetStats(0);
}

void make_measurement (TH1 * h_gen,
                       TH1 * h_rec,
                       double (* xsec)(double),
                       double minpt,
                       double maxpt,
                       TF1 * f_resolution,
                       unsigned long nevents,
                       double trigger)
{
    cout << ">>> truth and measurement" << endl;
    TRandom3 r; // TODO: seed
    double minresolution = f_resolution->GetXmin(),
           maxresolution = f_resolution->GetXmax();
    for (unsigned long i = 0 ; i < nevents ; i++)
    {
        // truth
        double pt_gen = minpt + r.Rndm()*(maxpt-minpt),
               weight = xsec(pt_gen);
        h_gen->Fill(pt_gen, weight);
        // measurement
        double x, pt_rec;
        if (f_resolution->GetParameter(1) == 0)
            x = f_resolution->GetParameter(0);
        else
            x = f_resolution->GetRandom(minresolution, maxresolution);
        pt_rec = pt_gen*(1-x);
        h_rec->Fill(pt_rec, weight);
    }
    h_gen->Scale(1./h_gen->Integral());
    h_rec->Scale(1./h_rec->Integral());

    short trigger_bin = h_rec->FindBin(trigger);
    double trigger_min_pt_content = h_rec->GetBinContent(trigger_bin),
           trigger_min_pt_error   = h_rec->GetBinError  (trigger_bin);
    for (short ibin = 0 ; ibin < trigger_bin ; ibin++)
    {
        h_rec->SetBinContent(ibin, trigger_min_pt_content);
        h_rec->SetBinError  (ibin, trigger_min_pt_error  );
    }

    // style
    h_rec->SetLineColor(kBlack);
    h_rec->SetMarkerColor(kBlack);
    h_rec->SetStats(0);
    h_rec->SetTitle("measurement");
    h_rec->SetMarkerStyle(20);
    h_rec->SetMarkerSize(0.5);
    h_rec->GetXaxis()->SetNoExponent();
    h_rec->GetXaxis()->SetMoreLogLabels();
    h_rec->GetYaxis()->SetRangeUser(1e-15, 1);
    h_gen->SetLineColor(kRed);
    h_gen->SetMarkerColor(kRed);
    h_gen->SetStats(0);
}

vector<double> find_binning (TH2 * h_RM, float minimal_stability, float minimal_purity)
{
    TH1D * px = h_RM->ProjectionX("px", 1, -1),
         * py = h_RM->ProjectionY("py", 1, -1);

    if (px->GetNbinsX() != py->GetNbinsX())
        throw TString("RM is not squared.");

    vector<double> new_edges;
    // determine the first bin (and if necessary, merge all non empty bins)
    int firstbin = 1;
    new_edges.push_back(h_RM->GetXaxis()->GetXmin());
    while (h_RM->GetBinContent(firstbin, firstbin) == 0 && firstbin <= h_RM->GetNbinsX()) firstbin++;
    // loop over bins
    // --> if purity and stability are good enough, then keep the up edge
    for (int ibin = firstbin ; ibin <= h_RM->GetNbinsX() ; ibin++)
    {
        cout << "ibin=" << ibin << "\tfirstbin=" << firstbin;
        double rec_content = h_RM->Integral(firstbin, ibin, 1, -1),
               gen_content = h_RM->Integral(1, -1, firstbin, ibin),
               matched_content = h_RM->Integral(firstbin, ibin, firstbin, ibin);
        cout << "\tmatched_content=" << matched_content << "\tgen_content=" << gen_content << "\trec_content=" << rec_content;
        const double stability = matched_content/gen_content,
                     purity    = matched_content/rec_content;
        cout << "\tstability=" << stability << "\tpurity=" << purity << endl;
        if (stability > minimal_stability && purity > minimal_purity)
        {
            cout << "------------------------" << endl;
            new_edges.push_back(px->GetXaxis()->GetBinUpEdge(ibin));
            firstbin = ibin+1;
        }
    }
    const int nbins = new_edges.size()-1;
    if (nbins < 2) throw TString("Invalid bin width");
    //if (nbins%2) new_edges.pop_back(); // TODO? see remark about rebinning in make_unfolding()
    new_edges.back() = px->GetXaxis()->GetXmax(); // last edge has to match /b/ old and new binnings
    // TODO: what about empty bins at high values?
    cout << "New binning (nbins=" << new_edges.size()-1 << "):";
    for (auto& new_edge: new_edges)
        cout << '\t' << new_edge;
    cout << endl;

    delete px;
    delete py;
    return new_edges;
}

vector<TH1 *> make_unfolding (TH2 * h_RM,
                              TH1 * h_rec)
{ // TODO: see whether the binning could not be different for the different cases
    vector<TH1 *> unfolded_spectra;

    TH1D * h_RM_projx = h_RM->ProjectionX("measurement", 0, -1, "oe"),
         * h_RM_projy = h_RM->ProjectionY("truth"      , 0, -1, "oe");
    TH2D * h_RM_noUnderflow = static_cast<TH2D *>(h_RM->Clone("RM_noUnderflow"));
    for (int xbin = 0 ; xbin < h_RM_noUnderflow->GetNbinsX() ; xbin++) h_RM_noUnderflow->SetBinContent(xbin, 0, 0);
    for (int ybin = 0 ; ybin < h_RM_noUnderflow->GetNbinsY() ; ybin++) h_RM_noUnderflow->SetBinContent(0, ybin, 0);
    //RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms
    // Set up from already-filled histograms.
    // "response" gives the response matrix, measured X truth.
    // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
    // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
    // in "truth" for unmeasured events (inefficiency).
    // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
    // to indicate, respectively, no fakes and/or no inefficiency.
    cout << "=== Unfolding:\n>>>Preparing RM" << endl;
    RooUnfoldResponse * RU_response = new RooUnfoldResponse(h_RM_projx, h_RM_projy, h_RM_noUnderflow);
    RooUnfold * RU_unfolding;
    // bin by bin
    cout << ">>> Bin-by-bin unfolding" << endl;
    RU_unfolding = new RooUnfoldBinByBin (RU_response, h_rec, "binbybin");
    unfolded_spectra.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("BinByBin")));
    unfolded_spectra.back()->SetTitle("BinByBin");
    // inversion
    cout << ">>> matrix-inversion unfolding" << endl;
    RU_unfolding = new RooUnfoldInvert (RU_response, h_rec, "inv");
    unfolded_spectra.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("Inverse")));
    unfolded_spectra.back()->SetTitle("Inverse");
    // Bayes
    cout << ">>> Bayes unfolding" << endl;
    RU_unfolding = new RooUnfoldBayes (RU_response, h_rec, 4, "bayes");
    unfolded_spectra.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("Bayes")));
    unfolded_spectra.back()->SetTitle("Bayes");
    // SVD
    cout << ">>> SVD unfolding" << endl;
    RU_unfolding = new RooUnfoldSvd (RU_response, h_rec, 20, "svd");
    unfolded_spectra.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("SVD")));
    unfolded_spectra.back()->SetTitle("SVD");
    //// TUnfold  TODO
    //cout << ">>> TUnfold unfolding" << endl;
    //RU_unfolding = new RooUnfoldTUnfold (RU_response, h_rec, TUnfold::kRegModeDerivative, "tu");
    //unfolded_spectra.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("TUnfold")));
    //unfolded_spectra.back()->SetTitle("TUnfold");

    vector<Color_t> colours = {kBlue, kMagenta+2, kOrange+7, kGreen+3, kCyan+2, kAzure+3};
    for (unsigned short i = 0 ; i < unfolded_spectra.size() ; i++)
    {
        unfolded_spectra[i]->SetLineColor(colours[i]);
        unfolded_spectra[i]->SetMarkerColor(colours[i]);
    }

    return unfolded_spectra;
}

pair<TH1 *, TH1 *> make_miss_fake (TH2 * h_RM, TH1 * h_gen, TH1 * h_rec)
{
    cout << ">>> miss and fake" << endl;
    TH1D * h_fake = h_RM->ProjectionX("fake", 0, 0, "oe"),
         * h_miss = h_RM->ProjectionY("miss", 0, 0, "oe");
    h_fake->Divide(h_RM->ProjectionX("fake_denominator", 0, -1, "oe"));
    h_miss->Divide(h_RM->ProjectionY("miss_denominator", 0, -1, "oe"));
    h_fake->Multiply(h_rec);
    h_miss->Multiply(h_gen);
    h_fake->SetTitle("fake");
    h_miss->SetTitle("miss");
    h_fake->SetMarkerColor(h_rec->GetMarkerColor());
    h_miss->SetMarkerColor(h_gen->GetMarkerColor());
    h_fake->SetLineColor(h_rec->GetLineColor());
    h_miss->SetLineColor(h_gen->GetLineColor());
    h_fake->SetLineStyle(2);
    h_miss->SetLineStyle(2);

    return {h_miss, h_fake};
}

vector<TH1 *> make_ABPS (TH2 * h_RM)
{
    cout << ">>> ABPS" << endl;
    // acceptance
    TH1 * acceptance = h_RM->ProjectionX("acceptance", 1, -1, "eo");
    acceptance->Divide(h_RM->ProjectionX("full_RM_px", 0, -1, "eo"));
    // background
    TH1 * background = h_RM->ProjectionY("background", 0,  0, "eo");
    background->Divide(h_RM->ProjectionY("full_RM_py", 0, -1, "eo"));
    // stability
    // purity
    TH1 * stability = static_cast<TH1 *>(h_RM->ProjectionX("RM_px", 1, -1, "oe")->Clone("stability"));
    TH1 * purity    = static_cast<TH1 *>(h_RM->ProjectionY("RM_py", 1, -1, "oe")->Clone("purity"   ));
    for (int xbin = 1 ; xbin <= stability->GetNbinsX() ; xbin++)
    {
        double diagonal_element = h_RM->GetBinContent(xbin,xbin);
        stability->SetBinContent(xbin, diagonal_element/stability->GetBinContent(xbin));
        stability->SetBinError  (xbin, 0);
        purity->SetBinContent(xbin, diagonal_element/purity->GetBinContent(xbin));
        purity->SetBinError  (xbin, 0);
    }

    acceptance->SetTitle("acceptance");  acceptance->SetLineColor(kRed );      acceptance->SetMarkerColor(kRed );
    background->SetTitle("background");  background->SetLineColor(kBlue);      background->SetMarkerColor(kBlue);
    stability ->SetTitle("stability" );  stability ->SetLineColor(kCyan  +2);  stability ->SetMarkerColor(kCyan  +2);
    purity    ->SetTitle("purity"    );  purity    ->SetLineColor(kOrange+2);  purity    ->SetMarkerColor(kOrange+2);

    acceptance->GetYaxis()->SetRangeUser(-0.05, 1.05);
    acceptance->GetXaxis()->SetMoreLogLabels();
    acceptance->GetXaxis()->SetNoExponent();
    acceptance->SetStats(0);

    return {acceptance, background, stability, purity};
}

TCanvas * make_canvas (TH1 * h_gen, 
                       TH1 * h_rec,
                       TH2 * h_RM,
                       TH2 * h_resolution,
                       TString true_xsec,
                       TString MC_xsec,
                       vector<TString> parameters)
{
    // defining canvas
    TCanvas * c = new TCanvas ("canvas");
    c->Divide(2,2); // RM, resolution, ABSP/pt spectra, ratios

    cout << "=== Plotting RM" << endl;
    c->cd(1);
    c->GetPad(1)->SetLogx();
    c->GetPad(1)->SetLogy();
    c->GetPad(1)->SetLogz();
    h_RM->GetXaxis()->SetNoExponent();
    h_RM->GetXaxis()->SetMoreLogLabels();
    h_RM->GetYaxis()->SetNoExponent();
    h_RM->GetYaxis()->SetMoreLogLabels();
    h_RM->DrawCopy("colz");
    h_RM->Write();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * RM_text = new TPaveText(0.7,0.2,0.89,0.3, "NBNDC");
    RM_text->SetTextFont(42);
    RM_text->SetFillColorAlpha(0,0);
    RM_text->AddText(MC_xsec);
    RM_text->Draw();

    cout << "=== Plotting resolution" << endl;
    c->cd(3);
    c->GetPad(3)->SetLogy();
    c->GetPad(3)->SetLogz();
    h_resolution->DrawCopy("colz");
    h_resolution->Write();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * resolution_text  = new TPaveText(0.7, 0.7 , 0.89, 0.89, "NBNDC");
    resolution_text->SetFillColorAlpha(0,0);
    resolution_text->SetTextFont(42);
    for (const TString& parameter: parameters)
        resolution_text->AddText(parameter);
    resolution_text->Draw();

    c->cd(2);
    c->GetPad(2)->Divide(1,2, 0);
    cout << "=== Plotting ABPS" << endl;
    c->GetPad(2)->cd(1);
    c->GetPad(2)->GetPad(1)->SetLogx();
    vector<TH1 *> ABPS = make_ABPS(h_RM);
    ABPS.front()->SetStats(0);
    ABPS.front()->Draw("hist");
    ABPS.front()->Write();
    for (unsigned short i = 1 ; i < ABPS.size() ; i++)
    {
        ABPS[i]->Draw("same hist");
        ABPS[i]->Write();
    }
    c->GetPad(2)->GetPad(1)->BuildLegend(0.4,0.35,0.6,0.6);
    c->GetPad(2)->GetPad(1)->RedrawAxis();

    pair<TH1 *, TH1 *> miss_and_fake = make_miss_fake(h_RM, h_gen, h_rec);
    //miss_and_fake.first ->Rebin(2);
    //miss_and_fake.second->Rebin(2);

    cout << "=== Plotting spectra" << endl;
    c->GetPad(2)->cd(2);
    c->GetPad(2)->GetPad(2)->SetLogx();
    c->GetPad(2)->GetPad(2)->SetLogy();
    //h_RM->Rebin2D(1,2);
    //h_gen->Rebin(2);

    // the following vector will contain all the spectra at hadron level: the generated spectrum and the unfolded spectra
    vector<TH1 *> unfolded_spectra = make_unfolding(h_RM, h_rec);

    //h_rec->Rebin(2); // TODO: keep fine binning at rec level for the unfolding??
    h_rec->Scale(1,"width");
    h_gen->Scale(1,"width");
    h_rec->DrawCopy();
    h_gen->DrawCopy("same");
    h_gen->Write();
    h_rec->Write();

    miss_and_fake.first ->Scale(1,"width");
    miss_and_fake.second->Scale(1,"width");
    miss_and_fake.first ->DrawCopy("same hist");
    miss_and_fake.second->DrawCopy("same hist");
    miss_and_fake.first ->Write();
    miss_and_fake.second->Write();

    for (auto& unfolded_spectrum: unfolded_spectra)
    {
        unfolded_spectrum->Scale(1,"width");
        unfolded_spectrum->DrawCopy("same");
        unfolded_spectrum->Write();
    }
    c->GetPad(2)->GetPad(2)->BuildLegend(0.3, 0.3, 0.5, 0.6);
    c->GetPad(2)->GetPad(2)->RedrawAxis();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * spectrum_text = new TPaveText(0.7,0.2,0.89,0.3, "NBNDC");
    spectrum_text->SetFillColorAlpha(0,0);
    spectrum_text->SetTextFont(42);
    spectrum_text->AddText(true_xsec);
    spectrum_text->Draw();

    cout << "=== Plotting ratios of unfolded spectra to truth" << endl;
    c->cd(4);
    vector<TH1 *> numerators = {h_gen};
    numerators.insert(end(numerators), begin(unfolded_spectra), end(unfolded_spectra));
    c->GetPad(4)->Divide(1,numerators.size(), 0,0);
    for (unsigned short i = 0 ; i < numerators.size() ; i++)
    {
        c->GetPad(4)->cd(i+1);
        c->GetPad(4)->GetPad(i+1)->SetLogx();
        TH1D * ratio = static_cast<TH1D*>(numerators[i]->Clone(TString::Format("%s/rec", numerators[i]->GetName())));
        ratio->Divide(h_rec);
        ratio->SetTitle(TString::Format(";;%s", numerators[i]->GetTitle()));
        ratio->GetYaxis()->SetTitleSize(0.2);
        ratio->GetYaxis()->SetTitleOffset(0.2);
        ratio->GetYaxis()->SetLabelSize(0.1);
        ratio->SetNdivisions(304, "Y");
        ratio->SetStats(0);
        ratio->GetYaxis()->SetRangeUser(0,2);
        ratio->DrawCopy();
        ratio->Write();
        delete ratio;
    }

    delete h_gen;
    delete h_rec;
    delete h_RM;
    delete h_resolution;
    delete miss_and_fake.first ;
    delete miss_and_fake.second;
    for (auto& unfolded_spectrum: unfolded_spectra)
        delete unfolded_spectrum;

    // esthetics
    cout << "=== Redrawing axes" << endl;
    for (unsigned short i = 1 ; i <= 3 ; i++)
        c->GetPad(i)->RedrawAxis();
    for (unsigned short i = 1 ; i <= numerators.size() ; i++)
        c->GetPad(4)->GetPad(i)->RedrawAxis();

    c->Write();

    cout << "=== The end" << endl;
    return c;
}

int main (int argc, char* argv[])
{
    gStyle->SetOptTitle(0);
    TH1D::SetDefaultSumw2(true); // correct error computation
    TColor::CreateColorWheel();

    // declaring TApplication (needed to use the full power of Root as a library)
    TApplication * rootapp = new TApplication ("unfolding", &argc, argv);

    // batch mode? (useful if bad connection or if submitted as a job)
    /*bool batch = (rootapp->Argc() > 1 && (TString(rootapp->Argv(1)) == "--batch" ||
                                          TString(rootapp->Argv(1)) ==   "batch" ||
                                          TString(rootapp->Argv(1)) ==   "b"     ||
                                          TString(rootapp->Argv(1)) ==  "-b"       ));
    if (batch)*/ gROOT->SetBatch();
    TFile * f = new TFile ("unfolding.root", "RECREATE");

    // parameters TODO
    const unsigned long nevents = 5e6;
    vector<TString> vsampling = {"perfect", "uniform", "core"};
    double N = 1, tau = 0, kL = -1, kR = 1, aL = -1, nL = 0, aR = 1, nR = 0; // TODO: play with this
    vector<double> vmu     ; for (unsigned short i = 0 ; i <= 0 ; i++) vmu     .push_back(i*0.01);
    vector<double> vsigma  ; for (unsigned short i = 2 ; i <= 2 ; i++) vsigma  .push_back(i*0.02);
    vector<double> vminSP  ; for (unsigned short i = 6 ; i <= 9 ; i++) vminSP  .push_back(i*0.10);
    vector<double> triggers; for (unsigned short i = 0 ; i <= 0 ; i++) triggers.push_back(i*  10);
    
    map<TString, double (*)(double)> MC_spectra = {{"flat spectrum", flat_spectrum},
        {"#frac{1}{p_{T}^{4}}", falling_spectrum},
        {"#frac{1}{p_{T}^{6}}", even_more_falling_spectrum}};

    // open PDF files for the different samplings
    for (const TString& sampling: vsampling)
    {
        (new TCanvas (sampling))->Print(sampling + ".pdf[");
        divider(sampling + " sampling")->Print(sampling + ".pdf"); // frontpage
    }
    vector<double> binning;
    for (int i = 18 ; i <= 153 ; i++)
        binning.push_back(i);
    //// Mikko's binning
    //vector<double> std_binning = {/*0,*/18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330/*,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1684,1784,1890,2000,2116,2238,2366,2500,2640,2787,2941,3103,3273,3450,3637,3832,4037,4252,4477,4713,4961,5220,5492,5777,6076,6389,6717,7000*/}; // std in SMP-j
    //// regular binning
    //vector<double> constant_binwidth_binning(binning.size()); // making binning in similar range but with constant binwidth
    //for (size_t i = 0 ; i < binning.size() ; i++) constant_binwidth_binning[i] = i*binning.back()/binning.size();

    for (auto& MC_spectrum: MC_spectra)
    {
        TCanvas * section = divider("MC spectrum ~ " + MC_spectrum.first);
        for (const TString& sampling: vsampling)
            section->Print(sampling + ".pdf");
        for (const double& trigger: triggers) for (const double& minSP: vminSP) for (const double& mu: vmu) for (const double& sigma: vsigma)
        {
            cout << "================================================================================"
                 << "\nParameters:"
                 << "\n\tminSP = " << minSP 
                 << "\n\tmu = " << mu
                 << "\n\tsigma = " << sigma 
                 << "\n\tnevents = " << nevents << endl;
            vector<TString> parameterisation = {TString::Format("#minSP=%f",minSP), TString::Format("#mu=%f",mu), TString::Format("#sigma=%f", sigma)};

            TString dirname = MC_spectrum.first + "_" + TString::Format("_minSP%f_mu%f_sigma%f", minSP, mu, sigma);
            dirname.ReplaceAll(".", "p");
            dirname.ReplaceAll(" ", "_");
            dirname.ReplaceAll("-", "_");
            dirname.ReplaceAll("#", "");
            dirname.ReplaceAll("{", "");
            dirname.ReplaceAll("}", "");
            dirname.ReplaceAll("^", "");
            f->mkdir(dirname)->cd();

            cout << "=== Defining resolution function (used only if sigma different from 0)" << endl;
            TF1 * f_resolution = new TF1 ("resolution", resolution, -1, 1, 2);
            f_resolution->SetParameters(N, mu, sigma, tau, kL, kR, aL, nL, aR, nR);

            cout << "=== Defining gen and rec histograms" << endl;
            TH1D * h_gen = new TH1D ("gen", "Truth"      , binning.size()-1, &binning[0]),
                 * h_rec = new TH1D ("rec", "Measurement", binning.size()-1, &binning[0]);
            const double minpt = 0, maxpt = h_gen->GetXaxis()->GetXmax()*1.1;
            make_measurement(h_gen, h_rec, falling_spectrum, minpt, maxpt, f_resolution, nevents, trigger);
            for (TString sampling: vsampling)
            {
                TCanvas * c;
                try
                {
                    cout << "=== Defining RM and resolution histograms" << endl;
                    TH2 * h_RM = new TH2D ("RM", "RM", binning.size()-1, &binning[0], binning.size()-1, &binning[0]),
                        * h_resolution = new TH2D("resolution", "resolution", 41, -1, 1, binning.size()-1, &binning[0]);

                    cout << "=== Filling histograms" << endl;
                    make_RM(h_RM, h_resolution, MC_spectrum.second, minpt, maxpt, f_resolution, nevents, sampling);

                    cout << "=== Rebinning" << endl;
                    vector<double> new_edges = find_binning(h_RM, minSP, minSP);

                    int nbins = new_edges.size()-1;
                    TH1 * rebinned_gen = h_gen->Rebin(nbins, sampling + TString("_rebinned_") + h_gen->GetName(), &new_edges[0]),
                        * rebinned_rec = h_rec->Rebin(nbins, sampling + TString("_rebinned_") + h_rec->GetName(), &new_edges[0]);
                    TH2 * rebinned_RM = new TH2D ("rebinned_RM", h_RM->GetTitle(), nbins, &new_edges[0], nbins, &new_edges[0]);
                    rebinned_RM->SetStats(0);
                    for (int ybin = 0 ; ybin <= h_RM->GetNbinsY()+1 ; ++ybin)
                        for (int xbin = 0 ; xbin <= h_RM->GetNbinsX()+1 ; ++xbin)
                            rebinned_RM->Fill(h_RM->GetXaxis()->GetBinCenter(xbin),
                                              h_RM->GetYaxis()->GetBinCenter(ybin),
                                              h_RM->GetBinContent(xbin, ybin));

                    cout << "=== Making canvas" << endl;
                    //make_canvas(static_cast<TH1 *>(h_gen->Clone("gen_" + sampling)), static_cast<TH1 *>(h_rec->Clone("rec_" + sampling)), h_RM, h_resolution, /*true*/ "#frac{1}{p_{T}^{4}}", /*MC*/ MC_spectrum.first, {TString::Format("#mu=%f",mu), TString::Format("#sigma=%f", sigma)} )->Print(sampling + ".pdf");

                    c = make_canvas(rebinned_gen, rebinned_rec, rebinned_RM, h_resolution, /*true*/ "#frac{1}{p_{T}^{4}}", /*MC*/ MC_spectrum.first, parameterisation);
                }
                catch (TString msg)
                {
                    vector<TString> pave = {msg, MC_spectrum.first};
                    pave.insert(end(pave), begin(parameterisation), end(parameterisation));
                    c = divider(pave, kWhite);
                    cout << "An error was thrown! Continuing." << endl;
                }
                if (c != nullptr) c->Print(sampling + ".pdf");
            }
        }
    }
    for (const TString& sampling: vsampling)
        (new TCanvas (sampling))->Print(sampling + ".pdf]");

    // ending
    f->Close();
    //if (batch)
        rootapp->Terminate();
    //else
    //    rootapp->Run();
    return EXIT_SUCCESS;
}

