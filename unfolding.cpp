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
// proper header
#include "unfolding.h"

using namespace std;

double get_nonnegative_minimum (const TH1 * h, const bool include_underflow = false)
{
    double minimum = h->GetMaximum();
    const unsigned short nbinsx = h->GetNbinsX(),
                         nbinsy = h->GetNbinsY(),
                         nbinsz = h->GetNbinsZ();
    switch (h->GetDimension())
    {
    case 1:
        for (unsigned short ibin = include_underflow ? 0 : 1 ; ibin <= nbinsx ; ibin++)
        {
            const double current_content = h->GetBinContent(ibin);
            if (current_content > 0 && current_content < minimum)
                minimum = current_content;
        }
    case 2:
        for (unsigned short ibin = include_underflow ? 0 : 1 ; ibin <= nbinsx ; ibin++)
            for (unsigned short jbin = include_underflow ? 0 : 1 ; jbin <= nbinsy ; jbin++)
            {
                const double current_content = h->GetBinContent(ibin,jbin);
                if (current_content > 0 && current_content < minimum)
                    minimum = current_content;
            }
    case 3:
        for (unsigned short ibin = include_underflow ? 0 : 1 ; ibin <= nbinsx ; ibin++)
            for (unsigned short jbin = include_underflow ? 0 : 1 ; jbin <= nbinsy ; jbin++)
                for (unsigned short kbin = include_underflow ? 0 : 1 ; kbin <= nbinsz ; kbin++)
                {
                    const double current_content = h->GetBinContent(ibin,jbin,kbin);
                    if (current_content > 0 && current_content < minimum)
                        minimum = current_content;
                }
    }
    return minimum;
}

//double flat_spectrum(double pt)
//{
//    return 1;
//}
//double falling_spectrum(double pt)
//{
//    if (pt == 0) return 0;
//    return pow(pt,-4);
//}
//double even_more_falling_spectrum(double pt)
//{
//    if (pt == 0) return 0;
//    return pow(pt,-6);
//}

void make_RM (TH2 * h_RM,
              TH2 * h_resolution,
              TF1 * xsec,
              //double (* xsec)(double),
              double minpt,
              double maxpt,
              TF1 * f_resolution,
              TF1 * f_sigma,
              unsigned long nevents,
              TString sampling)
{
    bool kDiagonalOnly = (f_resolution->GetParameter(2) == 0),
         kPerfectSampling = (sampling == "perfect"),
         kUniformSampling = (sampling == "uniform"),
         kCoreSampling = (sampling == "core");

    double minresolution = f_resolution->GetXmin(),
           maxresolution = f_resolution->GetXmax();

    TRandom3 r;
    cout << ">>> RM and resolution" << endl;
    const double step = 1;
    if (kPerfectSampling)
    {
        for (int ybin = 1 ; ybin <= h_RM->GetNbinsY() ; ybin++)
        {
            double ybinwidth = h_RM->GetYaxis()->GetBinWidth(ybin);
            for (double pt_gen = h_RM->GetYaxis()->GetBinLowEdge(ybin)+1 ; pt_gen <= h_RM->GetYaxis()->GetBinUpEdge(ybin) ; pt_gen += step)
            {
                double current_xsec = xsec->Eval(pt_gen);
                for (int xbin = 1 ; xbin <= h_RM->GetNbinsX() ; xbin++)
                {
                    double xbinwidth = h_RM->GetXaxis()->GetBinWidth(xbin);
                    for (double pt_rec = h_RM->GetXaxis()->GetBinLowEdge(xbin)+1 ; pt_rec <= h_RM->GetXaxis()->GetBinUpEdge(xbin) ; pt_rec += step)
                    {
                        f_resolution->SetParameter(2, f_sigma->Eval(pt_gen));
                        double x = (pt_gen-pt_rec)/pt_gen,
                               resolution = f_resolution->Eval(x);
                        h_RM->SetBinContent(xbin, ybin, h_RM->GetBinContent(xbin, ybin) + current_xsec*resolution/(xbinwidth*ybinwidth));
                        h_resolution->SetBinContent(h_resolution->GetXaxis()->FindBin(x), xbin, h_resolution->GetBinContent(h_resolution->GetXaxis()->FindBin(x), xbin) + current_xsec*resolution/(xbinwidth*ybinwidth));
                    }
                }
            }
        }
        // add an uncertainty TODO?
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
                   weight = xsec->Eval(pt_gen);
            // determine pt_rec and x according to sampling method
            double pt_rec, x;
            if (kDiagonalOnly)
            {
                pt_rec = pt_gen*(1-x);
                x = f_resolution->GetParameter(1);
            }
            else
            {
                f_resolution->SetParameter(2, f_sigma->Eval(pt_gen));
                if (kUniformSampling)
                {
                    pt_rec = minpt + r.Rndm()*(maxpt-minpt);
                    x = (pt_gen-pt_rec)/pt_gen;
                    if (x < minresolution)
                        continue; // to mimick core sampling
                    else if (x > maxresolution) // this should normally never happen
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
            }

            h_RM        ->Fill(pt_rec, pt_gen, weight);
            h_resolution->Fill(x     , pt_gen, weight);
        }
    }
    h_RM->Scale(1./h_RM->Integral());

    h_RM->SetMinimum(get_nonnegative_minimum(h_RM, true));
    h_RM->SetNdivisions(304, "Z");
    h_RM->GetXaxis()->SetNoExponent();
    h_RM->GetXaxis()->SetMoreLogLabels();
    h_RM->GetYaxis()->SetNoExponent();
    h_RM->GetYaxis()->SetMoreLogLabels();
    h_RM->SetStats(0);

    h_resolution->SetMinimum(get_nonnegative_minimum(h_resolution, true));
    h_resolution->SetNdivisions(304, "Z");
    h_resolution->GetYaxis()->SetNoExponent();
    h_resolution->GetYaxis()->SetMoreLogLabels();
    h_resolution->SetStats(0);
}

void make_measurement (TH1 * h_gen,
                       TH1 * h_rec,
                       TF1 * xsec,
                       TF1 * f_efficiency,
                       double minpt,
                       double maxpt,
                       TF1 * f_resolution,
                       TF1 * f_sigma,
                       unsigned long nevents)
{
    cout << ">>> truth and measurement" << endl;
    TRandom3 r;
    double minresolution = f_resolution->GetXmin(),
           maxresolution = f_resolution->GetXmax();
    for (unsigned long i = 0 ; i < nevents ; i++)
    {
        // truth
        double pt_gen = minpt + r.Rndm()*(maxpt-minpt),
               weight = xsec->Eval(pt_gen) * f_efficiency->Eval(pt_gen);
        h_gen->Fill(pt_gen, weight);
        // measurement
        double x, pt_rec;
        if (f_resolution->GetParameter(1) == 0)
            x = f_resolution->GetParameter(1);
        else
        {
            f_resolution->SetParameter(2, f_sigma->Eval(pt_gen));
            x = f_resolution->GetRandom(minresolution, maxresolution);
        }
        pt_rec = pt_gen*(1-x);
        h_rec->Fill(pt_rec, weight);
    }
    h_gen->Scale(1./h_gen->Integral());
    h_rec->Scale(1./h_rec->Integral());

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
        //cout << "ibin=" << ibin << "\tfirstbin=" << firstbin;
        double rec_content = h_RM->Integral(firstbin, ibin, 1, -1),
               gen_content = h_RM->Integral(1, -1, firstbin, ibin),
               matched_content = h_RM->Integral(firstbin, ibin, firstbin, ibin);
        //cout << "\tmatched_content=" << matched_content << "\tgen_content=" << gen_content << "\trec_content=" << rec_content;
        const double stability = matched_content/gen_content,
                     purity    = matched_content/rec_content;
        //cout << "\tstability=" << stability << "\tpurity=" << purity << endl;
        if (stability > minimal_stability && purity > minimal_purity)
        {
            //cout << "------------------------" << endl;
            new_edges.push_back(px->GetXaxis()->GetBinUpEdge(ibin));
            firstbin = ibin+1;
        }
    }
    const int nbins = new_edges.size()-1;
    if (nbins < 2) throw TString("Invalid bin width");
    //new_edges.back() = px->GetXaxis()->GetXmax(); // last edge has to match /b/ old and new binnings
    new_edges.push_back(px->GetXaxis()->GetXmax()); // in case of additional rebinning for TUnfold...
    cout << "New binning (nbins=" << new_edges.size()-1 << "):";
    for (auto& new_edge: new_edges)
        cout << '\t' << new_edge;
    cout << endl;

    delete px;
    delete py;
    return new_edges;
}

vector<TH1 *> make_unfolding (TH2 * rebinned_RM,
                              TH1 * rebinned_rec,
                              TH2 * semi_rebinned_RM,
                              TH1 * semi_rebinned_rec,
                              vector<UnfoldingParameters *> v_parameters)
{
    cout << "=== Unfolding:\n>>>Preparing RM" << endl;

    cout << ">>> Semi-rebinned RM" << endl;
    TH1D * semi_rebinned_RM_projx = semi_rebinned_RM->ProjectionX("semi_rebinned_measurement", 0, -1, "oe"),
         * semi_rebinned_RM_projy = semi_rebinned_RM->ProjectionY("semi_rebinned_truth"      , 0, -1, "oe");
    TH2D * semi_rebinned_RM_noUnderflow = static_cast<TH2D *>(semi_rebinned_RM->Clone("semi_rebinned_RM_noUnderflow"));
    for (int xbin = 0 ; xbin < semi_rebinned_RM_noUnderflow->GetNbinsX() ; xbin++) semi_rebinned_RM_noUnderflow->SetBinContent(xbin, 0, 0);
    for (int ybin = 0 ; ybin < semi_rebinned_RM_noUnderflow->GetNbinsY() ; ybin++) semi_rebinned_RM_noUnderflow->SetBinContent(0, ybin, 0);
    RooUnfoldResponse * semi_rebinned_RU_response = new RooUnfoldResponse(semi_rebinned_RM_projx, semi_rebinned_RM_projy, semi_rebinned_RM_noUnderflow);

    cout << ">>> Rebinned RM" << endl;
    TH1D * rebinned_RM_projx = rebinned_RM->ProjectionX("rebinned_measurement", 0, -1, "oe"),
         * rebinned_RM_projy = rebinned_RM->ProjectionY("rebinned_truth"      , 0, -1, "oe");
    TH2D * rebinned_RM_noUnderflow = static_cast<TH2D *>(rebinned_RM->Clone("rebinned_RM_noUnderflow"));
    for (int xbin = 0 ; xbin < rebinned_RM_noUnderflow->GetNbinsX() ; xbin++) rebinned_RM_noUnderflow->SetBinContent(xbin, 0, 0);
    for (int ybin = 0 ; ybin < rebinned_RM_noUnderflow->GetNbinsY() ; ybin++) rebinned_RM_noUnderflow->SetBinContent(0, ybin, 0);
    RooUnfoldResponse * rebinned_RU_response = new RooUnfoldResponse(rebinned_RM_projx, rebinned_RM_projy, rebinned_RM_noUnderflow);

    //RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms
    // Set up from already-filled histograms.
    // "response" gives the response matrix, measured X truth.
    // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
    // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
    // in "truth" for unmeasured events (inefficiency).
    // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
    // to indicate, respectively, no fakes and/or no inefficiency.

    cout << "=== Unfolding:\n>>>Applying unfolding" << endl;
    vector<TH1 *> unfolded_spectra;
    for (UnfoldingParameters * parameters: v_parameters)
    {
        TString title;
        RooUnfold * RU_unfolding;

        switch (parameters->type)
        {
            case UnfoldingParameters::kBinByBin:
            {
                cout << ">> Bin-by-bin unfolding" << endl;
                title = "Bin by bin";
                RU_unfolding = new RooUnfoldBinByBin (rebinned_RU_response, rebinned_rec, "binbybin");
                break;
            }
            case UnfoldingParameters::kInversion:
            {
                cout << ">> matrix-inversion unfolding" << endl;
                title = "Inversion";
                RU_unfolding = new RooUnfoldInvert (rebinned_RU_response, rebinned_rec, "inv");
                break;
            }
            case UnfoldingParameters::kSVD:
            {
                cout << ">> SVD unfolding" << endl;
                UnfoldingParametersSVD * svd_parameters = static_cast<UnfoldingParametersSVD *>(parameters);
                title = TString::Format("SVD (k=%d)", svd_parameters->kreg);
                RU_unfolding = new RooUnfoldSvd (rebinned_RU_response, rebinned_rec, svd_parameters->kreg, "svd");
                break;
            }
            case UnfoldingParameters::kBayes:
            {
                UnfoldingParametersBayes * bayes_parameters = static_cast<UnfoldingParametersBayes *>(parameters);
                title = TString::Format("Bayes (n=%d)", bayes_parameters->niterations); // TODO: include binning scheme?
                switch (bayes_parameters->rm_type)
                {
                    case UnfoldingParametersBayes::kSquare:
                        cout << ">> Bayes unfolding (squared RM)" << endl;
                        RU_unfolding = new RooUnfoldBayes (rebinned_RU_response, rebinned_rec, bayes_parameters->niterations, "bayes");
                        break;
                    case UnfoldingParametersBayes::kFine:
                        cout << ">> Bayes unfolding (fine binning)" << endl;
                        RU_unfolding = new RooUnfoldBayes (semi_rebinned_RU_response, semi_rebinned_rec, bayes_parameters->niterations, "bayes");
                        break;
                    default:
                        throw string("Unknown Bayes parameters.");
                }
                break;
            }
            case UnfoldingParameters::kTUnfold:
            {
                cout << ">>> TUnfold unfolding" << endl;
                UnfoldingParametersTUnfold * tunfold_parameters = static_cast<UnfoldingParametersTUnfold *>(parameters);
                title = "TUnfold"; // TODO: include regularisation scheme?
                RU_unfolding = new RooUnfoldTUnfold (semi_rebinned_RU_response, semi_rebinned_rec, tunfold_parameters->regularisation, "tu");
                break;
            }
            default:
                throw string("Unknown unfolding type.");
        }
        unfolded_spectra.push_back(RU_unfolding->Hreco());
        unfolded_spectra.back()->SetTitle(title);
        delete RU_unfolding;
    }

    vector<Color_t> colours = {kBlue, kMagenta+2, kOrange+7, kGreen+3, kCyan+2, kAzure+3, kOrange-5, kCyan-1, kYellow+4};
    for (unsigned short i = 0 ; i < unfolded_spectra.size() ; i++)
    {
        unfolded_spectra[i]->SetLineColor  (colours[i%colours.size()]);
        unfolded_spectra[i]->SetMarkerColor(colours[i%colours.size()]);
    }

    delete semi_rebinned_RU_response;
    delete semi_rebinned_RM_projx;
    delete semi_rebinned_RM_projy;
    delete semi_rebinned_RM_noUnderflow;
    delete rebinned_RU_response;
    delete rebinned_RM_projx;
    delete rebinned_RM_projy;
    delete rebinned_RM_noUnderflow;

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
    TH1 * stability = static_cast<TH1 *>(h_RM->ProjectionX("RM_px", 1, -1, "oe")->Clone("stability"));
    // purity
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

void draw_response_matrix (TVirtualPad * p,
                           TH2 * RM)
{
    const unsigned short nbinsx = RM->GetNbinsX(),
                         nbinsy = RM->GetNbinsY();
    RM->SetStats(0);
    RM->SetTitle(";p_{T}^{rec};p_{T}^{gen};N");
    RM->GetXaxis()->SetLabelSize(0);
    RM->GetYaxis()->SetLabelSize(0);

    // fake
    const string fake_name = string("fake_") + RM->GetName();
    TH2 * fake = static_cast<TH2 *>(RM->Clone(fake_name.c_str()));
    for (unsigned xbin = 0 ; xbin <= nbinsx ; xbin++)
    {
        for (unsigned ybin = 1 ; ybin <= nbinsy ; ybin++)
        {
            fake->SetBinContent(xbin, ybin, 0);
            fake->SetBinError  (xbin, ybin, 0);
        }
        fake->SetBinContent(xbin, 1, fake->GetBinContent(xbin, 0));
        fake->SetBinError  (xbin, 1, fake->GetBinContent(xbin, 0));
    }
    fake->RebinY(nbinsy);
    fake->GetYaxis()->SetLabelSize(0);
    fake->GetZaxis()->SetLabelSize(0);
    fake->SetTitle(";rec;;");
    // miss
    const string miss_name = string("miss_") + RM->GetName();
    TH2 * miss = static_cast<TH2 *>(RM->Clone(miss_name.c_str()));
    for (unsigned ybin = 0 ; ybin <= nbinsy ; ybin++)
    {
        for (unsigned xbin = 1 ; xbin <= nbinsx ; xbin++)
        {
            miss->SetBinContent(xbin, ybin, 0);
            miss->SetBinError  (xbin, ybin, 0);
        }
        miss->SetBinContent(1, ybin, miss->GetBinContent(0, ybin));
        miss->SetBinError  (1, ybin, miss->GetBinContent(0, ybin));
    }
    miss->RebinX(nbinsx);
    miss->GetXaxis()->SetLabelSize(0);
    miss->GetZaxis()->SetLabelSize(0);
    miss->SetTitle(";;gen;");

    // *********************************************
    // *      *                                    *
    // *      *          m                         *
    // *      *           a                        *
    // *   m  *            t                       *
    // *   i  *             c                      *
    // *   s  *              h                     *
    // *   s  *               e                    *
    // *      *                d                   *
    // *      *                                    *
    // *      *                                    *
    // *********************************************
    // *      *                                    *
    // *      *        fake                        *
    // *      *                                    *
    // *********************************************
    p->cd();
    p->SetRightMargin(0);
    p->SetLeftMargin(0);
    p->SetTopMargin(0.1);
    p->SetBottomMargin(0);
    
    const float lower_margin = 0.15,
                left_margin  = 0.15;
    TPad * p_fake = new TPad ("fake", "", left_margin, 0, 1, lower_margin),
         * p_miss = new TPad ("miss", "", 0, lower_margin, left_margin, 1),
         * p_main = new TPad ("main", "", left_margin, lower_margin, 1, 1),
         * p_labels = new TPad ("labels", "", 0, 0, left_margin, lower_margin);

    const double max = TMath::Max(TMath::Max(fake->GetMaximum(), miss->GetMaximum()), RM->GetMaximum()),
                 min = TMath::Min(TMath::Min(get_nonnegative_minimum(fake), get_nonnegative_minimum(miss)), get_nonnegative_minimum(RM));
    RM  ->SetMinimum(min);    RM  ->SetMaximum(max);
    miss->SetMinimum(min);    miss->SetMaximum(max);
    fake->SetMinimum(min);    fake->SetMaximum(max);

    const double pad_right_margin = 0.1,
                 pad_bottom_margin = 0.001,
                 pad_upper_margin = 0.03;
    p_main->SetTicks(1,1);
    p_main->SetRightMargin(pad_right_margin);
    p_main->SetLeftMargin(0);
    p_main->SetTopMargin(pad_upper_margin);
    p_main->SetBottomMargin(pad_bottom_margin);
    p_main->SetLogx();
    p_main->SetLogy();
    p_main->SetLogz();
    p_main->Draw();
    p_main->cd();
    RM->GetXaxis()->SetLabelSize(0);   RM->GetXaxis()->SetTitleSize(0);
    RM->GetYaxis()->SetLabelSize(0);   RM->GetYaxis()->SetTitleSize(0);
    RM->Draw("colz");

    p->cd();
    p_miss->SetTicks(1,1);
    p_miss->SetRightMargin(0.1);
    p_miss->SetLeftMargin(0.5);
    p_miss->SetTopMargin(pad_upper_margin);
    p_miss->SetBottomMargin(pad_bottom_margin);
    p_miss->SetLogy();
    p_miss->SetLogz();
    p_miss->Draw();
    p_miss->cd();
    miss->GetYaxis()->SetMoreLogLabels();
    miss->GetYaxis()->SetNoExponent();
    miss->GetYaxis()->SetLabelSize(0.15);
    miss->GetYaxis()->SetLabelFont(42);
    miss->GetYaxis()->SetTitleSize(0.2);
    miss->GetYaxis()->SetTitleOffset(0.5);
    miss->SetNdivisions(0, "X");
    miss->Draw("col");

    p->cd();
    p_fake->SetTicks(1,1);
    p_fake->SetRightMargin(pad_right_margin);
    p_fake->SetLeftMargin(0);
    p_fake->SetTopMargin(0.1);
    p_fake->SetBottomMargin(0.5);
    p_fake->SetLogx();
    p_fake->SetLogz();
    p_fake->Draw();
    p_fake->cd();
    fake->GetXaxis()->SetMoreLogLabels();
    fake->GetXaxis()->SetNoExponent();
    fake->GetXaxis()->SetLabelSize(0.3);
    fake->GetXaxis()->SetLabelFont(42);
    fake->GetXaxis()->SetTitleSize(0.3);
    fake->GetXaxis()->SetTitleOffset(0.5);
    fake->SetNdivisions(0, "Y");
    fake->Draw("col");

    p->cd();
    p_labels->Draw();
    p_labels->cd();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * miss_label = new TPaveText (0.2, 0.7, 0.5, 0.8, "NB"),
              * fake_label = new TPaveText (0.6, 0.2, 0.9, 0.3, "NB");
    const float text_size = 0.3;
    miss_label->SetFillColorAlpha(kWhite,0);   fake_label->SetFillColorAlpha(kWhite,0);   
    miss_label->SetTextSize(text_size);        fake_label->SetTextSize(text_size);
    miss_label->SetTextFont(42);               fake_label->SetTextFont(42);
    miss_label->AddText("miss");               fake_label->AddText("fake");
    miss_label->Draw();                        fake_label->Draw();
}

void draw_turnon_line (double turnon, double ymin, double ymax)
{
    double x[2] = {turnon, turnon},
           y[2] = {ymin,ymax};
    TGraph * line = new TGraph(2,x,y);
    line->Draw("sameL");
}

TCanvas * make_canvas (TH1 * h_gen, 
                       TH1 * h_rec,
                       TH2 * h_RM,
                       TH2 * h_resolution,
                       vector<double> new_edges,
                       TString true_xsec,
                       TString MC_xsec,
                       vector<UnfoldingParameters *> v_parameters,
                       vector<TString> pave_res_parameters,
                       vector<TString> pave_ABPS_eff,
                       double turnon)
{
    int nbins = new_edges.size()-1;
    TH1 * rebinned_gen = h_gen->Rebin(nbins, TString("rebinned_") + h_gen->GetName(), &new_edges[0]),
        * rebinned_rec = h_rec->Rebin(nbins, TString("rebinned_") + h_rec->GetName(), &new_edges[0]);
    TH2 * rebinned_RM = new TH2D ("rebinned_RM", h_RM->GetTitle(), nbins, &new_edges[0], nbins, &new_edges[0]),
        * semi_rebinned_RM = new TH2D ("semi_rebinned_RM", h_RM->GetTitle(), h_RM->GetNbinsX(), h_RM->GetXaxis()->GetXmin(), h_RM->GetXaxis()->GetXmax(), nbins, &new_edges[0]);
    rebinned_RM->SetStats(0);
    for (int ybin = 0 ; ybin <= h_RM->GetNbinsY()+1 ; ++ybin) for (int xbin = 0 ; xbin <= h_RM->GetNbinsX()+1 ; ++xbin)
    {
        rebinned_RM->Fill(h_RM->GetXaxis()->GetBinCenter(xbin),
                          h_RM->GetYaxis()->GetBinCenter(ybin),
                          h_RM->GetBinContent(xbin, ybin));
        semi_rebinned_RM->Fill(h_RM->GetXaxis()->GetBinCenter(xbin),
                               h_RM->GetYaxis()->GetBinCenter(ybin),
                               h_RM->GetBinContent(xbin, ybin));
    }

    // defining canvas
    TCanvas * c = new TCanvas ("canvas");
    c->Divide(2,2); // RM, resolution, ABSP/pt spectra, ratios

    cout << "=== Plotting RM" << endl;
    c->cd(1);
    h_RM->Write();
    draw_response_matrix(c->GetPad(1),h_RM);
    c->cd(1);
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
    c->GetPad(3)->SetTopMargin(0.05);
    h_resolution->DrawCopy("colz");
    h_resolution->Write();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * resolution_text  = new TPaveText(0.7, 0.6 , 0.89, 0.89, "NBNDC");
    resolution_text->SetFillColorAlpha(0,0);
    resolution_text->SetTextFont(42);
    for (const TString& parameter: pave_res_parameters)
        resolution_text->AddText(parameter);
    resolution_text->Draw();

    c->cd(2);
    c->GetPad(2)->Divide(1,2, 0);
    cout << "=== Plotting ABPS" << endl;
    c->GetPad(2)->cd(1);
    c->GetPad(2)->GetPad(1)->SetLogx();
    //c->GetPad(2)->GetPad(1)->SetGridx();
    c->GetPad(2)->GetPad(1)->SetTicks();
    vector<TH1 *> ABPS = make_ABPS(rebinned_RM);
    ABPS.front()->SetStats(0);
    ABPS.front()->Draw("hist");
    ABPS.front()->Write();
    for (unsigned short i = 1 ; i < ABPS.size() ; i++)
    {
        ABPS[i]->Draw("same hist");
        ABPS[i]->Write();
    }
    c->GetPad(2)->GetPad(1)->BuildLegend(0.4,0.2,0.6,0.5);
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br") 
    TPaveText * ABPS_text  = new TPaveText(0.6, 0.25, 0.89, 0.65, "NBNDC");
    ABPS_text->SetFillColorAlpha(0,0);
    ABPS_text->SetTextFont(42);
    for (const TString& line: pave_ABPS_eff)
        ABPS_text->AddText(line);
    ABPS_text->Draw();
    draw_turnon_line(turnon, ABPS.front()->GetYaxis()->GetXmin(), ABPS.front()->GetYaxis()->GetXmax());
    c->GetPad(2)->GetPad(1)->RedrawAxis();

    pair<TH1 *, TH1 *> miss_and_fake = make_miss_fake(rebinned_RM, rebinned_gen, rebinned_rec);

    cout << "=== Plotting spectra" << endl;
    c->GetPad(2)->cd(2);
    c->GetPad(2)->GetPad(2)->SetLogx();
    c->GetPad(2)->GetPad(2)->SetLogy();
    //c->GetPad(2)->GetPad(2)->SetGridx();
    c->GetPad(2)->GetPad(2)->SetTicks();

    // the following vector will contain all the spectra at hadron level: the generated spectrum and the unfolded spectra
    vector<TH1 *> unfolded_spectra = make_unfolding(rebinned_RM, rebinned_rec, semi_rebinned_RM, h_rec, v_parameters);

    rebinned_rec->GetYaxis()->SetMoreLogLabels();
    rebinned_rec->Scale(1,"width");
    rebinned_gen->Scale(1,"width");
    rebinned_rec->DrawCopy();
    rebinned_gen->DrawCopy("same");
    rebinned_rec->Write();
    rebinned_gen->Write();

    miss_and_fake.second->Scale(1,"width");
    miss_and_fake.first ->Scale(1,"width");
    miss_and_fake.second->DrawCopy("same hist");
    miss_and_fake.first ->DrawCopy("same hist");
    miss_and_fake.second->Write();
    miss_and_fake.first ->Write();

    for (auto& unfolded_spectrum: unfolded_spectra)
    {
        unfolded_spectrum->Scale(1,"width");
        unfolded_spectrum->DrawCopy("same");
        unfolded_spectrum->Write();
    }
    c->GetPad(2)->GetPad(2)->BuildLegend(0.3, 0.3, 0.5, 0.6);
    draw_turnon_line(turnon, rebinned_rec->GetYaxis()->GetXmin(), rebinned_rec->GetYaxis()->GetXmax());
    c->GetPad(2)->GetPad(2)->RedrawAxis();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * spectrum_text = new TPaveText(0.7,0.2,0.89,0.3, "NBNDC");
    spectrum_text->SetFillColorAlpha(0,0);
    spectrum_text->SetTextFont(42);
    spectrum_text->AddText(true_xsec);
    spectrum_text->Draw();

    cout << "=== Plotting ratios of unfolded spectra to truth" << endl;
    c->cd(4);
    vector<TH1 *> numerators = {rebinned_gen};
    numerators.insert(end(numerators), begin(unfolded_spectra), end(unfolded_spectra));
    c->GetPad(4)->Divide(1,numerators.size(), 0,0);
    for (unsigned short i = 0 ; i < numerators.size() ; i++)
    {
        c->GetPad(4)->cd(i+1);
        c->GetPad(4)->GetPad(i+1)->SetLogx();
        //c->GetPad(4)->GetPad(i+1)->SetGridx();
        c->GetPad(4)->GetPad(i+1)->SetTicks();
        TH1D * ratio = static_cast<TH1D*>(numerators[i]->Clone(TString::Format("%s/rec", numerators[i]->GetName())));
        ratio->Divide(rebinned_rec);
        ratio->SetTitle(TString::Format(";;%s", numerators[i]->GetTitle()));
        ratio->GetYaxis()->SetTitleSize(0.2);
        ratio->GetYaxis()->SetTitleOffset(0.2);
        ratio->GetYaxis()->SetLabelSize(0.1);
        ratio->SetNdivisions(304, "Y");
        ratio->SetTickLength (0.1, "X");
        ratio->SetStats(0);
        ratio->GetYaxis()->SetRangeUser(0,2);
        ratio->DrawCopy();
        ratio->Write();
        draw_turnon_line(turnon, ratio->GetYaxis()->GetXmin(), ratio->GetYaxis()->GetXmax());
        delete ratio;
    }

    //delete rebinned_gen;
    //delete rebinned_rec;
    //delete rebinned_RM;
    //delete h_resolution;
    //delete miss_and_fake.first ;
    //delete miss_and_fake.second;
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

