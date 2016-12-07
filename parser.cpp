#include <iostream>
#include <numeric>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <TF1.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>

#include "spectrum.h"
#include "efficiency.h"
#include "resolution.h"
#include "unfolding.h"

using namespace std;
namespace pt = boost::property_tree;


template<typename T, typename U> vector<T> convert_vector (vector<U> input)
{
    return vector<T>(input.begin(), input.end());
}
template<> vector<int> convert_vector<int, string> (vector<string> input)
{
    vector<int> output;
    transform(begin(input), end(input), back_inserter(output), [](string d) { return stoi(d); });
    return output;
}
template<> vector<double> convert_vector<double, string> (vector<string> input)
{
    vector<double> output;
    transform(begin(input), end(input), back_inserter(output), [](string d) { return stod(d); });
    return output;
}

#define DEFAULT_SEPARATOR ','
template<typename T> vector<T> get_elements (string s, const char & c = DEFAULT_SEPARATOR)
{
    vector<string> elements;
    size_t next_separator;
    do
    {
        while (s[0] == ' ')
            s = s.substr(1);
        next_separator = s.find(c);
        elements.push_back(s.substr(0,next_separator));
        if (next_separator < s.length())
            s = s.substr(next_separator+1); 
    }
    while (next_separator != string::npos);
    return convert_vector<T, string>(elements);
}

template<typename T> vector<T> get_vector_from_range (T minimum, T maximum, T step)
{
    if (step == 0)
        throw string("Step cannot be zero.");
    if (maximum < minimum)
        throw string("Minimum cannot be small than maximum.");

    int nelements = int((maximum-minimum)/step) + 1;
    vector<T> output(nelements, minimum);
    partial_sum(output.begin(), output.end(), output.begin(), [&](T x, T y) { return x+step; });
    
    return output;
}

template<typename T> vector<T> get_vector_from_range (string s, const char & c)
{
    vector<T> min_max_step = get_elements<T>(s, c);

    switch (min_max_step.size())
    {
    case 1:
        return {min_max_step[0]};
        break;
    case 2:
        return get_vector_from_range<T>(min_max_step[0], min_max_step[1], 1);
        break;
    case 3:
        return get_vector_from_range<T>(min_max_step[0], min_max_step[1], min_max_step[2]);
        break;
    default:
        throw string("Parser does not understand the given list: " + s + "\nPlease use the following syntax: min:max[:step]");
    }
}

template<typename T> vector<T> get_vector (string s)
{
    if (s.find(':') != string::npos)
        return get_vector_from_range<T>(s, ':');
    else
        return get_elements<T>(s, ',');
}

vector<TF1 *> get_functions (string s, double (* f)(double *, double *), double xmin, double xmax)
{
    cout << "Getting functions with following parameters:";
    vector<string> expressions = get_elements<string>(s, ';');
    vector<TF1 *> functions;
    for (string expression: expressions)
    {
        vector<double> parameters = get_elements<double>(expression, ',');
        replace(expression.begin(), expression.end(), ',', '_');
        functions.push_back(new TF1 (expression.c_str(), f, xmin, xmax, parameters.size()));
        cout << '\n';
        for (unsigned short iparam = 0 ; iparam < parameters.size() ; iparam++)
        {
            cout << '\t' << parameters[iparam];
            functions.back()->SetParameter(iparam+1, parameters[iparam]);
        }
    }
    cout << endl;
    return functions;
}


TCanvas * divider (const vector<string> lines, Color_t fillcolor = kGray)
{
    TCanvas * c = new TCanvas("divider", "Unfolding");
    c->Draw();
    c->cd();
    // TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
    TPaveText * pave  = new TPaveText(0.2, 0.2 , 0.8, 0.8, "NBNDC");
    pave->SetFillColor(fillcolor);
    pave->SetTextFont(42);
    for (const string line: lines)
        pave->AddText(line.c_str());
    pave->Draw();
    return c;
}
TCanvas * divider (string line, Color_t fillcolor = kGray)
{
    vector<string> lines = {line};
    return divider(lines, fillcolor);
}


int main (int argc, char * argv[])
{
    if (argc == 1)
    {
        cout << "Please specify a config file." << endl;
        return EXIT_FAILURE;
    }

    gROOT->SetBatch();
    gStyle->SetOptTitle(0);
    TH1D::SetDefaultSumw2(true); // correct error computation
    TColor::CreateColorWheel();

    for (int iarg = 1 ; iarg < argc ; ++iarg)
    {
        string ini_file = argv[iarg];
        cout << "Reading " << ini_file << endl;
        pt::ptree config;
        string rootfilename, filename, title, format;
        vector<double> v_eff_a, v_eff_mu, v_eff_sigma;
        vector<double> v_res_N, v_res_mu, v_res_tau, v_res_kL, v_res_kR, v_res_aL, v_res_nL, v_res_aR, v_res_nR;
        vector<TF1 *> v_res_sigma;
        vector<double> binning;
        vector<double> v_minS, v_minP;
        vector<string> v_sampling;
        TF1 * truth,
            * model;
        long model_nevents, truth_nevents;
        vector<UnfoldingParameters *> v_parameters;

        cout << "=== Parsing ini file" << endl;
        try
        {
            pt::ini_parser::read_ini(ini_file, config);

            // TODO: use default values if not found in the ini file
            // TODO: dump (either cout-ing or in the output file)

            {
                cout << ">> Extracting default binning and requirements on purity and stability" << endl;
                auto binning_section = config.get_child("binning");
                binning = get_vector<double>(binning_section.get<string>("default"));
                v_minS = get_vector<double>(binning_section.get<string>("min_stability"));
                v_minP = get_vector<double>(binning_section.get<string>("min_purity"));
            }
            double xmin = binning.front(), xmax = binning.back();

            {
                cout << ">> Extracting model" << endl;
                auto model_section = config.get_child("model");
                model = new TF1 ("model", Spectrum, xmin, xmax, NTERMS);
                model->SetTitle(model_section.get<string>("title").c_str());
                vector<double> parameters = get_elements<double>(model_section.get<string>("parameters"));
                for (unsigned short iparam = 0 ; iparam < NTERMS ; iparam++)
                    model->SetParameter(iparam, 0);
                for (unsigned short iparam = 0 ; iparam < parameters.size() ; iparam++)
                    model->SetParameter(iparam, parameters[iparam]);
                v_sampling = get_elements<string>(model_section.get<string>("sampling"), ',');
                model_nevents = model_section.get<long>("nevents");
            }

            {
                cout << ">> Extracting truth" << endl;
                auto truth_section = config.get_child("truth");
                truth = new TF1 ("truth", Spectrum, xmin, xmax, NTERMS);
                truth->SetTitle(truth_section.get<string>("title").c_str());
                vector<double> parameters = get_elements<double>(truth_section.get<string>("parameters"));
                for (unsigned short iparam = 0 ; iparam < NTERMS ; iparam++)
                    truth->SetParameter(iparam, 0);
                for (unsigned short iparam = 0 ; iparam < parameters.size() ; iparam++)
                    truth->SetParameter(iparam, parameters[iparam]);
                truth_nevents = truth_section.get<long>("nevents");
            }

            {
                cout << ">> Extracting output information" << endl;
                auto output_section = config.get_child("output");
                title = output_section.get<string>("title");
                format = output_section.get<string>("format");
                filename = output_section.get<string>("filename") + '.' + format;
                rootfilename = output_section.get<string>("filename") + ".root";
            }

            {
                cout << ">> Extracting resolution parameters" << endl;
                auto resolution_section = config.get_child("resolution");
                v_res_N = get_vector<double>(resolution_section.get<string>("N")); // TODO: function
                v_res_mu = get_vector<double>(resolution_section.get<string>("mu")); // TODO: function
                v_res_sigma = get_functions(resolution_section.get<string>("sigma"), Sigma, xmin, xmax);
                v_res_tau = get_vector<double>(resolution_section.get<string>("tau")); // TODO: function
                v_res_kL = get_vector<double>(resolution_section.get<string>("kL")); // ... ?
                v_res_aL = get_vector<double>(resolution_section.get<string>("aL"));
                v_res_nL = get_vector<double>(resolution_section.get<string>("nL"));
                v_res_kR = get_vector<double>(resolution_section.get<string>("kR"));
                v_res_aR = get_vector<double>(resolution_section.get<string>("aR"));
                v_res_nR = get_vector<double>(resolution_section.get<string>("nR"));
            }

            {
                cout << ">> Extracting efficiency function" << endl;
                auto efficiency_section = config.get_child("efficiency");
                v_eff_a = get_vector<double>(efficiency_section.get<string>("a"));
                v_eff_mu = get_vector<double>(efficiency_section.get<string>("mu"));
                v_eff_sigma = get_vector<double>(efficiency_section.get<string>("sigma"));
            }

            {
                cout << ">> Extracting unfolding parameters" << endl;
                auto unfolding_section = config.get_child("unfolding");
                for (auto& key: unfolding_section)
                {
                    if (key.first == "binbybin")
                    {
                        if (key.second.get<bool>(""))
                            v_parameters.push_back(new UnfoldingParameters(UnfoldingParameters::kBinByBin));
                    }
                    else if (key.first == "inversion")
                    {
                        if (key.second.get<bool>(""))
                            v_parameters.push_back(new UnfoldingParameters(UnfoldingParameters::kInversion));
                    }
                    else if (key.first == "SVD")
                    {
                        vector<int> v_kreg = get_elements<int>(key.second.get<string>(""), ';');
                        for (int& kreg: v_kreg)
                            v_parameters.push_back(new UnfoldingParametersSVD(UnfoldingParameters::kSVD, kreg));
                    }
                    else if (key.first == "Bayes")
                    {
                        vector<string> v_niterations_rm_type = get_elements<string>(key.second.get<string>(""), ';');
                        for (string& string_niterations_rm_type: v_niterations_rm_type)
                        {
                            vector<string> niterations_rm_type = get_elements<string>(string_niterations_rm_type, ',');
                            const int niterations = stoi(niterations_rm_type[0]);

                            UnfoldingParametersBayes::RMtype rm_type;
                            if (niterations_rm_type[1] == "square")
                                rm_type = UnfoldingParametersBayes::kSquare;
                            else if (niterations_rm_type[1] == "fine")
                                rm_type = UnfoldingParametersBayes::kFine;
                            else throw string("Unknown parameter for Bayesian unfolding: ") + niterations_rm_type[1];

                            v_parameters.push_back(new UnfoldingParametersBayes(UnfoldingParameters::kBayes, niterations, rm_type));
                        }
                    }
                    else if (key.first == "TUnfold")
                    {
                        vector<string> v_regularisation = get_elements<string>(key.second.get<string>(""), ';');
                        for (string& regularisation: v_regularisation)
                        {
                            TUnfold::ERegMode eregmode;
                                 if (regularisation == "none"      ) eregmode = TUnfold::kRegModeNone      ; // none
                            else if (regularisation == "size"      ) eregmode = TUnfold::kRegModeSize      ; // minimize the size of (x-x0)
                            else if (regularisation == "derivative") eregmode = TUnfold::kRegModeDerivative; // minimize the 1st derivative of (x-x0)
                            else if (regularisation == "curvature" ) eregmode = TUnfold::kRegModeCurvature ; // minimize the 2nd derivative of (x-x0) 
                            else throw string("Unknown parameter for TUnfold: ") + regularisation;

                            v_parameters.push_back(new UnfoldingParametersTUnfold(UnfoldingParameters::kTUnfold, eregmode));
                        }
                    }
                    else throw string("Unknown key in unfolding section: ") + key.first;
                }
            }
        }
        catch (const pt::ptree_error & e)
        {
            cerr << "Property tree: " << e.what() << "\nIgnoring " << ini_file << endl;
            continue;
        }
        catch (string s)
        {
            cerr << s << endl;
            continue;
        }


        cout << "=== Looping over parameters" << endl;
        TCanvas * c;
        TFile * f_output;
        try
        {
            c = divider(get_elements<string>(title, ';'));
            c->Draw();
            c->Print((filename + '(').c_str());
            f_output = new TFile (rootfilename.c_str(), "RECREATE");
            
            double xmin = binning.front(), xmax = binning.back();
            // loop over the resolution parameters
            for (const double res_N: v_res_N) for (const double res_mu: v_res_mu) for (TF1 * res_sigma: v_res_sigma) for (const double res_tau: v_res_tau) for (const double res_kL: v_res_kL) for (const double res_kR: v_res_kR) for (const double res_aL: v_res_aL) for (const double res_nL: v_res_nL) for (const double res_aR: v_res_aR) for (const double res_nR: v_res_nR)
            {
                cout << "-- resolution parameters:" << res_N << '\t' << res_mu << '\t' << res_sigma->Eval((xmax-xmin)/2) << '\t' << res_tau << '\t'
                                                    << res_kR << '\t' << res_kL << '\t' << res_aL << '\t' << res_nL << '\t' << res_aR << '\t' << res_nR << endl;

                // define new directory
                string dir_resolution_name = "resolution_" + to_string(res_N) + '_' + to_string(res_mu) + '_' + string(res_sigma->GetName()) + '_' + to_string(res_tau) + "__" + to_string(res_kL) + '_' + to_string(res_kR) + "__" + to_string(res_aL) + '_' + to_string(res_nL) + '_' + to_string(res_aR) + '_' + to_string(res_nR);
                TDirectory * dir_resolution = f_output->mkdir(dir_resolution_name.c_str());
                dir_resolution->cd();

                // define and save resolution from parameters
                TF1 * f_resolution = new TF1 ("resolution", Resolution, -1, 1, NRESOLUTION_PARAMETERS);
                f_resolution->SetParameters(res_N, res_mu, res_sigma->Eval((xmax-xmin)/2), res_tau, res_kL, res_kR, res_aL, res_nL, res_aR, res_nR);
                f_resolution->Write();

                // save parameter functions (if applicable)
                res_sigma->Write();

                // make and save measurement with default binning
                TH1D * h_gen = new TH1D ("gen", "Truth"      , binning.size()-1, &binning[0]),
                     * h_rec = new TH1D ("rec", "Measurement", binning.size()-1, &binning[0]);
                make_measurement(h_gen, h_rec, truth, /*f_efficiency,*/ xmin, xmax, f_resolution, res_sigma, truth_nevents);
                h_gen->Write();
                h_rec->Write();

                // loop on sampling methods 
                for (const string& sampling: v_sampling)
                {
                    cout << "-- sampling method: " << sampling << endl;
                    TDirectory * dir_sampling = dir_resolution->mkdir(sampling.c_str());
                    dir_sampling->cd();

                    // generate and save RM and differential resolution
                    TH2 * h_RM = new TH2D ("RM", "RM", binning.size()-1, &binning[0], binning.size()-1, &binning[0]),
                        * h_resolution = new TH2D("resolution", "resolution", 41, -1, 1, binning.size()-1, &binning[0]);
                    make_RM(h_RM, h_resolution, model, xmin, xmax, f_resolution, res_sigma, model_nevents, sampling.c_str());
                    h_RM->Write();
                    h_resolution->Write();

                    // loop on conditions on the binning
                    for (const double& minS: v_minS) for (const double& minP: v_minP)
                    {
                        cout << "-- requirement on binning: " << minS << '\t' << minP << endl;
                        string dir_binning_name = "binning_" + to_string(minS) + '_' + to_string(minP);
                        TDirectory * dir_binning = dir_sampling->mkdir(dir_binning_name.c_str());
                        dir_binning->cd();

                        // looping over efficiency parameters
                        for (const double& eff_a: v_eff_a) for (const double& eff_mu: v_eff_mu) for (const double& eff_sigma: v_eff_sigma)
                        {
                            cout << "-- efficiency parameters:" << eff_a << '\t' << eff_mu << '\t' << eff_sigma << endl;

                            // define new directory
                            string dir_efficiency_name = "efficiency_" + to_string(eff_a) + '_' + to_string(eff_mu) + '_' + to_string(eff_sigma);
                            replace(dir_efficiency_name.begin(), dir_efficiency_name.end(), '.', '_');
                            TDirectory * dir_efficiency = dir_binning->mkdir(dir_efficiency_name.c_str());
                            dir_efficiency->cd();

                            // define and save efficiency function
                            TF1 * f_efficiency = new TF1 ("efficiency", Efficiency, xmin, xmax, NEFFICIENCY_PARAMETERS);
                            f_efficiency->SetParameters(eff_a, eff_mu, eff_sigma);
                            f_efficiency->Write();
                            
                            // determine turn-on point
                            unsigned iturnon = binning.front();
                            if (eff_sigma > 0)
                                while (iturnon < binning.size() && f_efficiency->Eval(binning[iturnon]) < 0.99) iturnon++;
                            double turnon = binning[iturnon];

                            string new_name = dir_efficiency_name;
                            new_name.replace(new_name.find("efficiency"), 10, h_rec->GetName());
                            TH1 * h_rec_after_trigger = static_cast<TH1 *>(h_rec->Clone(new_name.c_str()));
                            if (eff_sigma > 0) h_rec_after_trigger->Multiply(f_efficiency);

                            // parameter information
                            vector<TString> pave_ABPS = {TString::Format("minS=%.2f", minS),
                                                         TString::Format("minP=%.2f", minP)},
                                            pave_trigger = {TString::Format("a=%.2f", eff_a),
                                                            TString::Format("#mu=%.2f", eff_mu),
                                                            TString::Format("#sigma=%.2f", eff_sigma),
                                                            TString::Format("turn-on=%.2f", turnon)},
                                            pave_resolution = {TString::Format("#sigma=%.3f+#frac{%.3f}{p_{T}^{%.3f}+%.3f*p_{T}}+%.3f*p_{T}", res_sigma->GetParameter(1), res_sigma->GetParameter(2), res_sigma->GetParameter(3), res_sigma->GetParameter(4), res_sigma->GetParameter(5))};
                            try
                            {
                                vector<double> new_edges = find_binning(h_RM, minS, minP);
                                c = make_canvas(h_gen, h_rec_after_trigger, h_RM, h_resolution, new_edges, truth->GetTitle(), model->GetTitle(), v_parameters, sampling.c_str(), pave_resolution, pave_ABPS, pave_trigger, turnon); // writing is done inside of the function
                            }
                            catch (TString s)
                            {
                                throw string(s.Data());
                            }
                            catch (string s)
                            {
                                c = divider(s, kWhite);
                            }
                            const string current_title = "Title:" + sampling;
                            c->Print(filename.c_str(), current_title.c_str());
                            c->Write();
                            dir_binning->cd();
                            dir_efficiency->Close();
                        }
                        dir_sampling->cd();
                        dir_binning->Close();
                    }
                    dir_resolution->cd();
                    dir_sampling->Close();
                }
                f_output->cd();
                dir_resolution->Close();
            } // end of loop over parameters

            // TODO: measure computation time

            c = new TCanvas ("back"); // should not be displayed
            c->Print((filename + ')').c_str()); 
            f_output->Close();
        }
        catch (TString s)
        {
            throw string(s.Data());
        }
        catch (string s)
        {
            if (c != nullptr)
                c->Print((filename + ')').c_str());
            if (f_output != nullptr)
                f_output->Close();
            cerr << s << '\n'
                 << "Continuing." << endl;
        }
    } // end of looping over the ini files

    return EXIT_SUCCESS;
}
