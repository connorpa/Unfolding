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

