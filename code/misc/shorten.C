// shorten histograms from 8500 to 7500 bins
//
{
    std::string name;
    std::cout << "Enter input filename: ";
    std::cin >> name;
    TFile infile(name.c_str(), "READ");
    std::cout << "Enter output filename: ";
    std::cin >> name;
    TFile outfile(name.c_str(), "RECREATE");
    std::vector<TH1F*> hist;
    std::vector<TH1F>  histOut;
    for ( int i = 0; i < 40; ++i ) {
        hist.push_back(dynamic_cast<TH1F*>(infile.Get(Form("h%i",i))));
        histOut.emplace_back(Form("h%i",i), Form("h%i",i), 7500, 0, 7500);
        for ( int j = 1; j <= 7500; ++j ) histOut[i].SetBinContent(j, hist[i]->GetBinContent(j));
        histOut[i].Write();
    }
    infile.Close();
    outfile.Close();
    gApplication->Terminate();
}

