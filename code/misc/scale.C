void scale( std::string name, int nscale = 1000 ) {

    TFile infile(name.c_str(), "READ");
    name += "_scaled";
    TFile outfile(name.c_str(), "RECREATE");
    std::vector<TH1F*> hist;
    for ( int i = 0; i < 40; ++i ) {
        hist.push_back(dynamic_cast<TH1F*>(infile.Get(Form("h%i",i))));
        hist[i]->Scale(nscale);
        hist[i]->Write();
    }
    infile.Close();
    outfile.Close();
}
