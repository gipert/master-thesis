void undo_normalization( std::string name ) {

    TFile f( name.c_str() );
    std::vector<TH1*> h;
    for ( int i = 0; i < 40; ++i ) h.push_back((TH1F*)f.Get(Form("h%i",i)));

    for ( int j = 0; j < 40; ++j ) {
        double min = 1000;
        int nbins = h[j]->GetNbinsX();
        double val;

        for ( int i = 1; i <= nbins; ++i ) {
            val = h[j]->GetBinContent(i);
            if ( val < min and val != 0 ) min = val;
        }
        h[j]->Scale(1./min);
    }

    name += "_new";
    TFile fout( name.c_str(), "RECREATE");
    for ( int j = 0; j < 40; ++j ) {
        h[j]->Write();
    }
    return;
}
