void foldposterior() {

    // build the asymmetric gaussian distribution
    double sup = 0.9;
    double sdown = 0.1;

    TF1 gausdown("gausdown", "[0]*TMath::Gaus(x,[1],[2])", 0, 2);
        gausdown.SetNpx(1000);
        gausdown.SetParameter(0, 1);
        gausdown.SetParameter(1, 1);
        gausdown.SetParameter(2, sdown);
    TF1 gausup("gausup", "[0]*TMath::Gaus(x,[1],[2])", 0, 2);
        gausup.SetNpx(1000);
        gausup.SetParameter(0, 1);
        gausup.SetParameter(1, 1);
        gausup.SetParameter(2, sup);
        gausup.SetLineColor(kBlue);

    // normalize everything to one
    double norm = 1./(gausdown.Integral(0,1) + gausup.Integral(1,2));
    gausdown.SetParameter(0, norm);
    gausup.SetParameter(0, norm);
    double lowarea = gausdown.Integral(0,1);

    TFile file("aof_post.root", "UPDATE");
    TH1D * hin = (TH1D*)file.Get("hist_0_2nbbLV");

    TH1D * hout = (TH1D*)hin->Clone();
    hout->Reset();
    int nentries;
    double bincent;
    double u;
    TRandom3 rand(0);
    for ( int i = 1; i < hin->GetNbinsX(); i++ ) {
        nentries = hin->GetBinContent(i);
        bincent = hin->GetBinCenter(i);
        for ( int j = 0; j < nentries; j++ ) {
            u = rand.Uniform(0,1);
            if ( u < lowarea ) hout->Fill(bincent*rand.Gaus(1,sdown));
            else hout->Fill(bincent*rand.Gaus(1,sup));
        }
    }
    TFile fout("aof_folded.root", "RECREATE");
    gausdown.Write();
    gausup.Write();
    hout->Write("hist_0_2nbbLV_out");
    return;
}
