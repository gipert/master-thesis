{
    TFile f1("aof_folded.root");
    TH1D * hist_0_2nbbLV = (TH1D*)f1.Get("hist_0_2nbbLV");
    if (!hist_0_2nbbLV) std::cout << "Zombie hist1!\n";
    TH1D * hist_0_2nbbLV_out = (TH1D*)f1.Get("hist_0_2nbbLV_out");
    if (!hist_0_2nbbLV_out) std::cout << "Zombie hist2!\n";
    std::ofstream f("aof.dat");
    f << "aof\tcounts\tcountsmear\n";

    for ( int i = 1; i < hist_0_2nbbLV->GetNbinsX(); i++ ) {
        f << hist_0_2nbbLV->GetBinLowEdge(i) << '\t' << hist_0_2nbbLV->GetBinContent(i) << '\t'
          << hist_0_2nbbLV_out->GetBinContent(i) << '\n';
    }
    f.close();
    gApplication->Terminate();
}
