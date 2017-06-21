{
    TChain c("MarkovChainTree_0");
    c.Add("markowChains_0.root");
    c.Add("markowChains_0.root/MarkovChainTree_1");
    c.Add("markowChains_0.root/MarkovChainTree_2");
    c.Add("markowChains_0.root/MarkovChainTree_3");
    c.Add("markowChains_0.root/MarkovChainTree_4");
    c.Add("markowChains_1.root");
    c.Add("markowChains_1.root/MarkovChainTree_1");
    c.Add("markowChains_1.root/MarkovChainTree_2");
    c.Add("markowChains_1.root/MarkovChainTree_3");
    c.Add("markowChains_1.root/MarkovChainTree_4");
    c.Add("markowChains_2.root");
    c.Add("markowChains_2.root/MarkovChainTree_1");
    c.Add("markowChains_2.root/MarkovChainTree_2");
    c.Add("markowChains_2.root/MarkovChainTree_3");
    c.Add("markowChains_2.root/MarkovChainTree_4");

    double a;
    TH1D * h = new TH1D("2nbbpost", "2nbbpost", 100, 0.49, 0.53);
    TH1D * hp= new TH1D("2nbbprior", "2nbbprior", 100, 0.49, 0.53);
    TH1D * hpp= new TH1D("2nbbprior2", "2nbbprior2", 100, 0.49, 0.53);
    TF1 invflat("inverse-flat", "1/x^2"  , 0.49, 0.53);
    TF1 invgaus("inverse-gaus", "exp(-((1/x-[0])^2)/(2*[1]^2))/(x^2)", 0.49, 0.53);
    invgaus.SetParameter(0, 1.926);
    invgaus.SetParameter(1, 0.095);
    hp->SetLineColor(kRed);
    hpp->SetLineColor(kBlue);
    c.SetBranchAddress("Parameter0", &a);
    int entries = c.GetEntries();

    for ( int i = 0; i < entries; i++ ) {
        c.GetEntry(i);
        h->Fill(a/417.483);
        hp->Fill(invgaus.GetRandom(0.49,0.53));
        hpp->Fill(invflat.GetRandom(0.49,0.53));
    }

    ofstream f("2nbb_post.dat");
    f << "point\tpost\tprior\n";
    for ( int i = 1; i < h->GetNbinsX(); ++i ) f << h->GetBinLowEdge(i) << '\t' << h->GetBinContent(i) << '\t' << hp->GetBinContent(i) << '\n';
    f.close();
    h->Draw();
    hp->Draw("same");
    hpp->Draw("same");
    //gPad->SetLogy();
    //gApplication->Terminate();
}
