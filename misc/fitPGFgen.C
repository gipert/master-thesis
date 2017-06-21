{
    std::cout << "Select BEGe or COAX [b/c]: ";
    char a;
    std::cin >> a;
    if ( a != 'b' and a != 'c' ) { std::cout << "bad input!\n"; gApplication->Terminate(); }
    std::cout << "Output filename: ";
    std::string name;
    std::cin >> name;
    std::ofstream fout(name.c_str());
    TFile::Open("outBATHist.root");

    if ( a == 'c' ) fout << "energy\tmidenergy\tdata\tsum\t2nbb\tK42homLAr\tK40fibers\tBi212Tl208fibers\tPb214Bi214fibers\tAc228holder\tCo60holder\tK40holder\tPb214Bi214holder\tK40cables\tBi212Tl208cables\tBi207minishroud\tPa234cables\talpha\tK42nPlus\n";
    else fout << "energy\tmidenergy\tdata\tsum\t2nbb\tK42homLAr\tK40fibers\tBi212Tl208fibers\tPb214Bi214fibers\tAc228holder\tCo60holder\tK40holder\tPb214Bi214holder\tK40cables\tBi212Tl208cables\tBi207minishroud\tPa234cables\talpha\n";

    std::vector<TH1D*> h;
    if ( a == 'b' ) {
        h.push_back(hDataBEGe);
        h.push_back(hsumBEGe);
        h.push_back(h2nbbBEGe);
        h.push_back(hK42homLArBEGe);
        h.push_back(hK40fibersBEGe);
            hBi212fibersBEGe->Add(hTl208fibersBEGe);
        h.push_back(hBi212fibersBEGe);
            hPb214fibersBEGe->Add(hBi214fibersBEGe);
        h.push_back(hPb214fibersBEGe);
        h.push_back(hAc228holderBEGe);
        h.push_back(hCo60holderBEGe);
        h.push_back(hK40holderBEGe);
            hPb214holderBEGe->Add(hBi214holderBEGe);
        h.push_back(hPb214holderBEGe);
        h.push_back(hK40cablesBEGe);
            hBi212cablesBEGe->Add(hTl208cablesBEGe);
        h.push_back(hBi212cablesBEGe);
        h.push_back(hBi207minishroudBEGe);
        h.push_back(Pa234cablesBEGe);
        h.push_back(hAlphaBEGe);
    }

    if ( a == 'c' ) {
        h.push_back(hDataCOAX);
        h.push_back(hsumCOAX);
        h.push_back(h2nbbCOAX);
        h.push_back(hK42homLArCOAX);
        h.push_back(hK40fibersCOAX);
            hBi212fibersCOAX->Add(hTl208fibersCOAX);
        h.push_back(hBi212fibersCOAX);
            hPb214fibersCOAX->Add(hBi214fibersCOAX);
        h.push_back(hPb214fibersCOAX);
        h.push_back(hAc228holderCOAX);
        h.push_back(hCo60holderCOAX);
        h.push_back(hK40holderCOAX);
            hPb214holderCOAX->Add(hBi214holderCOAX);
        h.push_back(hPb214holderCOAX);
        h.push_back(hK40cablesCOAX);
            hBi212cablesCOAX->Add(hTl208cablesCOAX);
        h.push_back(hBi212cablesCOAX);
        h.push_back(hBi207minishroudCOAX);
        h.push_back(Pa234cablesCOAX);
        h.push_back(hAlphaCOAX);
        h.push_back(hK42nPlusCOAX);
    }

    int dbin = hsumBEGe->GetXaxis()->FindBin(570);
    int ubin = hsumBEGe->GetXaxis()->FindBin(5300);
    int nbins = ubin - dbin;
    for ( int i = dbin-1; i <= ubin+1; ++i ) {
        fout << hsumBEGe->GetBinLowEdge(i) << '\t';
        fout << hsumBEGe->GetBinCenter(i) << '\t';
        for ( int j = 0; j < (int)h.size(); ++j ) fout << safelog(h[j]->GetBinContent(i)/h[j]->GetBinWidth(i)) << '\t';
        fout << '\n';
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i-dbin+1 << "/" << nbins << std::flush;
    }
    gApplication->Terminate();
}

double safelog( double a ) {
    if ( a <= 0 ) return -20;
    else return TMath::Log(a);
}
