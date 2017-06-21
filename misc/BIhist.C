// works only with 4 keV binning!
void BIhist4( TH1* h, double exposure ) {

    double intg = h->Integral(458, 523) - h->Integral(500, 502) - h->Integral(504, 506) - h->Integral(479, 491);
    cout << "Counts: " << intg << endl;
    int nbins = 47;
    cout << "range: " << nbins*4 << " keV\n";
    cout << "BI: " << intg/(exposure*nbins*4) << endl;
    return;
}
// works only with 1 keV binning!
void BIhist1( TH1* h, double exposure ) {

    double intg = h->Integral(1930, 2190) - h->Integral(2099, 2109) - h->Integral(2114, 2124) - h->Integral(2014, 2064);
    cout << "Counts: " << intg << endl;
    int nbins = 190;
    cout << "range: " << nbins << " keV\n";
    cout << "BI: " << intg/(exposure*nbins) << endl;
    return;
}
