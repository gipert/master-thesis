// works only with 4 keV binning!
void BI( TH1* h, double exposure ) {

    double intg = h->Integral(458, 523) - h->Integral(500, 502) - h->Integral(504, 506) - h->Integral(479, 491);
    cout << "Counts: " << intg << endl;
    int nbins = 47;
    cout << "range: " << nbins*4 << " keV\n";
    cout << "BI: " << intg/(exposure*nbins*4) << " cts/(keV*kg*yr)\n";
    return;
}
