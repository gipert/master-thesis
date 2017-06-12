void writePGF( TH1* h ) {

	h->Rebin(4);
	string s;
	cout << "Insert output filename: ";
	cin >> s;
	s = "../tex/data/" + s;
	ofstream f(s.c_str());
	for ( int i = 1; i <= 1875; ++i ) {
		if ( h->GetBinContent(i) == 0 ) {
			f << h->GetBinLowEdge(i) << '\t' << -10 << '\n';
		}
		else f << h->GetBinLowEdge(i) << '\t' << TMath::Log10(h->GetBinContent(i)) << '\n';
	}
}
