{
	TFile file("sumData.root","READ");
	auto energyBEGeAll = (TH1D*)file.Get("energyBEGeAll");
	auto energyEnrCoaxAll = (TH1D*)file.Get("energyEnrCoaxAll");
	auto energyNatCoaxAll = (TH1D*)file.Get("energyNatCoaxAll");
	energyBEGeAll->Rebin(4);
	energyEnrCoaxAll->Rebin(4);
	energyNatCoaxAll->Rebin(4);

	ofstream f1("BEGe.dat");
	ofstream f2("EnrCoax.dat");
	ofstream f3("NatCoax.dat");

	for ( int i = 1; i <= 1875; ++i ) {
		if ( energyBEGeAll->GetBinContent(i) == 0 ) {
			f1 << energyBEGeAll->GetBinLowEdge(i) << '\t' << 0 << '\n';
		}
		else if ( energyBEGeAll->GetBinContent(i) == 1 ) {
			f1 << energyBEGeAll->GetBinLowEdge(i) << '\t' << TMath::Log10(2) << '\n';
		}
		else f1 << energyBEGeAll->GetBinLowEdge(i) << '\t' << TMath::Log10(2*energyBEGeAll->GetBinContent(i)) << '\n';
	}
	for ( int i = 1; i <= 1875; ++i ) {
		if ( energyEnrCoaxAll->GetBinContent(i) == 0 ) {
			f2 << energyEnrCoaxAll->GetBinLowEdge(i) << '\t' << 0 << '\n';
		}
		else if ( energyEnrCoaxAll->GetBinContent(i) == 1 ) {
			f2 << energyEnrCoaxAll->GetBinLowEdge(i) << '\t' << TMath::Log10(2) << '\n';
		}
		else f2 << energyEnrCoaxAll->GetBinLowEdge(i) << '\t' << TMath::Log10(2*energyEnrCoaxAll->GetBinContent(i)) << '\n';
	}
	for ( int i = 1; i <= 1875; ++i ) {
		if ( energyNatCoaxAll->GetBinContent(i) == 0 ) {
			f3 << energyNatCoaxAll->GetBinLowEdge(i) << '\t' << 0 << '\n';
		}
		else if ( energyNatCoaxAll->GetBinContent(i) == 1 ) {
			f3 << energyNatCoaxAll->GetBinLowEdge(i) << '\t' << TMath::Log10(2) << '\n';
		}
		else f3 << energyNatCoaxAll->GetBinLowEdge(i) << '\t' << TMath::Log10(2*energyNatCoaxAll->GetBinContent(i)) << '\n';
	}
}
