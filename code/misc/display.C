void display( TH1* h ) {
   h->Rebin(4);
   h->GetXaxis()->SetRangeUser(570,2500);
   h->Draw("histo");
   gPad->SetLogy();
}
