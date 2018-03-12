{    
    // CPT conserving 2nbb
    TF1 bb   ( "bbspectra"  , "[0]*(pow(x,4)+10*pow(x,3)+40*x*x+60*x+30)*x*pow([1]-x,5)", 0, 4 );
    // total CPT breaking 2nbb
    TF1 bbLV ( "bbLVspectra", "[0]*(pow(x,4)+10*pow(x,3)+40*x*x+60*x+30)*x*pow([1]-x,4)", 0, 4 );
    // partial CPT breaking 2nbb
    // TF1 bbPlusLV ( "bbPlusLVspectra", "[0]*(pow(x,4)+10*pow(x,3)+40*x*x+60*x+30)*x*(pow([1]-x,5)+10*[2]*pow([1]-x,4))", 0, 4);
    float C1 = 1, C2 = 1;

    bb.SetTitle("Electron sum spectra");
    bb.SetParName(0, "C");
    bb.SetParName(1, "Qbb");

    bb.SetParameter("C", C1);
    bb.SetParameter("Qbb", 2039./511);

    bbLV.SetParName(0, "C");
    bbLV.SetParName(1, "Qbb");
    //bbLV.SetParName(2, "a3of");

    bbLV.SetParameter("C", C2);
    bbLV.SetParameter("Qbb", 2039./511);
    //bbLV.SetParameter("a3of", 0.5);

    // normalization
    C1 *= 1./bb.Integral(0,4);
    C2 *= 1./bbLV.Integral(0,4);

    bb.SetParameter("C", C1);
    bbLV.SetParameter("C", C2);
    
    // style
    bb.SetLineColor(kBlue);
    bbLV.SetLineColor(kRed);

    // draw
    TCanvas c("c", "Spectra", 1);
    c.cd();

    bb.GetXaxis()->SetTitle("K/m_{e}");
    bb.GetYaxis()->SetTitle("d#Gamma/dK - a.u.");

    bb.Draw("L");
    bbLV.Draw("LSAME");
}
