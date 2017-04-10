/* sumbbkgext.cxx
 *
 * Program to sum MaGe simulations for external bkg sources
 * into global energy spectra to be used in the fit
 *
 * supported sources:
 *   
 *   in LAr (homogeneous)
 *   K42
 *  
 *   on Fibers:
 *   K40 Bi212 Bi214 
 *   Pb214 Tl208
 *
 *   on contacts (p and n):
 *   K42
 *
 *   on Holder
 *   K40 Ac228 Bi212 
 *   Bi214 Co60 Pb214
 *   Tl208
 *
 * state-of-arts: enrCOAX are not summed. Different live times 
 * of the detectors (depending on the considered runs) are taken 
 * into account.
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 15/03/2017
 *
 */

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <map>

#include "TFile.h"
#include "TH1F.h"

#include "DataReader.h"

int main( int argc, char** argv ) {
    
    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; ++i ) args[i] = argv[i];
 
    if ( std::find(args.begin(), args.end(), "--help"   ) != args.end() ) {
        std::cout << "Available sources:\n\n"
                  << "homLAr:\n"
                  << "    K42\n\n"
                  << "on fibers:\n"
                  << "    K40, Bi212, Bi214, Pb214, Tl208\n\n"
                  << " on contacts (p and n):\n"
                  << "    K42\n\n"
                  << " on holder:\n"
                  << "    K40, Ac228, Bi212, Bi214, Co60, Pb214, Tl208\n";
        return 0;
    }
   
    std::string phys, place;
    if      ( std::find(args.begin(), args.end(), "--homLAr"   ) != args.end() ) place = "homLAr";
    else if ( std::find(args.begin(), args.end(), "--fibers"   ) != args.end() ) place = "fibers";
    else if ( std::find(args.begin(), args.end(), "--contacts" ) != args.end() ) place = "contacts";
    else if ( std::find(args.begin(), args.end(), "--holder"   ) != args.end() ) place = "holder";
    else { std::cout << "Please pecify place: --homLAr, --fibers, --contacts, --holder\n"; return -1; }
    
    if      ( std::find(args.begin(), args.end(), "--K42" ) != args.end() ) phys = "K42";
    else if ( std::find(args.begin(), args.end(), "--K40" ) != args.end() ) phys = "K40";
    else if ( std::find(args.begin(), args.end(), "--Bi212" ) != args.end() ) phys = "Bi212";
    else if ( std::find(args.begin(), args.end(), "--Tl208" ) != args.end() ) phys = "Tl208";
    else if ( std::find(args.begin(), args.end(), "--Bi214" ) != args.end() ) phys = "Bi214";
    else if ( std::find(args.begin(), args.end(), "--Pb214" ) != args.end() ) phys = "Pb214";
    else if ( std::find(args.begin(), args.end(), "--Ac228" ) != args.end() ) phys = "Ac228";
    else if ( std::find(args.begin(), args.end(), "--Co60" ) != args.end() ) phys = "Co60";
    else    { std::cout << "Please specify source: --K42, --K40, --Bi212, --Tl208, --Bi214, --Pb214, --Ac228, --Co60\n"; return -1; }
// ----------------------------------------------------------------------------------------------------------
    
    std::string rootpath = std::string(std::getenv("GERDACPTDIR"));

    // infos about runs
    GERDA::DataReader reader( rootpath + "/misc/paths.txt", false, "GELATIO");
    // detector types
    GERDA::DetectorSet set("GELATIO");

    // get detectorStatusMap for selected runs
    auto dsm = reader.GetDetectorStatusMap();

    // get live times saved in out/results.dat and store [s]
    std::string path = rootpath + "/data/sumData.dat";
    std::ifstream timeFile(path.c_str());
    std::map<int,int> timeMap;
    int runID, time;
    while ( timeFile >> runID >> time ) timeMap.insert(std::make_pair(runID,time));

    ////// calculate the total time each detector is ON --> detector status = 0
    std::vector<int> totalTime(40, 0);
    for ( const auto& i : dsm ) {
        for ( int j = 0; j < 40; ++j ) {
            if ( i.second[j] == 0 ) totalTime[j] += timeMap[i.first];
        }
    }
    
    long int Ngen;
    double M;
    if ( place == "homLAr" ) {
        Ngen = 5E09;
        M = 24563.385; // LAr mass [kg]
    }
    else if ( place == "fibers" ) {
        Ngen = 1E08;
        M = 1.3615078; // fibers' mass [kg]
    }

    else if ( place == "holder" ) {
        Ngen = 1E07;
        M = 0.658988; // holder's mass [kg]
    }

// -------------------------------------------------------------------------------------------------------------------    
    // SPECIAL CASE: contacts
    else if ( place == "contacts" ) {
        
        // simulation are already normalized to number of generated events
        //long int Ngen;
        double MBEGe_n = 0.07423; // [Kg]
        double MBEGe_p = 0.001; // [Kg] <----------- FAKE
        double MCOAX_n = 0.29962; // [Kg]
        double MCOAX_p = 0.0021176; // [Kg]
     
        // calculate total time BEGe and COAX are on
        int ltBEGe = 0;
        int ltCOAX = 0;
        for ( int i = 0; i < 40; ++i ) {
            if      ( set.GetDetectorTypes()[i] == 3 ) continue;
            else if ( set.GetDetectorTypes()[i] == 2 ) ltCOAX += totalTime[i];
            else                                       ltBEGe += totalTime[i];
        }
 
        if ( phys == "K42" ) {
            // retrieve simulations
            path = rootpath + "/out/processed/K42_nPlus.root";
            TFile file_n(path.c_str(), "READ");
            path = rootpath + "/out/processed/K42_pPlus.root";
            TFile file_p(path.c_str(), "READ");
        
            TH1F *hBEGe_n, *hBEGe_p, *hCOAX_n, *hCOAX_p;

            file_n.GetObject("energy_BEGe", hBEGe_n);
            file_n.GetObject("energy_COAX", hCOAX_n);
            file_p.GetObject("energy_BEGe", hBEGe_p);
            file_p.GetObject("energy_COAX", hCOAX_p);
    
            path = rootpath + "/data/sumMaGe_K42nPlus.root";
            TFile outfile_n(path.c_str(), "RECREATE");
            path = rootpath + "/data/sumMaGe_K42pPlus.root";
            TFile outfile_p(path.c_str(), "RECREATE");

            hBEGe_n->Scale(MBEGe_n*ltBEGe);
            hBEGe_p->Scale(MBEGe_p*ltBEGe);
            hCOAX_n->Scale(MCOAX_n*ltCOAX);
            hCOAX_p->Scale(MCOAX_p*ltCOAX);

            outfile_n.WriteObject(hBEGe_n, "energy_BEGe");
            outfile_n.WriteObject(hCOAX_n, "energy_COAX");
            outfile_p.WriteObject(hBEGe_p, "energy_BEGe");
            outfile_p.WriteObject(hCOAX_p, "energy_COAX");
        }

        return 0;
    }
// -------------------------------------------------------------------------------------------------------------------    
   
    else { std::cout << "Error: " << place << ": unknown place\n"; return -1; }

    std::vector<TH1F*> hist;

    // get original spectra
    path = rootpath + "/out/processed/" + phys + place + ".root";
    TFile infile(path.c_str(), "READ"); 
    if (!infile.IsOpen()) { std::cerr << "Zombie infile!\n"; return -1; }

    for ( int i = 0; i < 40; ++i ) { 
        hist.push_back(dynamic_cast<TH1F*>(infile.Get(Form("h%i", i))));
        if (hist[i]->IsZombie()) { std::cout << "Zombie h" << std::to_string(i) << "!\n"; return -1; }
    }
    
    // save
    path = rootpath + "/data/sumMaGe_" + phys + place + ".root";
    TFile outfile(path.c_str(), "RECREATE");
    for ( int i = 0; i < 40; ++i ) {
        hist[i]->Scale(totalTime[i]*M/Ngen);
        hist[i]->Write();
    }

    TH1F histBEGe("energy_BEGe", "BEGe global MaGe energy spectrum", 7500, 0, 7500);
    TH1F histCOAX("energy_COAX", "COAX global MaGe energy spectrum", 7500, 0, 7500);
    TH1F histTotAll("energy_total", "global MaGe energy spectrum", 7500, 0, 7500);
       
    for ( auto& h : hist ) histTotAll.Add(h);
    
    for ( int i = 0; i < 40; ++i ) {
        if      ( set.GetDetectorTypes()[i] == 3 ) continue;
        else if ( set.GetDetectorTypes()[i] == 2 ) histCOAX.Add(hist[i]);
        else                                       histBEGe.Add(hist[i]);
    }
    
    histCOAX.Write();
    histBEGe.Write();
    histTotAll.Write();

    outfile.Close();
    infile.Close();

    return 0;
}
