#include <iostream>
#include <fstream>
#include <ROOT/RVec.hxx>
#include <vector>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>

ROOT::RVec<double> convertToDouble(const ROOT::RVec<float> &column) {

  ROOT::RVec<double> col_value_double(column.size());
    
  for (size_t i = 0; i < column.size(); ++i) {
      col_value_double[i] = static_cast<double>(column[i]);
  }
    
  return col_value_double;
}

std::vector<fastjet::PseudoJet> createPseudoJets(ROOT::RVec<double>& px, 
                                                 ROOT::RVec<double>& py, 
                                                 ROOT::RVec<double>& pz, 
                                                 ROOT::RVec<double>& e) {
  std::vector<fastjet::PseudoJet> jets;
 
  for (size_t i = 0; i < px.size(); ++i) {

    fastjet::PseudoJet jet(px[i], py[i], pz[i], e[i]);
      
    jets.push_back(jet);
  }

  return jets;
}

std::vector<fastjet::PseudoJet> applyAntiKtAlgorithm(std::vector<fastjet::PseudoJet>& particles, double R) {
   
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  std::vector<fastjet::PseudoJet> jets = jet_def(particles);

  return sorted_by_pt(jets);
}

void file_checker()
{

  ROOT::RDataFrame df("events", "events_029930132.root");

  auto df0 = df.Define("Px", convertToDouble, {"Particle.momentum.x"});
  auto df1 = df0.Define("Py", convertToDouble, {"Particle.momentum.y"});
  auto df2 = df1.Define("Pz", convertToDouble, {"Particle.momentum.z"});
  auto df3 = df2.Define("M", "Particle.mass").Define("E", "sqrt(Px*Px + Py*Py + Pz*Pz + M*M)");
  auto df4 = df3.Define("Status", "Particle.generatorStatus");
  auto df5 = df4.Define("Stable_Particles", [](ROOT::RVec<int>& Status, ROOT::RVec<double>& Px, ROOT::RVec<double>& Py, ROOT::RVec<double>& Pz, ROOT::RVec<double>& E){
			  ROOT::RVec<double> Stable_Px, Stable_Py, Stable_Pz, Stable_E;
			  for (size_t i = 0; i < Status.size(); i++){
			    if (Status[i] == 1) 
			      {Stable_Px.push_back(Px[i]);
				Stable_Py.push_back(Py[i]);
				Stable_Pz.push_back(Pz[i]);
				Stable_E.push_back(E[i]);} } return std::make_tuple(Stable_Px, Stable_Py, Stable_Pz, Stable_E);}, {"Status", "Px", "Py", "Pz", "E"}).Define("Stable_Px", "std::get<0>(Stable_Particles)").Define("Stable_Py", "std::get<1>(Stable_Particles)").Define("Stable_Pz", "std::get<2>(Stable_Particles)").Define("Stable_E", "std::get<3>(Stable_Particles)");
  auto df6 = df5.Define("PJs", [](ROOT::RVec<double>& Px, 
                                  ROOT::RVec<double>& Py,
                                  ROOT::RVec<double>& Pz,
                                  ROOT::RVec<double>& E){
			  return createPseudoJets(Px, Py, Pz, E);}, {"Stable_Px", "Stable_Py", "Stable_Pz", "Stable_E"});
  auto df7 = df6.Define("Jets", [](std::vector<fastjet::PseudoJet>& PJs){
				   return applyAntiKtAlgorithm(PJs, 0.8);}, {"PJs"});
  auto df8 = df7.Define("Dijet_M", [](std::vector<fastjet::PseudoJet>& Jets){
			               return (Jets[0] + Jets[1]).m();}, {"Jets"});
  auto df9 = df8.Define("Dijet_Pt", [](std::vector<fastjet::PseudoJet>& Jets){
			  return (Jets[0] + Jets[1]).pt();}, {"Jets"});
  auto rdf = df9.Define("Dijet_Eta", [](std::vector<fastjet::PseudoJet>& Jets){
			  return (Jets[0] + Jets[1]).eta();}, {"Jets"});

auto h1 = rdf.Histo1D({"Mass", "Dijet [anti-kt, R = 0.8]", 250u, 0., 250.}, "Dijet_M");
auto h2 = rdf.Histo1D({"Pt", "Dijet [anti-kt, R = 0.8]", 120u, 0., 120.}, "Dijet_Pt");
auto h3 = rdf.Histo1D({"Eta", "Dijet [anti-kt, R = 0.8]", 10u, -5., 5.}, "Dijet_Eta");

h1->GetXaxis()->SetTitle("M [GeV]");
h1->GetYaxis()->SetTitle("Number of events");

auto c1 = new TCanvas("c1");
h1->DrawCopy();  

h2->GetXaxis()->SetTitle("pT [GeV]");
h2->GetYaxis()->SetTitle("Number of events");

auto c2 = new TCanvas("c2");
h2->DrawCopy();  

h3->GetXaxis()->SetTitle("#eta");
h3->GetYaxis()->SetTitle("Number of events");

auto c3 = new TCanvas("c3");
h3->DrawCopy();  

}










 

 
