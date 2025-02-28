#ifndef MatthewTransparency_H
#define MatthewTransparency_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "DataModel.h"

class MatthewTransparency: public Tool {
  
 public:

  MatthewTransparency();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  static constexpr int N{2068};
  
 private:


  bool ok{true};
  bool toolactive{true};
  std::array<double, N> wavelengths;
  std::vector<std::array<double, N>> transparency_traces;
  std::vector<std::array<int, 6>> transparency_times;
  
  struct Transparency {
    std::array<double, N> values;
    std::array<int, 6> date_time;
  };
  
  std::array<double, N> PopulateWavelength(TTree*) const ;
  std::array<double, N> RetrievePureValues(const int, const std::string) const;
  std::array<double, N> RetrievePureValuesFromFile(const std::string&, const std::string&) const;
  void SaveToMonthlyFile(const Transparency&);
  std::array<int, 6> GetDateTimeVec() const;
  bool Ready() const;
  double GetValAtWavelength(const Transparency&, const double) const;
  bool FileExists(std::string) const;
  void AmmendSetOfTraces(const Transparency&);
  std::array<double, N> DarkSubtract(TTree*, TTree*) const;
  void MakeQuickPlot(const std::array<double, N>&, const std::array<double, N>&, const std::string&) const;
  void MakeDoublePlot(const std::array<double, N>&, const std::array<double, N>&, const std::array<double, N>&, const std::string&) const;
  void MakeHeatMap(const std::vector< std::array<double,N> >&, const std::array<double, N>&, const std::string&) const;
  bool HandleErrors();
  std::map<std::string, std::pair<double, double>> CreateSamplesMap(const Transparency&) const;
  
  bool get_ok;
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
 
};


#endif
