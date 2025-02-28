#include "BenAnalysis.h"

/*double RosenBrock(const double *xx )
{
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
  }*/

BenAnalysis::BenAnalysis():Tool(){}


bool BenAnalysis::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;

  return true;
}


bool BenAnalysis::Execute(){
  
  if(m_data->state==Analyse){
    
    std::vector<double> value;
    std::vector<double> error;
    std::vector<double> wavelength;
    std::vector<double> darkvalue;
    std::vector<double> darkerror;
    std::vector<double> darkwavelength;
    
    short year;
    short month;
    short day;
    short hour;
    short min;
    
    std::vector<double>* a=&value;
    std::vector<double>* b=&error;
    std::vector<double>* c=&wavelength;
    std::vector<double>* darka=&darkvalue;
    std::vector<double>* darkb=&darkerror;
    std::vector<double>* darkc=&darkwavelength;
    
    
    
    for (std::map<std::string, TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); it++){
      
      std::cout<<"a1 "<<it->first<<std::endl;
      if(it->first.at(0)=='D') continue;
std::cout<<"a2 "<<std::endl;

      std::string dark="D"+(it->first);
std::cout<<"a3 "<<std::endl;
      
      m_data->m_trees[dark]->SetBranchAddress("value", &darka);
      m_data->m_trees[dark]->SetBranchAddress("error", &darkb);
      m_data->m_trees[dark]->SetBranchAddress("wavelength", &darkc);
      m_data->m_trees[dark]->SetBranchAddress("year", &year);
      m_data->m_trees[dark]->SetBranchAddress("month", &month);
      m_data->m_trees[dark]->SetBranchAddress("day", &day);
      m_data->m_trees[dark]->SetBranchAddress("hour", &hour);
      m_data->m_trees[dark]->SetBranchAddress("min", &min);

      std::cout<<"a4 "<<std::endl;
      m_data->m_trees[dark]->GetEntry(m_data->m_trees[dark]->GetEntriesFast()-1);
      std::cout<<"a5 "<<std::endl;
      
      it->second->SetBranchAddress("value", &a);
      it->second->SetBranchAddress("error", &b);
      it->second->SetBranchAddress("wavelength", &c);
      it->second->SetBranchAddress("year", &year);
      it->second->SetBranchAddress("month", &month);
      it->second->SetBranchAddress("day", &day);
      it->second->SetBranchAddress("hour", &hour);
      it->second->SetBranchAddress("min", &min);

      std::cout<<"a6 "<<std::endl;
      it->second->GetEntry(it->second->GetEntriesFast()-1);
      std::cout<<"a7 "<<std::endl;
      
      //TCanvas c1= new TCanvas();//"c1","A Simple Graph with error bars");//,200,10,700,500);
      TCanvas c1("c1","A Simple Graph with error bars",200,10,700,500);
      //c1.SetFillColor(42);
      //c1.SetGrid();
      //    c1.GetFrame()->SetFillColor(21);
      //c1.GetFrame()->SetBorderSize(12);
      //   const Int_t n = 10;
      // Double_t x[n]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
      //Double_t y[n]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
      // Double_t ex[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
      //Double_t ey[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
      
      /*     /////////////////////////////////////////////////// minimisation ///////////////////////
	     
	     ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
	     
	     minimum->SetMaxIterations(10000);  // for GSL
	     minimum->SetTolerance(0.001);
	     minimum->SetPrintLevel(1);
	     
	     ROOT::Math::Functor f(&RosenBrock,2);
	     double step[2] = {0.01,0.01};
	     
	     minimum->SetFunction(f);
	     
	     minimum->SetVariable(0,"x",variable[0], step[0]);
	     minimum->SetVariable(1,"y",variable[1], step[1]);
	     
	     minimum->Minimize();
	     
	     const double *xs = minimum->X();
	     std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
	     << minimum->MinValue()  << std::endl;
	     // expected minimum is 0
	     if ( minimum->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
	     std::cout << "Minimizer " << minName << " - " << algoName
	     << "   converged to the right minimum" << std::endl;
	     else {
	     std::cout << "Minimizer " << minName << " - " << algoName
	     << "   failed to converge !!!" << std::endl;
	     Error("NumericalMinimization","fail to converge");
	     }
	     
      //////////////////////////////*/
      /*      double sum=0;
      double darksum=0;
      double scale=0;
      
      for(int i=0; i<value.size(); i++){
	
	if(wavelength.at(i)>=330.0 && wavelength.at(i)<=355.0){
	  
	  darksum+=darkvalue.at(i);
	  sum+=value.at(i);
	}
      }

      
      scale=sum/darksum;
      
      std::cout<<"sum="<<sum<<" : "<<"darksum="<<darksum<<" : "<<"scale="<<scale<<std::endl;
      */
      std::cout<<"a8 "<<std::endl;
      for(int i=0; i<value.size(); i++){
     
	value.at(i)=value.at(i)-(darkvalue.at(i));//*scale);
	error.at(i)=sqrt((error.at(i)*error.at(i))+(darkerror.at(i)*darkerror.at(i)));
	
      }
    
      std::string newname="D"+(it->first)+"Corrected";
      TTree* tree=new TTree(newname.c_str(),newname.c_str());
      
      
      
      tree->Branch("value", &value);
      tree->Branch("error", &error);
      tree->Branch("wavelength", &wavelength);
      tree->Branch("year", &year);
      tree->Branch("month", &month);
      tree->Branch("day", &day);
      tree->Branch("hour", &hour);
      tree->Branch("min", &min);
      
      tree->Fill();

      m_data->AddTTree(newname,tree);
      
      std::cout<<"a9 "<<std::endl;
      
      std::cout<<"ggggggggggggg "<<value.size()<<":"<<value.at(0)<<std::endl;
      TGraphErrors gr(value.size(),&wavelength[0],&value[0],0,&error[0]);
      gr.SetTitle("TGraphErrors Example");
      //gr.SetMarkerColor(4);
      //gr.SetMarkerStyle(21);
      gr.Draw("AP");
      std::stringstream tmp;
      tmp<<it->first<<".png";
      c1.SaveAs(tmp.str().c_str());
      gr.Write(tmp.str().c_str(),TObject::kOverwrite);
      std::cout<<"a10 "<<std::endl;
    }
    std::cout<<"a11 "<<std::endl;
  }
  
  return true;
}


bool BenAnalysis::Finalise(){
  
  return true;
}
