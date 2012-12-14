#include <iostream>
#include <string>
#include <vector>
#include <sstream>

//using namespace std;

#include "bkgHistos.C"

//========================================================
// code to read input arguments 
template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
  const char delim = '=';
  // Skip if key or delimiter not there                                                                                                                                            
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;

  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key                                                                                                                       
  std::getline(sstr, sval, delim); // Get value                                                                                                                                    
  T tout;
  if (!(std::istringstream(sval) >> std::boolalpha >> tout))
    return;
  out = tout;
  std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}

// Special version for string without the conversion                                                                                                                                 
template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
  const char delim = '=';
  // Skip if key or delimiter not there                                                                                                                                            
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;
  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key                                                                                                                       
  std::getline(sstr, sval, delim); // Get value                                                                                                                                    
  out = sval;
  std::cout << skey << ": "  << out << std::endl;
}


//===================================================
int main(int argc, char* argv[]){

  // Set defaults
    int
        rapMin = 999,
        rapMax = 999,
        ptMin = 999,
        ptMax = 999,
        nState = 999;
    bool
        MC = false,
        f_BG_zero = false,
			 	doCtauUncer = false;

    // Loop over argument list
    for (int i=1; i < argc; i++)
      {
	std::string arg = argv[i];
        fromSplit("rapMin", arg, rapMin);
        fromSplit("rapMax", arg, rapMax);
        fromSplit("ptMin", arg, ptMin);
        fromSplit("ptMax", arg, ptMax);
        fromSplit("nState", arg, nState);
        fromSplit("MC", arg, MC);
        fromSplit("f_BG_zero", arg, f_BG_zero);      
        fromSplit("doCtauUncer", arg, doCtauUncer);      
      }

    std::cout << "-----------------------\n"
	      << "Creating background model for \n"
	      << "y bins " << rapMin << " - " << rapMax << "\n"
	      << "and pT bins "  << ptMin << " - " << ptMax << "\n"
	      << "-----------------------" << std::endl;

    for(int iRap = rapMin; iRap <= rapMax; iRap++){
      for(int iPT = ptMin; iPT <= ptMax; iPT++){

	std::stringstream temp;
	temp << "tmpFiles/fit_Psi" << nState-3 << "S_rap" << iRap << "_pt" << iPT << ".root";
	const std::string infilename = temp.str().c_str();
        
	bkgHistos(infilename.c_str(), iRap, iPT, nState, MC, f_BG_zero, doCtauUncer);

      }
    }

    return 0;
}