#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

void readp(const string& filename,
                 vector<double>& timeNumbers,
                 vector<double>& tele,
                 vector<vector<double>>& pfromNS,
                 int freq,
                 vector<double>& timeSteps,
                 vector<double>& timeStepsF,
                 vector<vector<double>>& pNS,
                 vector<vector<double>>& pNSF
                 ) {
    // Open the input file
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cerr << "Error opening file " << filename << endl;
        return;
    }
    
    // Read the file line by line
    vector<double> realNumbers;
    string line;
    while (getline(inputFile, line)) {
        if (line[0] == '#') {
           
            continue;
        }
        
        //cout << nColumns <<endl;
        istringstream iss(line);
        double num;
        char c;
        iss >> num;
        timeNumbers.push_back(num);
        while (iss >> num) {
                realNumbers.push_back(num);      
        }
    }
    
    inputFile.close();

    // Compute the matrices ufromNS, vfromNS, and wfromNS
    //ROOF 
    //int nColumns = 191;
    //int nColumns = 655;
    int nColumns = 256;
    int nRows = timeNumbers.size();
    pfromNS.resize(nRows, vector<double>(nColumns));

    for (int k = 0; k < nRows; k++) {
        for (int i = 0; i < nColumns; i++) {
            pfromNS[k][i] = realNumbers[i+nColumns*k];
        }
    }
  
     int ii=0;
    
   for (int k = 0; k < nRows; k++) {
          if (abs(pfromNS[k][3])<509){
            timeSteps.push_back(timeNumbers[k]);
              ii=ii+1;
          }                
   }

   int nRows2 = ii;  
   pNS.resize(nRows2, vector<double>(nColumns));
   int iii=0;
 for (int k = 0; k < nRows; k++) {
          if (abs(pfromNS[k][3])<509){
            for (int j = 0; j < nColumns;j++){
               pNS[iii][j]=pfromNS[k][j];
               }
               iii=iii+1;
          }
 }
    int jj=0;
   
for (int k = 0; k < nRows2; k++) {
          if ((timeSteps[k]-timeSteps[0])>=jj*0.002){
                     tele.push_back(k);
                     jj=jj+1;
                     }
          }  

      int nRows3 = jj;
       double sum1, sum0;
      timeStepsF.resize(nRows3, 1);
       pNSF.resize(nRows3, vector<double>(nColumns));
       int metr=0;
       for (int j = 0; j < nColumns;j++){
      for (int l = 1; l < nRows3; l++) {
      for (int k =tele[l-1]; k < tele[l]; k++) {
               sum0+=timeSteps[k];
               sum1+=pNS[k][j];
               metr+=1;
      }
           timeStepsF[l-1]=sum0 / metr;
           pNSF[l-1][j]=sum1 / metr;
           sum1=0;
           sum0=0;
           metr=0;
     
      }
      }
}


int main() {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//------------Steps4. Extract pressure time series by filtering the openfoam data-------------------------//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

    int freq = 1000;
    
    vector<double> timeNumbers,tele,timeStep,timeStepF;
    vector<vector<double>> pfromNS,pNS,pNSF;
    
    readp("/home/potsis/projects/def-statho/potsis/2024/TPUNI_12o12_r01_400/postProcessingwithWW/probesTPU/0/p",timeNumbers,tele, pfromNS,freq, timeStep,timeStepF,pNS,pNSF);
    cout << "Results from file p" << endl;
    cout << "-------------------" << endl;
     cout << timeNumbers.size()<<" became "<< timeStep.size()<<endl;
    int ntaps = 256;
    ofstream out_file;
    out_file.open("/home/potsis/projects/def-statho/potsis/Pressures/p_FINAL");
    out_file.precision(7);
     for (int k = 0; k < timeStepF.size()-1; k++) {
             out_file << timeStepF[k] <<" ";
             for (int j = 0; j < ntaps; j++) {
             out_file << pNSF[k][j]<<" ";
        }
        out_file << endl;
     }
    
   }
