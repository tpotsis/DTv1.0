#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;

namespace fs = std::filesystem;

struct Data {
    int ROW;
    int nocp;
    int tapU0;
    double infreq;
    double timeStep;
    double timeStart;
    double timeEnd;
    double quasiTS;
    std::string windTunnelDir;
    std::string writeturb;
    std::string caseDir;
    std::string includeV;
    std::string includeW;
    std::string udDir;
    std::string udrmsDir;
    std::string udrmsWDir;
    std::string taruDir;
    std::string tarurmsDir;
    std::string UNSDir;
    std::string UNSrmsDir;
    std::string modificationMean;
    std::string modificationRMS;
    std::string modificationRMSV;
    std::string modificationRMSW;
   int size;
   int size2;
   int nLines;
    std::vector<std::vector<double>> lmmu;
    std::vector<std::vector<double>> lmmurms;
    std::vector<std::vector<double>> lmmvrms;
    std::vector<std::vector<double>> lmmwrms;
};
void readU(const string& filename,
                int nColumns,
                 vector<double>& timeNumbers,
                 vector<vector<double>>& ufromNS,
                 vector<vector<double>>& vfromNS,
                 vector<vector<double>>& wfromNS,
                 vector<double>& avg,
                 vector<double>& rms) {
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
        istringstream iss(line);
        double num;
        char c;
        iss >> num;
        timeNumbers.push_back(num);
        while (iss >> c >> num) {
            if (c == '(') {
                realNumbers.push_back(num);
                iss >> num;
                realNumbers.push_back(num);
                iss >> num >> c;
                realNumbers.push_back(num);
            }
        }
    }
    inputFile.close();

    // Compute the matrices ufromNS, vfromNS, and wfromNS
   //int nColumns = 88;
    int nRows = timeNumbers.size();
    ufromNS.resize(nRows, vector<double>(nColumns));
    vfromNS.resize(nRows, vector<double>(nColumns));
    wfromNS.resize(nRows, vector<double>(nColumns));
    for (int k = 0; k < nRows; k++) {
        for (int i = 1; i < nColumns+1; i++) {
            ufromNS[k][i-1] = realNumbers[i*3-3+nColumns*3*k];
            vfromNS[k][i-1] = realNumbers[i*3-2+nColumns*3*k];
            wfromNS[k][i-1] = realNumbers[i*3-1+nColumns*3*k];
        }
    }

    // Compute the column averages and RMS values
    avg.resize(nColumns);
    rms.resize(nColumns);
    for (int i = 0; i < nColumns; i++) {
        double sum = 0.0;
        double sum_squares = 0.0;
        for (int k = 0; k < nRows; k++) {
            sum += ufromNS[k][i];
        }
        avg[i] = sum / nRows;
    }
    for (int i = 0; i < nColumns; i++) {
        double sum_squares = 0.0;
        for (int k = 0; k < nRows; k++) {
            sum_squares += pow(ufromNS[k][i]-avg[i], 2);
        } 
        rms[i] = sqrt(sum_squares / nRows);
    }

     
}

Data readInletData(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open file: " << filePath << '\n';
        exit(EXIT_FAILURE);
    }

    std::string line;
    Data data = {};
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string keyword;
        if (iss >> keyword) {
            if (keyword == "inlet_frequency") {
                iss >> data.infreq;
                data.timeStep=1/data.infreq;
                data.quasiTS=data.timeStep-0.00000001;
            }
            else if (keyword == "time_start") {
                iss >> data.timeStart;
            }
            else if (keyword == "time_end") {
                iss >> data.timeEnd;
                data.ROW=data.timeEnd/data.timeStep;
            }

            else if (keyword == "number_of_control_points") {
                iss >> data.nocp;
                
            }
           
           else if (keyword == "write") {
                iss >> data.writeturb;
                
            }
            else if (keyword == "wind_tunnel") {
                iss.ignore();
                iss >> data.windTunnelDir;
            }
            else if (keyword == "case" ) {
                iss.ignore();
                iss >> data.caseDir;
            }
            else if (keyword == "include_v") {
                iss.ignore();
                iss >> data.includeV;
                            }
            else if (keyword == "include_w") {
                iss.ignore();
                iss >> data.includeW;
            }
             else if (keyword == "inlet_nLines") {
                iss.ignore();
                iss >> data.nLines;
            }

            else if (keyword == "modification_mean") {
                iss.ignore();
                    iss >> data.modificationMean;
            }

            else if (keyword == "user_defined_directory_mean") {
                iss.ignore();
                    iss >> data.udDir;
            } 

            else if (keyword == "targetProfileU_directory") {
                iss.ignore();
                    iss >> data.taruDir;
            } 
  
            else if (keyword == "results_U0_directory") {
                iss.ignore();
                    iss >> data.UNSDir;
            } 

            else if (keyword == "number_of_probes_in_U0") {
                iss.ignore();
                    iss >> data.tapU0;
                    if(data.tapU0!=22){
                   cout<<"Please check Incident Profile Data in " <<data.UNSDir<< endl;
                    }
            } 
            
             else if (keyword == "modification_rms") {
                 iss.ignore();
                iss >> data.modificationRMS;
             }
             else if (keyword == "modification_rmsV") {
                 iss.ignore();
                iss >> data.modificationRMSV;
             }
             else if (keyword == "modification_rmsW") {
                 iss.ignore();
                iss >> data.modificationRMSW;
             }

             else if (keyword == "user_defined_directory_rmsW") {
                 iss.ignore();
                iss >> data.udrmsWDir;
             }

             else if (keyword == "user_defined_directory_rms") {
                 iss.ignore();
                iss >> data.udrmsDir;
             }
    
            else if (keyword == "targetProfileUrms_directory") {
                 iss.ignore();
                iss >> data.tarurmsDir;
             }
             else if (keyword == "results_U0rms_directory") {
                 iss.ignore();
                iss >> data.UNSrmsDir;
             }
  }
         }
       
             if (data.modificationMean == "none") {
                data.lmmu.resize(1, std::vector<double>(data.nocp)); // create a 1xsize matrix
                for (int i = 0; i < data.nocp; i++) {
                data.lmmu[0][i] = 1; // read the values into a vector
                                                      }
               
                                                    }

                 

             if (data.modificationMean == "user_defined") {
                
                std::ifstream infile(data.udDir);
                infile >> data.size; // read the first number, which is the size of the array
                data.lmmu.resize(1, std::vector<double>(data.size));
                for (int i = 0; i < data.size; i++) {
                 infile >> data.lmmu[0][i]; // read the values into a vector
                 }
                 infile.close();
                                                           }
               
             if (data.modificationMean == "linear") {
                 vector<double> timeNumbers;
                vector<vector<double>> ufromNS, vfromNS, wfromNS;
                vector<double> avgNS, rmsNS;
                readU(data.UNSDir, data.tapU0,timeNumbers, ufromNS, vfromNS, wfromNS, avgNS, rmsNS);

                std::ifstream infile(data.taruDir);
                infile >> data.size; // read the first number, which is the size of the array
                data.lmmu.resize(1, std::vector<double>(data.size));
                for (int i = 0; i < data.size; i++) {
                 infile >> data.lmmu[0][i]; // read the values into a vector
                 data.lmmu[0][i] /=avgNS[i];
                 }
                 infile.close();
                                                         }
                 
                if (data.modificationRMS == "none") {
                data.lmmurms.resize(1, std::vector<double>(data.nocp)); // create a 1xsize matrix
                for (int i = 0; i < data.nocp; i++) {
                data.lmmurms[0][i] = 1; // read the values into a vector
                                                      }
                                                     }                                            
    
                if (data.modificationRMS == "user_defined") {
                
                std::ifstream infile(data.udrmsDir);
                infile >> data.size; // read the first number, which is the size of the array
                data.lmmurms.resize(1, std::vector<double>(data.size));
                for (int i = 0; i < data.size; i++) {
                 infile >> data.lmmurms[0][i]; // read the values into a vector
                 }
                 infile.close();
                                                           }
               if (data.modificationRMS == "linear") {
                 vector<double> timeNumbers;
                vector<vector<double>> ufromNS, vfromNS, wfromNS;
                vector<double> avgNS, rmsNS;
                readU(data.UNSrmsDir,data.tapU0,timeNumbers, ufromNS, vfromNS, wfromNS, avgNS, rmsNS);

                std::ifstream infile(data.tarurmsDir);
                infile >> data.size; // read the first number, which is the size of the array
                data.lmmurms.resize(1, std::vector<double>(data.size));
                for (int i = 0; i < data.size; i++) {
                 infile >> data.lmmurms[0][i]; // read the values into a vector
                 data.lmmurms[0][i] /=rmsNS[i];
                 }
                 infile.close();
                                                         }

                 if (data.modificationRMS == "nonlinear") {
                 vector<double> timeNumbers;
                vector<vector<double>> ufromNS, vfromNS, wfromNS;
                vector<double> avgNS, rmsNS;
                readU(data.UNSrmsDir, data.tapU0, timeNumbers, ufromNS, vfromNS, wfromNS, avgNS, rmsNS);

                std::ifstream infile(data.tarurmsDir);
                infile >> data.size; // read the first number, which is the size of the array
                data.lmmurms.resize(1, std::vector<double>(data.size));
                for (int i = 0; i < data.size; i++) {
                 infile >> data.lmmurms[0][i]; // read the values into a vector
                 data.lmmurms[0][i] =pow(data.lmmurms[0][i]/rmsNS[i],2);
                 }
                 infile.close();
                                                         }  

                              
                          if (data.modificationRMSV == "yes") {
                             data.lmmvrms.resize(1, std::vector<double>(data.nocp)); // create a 1xsize matrix
                for (int i = 0; i < data.nocp; i++) {
                 data.lmmvrms[0][i]=data.lmmurms[0][i]; // read the values into a vector
                 }
                          }   
                          if (data.modificationRMSV == "no") {
                             data.lmmvrms.resize(1, std::vector<double>(data.nocp)); // create a 1xsize matrix
                for (int i = 0; i < data.nocp; i++) {
                 data.lmmvrms[0][i]=1; // read the values into a vector
                 }
                          }   

                          if (data.modificationRMSW == "yes") {
                             data.lmmwrms.resize(1, std::vector<double>(data.nocp)); // create a 1xsize matrix
                for (int i = 0; i < data.nocp; i++) {
                 data.lmmwrms[0][i]=data.lmmurms[0][i]; // read the values into a vector
                 }
                          }     
                          if (data.modificationRMSW == "no") {
                             data.lmmwrms.resize(1, std::vector<double>(data.nocp)); // create a 1xsize matrix
                for (int i = 0; i < data.nocp; i++) {
                 data.lmmwrms[0][i]=1; // read the values into a vector
                 }
                          }         

                          if (data.modificationRMSW == "user_defined") {
                std::ifstream infile(data.udrmsWDir);
                infile >> data.size; // read the first number, which is the size of the array
                data.lmmwrms.resize(1, std::vector<double>(data.size));
                for (int i = 0; i < data.size; i++) {
                 infile >> data.lmmwrms[0][i]; // read the values into a vector
                 }
                 infile.close();

                
                          }                                                 
                                                  
    file.close();
    return data;
}

void writeHeader() {
    cout << "/*---------------------------------------------------------------------------*\\" << endl;
    cout << "|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|" << endl;
    cout << "| ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ |" << endl;
    cout << "| ~                             Dynamic Terrain                             ~ |" << endl;
    cout << "| ~     Modelling wind flow and wind-induced structural loads in the ABL    ~ |" << endl;
    cout << "| ~                                  v1.0                                   ~ |" << endl;
    cout << "| ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ |" << endl;
    cout << "|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|" << endl;
    cout << "\\*---------------------------------------------------------------------------*/" << endl;
    cout << endl << endl;
}

void readDataFromFile(string filename, int ROWS, int COLS, vector<vector<double>>& matrix1, vector<vector<double>>& matrix2, vector<vector<double>>& matrix3) {
    int row = 0; // current row being read

    ifstream infile(filename);

    if (!infile) { // check if file was opened successfully
        cerr << "Error opening file " << filename << endl;
        return;
    }

    // read data from file
    while (!infile.eof() && row < ROWS) {
        double num1, num2, num3;
        infile >> num1 >> num2 >> num3;

        // store data in matrices
        matrix1[row][0] = num1;
        matrix2[row][0] = num2;
        matrix3[row][0] = num3;

        row++;
    }

    infile.close(); // close file
}

// Function to check if a string contains only digits
bool is_digits(const std::string& str) {
    return std::all_of(str.begin(), str.end(), ::isdigit);
}

// Function to sort filenames in ascending order based on the numerical part of their name
bool sort_filenames(const std::string& a, const std::string& b) {
    std::string a_num = a.substr(0, a.find(".txt"));
    std::string b_num = b.substr(0, b.find(".txt"));
    return std::stoi(a_num) < std::stoi(b_num);
}


void createFile(string dir_path, vector<vector<double>> matrix_values, int nLines) {
    // Create the directory if it does not exist
    //mkdir(dir_path.c_str());
    fs::create_directory(dir_path.c_str());
    //int nLines;
    // Open the file for writing
    ofstream out_file;
    out_file.open(dir_path + "/U");
    out_file << "(" << endl;
    // Write the matrix values to the file 11 times, one under the other
    out_file << fixed;
    out_file.precision(4);
    for (int k = 0; k < nLines; k++) {
        for (int i = 0; i < matrix_values[0].size(); i++) {
            out_file << "(";
            for (int j = 0; j < matrix_values.size(); j++) {
                out_file << matrix_values[j][i];
                if (j != matrix_values.size() - 1) {
                    out_file << " ";
                }
            }
            out_file << ")" << endl;
        }
    }

    // t lose the file
    out_file << ")" << endl;
    out_file.close();

}




int main() {
    /////////////////////////////////////////////////////////////////////////////////////
     //-------------------------Path of Data file-------------------------------------//
     ////////////////////////////////////////////////////////////////////////////////////

     std::string DataPath = "/home/potsis/projects/def-statho/potsis/GitHub/Data";

    //////////////////////////////////////////////////////////////////////////////////
     //--------------------------Read Data file-------------------------------------//
     /////////////////////////////////////////////////////////////////////////////////
    Data data = readInletData(DataPath);
    int ROWS = data.ROW;
    int COLS = 3;
    writeHeader();
    std::string path = data.windTunnelDir;
    std::string path1 = data.caseDir;
    //////////////////////////////////////////////////////////////////////////////////////////
    //--------------------------Step 1. Wind Tunnel Data-----------------------------------//
    /////////////////////////////////////////////////////////////////////////////////////////
    
    std::vector<std::string> filenames;
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        
        if (entry.is_regular_file() && entry.path().extension() == ".txt" && is_digits(entry.path().stem().string())) {
        filenames.push_back(entry.path().filename().string());
        }
    }
    // Sort filenames in ascending order based on the numerical part of their name
    std::sort(filenames.begin(), filenames.end(), sort_filenames);
    // Store filenames in increasing order matrix named heights
    int n_files = filenames.size();
    double nheights[n_files][1];
    std::string heights[n_files][1];
    for (int i = 0; i < n_files; i++) {
    heights[i][0] = filenames[i];
    nheights[i][0]=stod(filenames[i]);
    }
    
    std::cout << "Reading Wind Tunnel Data from " <<  n_files << " Locations" <<std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
    /*
    for (int i = 0; i < n_files; i++) {
    std::cout << heights[i][0] << std::endl;
    std::cout << nheights[i][0]<< " mm"<< std::endl;
    }
    */
    double udata[ROWS][n_files];
    double vdata[ROWS][n_files];
    double wdata[ROWS][n_files];
    for (int i = 0; i < n_files; i++) {
    string filename1 = path + "/" + heights[i][0];
    vector<vector<double>> matrix1(ROWS, vector<double>(COLS));
    vector<vector<double>> matrix2(ROWS, vector<double>(COLS));
    vector<vector<double>> matrix3(ROWS, vector<double>(COLS));
    readDataFromFile(filename1,ROWS, COLS,matrix1, matrix2, matrix3);
   
   for (int k = 0; k < ROWS; k++) {
            udata[k][i]=matrix1[k][0];
            vdata[k][i]=matrix2[k][0];
            wdata[k][i]=matrix3[k][0];
    }

    }
      if (data.includeV == "no"){
   for (int i = 0; i < n_files; i++) {
   for (int k = 0; k < ROWS; k++) {
            vdata[k][i]=0;
    }

      }
      }
if (data.includeW == "no"){
   for (int i = 0; i < n_files; i++) {
   for (int k = 0; k < ROWS; k++) {
            wdata[k][i]=0;
    }

      }
      }

    vector<double> avgu, avgv, avgw, rmsu;   
    avgu.resize(n_files);
    avgv.resize(n_files);
    avgw.resize(n_files);
    rmsu.resize(n_files);

    for (int i = 0; i < n_files; i++) {
     
        double sumu = 0.0;
        double sumv = 0.0;
        double sumw = 0.0;
      for (int f = 0; f < ROWS; f++) {
            sumu += udata[f][i];
            sumv += vdata[f][i];
            sumw += wdata[f][i];
               }
      avgu[i] = sumu / ROWS;
      avgv[i] = sumv / ROWS;
      avgw[i] = sumw / ROWS;
     }

     for (int i = 0; i < n_files; i++) {
             double sum_squares = 0.0;
      for (int f = 0; f < ROWS; f++) {

             sum_squares += pow(udata[f][i]-avgu[i], 2);
    }
     rmsu[i] = sqrt(sum_squares/ ROWS);
     }
cout << endl << endl;
    cout << "Statistics of Wind Tunnel Data " << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Z (m) Umean (m/s) Urms (m/s) Turbulence Intensity (%):" << endl;
     for (int k = 0; k < avgu.size(); k++) {
             cout <<nheights[k][0]/1000<<"  "<< avgu[k]<<" "<< rmsu[k]<< " "<<(rmsu[k]/avgu[k])*100<< " %"<< endl;
        } 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//------------Steps 2 & 3. Creating Turbuelence based on Inlet frequency and mean and std modifications----------------------//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int f = 0; f < ROWS; f++) {
      vector<vector<double>> uinlet(3, vector<double>(n_files));
 
     // double udataf[ROWS][n_files];
      for (int i = 0; i < n_files; i++) {
            
            uinlet[0][i]=(udata[f][i]-avgu[i])*data.lmmurms[0][i]+avgu[i]*data.lmmu[0][i];
            uinlet[1][i]=(vdata[f][i]-avgv[i])*data.lmmvrms[0][i]+avgv[i];
            uinlet[2][i]=(wdata[f][i]-avgw[i])*data.lmmwrms[0][i]+avgw[i];
    }
     
      double time=f*data.timeStep;
      double timeq=f*data.timeStep+data.quasiTS;

      std::stringstream ss;
     ss << std::fixed << std::setprecision(8) << time; // set precision to 8 decimal places
       std::string str = ss.str();

      std::stringstream ssq;
      ssq << std::fixed << std::setprecision(8) << timeq; // set precision to 8 decimal places
      std::string strq = ssq.str();

  
    std::string dir_path = path1 + "/"+ str; 
    std::string dir_pathq = path1 + "/"+ strq; 
    // print data from matrices
     if (data.writeturb == "yes"){
    createFile(dir_path, uinlet,data.nLines);
    createFile(dir_pathq, uinlet,data.nLines);
     }
    }
     vector<double> avguf, avgvf, avgwf, rmsuf,rmsvf,rmswf;   
    avguf.resize(n_files);
    avgvf.resize(n_files);
    avgwf.resize(n_files);
    rmsuf.resize(n_files);
    rmsvf.resize(n_files);
    rmswf.resize(n_files);

    for (int i = 0; i < n_files; i++) {
     
        double sumuf = 0.0;
        double sumvf = 0.0;
        double sumwf = 0.0;
      for (int f = 0; f < ROWS; f++) {
            sumuf += (udata[f][i]-avgu[i])*data.lmmurms[0][i]+avgu[i]*data.lmmu[0][i];
            sumvf += (vdata[f][i]-avgv[i])*data.lmmvrms[0][i]+avgv[i];
            sumwf += (wdata[f][i]-avgw[i])*data.lmmwrms[0][i]+avgw[i];
               }
      avguf[i] = sumuf / ROWS;
      avgvf[i] = sumvf / ROWS;
      avgwf[i] = sumwf / ROWS;
     }
     for (int i = 0; i < n_files; i++) {
             double sum_squaresf = 0.0;
             double sum_squaresfv = 0.0;
             double sum_squaresfw = 0.0;
      for (int f = 0; f < ROWS; f++) {

             sum_squaresf += pow((udata[f][i]-avgu[i])*data.lmmurms[0][i]+avgu[i]*data.lmmu[0][i]-avguf[i], 2);
             sum_squaresfv += pow((vdata[f][i]-avgv[i])*data.lmmvrms[0][i]+avgv[i], 2);
             sum_squaresfw += pow((wdata[f][i]-avgw[i])*data.lmmwrms[0][i]+avgw[i], 2);
    }
     rmsuf[i] = sqrt(sum_squaresf/ ROWS);
     rmsvf[i] = sqrt(sum_squaresfv/ ROWS);
     rmswf[i] = sqrt(sum_squaresfw/ ROWS);
     }
    cout << endl << endl;
    cout << "Statistics of Imposed Turbulence at Inlet" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout <<"Z(m) MeanU(m/s) MeanV(m/s) MeanW(m/s) rmsU(m/s) rmsV(m/s) rmsW(m/s) IntensityU(%) IntensityV(%) IntensityW(%)" << endl;
     for (int k = 0; k < avgu.size(); k++) {
        
             cout <<nheights[k][0]/1000<<"  "<< avguf[k]<<"  "<< avgvf[k]<<"  "<< avgwf[k]<<" "<< rmsuf[k]<< " "<< rmsvf[k]<< " "<< rmswf[k]<< " "<<(rmsuf[k]/avguf[k])*100<< " % "<<(rmsvf[k]/avguf[k])*100<< " % "<<(rmswf[k]/avguf[k])*100<< " % "<< endl;
        } 

    cout <<data.writeturb<< endl;

    return 0;
}
