//
//  readPlink.cpp
//  bivas_Xcode
//
//  Created by CAI Mingxuan on 2017/3/28.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#include <RcppArmadillo.h>
#include<iostream>
#include<string>
using namespace Rcpp;
using namespace std;
using namespace arma;
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "plinkfun.hpp"

void read_phenotypes(string filename, int N, double* y){
    std::ifstream ifs(filename.c_str());

    std::string line;

    string snpname;
    float pheno = 0;
    string tmp_str;
    int idx = 0;
    int phenopos = 5;
    while(std::getline(ifs, line)) // read one line from ifs
    {
        std::istringstream iss(line); // access line as a stream
        for(int i = 0; i < phenopos; i++){
            iss >> tmp_str;
        }
        iss >> pheno;
        y[idx] = pheno;
        idx ++;

    }
    ifs.close();
}

vector<string> read_snpnames(string filename, int P){
    vector<string> snpnames;
    std::ifstream ifs(filename.c_str());
    std::string line;
    int chromsome;
    string snpname;
    vector <string> fields;
    while(std::getline(ifs, line)) // read one line from ifs
    {
        std::istringstream iss(line); // access line as a stream
        boost::split( fields, line, boost::is_any_of(" \t *"));
        chromsome = (int)atoi(fields[0].c_str());
        snpname = fields[1];
        snpnames.push_back(snpname);
    }
    ifs.close();
    return snpnames;
}

// [[Rcpp::export]]
RcppExport SEXP read_data(Rcpp::String stringname, std::string fillMiss="none") {
    String famfile = stringname;
    famfile += ".fam";
    String bimfile = stringname;
    bimfile += ".bim";
    long long N =  getLineNum(famfile);
    long long P =  getLineNum(bimfile);
    int* X = new int[N*P];

    readPlink(stringname,N, P, X);
    arma::Mat<int> X1(X,N,P,false);
    if(fillMiss.compare("zero")==0)
      // X1.elem(find(X1==3)).zeros();
      X1.replace(3,0);

    double* y = new double[N];
    arma::Col<double> y1(y,N,false);
    read_phenotypes(famfile, N, y);
    vector<string> snps = read_snpnames(bimfile, P);

    List ret;
    ret["X"] = X1;
    ret["y"] = y1;
    ret["snps"] = snps;
    delete [] y;
    delete [] X;
    return ret;
}
