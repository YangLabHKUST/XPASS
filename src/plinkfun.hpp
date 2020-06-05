//
//  plinkfun.hpp
//  bivas_Xcode
//
//  Created by CAI Mingxuan on 2017/3/28.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#ifndef plinkfun_hpp
#define plinkfun_hpp

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <bitset>
using namespace std;


void readPlink(string stringname,int N, int P, int* X);
int getLineNum(string filename);

#endif /* plinkfun_hpp */
