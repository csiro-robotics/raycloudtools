// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloud.h"
#include "rayalignment.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <complex>
using namespace std;
using namespace Eigen;
using namespace RAY;

void usage(bool error=false)
{
  cout << "Align raycloudA onto raycloudB, rigidly. Outputs the transformed version of raycloudA." << endl;
  cout << "usage:" << endl;
  cout << "rayalign raycloudA raycloudB." << endl;
  exit(error);
}

int main(int argc, char *argv[])
{
  // TODO: This method works when there is more than 30% overlap. 
  // For less, we can additionally repeat this procedure for small scale cubic sections (perhaps 20% of map width)
  // then use the fft power spectrum (low res) as a descriptor into the knn to find matching locations, 
  // then pick the best match.
  
  if (argc != 3)
    usage();

  string fileA = argv[1];
  string fileB = argv[2];

  AlignTranslationYaw aligner;
  aligner.clouds[0].load(fileA);
  aligner.clouds[1].load(fileB);
  
  aligner.alignCloud0ToCloud1(0.5);

  string fileStub = fileA;
  if (fileStub.substr(fileStub.length()-4)==".ply")
    fileStub = fileStub.substr(0,fileStub.length()-4);
  aligner.clouds[0].save(fileStub + "_aligned.ply");  

  return true;
}
