#include "rayreader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void usage(bool error=false, bool wait=false)
{
  cout << "usage:" << endl;
  cout << "raydiff file.ply ... " << endl;
  // ....
  if (wait)
  {
    cout << "<press ENTER>" << endl;
    getc(stdin);
  }
  exit(error);
}

int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    usage();
  }
  else
  {
  }
}
