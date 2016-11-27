// ex1.cpp : Example to writing to a binary file
//

#include <iostream>
#include <string>

using namespace std;

enum writedir { NOTOPEN, READ, WRITE, READWRITE };

class RandomAccess {
private:
  FILE * file;
  string filename;
  writedir dir;
  
public:
  RandomAccess(const string FileName) : filename(FileName),dir(NOTOPEN) {};
  bool OpenWrite();
  bool FileSuccess(string reason);
  bool Write(string text);
  bool Close();
};

bool RandomAccess::OpenWrite() {
  if (dir != NOTOPEN )
    {
      file =0;
      return(FileSuccess("File already Opened"));
    }

  file = fopen(filename.c_str(), "wb");
  dir = WRITE;
  return (FileSuccess("Opening File: "));
}

bool RandomAccess::FileSuccess(string reason) {
  reason += filename + " Result :" +  ( file ? "Succeeded" : "Failed" );
  //OutputDebugString( reason.c_str() );
  return (file != 0 );
}

bool RandomAccess::Close() {
  int errornum = fclose( file );
  return (errornum == 0);
}


bool RandomAccess::Write(string text) {
  size_t numwritten= fwrite(text.c_str(),sizeof(char),text.size(), file);
  return (numwritten >0); 
}


int main(int argc, char * argv[])
{
  const string filename="test.txt";
  const string mytext="Once upon a time there were three bears.";

  RandomAccess ra( filename);
  if ( ra.OpenWrite() ) 
    {
      if (!ra.Write( mytext )) 
	cout << "Failed to write to file " << filename << endl;
      ra.Close();
    }
  else
    cout << "Failed to open " << filename << " fior writing " << endl;

  return 0;
}
