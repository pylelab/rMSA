const char* docstring=""
"fasta2pfam metaclust.fasta metaclust.pfam\n"
"    convert FASTA format alignment metaclust.fasta\n"
"    to tab-eliminated table metaclust.pfam\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

int fasta2pfam(const string infile="-", const string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,line,header;
    int nseqs=0;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (sequence.length()>0)
            {
                if (outfile!="-") fp_out<<header<<'\t'<<sequence<<endl;
                else                cout<<header<<'\t'<<sequence<<endl;
            }
            sequence.clear();
            header=line.substr(1,line.size()-1);
            nseqs++;
        }
        else
            sequence+=line;
    }
    fp_in.close();
    if (outfile!="-") fp_out<<header<<'\t'<<sequence<<endl;
    else                cout<<header<<'\t'<<sequence<<endl;
    fp_out.close();
    sequence.clear();
    header.clear();
    return nseqs;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string infile=argv[1];
    string outfile=(argc<=2)?"-":argv[2];
    fasta2pfam(infile,outfile);
    return 0;
}
