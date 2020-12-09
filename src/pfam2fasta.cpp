const char* docstring=""
"pfam2fasta metaclust.pfam metaclust.fasta\n"
"    convert tab-eliminated table metaclust.pfam\n"
"    to FASTA format alignment metaclust.fasta\n"
"\n"
"pfam2fasta metaclust.pfam metaclust.fasta 60\n"
"    convert tab-eliminated table metaclust.pfam\n"
"    to FASTA format alignment metaclust.fasta\n"
"    where the sequence is wrapped into lines with up to 60 characters\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

int pfam2fasta(const string infile="-", const string outfile="-",
    const int maxLineAAnum=0)
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string txt,line;
    int nseqs=0;
    int i,j;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        for (i=0;i<line.size();i++)
        {
            if (line[i]=='\t')
            {
                nseqs++;
                txt='>'+line.substr(0,i)+'\n';
                line=line.substr(i+1);
                if (maxLineAAnum==0) txt+=line;
                else
                {
                    for (j=0;j<line.size();j+=60)
                    {
                        if (j+60>line.size())
                            txt+=line.substr(j,line.size()-j)+'\n';
                        else
                            txt+=line.substr(j,60)+'\n';
                    }
                }
                if (outfile!="-") fp_out<<txt.substr(0,txt.size()-1)<<endl;
                else                cout<<txt.substr(0,txt.size()-1)<<endl;
                txt.clear();
                break;
            }
        }
    }
    fp_in.close();
    fp_out.close();
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
    int maxLineAAnum=(argc<=3)?0:atoi(argv[3]);
    pfam2fasta(infile,outfile,maxLineAAnum);
    return 0;
}
