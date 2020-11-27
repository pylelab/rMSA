const char* docstring=""
"fastaNA rnacentral_species_specific_ids.fasta > rnacentral_species_specific_ids.fasta\n"
"    convert non-standard nucleotide sequence to -.*ATCGN\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

size_t fastaNA(const string infile="-", const string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,line;
    size_t nseqs=0;
    size_t i;
    char na;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        else if (line[0]=='>')
        {
            if (outfile!="-") fp_out<<line<<endl;
            else                cout<<line<<endl;
            nseqs++;
        }
        else
        {
            for (i=0;i<line.length();i++)
            {
                na=line[i];
                if ('a'<=na && na<='z') na-=32;
                if      (na=='I') na='A';
                else if (na=='U') na='T';
                if ('A'<=na && na<='Z' && (na!='A' && na!='T'
                                       &&  na!='C' && na!='G')) na='N';
                sequence+=na;
            }
            if (outfile!="-") fp_out<<sequence<<endl;
            else                cout<<sequence<<endl;
            sequence.clear();
        }
    }
    fp_in.close();
    fp_out.close();
    sequence.clear();
    line.clear();
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
    fastaNA(infile,outfile);
    return 0;
}
