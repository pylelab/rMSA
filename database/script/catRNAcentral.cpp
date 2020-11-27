const char* docstring=""
"catRNAcentral input.fasta output.fasta output.tsv\n"
"    Combine different RNAcentral entries from different specieis in\n"
"    input.fasta for the same sequence into a single entry in output.fasta.\n"
"    Also output the table for entry name vs species to output.tsv.\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

void catRNAcentral(const string infile="-", 
    const string outfasta="-",const string outtsv="-")
{
    ifstream fp_in;
    ofstream fp_fasta;
    if (infile!="-")   fp_in.open(infile.c_str(),ios::in);
    if (outfasta!="-") fp_fasta.open(outfasta.c_str(),ofstream::out);
    string name,line;
    size_t i;
    bool readseq;
    vector<string> name_vec;
    map<string,string> species_map;
    map<string,size_t> len_map;

    /* write fasta */
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            for (i=1;i<line.size();i++)
                if (line[i]=='_') break;
            name=line.substr(1,i-1);
            if (species_map.count(name))
            {
                readseq=false;
                species_map[name]+=","+line.substr(i+1);
                continue;
            }
            readseq=true;
            species_map[name]=line.substr(i+1);
            name_vec.push_back(name);
            len_map[name]=0;
            if (outfasta=="-") cout<<line.substr(0,i)<<endl;
            else           fp_fasta<<line.substr(0,i)<<endl;
        }
        else if (readseq)
        {
            if (outfasta=="-") cout<<line<<endl;
            else           fp_fasta<<line<<endl;
            len_map[name]+=line.size();
            line.clear();
        }
    }
    fp_in.close();
    fp_fasta.close();
    line.clear();
    
    /* write species tsv */
    ofstream fp_tsv;
    if (outtsv!="-") fp_tsv.open(outtsv.c_str(),ofstream::out);
    size_t L;
    for (i=0;i<name_vec.size();i++)
    {
        name=name_vec[i];
        L=len_map[name];
        if (outtsv=="-") cout<<name<<'\t'<<L<<'\t'<<species_map[name]<<endl;
        else           fp_tsv<<name<<'\t'<<L<<'\t'<<species_map[name]<<endl;
    }
    fp_tsv.close();
    name.clear();
    vector<string>().swap(name_vec);
    map<string,string>().swap(species_map);
    return;
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
    string outfasta=(argc<=2)?"-":argv[2];
    string outtsv  =(argc<=3)?"-":argv[3];
    catRNAcentral(infile,outfasta,outtsv);
    return 0;
}
