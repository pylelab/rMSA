const char* docstring=""
"RemoveNonQueryPosition clustalo.fasta > clustalo.noquerygap.fasta\n"
"    delete any position corresponding to gap in query\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

int RemoveNonQueryPosition(const string infile="-", const string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,line;
    int nseqs=0;

    int i;
    vector <int> nongap_pos; // position not corresponding to gap in query
    string no_query_gap_sequence;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (sequence.length()>0)
            {
                if (nongap_pos.size()==0)
                    for (i=0;i<sequence.size();i++)
                        if (sequence[i]!='-') 
                            nongap_pos.push_back(i);

                no_query_gap_sequence.clear();
                for (i=0;i<nongap_pos.size();i++)
                    no_query_gap_sequence+=sequence[nongap_pos[i]];

                if (outfile!="-") fp_out<<no_query_gap_sequence<<endl;
                else                cout<<no_query_gap_sequence<<endl;
            }
            sequence.clear();
            nseqs++;
            if (outfile!="-") fp_out<<line<<endl;
            else                cout<<line<<endl;
        }
        else
            sequence+=line;
    }
    fp_in.close();

    if (nongap_pos.size()==0)
        for (i=0;i<sequence.size();i++)
            if (sequence[i]!='-') 
                nongap_pos.push_back(i);

    no_query_gap_sequence.clear();
    for (i=0;i<nongap_pos.size();i++)
        no_query_gap_sequence+=sequence[nongap_pos[i]];

    if (outfile!="-") fp_out<<no_query_gap_sequence<<endl;
    else                cout<<no_query_gap_sequence<<endl;

    fp_out.close();
    sequence.clear();
    nongap_pos.clear();
    no_query_gap_sequence.clear();
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
    RemoveNonQueryPosition(infile,outfile);
    return 0;
}
