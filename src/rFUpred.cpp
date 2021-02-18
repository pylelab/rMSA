const char* docstring=""
"rFUpred seq.ct seq.FU\n"
"    domain partition of RNA secondary structure in ct format\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <algorithm>

using namespace std;

void split(const string &line, vector<string> &line_vec)
{
    bool within_word = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (line[pos]==' ' || line[pos]=='\t')
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

void rFUpred(const string infile="-", const string outfile="-")
{
    /* parse input file */
    ifstream fp_in;
    ofstream fp_out;
    string line;
    vector<string> line_vec;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    vector<long int> resi1_vec;
    vector<long int> resi2_vec;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);
        if (line.size()==0 || line[0]=='#') continue;
        split(line,line_vec);
        if (line_vec.size()>=6)
        {
            resi1_vec.push_back(atol(line_vec[0].c_str()));
            resi2_vec.push_back(atol(line_vec[4].c_str()));
        }
        line_vec.clear();
        line.clear();
    }
    vector<string>().swap(line_vec);
    string().swap(line);
    fp_in.close();


    /* calculate FU score */
    long int L=resi1_vec.size();
    vector<bool> tmp_bool_vec(L,0);
    vector<vector<bool> > contact_map(L,tmp_bool_vec);
    vector<bool>().swap(tmp_bool_vec);
    long int i,j,pos;
    for (pos=0;pos<L;pos++)
    {
        i=resi1_vec[pos];
        j=resi2_vec[pos];
        if (j==0) continue;
        contact_map[i-1][j-1]=1;
    }
    
    vector<double>FUscore_list(L,0);
    double N12,N1,N2;
    for (pos=1;pos<L;pos++)
    {
        N12=N1=N2=1; // 1 is pseudo count to avoid division of 0
        for (i=0;i<pos;i++) for (j=pos;j<L;j++) N12+=contact_map[i][j];
        for (i=0;i<pos;i++) for (j=0;j<pos;j++) N1 +=contact_map[i][j];
        for (i=pos;i<L;i++) for (j=pos;j<L;j++) N2 +=contact_map[i][j];
        FUscore_list[pos]=2*N12*(1./N1+1./N2);
    }
    FUscore_list[0]=FUscore_list[1];

    vector<pair<double,long int> > score_list;
    for (pos=2;pos<L-1;pos++)
        if (FUscore_list[pos]<=FUscore_list[pos-1] &&
            FUscore_list[pos]<=FUscore_list[pos+1])
            score_list.push_back(make_pair(FUscore_list[pos],pos));

    vector<pair<double,pair<long int,long int> > > linker_list;
    vector<pair<double,long int> >linker;
    double score;
    for (pos=0;pos<score_list.size();pos++)
    {
        score=score_list[pos].first;
        i    =score_list[pos].second;
        if (linker.size() && linker.back().second!=i-1)
        {
            linker_list.push_back(make_pair(linker[0].first,
                make_pair(linker[0].second,linker.back().second)));
            linker.clear();
        }
        linker.push_back(make_pair(score,i));
    }
    if (linker.size())
        linker_list.push_back(make_pair(score,
            make_pair(linker[0].second,linker.back().second)));
    vector<pair<double,long int> >().swap(linker);
    vector<pair<double,long int> >().swap(score_list);
    sort(linker_list.begin(),linker_list.end());

    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    long int mid;
    if (outfile=="-") cout<<"#(domain1)(domain2)\t"
        <<"(domain1)linker(domain2)\tFUscore"<<endl;
    else            fp_out<<"#(domain1)(domain2)\t"
        <<"(domain1)linker(domain2)\tFUscore"<<endl;
    for (pos=0;pos<linker_list.size();pos++)
    {
        score=linker_list[pos].first;
        i=linker_list[pos].second.first;
        j=linker_list[pos].second.second;
        mid=(long int)((i+j)/2);
        if (outfile=="-") cout<<"(1-"<<mid<<")("<<mid+1<<"-"<<L<<")\t(1-"
            <<i<<")"<<i+1<<"-"<<j<<"("<<j+1<<"-"<<L<<")\t"<<score<<endl;
        else            fp_out<<"(1-"<<mid<<")("<<mid+1<<"-"<<L<<")\t(1-"
            <<i<<")"<<i+1<<"-"<<j<<"("<<j+1<<"-"<<L<<")\t"<<score<<endl;
    }
    vector<pair<double,pair<long int,long int> > >().swap(linker_list);
    fp_out.close();
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
    string outfile=(argc<=2)?"-":argv[2];
    rFUpred(infile,outfile);
    return 0;
}
