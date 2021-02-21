const char* docstring=""
"rFUpred seq.ct seq.FU\n"
"    Split input RNA secondary structure file seq.ct into domains by\n"
"    the FUpred algorithm. Output likely domain boundaries to seq.FU.\n"
"    Each domain must have at least one base pair.\n"
"\n"
"Output format:\n"
"    FUscore - a smaller Folding Unit score means a more likely boundary.\n"
"              FUscore>=1 means the predicted boundary is unreliable\n"
"    DC      - discontinuity. C for two continuous domains;\n"
"              D for a discontinous domain inserted by another domain.\n"
"    (domain1)(domain2) - boundaries of two domains without a linker\n"
"    (domain1)linker(domain2) - boundaries of two domains with linker.\n"
"             linker is in the format is i-j. if i<=j, i and j are the\n"
"             linker boundaries; if i=j-1, there is no linker.\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>

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
    long int i,j;
    long int L=0;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);
        if (line.size()==0 || line[0]=='#') continue;
        split(line,line_vec);
        if (line_vec.size()>=6)
        {
            L++;
            i=atol(line_vec[0].c_str());
            j=atol(line_vec[4].c_str());
            if (i<j)
            {
                resi1_vec.push_back(i-1);
                resi2_vec.push_back(j-1);
            }
        }
        line_vec.clear();
        line.clear();
    }
    vector<string>().swap(line_vec);
    fp_in.close();

    /* calculate FU score */
    vector<double>FUscore_list(L,0);
    vector<vector<double> >FUscore2d_mat(L,FUscore_list);
    vector<bool> FUsele_list(L,false);
    vector<vector<bool> >FUsele2d_mat(L,FUsele_list);
    vector<bool>().swap(FUsele_list);
    long int pos, bp;
    double N12,N1,N2;
    for (pos=1;pos<L;pos++)
    {
        N12=N1=N2=1; // 1 is pseudo count to avoid division of 0
        for (bp=0;bp<resi1_vec.size();bp++)
        {
            i=resi1_vec[bp];
            j=resi2_vec[bp];
            if      (j<pos) N1++;
            else if (i<pos  && pos<=j) N12++;
            else if (pos<=i) N2++;
        }
        FUscore_list[pos]=2*N12*(1./N1+1./N2);
        FUsele2d_mat[pos][pos]=(N1>1 && N2>1);
    }
    if (L>1) FUscore_list[0]=FUscore_list[1];
    
    /* calculate FU score 2d */
    long int pos1,pos2;
    double Nn2,N2c,Nn,Nc,Nnc;
    for (pos1=1;pos1<L;pos1++)
    {
        FUscore2d_mat[pos1][pos1]=FUscore_list[pos1];
        for (pos2=pos1+1;pos2<L;pos2++)
        {
            Nn=Nnc=N2=Nc=1;
            Nn2=N2c=1;
            for (bp=0;bp<resi1_vec.size();bp++)
            {
                i=resi1_vec[bp];
                j=resi2_vec[bp];
                if      (j<pos1) Nn++;
                else if (i<pos1 && pos1<=j && j<pos2) Nn2++;
                else if (i<pos1 && pos2<=j) Nnc++;
                else if (pos1<=i && j<pos2) N2++;
                else if (pos1<=i && i<pos2 && pos2<=j) N2c++;
                else if (pos2<=i) Nc++;
            }
            FUscore2d_mat[pos1][pos2]=FUscore2d_mat[pos2][pos1]=
                //2*(Nn2+N2c)*(1./(Nn+Nc)+1./(2*Nnc)+1./N2); // in RNA
                2*(Nn2+N2c)*(1./(Nn+Nc+2*Nnc)+1./N2); // in protein
            FUsele2d_mat[pos1][pos2]=(Nn>1 && Nc>1 && N2>1 && Nnc>1);
        }
        if (L>1) FUscore2d_mat[pos1][0]=
                 FUscore2d_mat[0][pos1]=FUscore2d_mat[pos1][1];
    }
    if (L>1) FUscore2d_mat[0][0]=FUscore2d_mat[1][1];

    /* sort position by FU score */
    vector<long int>().swap(resi1_vec);
    vector<long int>().swap(resi2_vec);
    vector<pair<double,string> > linker_list;
    vector<char> type_list;
    //vector<pair<double,vector<long int> > > resi_list;
    //vector<long int> resi_vec;
    long int mid,resi1,resi2;
    double score;
    bool accept;
    stringstream ss;
    for (pos=2;pos<L-1;pos++)
    {
        accept=true;
        score=FUscore_list[pos];
        if (score>=FUscore_list[pos-1] ||
            score> FUscore_list[pos+1]) continue;
        resi1=pos;
        for (resi2=pos;resi2<L-1;resi2++)
        {
            if (score>FUscore_list[resi2+1])
            {
                accept=false;
                break;
            }
            else if (score<FUscore_list[resi2+1]) break;
        }
        if (accept==false) continue;
        mid=(resi1+resi2)/2;
        if (!FUsele2d_mat[mid][mid]) continue;
        ss<<setiosflags(ios::fixed)<<setprecision(6)<<score
          <<"\tC\t(1-"<<mid<<")("<<mid+1<<"-"<<L<<")\t(1-"<<resi1
          <<")"<<resi1+1<<"-"<<resi2<<"("<<resi2+1<<"-"<<L<<")\t";
        ss.flush();
        linker_list.push_back(make_pair(score,ss.str()));
        type_list.push_back('c');
        ss.str("");
        //for (i=resi1;i<=resi2;i++) resi_vec.push_back(i);
        //resi_list.push_back(make_pair(score,resi_vec));
        //resi_vec.clear();
    }

    /* sort position by FU score 2d */
    long int midn,midc,resi1n,resi2n,resi1c,resi2c;
    for (pos1=2;pos1<L-1;pos1++)
    {
        resi1n=resi2n=pos1;
        for (pos2=pos1+3;pos2<L-1;pos2++)
        {
            resi1c=resi2c=pos2;
            score=FUscore2d_mat[pos1][pos2];
            accept=true;
            if (score>=FUscore_list[pos1] ||
                score>=FUscore_list[pos2] ||
                score>=FUscore2d_mat[pos1-1][pos2-1] ||
                score>=FUscore2d_mat[pos1-1][pos2]   ||
                score>=FUscore2d_mat[pos1][pos2-1]   ||
                score> FUscore2d_mat[pos1+1][pos2]   ||
                score> FUscore2d_mat[pos1][pos2+1]   ||
                score> FUscore2d_mat[pos1+1][pos2+1]) continue;
            for (resi2n=pos1;resi2n<=pos2;resi2n++)
            {
                if (score>=FUscore_list[resi2n] ||
                    score> FUscore2d_mat[resi2n+1][pos2])
                {
                    accept=false;
                    break;
                }
                else if (score<FUscore2d_mat[resi2n+1][pos2]) break;
            }
            if (accept==false) continue;
            for (resi2c=pos2;resi2c<L-1;resi2c++)
            {
                if (score>=FUscore_list[resi2c] ||
                    score> FUscore2d_mat[resi2n][resi2c+1])
                {
                    accept=false;
                    break;
                }
                if (score<FUscore2d_mat[resi2n][resi2c+1]) break;
            }
            if (accept==false) continue;
            if (score>=FUscore2d_mat[resi1n-1][resi2c+1] ||
                score>=FUscore2d_mat[resi2n+1][resi1c-1] ||
                score>=FUscore2d_mat[resi2n+1][resi2c+1]) continue;
            for (i=resi1n;i<=resi2n;i++)
            {
                if (score<FUscore2d_mat[i][resi1c-1] ||
                    score<FUscore2d_mat[i][resi2c+1]) continue;
                accept=true;
                break;
            }
            if (accept==false) continue;
            for (j=resi1c;j<=resi2c;j++)
            {
                if (score<FUscore2d_mat[resi1n-1][j] ||
                    score<FUscore2d_mat[resi2n+1][j]) continue;
                accept=true;
                break;
            }
            if (accept==false) continue;
            midn=(resi1n+resi2n)/2;
            midc=(resi1c+resi2c)/2;
            if (!FUsele2d_mat[midn][midc]) continue;
            ss<<setiosflags(ios::fixed)<<setprecision(6)<<score
              <<"\tD\t(1-"<<midn<<","<<midc+1<<"-"<<L<<")("<<midn+1
              <<"-"<<midc<<")\t(1-"<<resi1n<<","<<resi2c+1<<"-"<<L<<")"
              <<resi1n+1<<"-"<<resi2n<<","<<resi1c+1<<"-"<<resi2c
              <<"("<<resi2n+1<<"-"<<resi1c<<")\t";
            ss.flush();
            linker_list.push_back(make_pair(score,ss.str()));
            type_list.push_back('D');
            ss.str("");
            //for (i=resi1n;i<=resi2n;i++) resi_vec.push_back(i);
            //for (j=resi1c;j<=resi2c;j++) resi_vec.push_back(j);
            //resi_list.push_back(make_pair(score,resi_vec));
            //resi_vec.clear();
        }
    }
    vector<double>().swap(FUscore_list);
    vector<vector<double> >().swap(FUscore2d_mat);
    vector<vector<bool> >().swap(FUsele2d_mat);

    /* output result */
    sort(linker_list.begin(),linker_list.end());
    //sort(resi_list.begin(),resi_list.end());
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    line="#FUscore\tDC\t(domain1)(domain2)\t(domain1)linker(domain2)";
    if (outfile=="-") cout<<line<<endl;
    else            fp_out<<line<<endl;
    vector<bool> sele_list(L,false);
    int total_accepted=0;
    for (pos=0;pos<linker_list.size();pos++)
    {
        //accept=true;
        //for (i=0;i<resi_list[pos].second.size();i++)
        //{
            //if (sele_list[resi_list[pos].second[i]]==false) continue;
            //accept=false;
            //break;
        //}
        //if (accept==false) continue;

        if (total_accepted<10) total_accepted++;
        else if (linker_list[pos].first>=1) break;

        if (outfile=="-") cout<<linker_list[pos].second<<endl;
        else            fp_out<<linker_list[pos].second<<endl;
        
        //for (i=0;i<resi_list[pos].second.size();i++)
            //sele_list[resi_list[pos].second[i]]=true;
    }
    vector<bool>().swap(sele_list);
    vector<pair<double,string> >().swap(linker_list);
    //vector<pair<double,vector<long int> > >().swap(resi_list);
    fp_out.close();
    line.clear();
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
