const char* docstring=""
"fastNf seq.aln 0.8\n"
"    calculate Nf, number of effective sequence defined in gremlin,\n"
"    using 0.8 (default) sequence identity cutoff.\n"
"\n"
"fastNf seq.aln 0.8 0\n"
"    The third optional argument is normalization methods:\n"
"       0 - normalized by L^0.5 (default)\n"
"       1 - normalized by L\n"
"       2 - do not normalize\n"
"\n"
"fastNf seq.aln 0.8 0 128\n"
"    The fourth optional argument is target Nf. Stop Nf calculation if\n"
"    it is already greater than target Nf. default is 0 (no target Nf)\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <map>
#include <bits/stdc++.h> 

using namespace std;

const char *aa_list="-ACDEFGHIKLMNPQRSTVWY";

template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
}

template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
}

inline bool iverson_bracket(char *aln_n,char  *aln_m,const int L,const int maxLdiff)
{
    static int i=0;
    int Ldiff=0;
    for (i=0;i<L;i++)
    {
        Ldiff+=(aln_n[i]!=aln_m[i]);
        if (Ldiff>maxLdiff) return 0; // I[S_{m.n} >= Scut]
    }
    return 1;
}

int aa2int(const string&sequence,char *aln_n)
{
    int i=0;
    int j=0;
    int a=0;
    for (;i<sequence.size();i++)
    {
        if (sequence[i]=='-' || ('A'<=sequence[i] && sequence[i]<='Z'))
        {
            aln_n[j]=0;
            for (a=0;a<21;a++)
            {
                if (sequence[i]==aa_list[a])
                {
                    aln_n[j]=a;
                    break;
                }
            }
            j++;
        }
    }
    return j; // sequence length
}

/* return true if length matches L */
bool getupperseq(const string&sequence,string&upperseq,const int L)
{
    int i=0;
    int j=0;
    int a=0;
    upperseq.assign(L,'-');
    for (;i<sequence.size();i++)
    {
        if (sequence[i]=='-' || ('A'<=sequence[i] && sequence[i]<='Z'))
        {
            if (j==L) return false;
            upperseq[j++]=sequence[i];
        }
    }
    return (j==L);
}

double fastNf(const string infile, const double id_cut=0.8, const int norm=0,
    double target_Nf=0)
{
    /* read alignment */
    size_t i,j; // index of residue
    size_t m,n; // index of sequence
    size_t Nseq=0;
    size_t L=0;
    vector<string>aln;
    string sequence,upperseq;
    ifstream fp;
    if (infile!="-") fp.open(infile.c_str(),ios::in);
    L=0;
    while ((infile!="-")?fp.good():cin.good())
    {
        if (infile!="-") getline(fp,sequence);
        else getline(cin,sequence);

        if (sequence.length()==0||sequence[0]=='>') continue;
        if (L==0)
        {
            for (i=0;i<sequence.size();i++) 
                L+=(sequence[i]=='-' || ('A'<=sequence[i] && sequence[i]<='Z'));
        }
        if (!getupperseq(sequence,upperseq,L))
        {
            cerr<<"ERROR! length (L="<<L
                <<" mismatch for sequence "<<aln.size()<<endl;
            exit(0);
        }
        aln.push_back(upperseq);
    }
    if (infile!="-") fp.close();

    /* convert to int */
    Nseq=aln.size();
    char **msa;
    NewArray(&msa, Nseq, L);
    for (n=0;n<Nseq;n++) aa2int(aln[n],msa[n]);
    vector<string>().swap(aln);

    /* scale target_Nf by L */
    if (target_Nf>0)
    {
        if (norm==0) target_Nf*=sqrt(L);
        else if (norm==1) target_Nf*=L;
    }

    /* calculate weight */
    size_t *weight_list=new size_t[Nseq];
    for (n=0;n<Nseq;n++) weight_list[n]=1;
    double Nf=0;
    size_t maxLdiff=(1-id_cut)*L;
    bool geScut=false; // greater than or equal to seqID cut?
    for (n=0;n<Nseq;n++)
    {
        for (m=n+1;m<Nseq;m++)
        {
            geScut=iverson_bracket(msa[n],msa[m],L,maxLdiff);
            weight_list[n]+=geScut;
            weight_list[m]+=geScut;
        }
        Nf+=1./weight_list[n];
        if (target_Nf>0 && Nf>target_Nf) break;
    }

    /* normalize Nf */
    if (norm==0) Nf/=sqrt(L);
    else if (norm==1) Nf/=L;

    /* clean up */
    DeleteArray(&msa,Nseq);
    delete[]weight_list;
    return Nf;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    double id_cut=0.8; // defined by gremlin
    int norm=0; // 0 - L^0.5, 1 - L, 2 - no normalize
    double target_Nf=0;
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string infile=argv[1];
    if (argc>2) id_cut=atof(argv[2]);
    if (id_cut>1) id_cut/=100.;
    if (argc>3) norm=atoi(argv[3]);
    if (argc>4) target_Nf=atof(argv[4]);
    
    /* calculate Nf*/
    double Nf=fastNf(infile,id_cut,norm,target_Nf);
    cout<<Nf<<endl;
    return 0;
}
