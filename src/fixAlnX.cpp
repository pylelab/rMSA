const char* docstring=""
"fixAlnX seq.aln N fix.aln\n"
"    replace residue type 'N' in alignment file seq.afa by the most frequent\n"
"    residue type in that position\n"
"\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <map>
#include <climits>
#include <iomanip>

using namespace std;

size_t fixAlnX(const string infile, const char replace, const string outfile)
{
    /* read aln file */
    vector<string> aln;
    vector<string> header_list;
    string sequence;
    size_t L=0;
    ifstream fp;
    if (infile!="-") fp.open(infile.c_str(),ios::in);
    while ((infile!="-")?fp.good():cin.good())
    {
        if (infile!="-") getline(fp,sequence);
        else getline(cin,sequence);

        if (sequence.length()==0) continue;
        else if (sequence[0]=='>')
        {
            header_list.push_back(sequence);
            continue;
        }
        if (L==0) L=sequence.length();
        else if (sequence.length()!=L)
        {
            cerr<<"ERROR! length not match for sequence\n"<<aln.size();
            exit(0);
        }

        aln.push_back(sequence);
        if (aln.size()==INT_MAX)
        {
            cerr<<"WARNING! Cannot read beyond sequence number"<<INT_MAX<<endl;
            break;
        }
    }
    fp.close();

    /* get residue type */
    string aa_list="";
    vector<vector<size_t> >aa_count_mat;
    vector<size_t> aa_count_vec(L,0);
    size_t seq_num=aln.size();
    size_t a,i,j;
    char aa;
    bool new_aa;
    for (i=0;i<seq_num;i++)
    {
        for (j=0;j<L;j++)
        {
            aa=aln[i][j];
            if (aa=='-' || aa=='.' || aa==replace) continue;
            new_aa=true;
            for (a=0;a<aa_list.size();a++)
            {
                if (aa==aa_list[a])
                {
                    new_aa=false;
                    break;
                }
            }
            if (new_aa)
            {
                aa_list+=aa;
                aa_count_mat.push_back(aa_count_vec);
            }
            aa_count_mat[a][j]++;
        }
    }

    /* get the most frequent aa type per column */
    for (j=0;j<L;j++)
        for (a=0;a<aa_list.size();a++)
            if (aa_count_mat[a][j]>aa_count_mat[aa_count_vec[j]][j])
                aa_count_vec[j]=a;
    
    vector<vector<size_t> >().swap(aa_count_mat);
    
    /* output */
    string txt;
    for (i=0;i<seq_num;i++)
    {
        if (header_list.size()>i) txt+=header_list[i]+'\n';
        for (j=0;j<L;j++)
        {
            if (aln[i][j]!=replace) txt+=aln[i][j];
            else txt+=aa_list[aa_count_vec[j]];
        }
        txt+='\n';
    }

    if (outfile=="-") cout<<txt<<flush;
    else
    {
        ofstream fp_out;
        fp_out.open(outfile.c_str(),ofstream::out);
        fp_out<<txt;
        fp_out.close();
    }

    /* clean up */
    vector<string> ().swap(aln);
    vector<string> ().swap(header_list);
    sequence.clear();
    aa_count_vec.clear();
    txt.clear();
    return aln.size();
}


int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string infile =argv[1];
    char   replace=argv[2][0];
    string outfile=(argc<=3)?"-":argv[3];
    fixAlnX(infile,replace,outfile);
    return 0;
}
