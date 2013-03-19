/*
 * =====================================================================================
 *
 *       Filename:  tannbc.cc
 *
      Tree augmented naive bayes classifier
 *    Description:  Input to the program is the data matrix with the class information 
 *    in the first column and the rest of columns has the feature values. One line of the 
 *    inputfile looks like below:
 *    r,n,y,n,y,y,y,n,n,n,y,?,y,y,y,n,y
 *    All variable values are single characters including class variable
 *
 *    Program deals with two class classification problem. Attribute values should be
 *    categorical.
 *
 *    The program takes advantage of the fact that domain from where values of attributes 
 *    is taken is the same('y','n',?') for all the attributes. 
 *
 *    Extension to attribute values with strings or nominal values (numbers representing
 *    categories) is straightforward. 
 *    This program is economized to take character inputs to avoid string comparison 
 *    operations.
 *
 *        Version:  1.0
 *        Created:  11/04/2011 01:28:19 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jeyanthi S N (jey), jsalemna@eecs.wsu.edu
 *        Company:  Washington State University, Pullman
 *
 * =====================================================================================
 */
#include <iostream>
#include <assert.h>
#include <fstream>
#include <cstdlib>//for rand
#include <algorithm>
#include <iterator>
#include <utility>
#include <cmath>//for log function
#include <queue>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <stack>
#include <sys/resource.h>//for getrusage

#define SPACES " \t\n\r"
using namespace std;
typedef vector <char> v_char;
/*
 * data structure used for MST construction
 */
class pqds{
        public:
                int sv; int ev; double cmival;

                bool operator< (const pqds &p) const
                {
                        return cmival < p.cmival;
                }
};

int main(int argc, char *argv[])
{
        int numline=0, i, ccnt[2]={0,0}, numofsamp; 
        double priorc1, priorc2;
        ifstream ifs, ifs1;
        string token, testsample;
        vector <v_char> inputmat;
        vector <v_char> bk_inputmat;
        vector<string>attrnames;
        vector<double> prod;
        v_char tempvchar;
        v_char classvec;
        vector<pair<int, char> > attidxvals;//used for passing values to find their count
        void remove_space(string &);
        void cmi(vector<v_char>, string, vector<string>, int[2], vector<double>&);
        int nbc(vector<v_char>, string, int[2], vector<double>&);
        int calc_freq(vector<pair<int, char> >, vector<v_char>);
        void rocplot(vector <pair <double, double> >, vector<pair <double, double> >,
                        v_char, int[2]);
#ifdef TREEDRAW
        if (argc != 3)
        {
                cerr<<"Usage: Input file, file with attr names"<<endl;
                exit(0);
        }
        ifs1.open(argv[2]);
        getline(ifs1,token);
        while(!ifs1.eof())
        {
                remove_space(token);
                attrnames.push_back(token);
                getline(ifs1,token);
        }
#else
        int j, k, sum_cmi = 0, sum_nbc = 0;
        double normconst, secelap_tan=0, secelap_nbc=0;
        struct rusage resuse[2];
        vector <pair <double, double> >cmi_probval;
        vector <pair <double, double> >nbc_probval;
        if (argc != 2)
        {
                cerr<<"Usage: Input file"<<endl;
                exit(0);
        }
#endif
        ifs.open(argv[1]);

        getline(ifs, token);
        numline++;
        while (!ifs.eof())
        {
                remove_space(token);
                if (token.size() == 0)
                {
                        cerr<<"cannot have blank lines in the input file"<<endl;
                        exit(0);
                }
                if (token[0] != ',')
                        classvec.push_back(token[0]);
                else
                {
                        cerr<<"First character of a line is class"<<endl;
                        exit(0);
                }
                for (i = 0; i < (int)token.length(); i++)
                {
                        if (token[i]!=',') 
                                tempvchar.push_back(token[i]);

                }
                inputmat.push_back(tempvchar);
                tempvchar.clear();
                getline(ifs, token);
                numline++;
        }
        numline--;
        numofsamp = numline;
        assert(numofsamp == (int)inputmat.size());
        //verify if the number of characters in each line is the same
        for (i = 0; i < (int)inputmat.size()-1; i++)
        {
                if(inputmat[i].size()-inputmat[i+1].size() != 0)
                {
                        cerr<<"The number of columns do not match "
                                <<i<<" "<<inputmat[i].size()
                                <<" "<<i+1<<" "<<inputmat[i+1].size()<<endl;
                        exit(0);
                } 
        } 
        attidxvals.push_back(make_pair(0, 'd'));
        ccnt[0] = calc_freq(attidxvals, inputmat);
        attidxvals.clear();
        attidxvals.push_back(make_pair(0, 'r'));
        ccnt[1] = calc_freq(attidxvals, inputmat);
        attidxvals.clear();
#ifdef LAPLACE
        priorc1 = (double)(ccnt[0]+1)/(double)(numofsamp+2);//laplace estimate
        priorc2 = (double)(ccnt[1]+1)/(double)(numofsamp+2);
#else
        priorc1 = (double)(ccnt[0])/(double)(numofsamp);
        priorc2 = (double)(ccnt[1])/(double)(numofsamp);
#endif
#ifdef TREEDRAW 
        cmi(inputmat, "", attrnames, ccnt, prod);
#else
        /*
         * Leave one out testing 
         */

        for (i = 0; i < (int)inputmat.size(); i++)
        {
                //form the current dataset 
                for(j = 0, k = 0;j < (int)inputmat.size();j++)
                {
                        if(j!=i)
                        {
                                tempvchar = inputmat[j];
                                bk_inputmat.push_back(tempvchar);   
                                k++;
                        }
                }
                //check if the reduced inputmat size is correct
                if (k != (int) inputmat.size()-1)
                {
                        cerr<<"The k value: "<<k<<" The j value: "<<j<<endl;
                        exit(0);
                }
                //construct the test sample
                for (j = 0; j < (int)inputmat[i].size();j++)
                        testsample+=inputmat[i][j];

                prod.push_back(priorc1);
                prod.push_back(priorc2);
                //attrnames will not be used 
                getrusage(RUSAGE_SELF, resuse+0);
                cmi(bk_inputmat, testsample, attrnames, ccnt, prod);
                getrusage(RUSAGE_SELF, resuse+1);
                secelap_tan += (double)(
                (resuse[1].ru_utime.tv_usec*1e-6+resuse[1].ru_stime.tv_usec*1e-6+
                resuse[1].ru_utime.tv_sec+resuse[1].ru_stime.tv_sec)-
                (resuse[0].ru_utime.tv_usec*1e-6+resuse[0].ru_stime.tv_usec*1e-6+
                resuse[0].ru_utime.tv_sec+resuse[0].ru_stime.tv_sec));
                //marginal likelihood of the evidence
                normconst = prod[0] + prod[1];
                prod[0] = prod[0]/normconst;
                prod[1] = prod[1]/normconst;
                if ((prod[0] > prod[1] && testsample[0] == 'd') ||
                                (prod[0] < prod[1] && testsample[0] == 'r'))
                        sum_cmi += 1;
                cmi_probval.push_back(make_pair(prod[0], prod[1]));
                prod.clear();
                prod.push_back(priorc1);
                prod.push_back(priorc2);
                getrusage(RUSAGE_SELF, resuse+0);
                nbc(bk_inputmat, testsample, ccnt, prod);
                getrusage(RUSAGE_SELF, resuse+1);
                secelap_nbc += (double)(
                (resuse[1].ru_utime.tv_usec*1e-6+resuse[1].ru_stime.tv_usec*1e-6+
                resuse[1].ru_utime.tv_sec+resuse[1].ru_stime.tv_sec)-
                (resuse[0].ru_utime.tv_usec*1e-6+resuse[0].ru_stime.tv_usec*1e-6+
                resuse[0].ru_utime.tv_sec+resuse[0].ru_stime.tv_sec));

                //marginal likelihood of the evidence
                normconst = prod[0] + prod[1];
                prod[0] = prod[0]/normconst;
                prod[1] = prod[1]/normconst;
                if ((prod[0] > prod[1] && testsample[0] == 'd') ||
                                (prod[0] < prod[1] && testsample[0] == 'r'))
                        sum_nbc += 1;
                nbc_probval.push_back(make_pair(prod[0], prod[1]));
                bk_inputmat.clear();
                testsample.clear();
                tempvchar.clear();
                prod.clear();
        }//end for i
        cout<<"The l-o-o acc for TAN: "<< (double)sum_cmi/(double)(inputmat.size())<<endl;
        cout<<"The l-o-o acc for NBC: "<< (double)sum_nbc/(double)(inputmat.size())<<endl;
        cout<<"Average time for one call to TAN:"<<secelap_tan/(double)numofsamp<<" seconds"<<endl;
        cout<<"Average time for one call to NBC:"<<secelap_nbc/(double)numofsamp<<" seconds"<<endl;
        rocplot(cmi_probval, nbc_probval, classvec, ccnt);
#endif

}//end main

/*
 * Roc plot for the classifiers for both target classes
 */
void rocplot(vector< pair<double, double> > cmi_probval, 
                vector< pair<double, double> >nbc_probval, v_char classvec, int ccnt[2])
{
        vector< pair<double, double> >outcmi_d, outnbc_d, outcmi_r, outnbc_r; 
        double t=0.0, inc=0.1, tmax=1.0;
        int i;
        int tp_cmi_d, tp_nbc_d, fp_cmi_d, fp_nbc_d, tp_cmi_r, tp_nbc_r, fp_cmi_r, fp_nbc_r;
        //initialize the TP and FP rates
        while (t <= tmax)
        {
                tp_cmi_d= tp_nbc_d= fp_cmi_d= fp_nbc_d= tp_cmi_r= 
                        tp_nbc_r= fp_cmi_r= fp_nbc_r=0;
                for ( i = 0; i < (int)classvec.size(); i++)
                {
                        if (cmi_probval[i].first >= t)
                        {
                                if(classvec[i] == 'd')
                                        tp_cmi_d++;
                                else
                                        fp_cmi_d++;
                        }
                        if(cmi_probval[i].second >= t)
                        {
                                if(classvec[i] == 'r')
                                        tp_cmi_r++;
                                else
                                        fp_cmi_r++;

                        }
                        if (nbc_probval[i].first >= t)
                        {
                                if(classvec[i] == 'd')
                                        tp_nbc_d++;
                                else
                                        fp_nbc_d++;
                        }
                        if (nbc_probval[i].second >= t)
                        {
                                if(classvec[i] == 'r')
                                        tp_nbc_r++;
                                else
                                        fp_nbc_r++;
                        }
                }//end for
                //save the TP rate and FP rate
                outcmi_d.push_back(make_pair((double)fp_cmi_d/(double)ccnt[1], (double)tp_cmi_d/(double)ccnt[0]));
                outnbc_d.push_back(make_pair((double)fp_nbc_d/(double)ccnt[1], (double)tp_nbc_d/(double)ccnt[0]));
                outcmi_r.push_back(make_pair((double)fp_cmi_r/(double)ccnt[0], (double)tp_cmi_r/(double)ccnt[1]));
                outnbc_r.push_back(make_pair((double)fp_nbc_r/(double)ccnt[0], (double)tp_nbc_r/(double)ccnt[1]));
                t = t + inc;
        }//end while
        ofstream ofs1("cmiroc_d.txt");
        ofstream ofs2("nbcroc_d.txt");
        ofstream ofs3("cmiroc_r.txt");
        ofstream ofs4("nbcroc_r.txt");
        //write to files
        for(i = 0; i < (int) outcmi_d.size(); i++)
        {
                ofs1<<outcmi_d[i].first<<" "<<outcmi_d[i].second<<endl;
                ofs2<<outnbc_d[i].first<<" "<<outnbc_d[i].second<<endl;
                ofs3<<outcmi_r[i].first<<" "<<outcmi_r[i].second<<endl;
                ofs4<<outnbc_r[i].first<<" "<<outnbc_r[i].second<<endl;
        }
        ofs1.close();
        ofs2.close();
        ofs3.close();
        ofs4.close();
}
/*
 * Function to calculate the CMI for each pair of attributes in the input matrix
 */
void cmi(vector<v_char>inputmat, string testsample, vector<string> attrnames, int ccnt[2], 
                vector<double>& prod)
{
        void write_dotout(map<int, list<int> >, int, vector<string>);
        void dfs ( map<int, list<int> >, map<int, list<int> >&, vector<pair<int, int> >&);
        priority_queue<pqds> pq;
        typedef list <int> intlisttype;
        intlisttype neighbors;
        map <int, list<int> > adjlist, unsymmadjlist; 
        map <int, list<int> >::iterator adjlist_it, adjlist_it1;
        map<int, map<string,vector<double> > > cpts;
        map<string, vector<double> > strmap;
        vector <pair<int, int> > parentlist;
        vector <double> problist;
        vector <double>::iterator problist_it;
        map<string, vector<double> >::iterator strmap_it;
        int numofattr, i, j, k, l, m, numofsamp, ct;
        int t1, t2, t3;//terms under one summation of cmi formula
        int prsntflag_sv, prsntflag_ev; 
        double log2base10 = log10((double)2.0);
        double temp1, temp2, temp3, temp4, sum=0.0; 
        pqds temppqds;
        vector<pair<int, char> > attidxvals;//used for passing values to find their count
        string attrvals="yn?";
        string cvals ="dr";//this order matters
        //variable tempstr to construct the parent vals while inferencing
        string tempstr;
        int calc_freq(vector<pair<int, char> >, vector<v_char>);
        void calc_cpts(vector<pair<int, int> >, vector<v_char >, map<int, map<string, 
                        vector<double> > >&, string, string, int [2]);

        assert(!inputmat.empty());
        numofsamp = inputmat.size();
        numofattr = (int) (inputmat[0].size()-1);//exclude the class info 
        //first two for loops correspond to taking pair of attributes
        for ( i = 1; i < numofattr;i++)//inputmat has features from 2nd column
        {
                for (j = i+1; j <= numofattr; j++)
                {
                        sum = 0.0;//var to accumulate individual terms of cmi formula
                        //for all the class, and other feature values from the domain
                        for ( k = 0; k < (int)cvals.size();k++)
                        {
                                ct = ccnt[k];//common term involving count of class
                                for (l = 0; l < (int)attrvals.size();l++)
                                {
                                        attidxvals.clear();
                                        attidxvals.push_back(make_pair(0, cvals[k]));
                                        attidxvals.push_back(make_pair(i, attrvals[l]));
                                        t2 = calc_freq(attidxvals, inputmat);
#ifdef LAPLACE
                                        t2 = t2 + 1;//laplace estimate
#endif
                                        for (m= 0; m < (int)attrvals.size();m++)
                                        {
                                                attidxvals.clear();
                                                attidxvals.push_back(make_pair(0, cvals[k]));
                                                attidxvals.push_back(make_pair(i, attrvals[l]));
                                                attidxvals.push_back(make_pair(j, attrvals[m]));
                                                t1 = calc_freq(attidxvals, inputmat);
#ifdef LAPLACE
                                                t1 = t1 + 1;//laplace estimate
#endif
                                                attidxvals.clear();
                                                attidxvals.push_back(make_pair(0, cvals[k]));
                                                attidxvals.push_back(make_pair(j, attrvals[m]));
                                                t3 = calc_freq(attidxvals, inputmat);
#ifdef LAPLACE
                                                t3 = t3 + 1;
                                                temp1 = ((double)t1/(double)(numofsamp+18));
                                                temp2 = (double)t1/(double)(ct+9);
                                                temp3 = (double)t2/(double)(ct+3);
                                                temp4 = (double)t3/(double)(ct+3);
#else
                                                temp1 = ((double)t1/(double)(numofsamp));
                                                temp2 = (double)t1/(double)(ct);
                                                temp3 = (double)t2/(double)(ct);
                                                temp4 = (double)t3/(double)(ct);

#endif
                                                temp2 = temp2/(temp3 * temp4);
                                                temp3 = log10(temp2)/log2base10;
                                                sum = sum + temp1*temp3;
                                        }//end for m
                                }//end for l
                        }//end for k
                        temppqds.sv = i;
                        temppqds.ev = j;
                        temppqds.cmival = sum;
                        pq.push(temppqds);
                }//end for j
        }// end for i
        temppqds = pq.top();

        /*
         * keep removing edges from the pq till n-1 edges are added or n unique vertices  
         * are added - MST construction
         */
        while (1)
        {
                prsntflag_sv =  prsntflag_ev = 0;
                //undirected graph construction
                temppqds = pq.top();
                pq.pop();
                adjlist_it = adjlist.find(temppqds.sv);
                adjlist_it1 = adjlist.find(temppqds.ev);
                if (adjlist_it != adjlist.end())
                        prsntflag_sv = 1;
                if (adjlist_it1 != adjlist.end())
                        prsntflag_ev = 1;
                if (prsntflag_sv == 1 && prsntflag_ev == 1)
                        continue;//both vertices are already considered
                if (prsntflag_sv == 0) 
                {
                        neighbors.clear();
                        neighbors.push_back(temppqds.ev);
                        adjlist[temppqds.sv]=neighbors;
                }
                else
                        (adjlist_it->second).push_back(temppqds.ev);
                if (prsntflag_ev == 0) 
                {
                        neighbors.clear();
                        neighbors.push_back(temppqds.sv);
                        adjlist[temppqds.ev]=neighbors;
                }
                else
                        (adjlist_it1->second).push_back(temppqds.sv);
                if ( (int) adjlist.size() == numofattr)
                        break;
        }//end while

        /*Depth first search - parent list is tuned for TAN - just one parent from
          features for each attribute, class parent would be added while inferencing
          */
        dfs(adjlist, unsymmadjlist, parentlist);
#ifdef TREEDRAW   
        write_dotout(unsymmadjlist, numofattr, attrnames);
#else
        int findflag;
        /*
         * calculate the conditional probability table for each of the dependencies of the DAG
         */
        calc_cpts(parentlist, inputmat, cpts, cvals, attrvals, ccnt); 

        /*
         * inferencing procedure follows for the passed testsample
         */
        assert((int)testsample.length()==numofattr+1);
        for (i  = 1; i < (int)testsample.length(); i++)
        {
                strmap = cpts[i];
                if(strmap.empty())
                {
                        cerr<<"The CPT for "<<i<<" is empty"<<endl;
                        exit(0);
                }
                for (j = 0; j <(int) parentlist.size(); j++)
                {
                        if (parentlist[j].first ==  i)
                                break;
                }

                for ( k = 0; k < (int)cvals.size();k++)
                {
                        tempstr.clear();

                        if (parentlist[j].second == 0)
                        {
                                tempstr += cvals[k];
                        }
                        else
                        {
                                tempstr += cvals[k];
                                tempstr += testsample[parentlist[j].second];
                        }
                        strmap_it = strmap.find(tempstr);
                        if(strmap_it == strmap.end())
                        {
                                cerr<<"tempstr: "<<tempstr<<endl;
                                cerr<<"string not found"<<endl;
                                exit(0);
                        }

                        problist = strmap[tempstr];
                        findflag = (int)attrvals.find(testsample[i]);
                        assert(findflag != (int)string::npos);
                        prod[k] *= problist[findflag];
                        tempstr.clear();
                }//end for k
        }//end for i
#endif
}//end of cmi

#ifndef TREEDRAW
void nbc(vector<v_char> inputmat, string testsample, int ccnt[2], vector<double>& prod)
{
        int i, j, k, findflag, numofsamp, numofattr;
        vector <pair<int, int> > parentlist;
        map<int, map<string, vector<double> > > cpts;
        map<string, vector<double> >strmap;
        map<string, vector<double> >::iterator strmap_it;
        string tempstr;
        vector <double> problist;
        string attrvals="yn?";
        string cvals ="dr";//this order matters
        vector<pair<int, char> > attidxvals;//used for passing values to find their count
        void calc_cpts(vector<pair<int, int> >, vector<v_char >, 
                        map<int, map<string, vector<double> > >&,
                        string, string, int [2]);
        int calc_freq(vector<pair<int, char> >, vector<v_char>);

        numofsamp = inputmat.size();
        numofattr = (int) (inputmat[0].size()-1);//exclude the class info

        /*class is the only parent for all the attributes - manually set the parent
         * to class id: 0
         */
        for (i  = 1; i <= numofattr; i++)
        {
                parentlist.push_back(make_pair(i, 0));
        }
        calc_cpts(parentlist, inputmat, cpts, cvals, attrvals, ccnt); 

        /*
         * inferencing procedure follows for the passed testsample
         */
        for (i  = 1; i < (int)testsample.length(); i++)
        {
                strmap = cpts[i];
                if(strmap.empty())
                {
                        cerr<<"The CPT of "<<i<<" is empty"<<endl;
                        exit(0);
                }
                for (j = 0; j <(int) parentlist.size(); j++)
                {
                        if (parentlist[j].first ==  i)
                                break;
                }

                for ( k = 0; k < (int)cvals.size();k++)
                {
                        tempstr.clear();
                        tempstr += cvals[k];
                        strmap_it = strmap.find(tempstr);
                        if(strmap_it == strmap.end())
                        {
                                cerr<<"tempstr: "<<tempstr<<endl;
                                cerr<<"string not found in CPT"<<endl;
                                exit(0);
                        }

                        problist = strmap[tempstr];
                        findflag = (int)attrvals.find(testsample[i]);
                        assert(findflag != (int)string::npos);
                        prod[k] *= problist[findflag];
                        tempstr.clear();
                }//end for k
        }//end for i

}//end nbc

/*
 * function calc_cpts to calculate and store the conditional probability tables based on the
 * dependencies from the graph (stored in parentlist) and the input dataset. Returns 'cpts'
 * - datastructure that holds all the cpts for all the features except for class feature. cpts
 * returned will be used for inferencing based on TAN and NBC
 *
 */
void calc_cpts(vector<pair<int, int> >parentlist, vector<v_char >inputmat, map<int, 
                map<string, vector<double> > >&cpts, string cvals, string attrvals, int ccnt[2])
{
        vector<pair<int, char> > attidxvals;//used for passing values to find their count
        vector <double> problist;
        map<string, vector<double> > strmap;
        vector <double>::iterator problist_it;
        map<string, vector<double> >::iterator strmap_it;
        string strparents;
        int i, j, k, t1, t2, l;
        int calc_freq(vector<pair<int, char> >, vector<v_char>);

        for ( i = 0; i < (int)parentlist.size();i++)//for attributes in parentlist
        {
                if (parentlist[i].second != 0)
                {
                        for(j = 0; j < (int)cvals.length();j++)//first parent - class
                        {
                                attidxvals.clear();
                                attidxvals.push_back(make_pair(0, cvals[j]));
                                for (k = 0;k < (int)attrvals.length();k++)//second parent from parentlist
                                {
                                        attidxvals.push_back(make_pair
                                                        (parentlist[i].second,attrvals[k]));
#ifdef LAPLACE
                                        t2 = calc_freq(attidxvals, inputmat)+3;
#else
                                        t2 = calc_freq(attidxvals, inputmat);
#endif
                                        problist.clear();
                                        for(l = 0; l < (int) attrvals.length();l++)//for each value of attribute
                                        {
                                                attidxvals.push_back(make_pair(parentlist[i].first,
                                                                        attrvals[l]));
#ifdef LAPLACE
                                                t1 = calc_freq(attidxvals, inputmat)+1;
#else
                                                t1 = calc_freq(attidxvals, inputmat);
#endif
                                                problist.push_back((double)t1/(double)t2);
                                                attidxvals.pop_back();
                                        }
                                        strparents += cvals[j];
                                        strparents += attrvals[k]; 
                                        assert(!problist.empty());
                                        strmap[strparents]=problist;
                                        //keep the datastructures ready for next iteration
                                        strparents.clear();
                                        problist.clear();
                                        //last entry correspond to one parent 
                                        attidxvals.pop_back();
                                }//end for k
                                attidxvals.clear();
                        }//end for j
                }//end if parentlist.second
                //case for root of the DFS tree - has class attribute as single parent 
                else
                {                
                        strmap.clear();
                        attidxvals.clear();
                        for(j = 0; j < (int)cvals.length();j++)//first parent - class
                        {
                                attidxvals.push_back(make_pair(0, cvals[j]));        
#ifdef LAPLACE
                                t2 = ccnt[j] + 3;
#else
                                t2 = ccnt[j];
#endif
                                problist.clear();
                                for (k = 0;k < (int)attrvals.length();k++)//for each value of the attr
                                {
                                        attidxvals.push_back(make_pair(parentlist[i].first,
                                                                attrvals[k]));        
#ifdef LAPLACE
                                        t1 = calc_freq(attidxvals, inputmat) + 1;
#else
                                        t1 = calc_freq(attidxvals, inputmat);
#endif
                                        problist.push_back((double)t1/(double)t2);
                                        attidxvals.pop_back();
                                }//end for k
                                strparents+=cvals[j];
                                assert(!problist.empty());
                                strmap[strparents]=problist;
                                strparents.clear();
                                problist.clear();
                                attidxvals.pop_back();
                        }//end for j
                }//end else
                assert(!strmap.empty());
                cpts[parentlist[i].first] = strmap;
                strmap.clear();
        }//end for i

}//end calc_cpts
#endif

/*
 * function dfs to find the directed traversal of the Max weighted tree 
 */
void dfs (map<int, list<int> >adjlist, map<int, list<int> >&outadjlist,
                vector<pair<int, int> > &parentlist)
{
        int rannum, numnodes, i;
        pair <int, int> tpair;
        list <int> tlist, tlisto;
        list <int>::iterator tlist_it, tlist_ito;
        map<int, list<int> >::iterator adjlist_it, adjlist_ito;
        //stack stores id of node and its parent
        stack<pair<int, int> > np_stack; 
        vector<int> visitedarr;
        vector<int>::iterator visitedarr_it;
        numnodes = (int)adjlist.size();
        //let a random node be the root of the tree
        rannum = (int)((rand()/(RAND_MAX+1.0))*(double)numnodes);
        //for the root of the tree, Class (0) is the parent
        adjlist_it = adjlist.begin();
        for (i = 0; i < rannum; i++)
                adjlist_it++;
        np_stack.push(make_pair(adjlist_it->first, 0));
        /*
         * stopping condition is the number of nodes added so far to the visited
         * array
         */
        while ((int)visitedarr.size() < numnodes)
        {
                if(np_stack.empty())
                {
                        cerr<<"Stack is empty: disconnected graph"<<endl;
                        exit(0);
                }

                tpair = np_stack.top();
                np_stack.pop();
                visitedarr_it = find(visitedarr.begin(), visitedarr.end(), tpair.first);
                //continue to next node if curr node is already visited
                if (visitedarr_it != visitedarr.end())
                        continue;
                visitedarr.push_back(tpair.first);
                adjlist_it = adjlist.find(tpair.first);
                assert(adjlist_it != adjlist.end());
                tlist = adjlist_it->second;
                //push child of the curr node along with its parent
                for(tlist_it=tlist.begin();tlist_it!=tlist.end();tlist_it++)
                {
                        np_stack.push(make_pair(*tlist_it, tpair.first));
                }
                //parent is a non-zero number, make an entry in outadj
                if (tpair.second !=0)
                {
                        adjlist_ito = outadjlist.find(tpair.second); 
                        if (adjlist_ito == outadjlist.end())
                        {
                                //create a temp list and insert back
                                tlisto.clear();
                                tlisto.push_back(tpair.first);
                                outadjlist[tpair.second] = tlisto;
                        }
                        else
                        {
                                outadjlist[tpair.second].push_back(tpair.first); 
                        }
                }//end if

                /*only the root attribtue will have a zero entry for parent 
                 * fine tuned for TAN based tree structure
                 */
                parentlist.push_back(make_pair(tpair.first, tpair.second));//parent in second entry
        }//end while
}//end dfs

/*
 * Function to count the occurrence of given set of attributes in the vector 
 * aiv in the inputmat
 * first attribute is the class
 */
int calc_freq(vector<pair<int, char> > aiv, vector<v_char>inputmat)
{
        int i, j, flag = 0, sum = 0;
        for (i = 0; i < (int)inputmat.size();i++)
        {
                for (j = 0; j < (int) aiv.size();j++)
                {
                        if (inputmat[i][aiv[j].first] == aiv[j].second)
                                flag++;
                }
                if (flag == (int)aiv.size())
                {
                        sum++;
                }
                flag = 0;
        }
        return sum;
}

/*function to remove extra leading and trailing spaces around
* a string
*/
void remove_space(string &node)
{
        size_t spacestart, spaceend;

        spacestart = node.find_first_not_of(SPACES);
        spaceend = node.find_last_not_of(SPACES);
        if (spacestart == string::npos)
                node="";
        else
                node = node.substr(spacestart, spaceend+1-spacestart);
}


/*
 * function to take an adjacency list as input and create a dot file.
 * A file with names of the attributes is to be given to see the names of the
 * attributes in the produced graph
 */
void write_dotout(map<int, list<int> >adjlist, int numofattr, vector<string>attrnames)
{
        map<int, list<int> >::iterator adjlist_it;
        list <int> templist;        
        list <int>::iterator list_it;
        int tempfirst, i;
        assert(numofattr == (int)attrnames.size()-1);
        ofstream outfile("tree.dot");
        outfile<<"digraph TAN{"<<endl;
        for (i = 1; i <= numofattr; i++)
                templist.push_back(i);
        adjlist[0]=templist;
        templist.clear();
        for (i = 0; i <=numofattr; i++)
        {
                outfile<<i<<" [label=\""<<attrnames[i]<<"\"];"<<endl;
        }
        for ( adjlist_it=adjlist.begin(); adjlist_it!= adjlist.end();adjlist_it++)
        {
                tempfirst = adjlist_it->first;
                templist = adjlist_it->second;
                for ( list_it = templist.begin();list_it != templist.end(); list_it++)
                {
                        outfile<<tempfirst<<"->"<<*list_it<<";"<<endl;
                }
                templist.clear();
        }
        outfile<<"}"<<endl;
        outfile.close();
}
