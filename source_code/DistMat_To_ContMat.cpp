#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
using namespace std;


//======================= I/O related ==========================//
//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//---------- input a string, output a vector -----//
int String_To_Vector(string &input,vector <double> &output)
{
	istringstream www(input);
	output.clear();
	int count=0;
	double value;
	for(;;)
	{
		if(! (www>>value) )break;
		output.push_back(value);
		count++;
	}
	return count;
}

//---------- input a matrix --------//
int Input_Matrix(string &input_file, vector < vector < double > > &output_mat)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		return -1;
	}
	//load
	int count=0;
	int colnum;
	int colcur;
	int first=1;
	output_mat.clear();
	vector <double> tmp_rec;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		colcur=String_To_Vector(buf,tmp_rec);
		if(first==1)
		{
			first=0;
			colnum=colcur;
		}
		else
		{
			if(colcur!=colnum)
			{
				fprintf(stderr,"current column number %d not equal to the first column number %d \n",
					colcur,colnum);
				return -1;
			}
		}
		output_mat.push_back(tmp_rec);
		count++;
	}
	//return
	return count;
}

//------------ filter DistMat -----------//
void DistMat_To_ContMat(vector < vector < double > > &input_mat, double thres,
	vector < vector < double > > &output_mat)
{
	//init
	output_mat=input_mat;
	//get start
	int i,j;
	int start=0;
	for(i=0;i<(int)input_mat.size();i++)
		for(j=0;j<(int)input_mat[i].size();j++)
		{
			//-> '-1' check
			if(input_mat[i][j]<0)
			{
				output_mat[i][j]=0;
				continue;
			}
			//-> contact matrix
			if(input_mat[i][j]<thres)output_mat[i][j]=1;
			else output_mat[i][j]=0;
		}
}

//------------ main -------------//
int main(int argc, char** argv)
{
	//------- DistMat_To_ContMat -----//
	{
		if(argc<3)
		{
			fprintf(stderr,"DistMat_To_ContMat <input_symm_distmat> <dist_thres>  \n");
			exit(-1);
		}
		string input_symm_distmat=argv[1];
		double dist_thres=atof(argv[2]);
		//load matrix
		vector < vector < double > > input_mat;
		int retv=Input_Matrix(input_symm_distmat, input_mat);
		if(retv<=0)exit(-1);
		//process
		vector < vector < double > > output_mat;
		DistMat_To_ContMat(input_mat, dist_thres, output_mat);
		//output
		for(int i=0;i<(int)input_mat.size();i++)
		{
			for(int j=0;j<(int)input_mat[i].size();j++)printf("%d ",(int)(output_mat[i][j]));
			printf("\n");
		}
		//exit
		exit(0);
	}
}
