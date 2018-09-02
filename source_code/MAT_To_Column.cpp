#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include "Fast_Sort.h"
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

//---------- input a symmetric matrix --------//
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
int Input_Symmetric_Matrix(string &input_file, vector < vector < double > > &output_mat)
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
	//final check
	if(count!=colnum)
	{
		fprintf(stderr,"row number %d not equal to column number %d \n",
			count,colnum);
		return -1;
	}
	//symmetric check
	double epsilu=1.e-9;
	for(int i=0;i<count;i++)
	{
		for(int j=i+1;j<count;j++)
		{
			double value=fabs(output_mat[i][j]-output_mat[j][i]);
			if(value>epsilu)
			{
				fprintf(stderr,"WARNING: symmetric condition broken at [%d,%d] -> %lf not equal to %lf \n",
					i,j,output_mat[i][j],output_mat[j][i]);
			}
		}
	}
	//return
	return count;
}
int Input_NonSymmetric_Matrix(string &input_file, 
	vector < vector < double > > &output_mat,int &len1,int &len2)
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
	len1=count;
	len2=colnum;
	return count;
}

//---- calculate matrix Z_score ------//
void Calc_Mat_Zscore(vector < vector < double > > &input_mat, double &mean, double &vari)
{
	int i,j;
	int size=(int)input_mat.size();
	int count;
	//mean
	mean=0;
	count=0;
	for(i=0;i<size;i++)
	{
		for(j=i+6;j<size;j++)
		{
			mean+=input_mat[i][j];
			count++;
		}
	}
	mean/=count;
	//vari
	vari=0;
	count=0;
	for(i=0;i<size;i++)
	{
		for(j=i+6;j<size;j++)
		{
			vari+=(input_mat[i][j]-mean)*(input_mat[i][j]-mean);
			count++;
		}
	}		
	vari=1.0*sqrt(1.0*vari/count);
}

//----- Calculate vector Z-score ----//
void Calc_Vec_Zscore(vector <double> &in,double &mean,double &vari, vector <double> &out)
{
	//init
	mean=0;
	vari=1;
	//proc
	int i;
	int size=(int)in.size();
	if(size==0)return;
	//-> calculate mean
	for(i=0;i<size;i++)mean+=in[i];
	mean/=size;
	//-> calculate vari
	vari=0;
	for(i=0;i<size;i++)vari+=(in[i]-mean)*(in[i]-mean);
	vari=1.0*sqrt(vari/size);
	//-> calculate Z-score
	out.resize(size);
	for(i=0;i<size;i++)out[i]=1.0*(in[i]-mean)/vari;
}

//--------- symmetric matrix to column stype --------//
//-> note: posi starts from 0 !!
int SymmMat_To_Column_Interface(vector < vector < double > > &input_mat,
	vector<pair<int, int> > &posi, vector <double> &value,
	int len1,int len2)
{
	int i,j;
	int size=(int)input_mat.size();
	//length check
	if(size!=(len1+len2))
	{
		fprintf(stderr,"size %d not equal to len1 %d + len2 %d \n",
			size,len1,len2);
		exit(-1);
	}
	//process
	posi.clear();
	value.clear();
	int count=0;
	for(i=0;i<len1;i++)
	{
		for(j=len1;j<size;j++)
		{
			posi.push_back(pair<int,int>(i,j-len1));
			value.push_back(0.5*(input_mat[i][j]+input_mat[j][i]));
			count++;
		}
	}
	//return
	return count;
}
int NonSymmMat_To_Column_Interface(
	vector < vector < double > > &input_mat,
	vector<pair<int, int> > &posi, vector <double> &value,
	int len1,int len2)
{
	int i,j;
	int size=(int)input_mat.size();
	//process
	posi.clear();
	value.clear();
	int count=0;
	for(i=0;i<len1;i++)
	{
		for(j=0;j<len2;j++)
		{
			posi.push_back(pair<int,int>(i,j));
			value.push_back(input_mat[i][j]);
			count++;
		}
	}
	//return
	return count;
}



//----------- vector sort ----------//
//-> UPorDOWN: UP for ascending order sort, DOWN for descending order sort
void Vector_Sort(vector <double> &input, vector <int> &order,int UPorDOWN=0)
{
	order.clear();
	//prepare Fast_Sort
	Fast_Sort <double> fast_sort;
	int size=(int)input.size();
	double *temp_score=new double[size];
	int *temp_index=new int[size];
	for(int i=0;i<size;i++)temp_score[i]=input[i];
	//Fast_Sort
	if(UPorDOWN==1)fast_sort.fast_sort_1up(temp_score,temp_index,size);
	else fast_sort.fast_sort_1(temp_score,temp_index,size);
	//output
	order.resize(size);
	for(int i=0;i<size;i++)order[i]=temp_index[i];
	//delete
	delete [] temp_score;
	delete [] temp_index;
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//------- Matrix_Sort -----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Version 1.02\n");
			fprintf(stderr,"MAT_To_Column <input_mat> <out_file> <sort_order> \n");
			fprintf(stderr,"[note]: sort matrix to column in DESCEND or ASCEND order \n");
			fprintf(stderr,"        sort_order should be [0]: DESCEND, or 1: ASCEND \n");
			exit(-1);
		}

		//----- part 0: input arguments -------//
		//-> input files
		string input_file=argv[1];
		//-> output files
		string out_file=argv[2];
		//-> parameters
		int sort_order=atoi(argv[3]);    //should be [0,1]

		//------ part 2: input matrix --------//
		vector < vector < double > > nonsymm_mat;
		int len1,len2;
		int retv=Input_NonSymmetric_Matrix(input_file,nonsymm_mat,len1,len2);
		if(retv<=0)exit(-1);
		vector<pair<int, int> > posi;
		vector <double> value;
		int totnum=NonSymmMat_To_Column_Interface(nonsymm_mat,posi,value,len1,len2);

		//------ part 3: main process --------//
		//calculate zscore
		vector <double> value_zsco;
		double mean,vari;
		Calc_Vec_Zscore(value, mean, vari,value_zsco);
		//process sort
		vector <int> order;
		Vector_Sort(value, order, sort_order);  //0 for DESCEND, 1 for ASCEND

		//------ part 4: build modeller script ------------//
		FILE *fp;
		fp=fopen(out_file.c_str(),"wb");
		for(int i=0;i<totnum;i++)
		{
			//-> get index
			int index=order[i];
			//-> get score
			int ii=posi[index].first;
			int jj=posi[index].second;
			double zsco=1.0*(value[index]-mean)/vari;
			//-> output
			fprintf(fp,"%5d %5d %.3f %.3f \n",
				ii+1,jj+1,zsco,value[index]);
		}
		fclose(fp);
		//exit
		exit(0);
	}
}
