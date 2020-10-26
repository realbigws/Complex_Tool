#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <getopt.h>
#include "Confo_Back.h"
#include "Confo_Beta.h"
#include "Confo_Lett.h"
#include "Acc_Surface.h"
#include "PDB_Chain_Fold.h"
#include "Mol_File.h"
using namespace std;


//-------- utility ------//
/*
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
*/


//========================= output ".ami", ".sse", ".cle", ".acc", ".feature" ============================//start

//------ output protein feature files -------//2015_02_20//
//--> output AMI_SSE_CLE
void Output_Protein_Features_AMI_SSE_CLE(
	string &outroot,string &outname,int moln,PDB_Residue *pdb,
	string &amino, vector <int> &mapping_data)
{
	//init
	int i;
	FILE *fpp;
	string file;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
	int retv;
	//calc
	retv=chain_fold.calculate_CLE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]CLE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	retv=chain_fold.calculate_SSE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]SSE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	//------ output AMI,CLE,SSE ------//
	string CLE,SSE;
	CLE=chain_fold.get_CLE();
	SSE=chain_fold.get_SSE();
	//output AMI
	file="";
	file=file+outroot+"/"+outname+".ami";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		fprintf(fpp,"%s\n",amino.c_str());
		fclose(fpp);
	}
	//output CLE
	file="";
	file=file+outroot+"/"+outname+".cle";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		for(i=0;i<(int)mapping_data.size();i++)
		{
			int index=mapping_data[i];
			if(index<0)fprintf(fpp,"R");
			else fprintf(fpp,"%c",CLE[index]);
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
	//output SSE
	file="";
	file=file+outroot+"/"+outname+".sse";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		for(i=0;i<(int)mapping_data.size();i++)
		{
			int index=mapping_data[i];
			if(index<0)fprintf(fpp,"L");
			else fprintf(fpp,"%c",SSE[index]);
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
	//output 0/1
	file="";
	file=file+outroot+"/"+outname+".miss";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		for(i=0;i<(int)mapping_data.size();i++)
		{
			int index=mapping_data[i];
			if(index<0)fprintf(fpp,"1");
			else fprintf(fpp,"0");
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
}

//--> output ACC and ACC_Value
void Output_Protein_Features_ACC(
	string &outroot,string &outname,Acc_Surface *acc_surface,
	int moln,PDB_Residue *pdb,XYZ **mcc,int *mcc_side,char *ami,char *acc,
	vector <int> &mapping_data)
{
	//init
	int i;
	FILE *fpp;
	string file;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
	//------ output ACC_Code and ACC_Value ------//
	for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
	acc_surface->AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);
	//output ACC_Code
	file="";
	file=file+outroot+"/"+outname+".acc";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		for(i=0;i<(int)mapping_data.size();i++)
		{
			int index=mapping_data[i];
			if(index<0)fprintf(fpp,"E");
			else fprintf(fpp,"%c",acc[index]);
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
	//output ACC_Value
	file="";
	file=file+outroot+"/"+outname+".acc_value";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		for(i=0;i<(int)mapping_data.size();i++)
		{
			int minus=-1;
			int index=mapping_data[i];
			if(index<0)fprintf(fpp,"%3d %3d\n",minus,minus);
			else fprintf(fpp,"%3d %3d\n",acc_surface->AC_normal[index],acc_surface->AC_output[index]);
		}
		fclose(fpp);
	}
}

//------ output protein feature in one file -------//2015_02_20//
// file format should be consistent with TPL file
/*
//////////// Features
  Num Res  Missing   SSE    CLE   ACC   pACC  CNa CNb   Xca       Yca       Zca       Xcb       Ycb       Zcb
   1   E      0       L      1     2     87   2   1     19.400     4.600    31.600    17.889     4.542    31.872
   2   N      0       E      5     2     48   4   3     20.500     7.600    33.700    21.326     8.630    32.887
   3   I      0       E      5     1     18   3   5     18.900     9.100    36.800    18.496     8.208    37.931
   4   E      0       E      5     2     53   3   2     19.900    12.700    37.700    19.630    13.781    36.734
   5   V      0       E      5     0      0   6   9     20.200    13.600    41.400    20.814    12.502    42.276
   6   H      0       E      5     1     32   6   3     20.700    17.200    42.800    19.590    18.263    42.565
   7   M      0       E      5     0      0   9   8     22.700    17.900    45.900    24.187    17.551    45.967
   8   L      0       E      5     1     12   5   4     21.100    20.900    47.700    19.558    20.886    47.424
   9   N      0       E      5     1     33   4   3     21.400    23.000    50.900    22.004    24.410    50.855
*/
//-------- ACC<->Int -------//
int ACC_To_Int(char c)
{
	switch(c)
	{
		case 'B': return 0;
		case 'M': return 1;
		case 'E': return 2;
		default: return 1;
	}
}
//----- output protein features -----//
void Output_Protein_Features(
	string &outroot,string &outname,Acc_Surface *acc_surface,XYZ *mol,XYZ *mcb,
	int moln,PDB_Residue *pdb,XYZ **mcc,int *mcc_side,char *ami,char *acc,
	string &amino, vector <int> &mapping_data)
{
	//init
	int i,j;
	int retv;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);

	//------ calculate AMI,CLE,SSE ------//
	retv=chain_fold.calculate_CLE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]CLE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	retv=chain_fold.calculate_SSE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]SSE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	string CLE,SSE;
	CLE=chain_fold.get_CLE();
	SSE=chain_fold.get_SSE();

	//------ calculate ACC_Code and ACC_Value ------//
	for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
	acc_surface->AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);

	//------ calculate contact number for CA and CB ---//
	double DIST_CUTOFF = 64; //-> 8.0A (note that in older version of TPL, the CA/CB contact cutoff is 7.0A
	vector <int> cn_ca(moln,0);
	vector <int> cn_cb(moln,0);
	for(i=0;i<moln;i++)
	{
		for(j=0;j<moln;j++)
		{
			double dist_ca=mol[i].distance_square(mol[j]);
			double dist_cb=mcb[i].distance_square(mcb[j]);
			//-> normal condition
			if( abs(i-j)>=4 )
			{
				if(dist_ca <= DIST_CUTOFF)cn_ca[i]++;
				if(dist_cb <= DIST_CUTOFF)cn_cb[i]++;
			}
			//-> chain broken
			else if( abs(i-j)>=1 )
			{
				//check PDB_residue
				string str1,str2;
				pdb[i].get_PDB_residue_number(str1);
				pdb[j].get_PDB_residue_number(str2);
				str1=str1.substr(1,4);
				str2=str2.substr(1,4);
				int pos1=atoi(str1.c_str());
				int pos2=atoi(str2.c_str());
				if( abs(pos1-pos2)>=4 )
				{
					if(dist_ca <= DIST_CUTOFF)cn_ca[i]++;
					if(dist_cb <= DIST_CUTOFF)cn_cb[i]++;
				}
			}
		}
	}

/*
	//------ calculate core region ------//
	int MinHelixLen = 4;
	int MinBetaLen = 3;
	int MinContact = 1;
	vector <int> Core(moln,1);
	for(i=0;i<moln;i++)
	{
		//calculate core
		Core[i] = 1;
		if(SSE[i]!='H' && SSE[i]!='E') continue;
		Core[i] = 2;
		int sslen = 1;
		for(j=i-1;j>0;j--)
		{
			if(SSE[j]!=SSE[i]) break;
			sslen++;
		}
		for(j=i+1;j<moln;j++)
		{
			if(SSE[j]!=SSE[i]) break;
			sslen++;
		}
		if(SSE[i]=='H' && sslen >= MinHelixLen && cn_ca[i]>MinContact)Core[i] = 5;
		if(SSE[i]=='E' && sslen >= MinBetaLen  && cn_ca[i]>MinContact)Core[i] = 5;
	}
*/

	//------ output feature files -------//
	string file=outroot+"/"+outname+".feature";
	FILE *fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,"  Num Res  Missing   SSE    CLE   ACC   pACC  CNa CNb   Xca       Yca       Zca       Xcb       Ycb       Zcb\n");
		for(i=0;i<(int)mapping_data.size();i++)
		{
			int index=mapping_data[i];
			if(index<0)fprintf(fpp,"%4d   %c      1       L      R     2     45   0   0   \n",i+1,amino[i]);
			else
			{
				char c=ACC_To_Int(acc[index]);
				fprintf(fpp,"%4d   %c      0       %c      %c     %1d    %3d  %2d  %2d   %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
					i+1,amino[i],SSE[index],CLE[index],c,acc_surface->AC_normal[index],
					cn_ca[index],cn_cb[index],mol[index].X,mol[index].Y,mol[index].Z,mcb[index].X,mcb[index].Y,mcb[index].Z);
			}
		}
	}
	fclose(fpp);
}

//--------- main process ----------//
void Feat_Main_Process(vector <PDB_Residue> &pdbmol,
	string &output,int OutFifi)
{
	//init check
	if(OutFifi==0)return;
	
	int i,j;
	//---- get maximal length -----//
	int totlen=pdbmol.size();
	//---- create data structure --//
	PDB_Residue *pdb=new PDB_Residue[totlen];
	XYZ *mol=new XYZ[totlen];   //CA
	XYZ *mcb=new XYZ[totlen];   //CB
	XYZ **mbb;                  //BackBone (N,CA,C,O,CB)
	XYZ **mcc;                  //BackBone (N,CA,C,O,CB,...,)
	NewArray2D(&mbb,totlen,5);
	NewArray2D(&mcc,totlen,15);
	int *mcc_side=new int[totlen];
	char *ami=new char[totlen+1];
	char *cle=new char[totlen+1];
	char *acc=new char[totlen+1];

	//----- get NO_MISS structure and mapping data -----------//
	string amino;
	vector <int> mapping_data(totlen,-1);
	int moln=0;
	for(i=0;i<totlen;i++)
	{
		//init
		PDB_Residue PDB=pdbmol[i];
		char amino_char=PDB.get_AA();
		amino.push_back(amino_char);
		string pdbind_;
		PDB.get_PDB_residue_number(pdbind_);
		int miss_or_not;
		if(pdbind_[5]=='!')miss_or_not=1;
		else miss_or_not=0;
		//record NO_MISS structure
		if(miss_or_not==0)
		{
			mapping_data[i]=moln;
			pdb[moln]=PDB;
			ami[moln]=amino_char;
			PDB.get_backbone_atom(1,mol[moln]);
			moln++;
		}
	}
	ami[moln]='\0';

	//----- reconstruct CB and HeavyAtom -----//
	Confo_Lett confo_lett;
	Confo_Beta *confo_beta=new Confo_Beta(totlen);
	Confo_Back *confo_back=new Confo_Back(totlen);
	Acc_Surface *acc_surface=new Acc_Surface(totlen);

	//[1]recon
	confo_lett.btb_ori(0,0,0,moln,mol,cle);
	confo_back->Recon_Back_Main(mol,cle,moln,mbb);      //given CA, recon BackBone (N,CA,C,O,CB)
	confo_beta->Recon_Beta_21(mol,mcb,moln,ami,cle);    //given CA, recon CB
	//[2]assign
	for(i=0;i<moln;i++)
	{
		//CA missing!!
		if(pdb[i].get_backbone_part_index(1)==0)
		{
			//backbone (N,CA,C,O)
			for(j=0;j<4;j++)pdb[i].set_backbone_atom(j,mbb[i][j]);
			//sidechain (CB)
			if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
			else pdb[i].set_sidechain_atom(0,mcb[i]);
		}
		else
		{
			//backbone (N,CA,C,O)
			for(j=0;j<4;j++)
			{
				if(pdb[i].get_backbone_part_index(j)==0)
				{
					pdb[i].set_backbone_atom(j,mbb[i][j]);
				}
			}
			//sidechain (CB)
			if(pdb[i].get_sidechain_part_index(0)==0)
			{
				if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
				else pdb[i].set_sidechain_atom(0,mcb[i]);
			}
		}
	}
	//get CB
	for(i=0;i<moln;i++) pdb[i].get_sidechain_atom( "CB ",mcb[i] );

	//----- calculate ALL features ------//
	//output others
	if(OutFifi!=0)
	{
		//-> get name and root
		string outroot,outname;
		getBaseName(output,outname,'/','.');
		getRootName(output,outroot,'/');
		//-> output files
		if(OutFifi==1 || OutFifi==3 || OutFifi==5 || OutFifi==7)
			Output_Protein_Features_AMI_SSE_CLE(outroot,outname,moln,pdb,amino,mapping_data);
		if(OutFifi==2 || OutFifi==3 || OutFifi==6 || OutFifi==7)
			Output_Protein_Features_ACC(outroot,outname,acc_surface,moln,pdb,mcc,mcc_side,ami,acc,mapping_data);
		if(OutFifi==4 || OutFifi==5 || OutFifi==6 || OutFifi==7)
			Output_Protein_Features(outroot,outname,acc_surface,mol,mcb,moln,pdb,mcc,mcc_side,ami,acc,amino,mapping_data);
	}

	//terminal
	delete [] pdb;
	delete [] ami;
	delete [] cle;
	delete [] acc;
	delete [] mol;
	delete [] mcb;
	DeleteArray2D(&mbb,totlen);
	DeleteArray2D(&mcc,totlen);
	delete [] mcc_side;
	delete confo_beta;
	delete confo_back;
	delete acc_surface;
}


//========================= output ".ami", ".sse", ".cle", ".acc", ".feature" ============================//over



//---------- dynamic programming ----------//
int WWW_Advance_Align_Dyna_Prog_Double(int n1,int n2,const vector<double> &score,
								   double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
								   double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
								   vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;

	//create D and M
	vector <int> D[3];      // the path (directions) matrix
	vector <double> M[3];   // the current scores (values) matrix
	//resize(m,n)
	for (i = 0; i < 3; ++i) 
	{
		D[i].resize(m*n);
		M[i].resize(m*n);
	}
	//init()
	double WS_MIN=-1000000;
	D[_S_][0*DP_maximal+ 0] = -1;
	D[_H_][0*DP_maximal+ 0] = -1;
	D[_V_][0*DP_maximal+ 0] = -1;
	M[_S_][0*DP_maximal+ 0] = 0;
	M[_H_][0*DP_maximal+ 0] = WS_MIN;
	M[_V_][0*DP_maximal+ 0] = WS_MIN;
	for (i = 1; i < m; i++) 
	{
		D[_S_][i*DP_maximal+ 0] = _V_;
		D[_H_][i*DP_maximal+ 0] = _V_;
		D[_V_][i*DP_maximal+ 0] = _V_;
		M[_S_][i*DP_maximal+ 0] = WS_MIN;
		M[_H_][i*DP_maximal+ 0] = WS_MIN;
		M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		D[_S_][0*DP_maximal+ j] = _H_;
		D[_H_][0*DP_maximal+ j] = _H_;
		D[_V_][0*DP_maximal+ j] = _H_;
		M[_S_][0*DP_maximal+ j] = WS_MIN;
		M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][0*DP_maximal+ j] = WS_MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		for (j = 1; j < n; j++) 
		{
			//condition upper
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			v1 = M[_V_][(i-1)*DP_maximal+ j] + gap_ext;
			v2 = M[_S_][(i-1)*DP_maximal+ j] + gap_open;
			v3 = M[_H_][(i-1)*DP_maximal+ j] + gap_open;
			M[_V_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_V_][i*DP_maximal+ j] == v1) D[_V_][i*DP_maximal+ j] = _V_;
			else if(M[_V_][i*DP_maximal+ j] == v2) D[_V_][i*DP_maximal+ j] = _S_;
			else D[_V_][i*DP_maximal+ j] = _H_;
			//condition left
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			v1 = M[_H_][i*DP_maximal+ j-1] + gap_ext;
			v2 = M[_S_][i*DP_maximal+ j-1] + gap_open;
			v3 = M[_V_][i*DP_maximal+ j-1] + gap_open;
			M[_H_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_H_][i*DP_maximal+ j] == v1) D[_H_][i*DP_maximal+ j] = _H_;
			else if(M[_H_][i*DP_maximal+ j] == v2) D[_H_][i*DP_maximal+ j] = _S_;
			else D[_H_][i*DP_maximal+ j] = _V_;
			//condition diag
			dist = score.at((i-1)*IN_maximal+ j-1);  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = M[_V_][(i-1)*DP_maximal+ j-1] + dist;
			v2 = M[_H_][(i-1)*DP_maximal+ j-1] + dist;
			v3 = M[_S_][(i-1)*DP_maximal+ j-1] + dist;
			M[_S_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_S_][i*DP_maximal+ j] == v3) D[_S_][i*DP_maximal+ j] = _S_;
			else if (M[_S_][i*DP_maximal+ j] == v1) D[_S_][i*DP_maximal+ j] = _V_;
			else D[_S_][i*DP_maximal+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=M[_V_][i*DP_maximal+ j];
	v2=M[_H_][i*DP_maximal+ j];
	v3=M[_S_][i*DP_maximal+ j];
	double maximal = std::max(v1, std::max(v2, v3));
	int k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	alignment.clear();
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	for(;;)
	{
		if(i==0||j==0)break;
		pre_case=D[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				alignment.push_back(pair<int,int>(i,j)); 
				i--;
				j--;
				++matches;
				break;
			case _V_:
				alignment.push_back(pair<int,int>(i,-j)); 
				i--;
				break;
			case _H_:
				alignment.push_back(pair<int,int>(-i,j)); 
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< D[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) alignment.push_back(pair<int,int>(-i,j)),j--;
	while (i> 0) alignment.push_back(pair<int,int>(i,0)), i--;
	reverse(alignment.begin(), alignment.end());
	ali_sco=maximal;
	return matches;
}

//---------- old mapping ------------//
int analyse_seperate(int i1,char c1,int i2,char c2,int &sep) 
{
	if(c1==' ')
	{
		if(c2==' ')sep=i2-i1;
		else
		{
			if(i1==i2)
			{
				sep=c2-'A';
				return -1;
			}
			else
			{
				sep=0;
				return -1;
			}
		}
	}
	else
	{
		if(c2!=' ')
		{
			if(c1==c2)sep=i2-i1;
			else
			{
				if(i1==i2)sep=abs(c2-c1);
				else
				{
					sep=0;
					return -1;
				}
			}
		}
		else
		{
			sep=0;
			return -1;
		}
	}
	if(sep<0)return -1;
	if(i1*i2<0)
	{
		sep--;
		return 2;  // represent:"0"	
	}
	else return 1;  // normal_return
}
int process_oriami_record(char *seqres,char *ami_,int *int_,char *ins_,char *tag_,
	vector <int> &mapping,string &out1,string &out2)
{
	int i,j;
	int head=0;
	int len;
	int totnum;
	int ii,jj;
	int seperate;
	int ret_val;
	int ws_rec_num;
	int n1,n2;

	//--[0]check
	len=(int)strlen(seqres);
	totnum=(int)strlen(ami_);

	//--[1]dynamic_programming	
	n1=len;    //SEQRES
	n2=totnum; //ATOM
	vector <double> score;
	score.resize(len*totnum);
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			if(seqres[i]==ami_[j+head])score.at(i*n2+j)=10;
			else
			{
				if(seqres[i]=='X'||seqres[i]=='Z'||seqres[i]=='.')score.at(i*n2+j)=0;
				else if(ami_[j+head]=='X'||ami_[j+head]=='Z'||ami_[j+head]=='.')score.at(i*n2+j)=0;
				else score.at(i*n2+j)=-15;
			}
		}
	}
	double sco;
	int matchs;
	vector<pair<int,int> > alignment;
	matchs=WWW_Advance_Align_Dyna_Prog_Double(n1,n2,score,-11,-1,-110,-10,-22,-22,-110,-110,alignment,sco);
	int lcmp=(int)alignment.size();

	//--Output_DP_Result
	vector <int> wali1;
	vector <int> wali2;
	wali1.resize(n1);
	wali2.resize(n2);
	for(j=0;j<n1;j++)wali1[j]=-1;
	for(j=0;j<n2;j++)wali2[j]=-1;

	//record_neo
	i=0;
	ii=0;
	int first,second; //first SEQUES, second ATOM
	int IsInsert=1;
	for(j=0;j<lcmp;j++)
	{
		first=alignment.at(j).first;
		second=alignment.at(j).second;
		if(first<=0)
		{
			if(second>0)
			{
				IsInsert=-1;
				tag_[i+head]='i';
				i++;
			}
		}
		else
		{
			if(second>0)
			{
				wali1.at(ii)=i;  //seqres_ami
				wali2.at(i)=ii;  //atom_ami
				if(seqres[ii]!=ami_[i+head])tag_[i+head]='/';
				i++;
			}
			ii++;
		}
	}
	//bad process
	if(IsInsert==-1)return IsInsert;
	if(matchs<2)return -1;

	//head_tail_tag	
	ii=wali2.at(totnum-1);
	if(ii!=-1)for(i=ii+1;i<len;i++)wali1.at(i)=-2; //tail_tag
	ii=wali2.at(0);
	if(ii!=-1)for(i=0;i<ii;i++)wali1.at(i)=-2; //head_tag

	//analyse_main_backword
	ws_rec_num=0;
	for(i=totnum-1;i>=1;i--)
	{
		if(tag_[i+head]=='i')continue; //__Found_Insert__//__080326__//
		ret_val=analyse_seperate(int_[i-1+head],ins_[i-1+head],int_[i+head],ins_[i+head],seperate);
		ii=wali2.at(i)-seperate; //expected_position
		jj=wali2.at(i-1);        //current_position
		if(ii!=jj && ret_val!=-1 && ii>=0 && ii<len) //error
		{
			if(ws_rec_num>=8)tag_[i+head]*=-1;  //solid!!
			ws_rec_num=0;
			for(j=0;j<ret_val;j++)
			{
				if(ii-j<0)break;
				if(wali1.at(ii-j)==-1)
				{
					if(ami_[i-1+head]==seqres[ii-j])
					{
						if(tag_[i-1+head]=='/')tag_[i-1+head]=' ';
					}
					else if(tag_[i-1+head]!='/')continue;
					if(jj>=0 && jj<n1)wali1.at(jj)=-1;
					wali1.at(ii-j)=i-1;
					wali2.at(i-1)=ii-j;				
					if(seperate==1)tag_[i-1+head]*=-1;  //solid!!
					break;
				}
			}
		}
		else ws_rec_num++;
	}

	//analyse_main_forward
	for(i=0;i<totnum-1;i++)
	{
		if(tag_[i+head]=='i')continue; //__Found_Insert__//__080326__//
		ret_val=analyse_seperate(int_[i+head],ins_[i+head],int_[i+1+head],ins_[i+1+head],seperate);
		ii=wali2.at(i)+seperate; //expected_position
		jj=wali2.at(i+1);        //current_position
		if(ii!=jj && ret_val!=-1 && ii>=0 && ii<len) //error
		{
			if(seperate!=1 && tag_[i+1+head]<0)continue;
			for(j=0;j<ret_val;j++)
			{
				if(ii+j>=len)break;
				if(wali1.at(ii+j)==-1)
				{
					if(ami_[i+1+head]==seqres[ii+j])
					{
						if(tag_[i+1+head]=='/')tag_[i+1+head]=' ';
					}
					else if(tag_[i+1+head]!='/')continue;
					if(jj>=0 && jj<n1)wali1.at(jj)=-1;
					wali1.at(ii+j)=i+1;
					wali2.at(i+1)=ii+j;
					break;
				}
			}
		}
	}

	//[final correction]
	int cur;
	//head_correct
	cur=0;
	ii=wali2.at(cur);     //current
	jj=wali2.at(cur+1)-1; //mapping
	if(ii!=jj && jj>=0 && jj<len)
	{
		if(wali1.at(jj)==-1)
		{
			if(ami_[cur+head]==seqres[jj])
			{
				if(tag_[cur+head]=='/')tag_[cur+head]=' ';
				wali1.at(ii)=-1;
				wali1.at(jj)=cur;
				wali2.at(cur)=jj;
			}
		}
	}
	//tail_correct
	cur=n2-1;
	ii=wali2.at(cur);     //current
	jj=wali2.at(cur-1)+1; //mapping
	if(ii!=jj && jj>=0 && jj<len)
	{
		if(wali1.at(jj)==-1)
		{
			if(ami_[cur+head]==seqres[jj])
			{
				if(tag_[cur+head]=='/')tag_[cur+head]=' ';
				wali1.at(ii)=-1;
				wali1.at(jj)=cur;
				wali2.at(cur)=jj;
			}
		}
	}

	//---- output mapping ------//
	//extract
	mapping.resize(len);
	for(i=0;i<len;i++)mapping[i]=-1; //default: NO
	int lali=0;
	out1="";
	out2="";
	for(i=0;i<lcmp;i++)
	{
		first=alignment[i].first;
		second=alignment[i].second;
		if(first<=0)continue;
		ii=first-1;
		out1+=seqres[ii];
		if(wali1.at(ii)<=-1)
		{
			out2+='-';
			continue;
		}
		jj=wali1.at(ii);
		out2+=ami_[jj];
		//get mapping
		mapping[ii]=jj;
		lali++;
	}

	//return
	return lali;
}

//--------- Extract Mapping --------//
void Extract_Mapping_XYZ(vector <int> &mapping, vector <XYZ> &input, vector <XYZ> &output)
{
	int i;
	int size=(int)mapping.size();
	output.resize(size);
	XYZ xyz;
	for(i=0;i<size;i++)
	{
		if(mapping[i]<0)
		{
			xyz=-99999.0;
			output[i]=xyz;
		}
		else
		{
			output[i]=input.at(mapping[i]);
		}
	}
}
void Extract_Mapping_PDB(string &seqres,vector <int> &mapping, 
	vector <PDB_Residue> &input, vector <PDB_Residue> &output,char Chain_ID)
{
	int i;
	int size=(int)mapping.size();
	output.resize(size);
	PDB_Residue residue;
	char tmp_rec[7];
	for(i=0;i<size;i++)
	{
		if(mapping[i]<0)
		{
			//get pdbind
			sprintf(tmp_rec,"%c%4d!",Chain_ID,i+1); //-> we start from zero !!
			//assign missing
			XYZ xyz=-99999.0;
			residue.PDB_residue_backbone_initialize('G');
			residue.set_backbone_atom( 1,xyz );
			//assign residue name
			residue.set_AA(seqres[i]);
			string tmp_rec_str=tmp_rec;
			residue.set_PDB_residue_number(tmp_rec_str);
			output[i]=residue;
		}
		else
		{
			//get pdbind
			sprintf(tmp_rec,"%c%4d ",Chain_ID,i+1); //-> we start from zero !!
			//record
			residue=input.at(mapping[i]);
			string tmp_rec_str=tmp_rec;
			residue.set_PDB_residue_number(tmp_rec_str);
			output[i]=residue;
		}
	}
}


//------------------ min distance square for PDB version ------------//
//--> collect all atoms from PDB_Residue
void Cllect_AllAtoms_from_PDB(PDB_Residue &PDB, vector <XYZ> &output_xyz)
{
	int k;
	int number;
	XYZ xyz;
	char amino=PDB.get_AA();
	output_xyz.clear();
	//hydrogen_out
	number=PDB.hydro_rec.size();
	for(k=0;k<number;k++)
	{
		//get_hydrogen
		xyz=PDB.hydro_rec[k];
		output_xyz.push_back(xyz);
	}
	//backbone_out
	number=PDB.get_backbone_totnum();  //this value should be 4
	for(k=0;k<number;k++)
	{
		//get_backbone
		if(PDB.get_backbone_part_index(k)==0)continue;
		PDB.get_backbone_atom(k,xyz);
		output_xyz.push_back(xyz);
	}
	//sidechain_out
	if(amino=='G')return;
	number=PDB.get_sidechain_totnum();
	for(k=0;k<number;k++)
	{
		//get_sidechain
		if(PDB.get_sidechain_part_index(k)==0)continue;
		PDB.get_sidechain_atom(k,xyz);
		output_xyz.push_back(xyz);
	}
}
//---> main calculate
double PDB_Distance_Square(PDB_Residue &r1, PDB_Residue &r2)
{
	//-> collect all atoms from PDB_Residue
	vector <XYZ> out_xyz_1;
	vector <XYZ> out_xyz_2;
	Cllect_AllAtoms_from_PDB(r1, out_xyz_1);
	Cllect_AllAtoms_from_PDB(r2, out_xyz_2);
	//-> calculate minimal distance square
	int i,j;
	int n1=(int)out_xyz_1.size();
	int n2=(int)out_xyz_2.size();
	double min_dist2=999999;
	double dist2;
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			dist2=out_xyz_1[i].distance_square(out_xyz_2[j]);
			if(dist2<min_dist2)min_dist2=dist2;
		}
	}
	//return
	return min_dist2;
}


//--------- output inter/intra molecule distance map --------------//
void Output_XYZ_Miss_contact(FILE *fp,
	vector <XYZ> & chain1, vector <XYZ> & chain2)
{
	int i,j;
	int ll=(int)chain1.size();
	int pl=(int)chain2.size();
	//distance output
	for(i=0;i<ll;i++)
	{
		double distance;
		for(j=0;j<pl;j++)
		{
			//calc distance
			if(chain1[i].X==-99999.0 || chain2[j].X==-99999.0)
			{
				distance=-1.0;
			}
			else
			{
				distance=chain1[i].distance_square(chain2[j]);
				distance=sqrt(distance);
			}
			//output
			fprintf(fp,"%8.2f\t",distance);
		}
		fprintf(fp,"\n");
	}
}
void Output_PDB_Miss_contact(FILE *fp,
	vector <PDB_Residue> & chain1, vector <PDB_Residue> & chain2)
{
	int i,j;
	int ll=(int)chain1.size();
	int pl=(int)chain2.size();
	//distance output
	for(i=0;i<ll;i++)
	{
		double distance;
		for(j=0;j<pl;j++)
		{
			//calc distance
			XYZ xyz1,xyz2;
			chain1[i].get_backbone_atom(1,xyz1);
			chain2[j].get_backbone_atom(1,xyz2);
			if(xyz1.X==-99999.0 || xyz2.X==-99999.0)
			{
				distance=-1.0;
			}
			else
			{
				distance=PDB_Distance_Square(chain1[i],chain2[j]);
				distance=sqrt(distance);
			}
			//output
			fprintf(fp,"%8.2f\t",distance);
		}
		fprintf(fp,"\n");
	}
}


//--------- output PDB format with MISS --------------//__2014_09_15__//
int Output_PDB_Miss(FILE *fp,vector <PDB_Residue> &mol,char Chain_ID, int Hydro) //output full-atom PDB file
{
	int i,k;
	string TER="TER                                                                             ";
	string buf="              ";
	int number;
	char amino;
	const char *dummyaa;
	const char *atomname;
	string pdbind_;
	string pdbind;
	double x,y,z;
	double rfactor,temperature;
	int numb;
	char Ini_Chain;
	char Real_Chain;
	PDB_Residue PDB;
	XYZ xyz;
	int cur_pos;
	char cur_chain;
	char nxt_chain;
	//ws_new//__110430__//
	char output[100];
	int miss_or_not;

	//judge
	int len=mol.size();
	if(len<=0)return -1;
	Ini_Chain=Chain_ID;
	if(Ini_Chain=='!'||Ini_Chain==-1||Ini_Chain=='_')Ini_Chain='_'; //put the chain from pdbind//__110230__//

	//------ get HEAD and TAIL missing regions --------//__180827__//
	vector <int> miss_region(len,0);
	for(i=0;i<len;i++)
	{
		//init
		PDB=mol[i];
		PDB.get_PDB_residue_number(pdbind_);
		if(pdbind_[5]=='!')miss_or_not=1;
		else miss_or_not=0;
		//record miss
		if(miss_or_not==1)miss_region[i]=1;
		else break;
	}
	for(i=len-1;i>=0;i--)
	{
		//init
		PDB=mol[i];
		PDB.get_PDB_residue_number(pdbind_);
		if(pdbind_[5]=='!')miss_or_not=1;
		else miss_or_not=0;
		//record miss
		if(miss_or_not==1)miss_region[i]=-1;
		else break;
	}
	string HEAD="HEAD";
	string TAIL="TAIL";
	string BODY="MISS";

	//process
	cur_pos=0;
	for(i=0;i<len;i++)
	{
		//init
		PDB=mol[i];
		amino=PDB.get_AA();
		dummyaa=One2Three_III(amino);
		PDB.get_PDB_residue_number(pdbind_);
		pdbind=pdbind_.substr(1,4);
		if(Ini_Chain=='_')Real_Chain=pdbind_[0];
		else Real_Chain=Ini_Chain;
		cur_chain=pdbind_[0];
		if(pdbind_[5]=='!')miss_or_not=1;
		else miss_or_not=0;
		//output miss
		if(miss_or_not==1)
		{
			string miss=BODY;
			if(miss_region[i]==1)miss=HEAD;
			else if(miss_region[i]==-1)miss=TAIL;
			sprintf(output,"%s             %3s %c%4s       0.000   0.000   0.000  0.00  0.00%s\n",
				miss.c_str(),dummyaa,Real_Chain,pdbind.c_str(),buf.c_str());
			fprintf(fp,"%s",output);
			continue;
		}
		//backbone_out
		number=PDB.get_backbone_totnum();  //this value should be 4
		for(k=0;k<number;k++)
		{
			//get_backbone
			if(PDB.get_backbone_part_index(k)==0)
			{
				continue;
			}
			atomname=backbone_atom_name_decode(k);
			PDB.get_backbone_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			{
				sprintf(output,"ATOM  %5d %4s %3s %c%4s    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			//real_output//__110430__//
			{
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}
		//sidechain_out
		number=PDB.get_sidechain_totnum();
		for(k=0;k<number;k++)
		{
			if(amino=='G')continue;
			//get_sidechain
			if(PDB.get_sidechain_part_index(k)==0)
			{
				continue;
			}
			atomname=sidechain_atom_name_decode(k,amino);
			PDB.get_sidechain_atom(k,xyz, numb, rfactor, temperature);
			x=xyz.X;
			y=xyz.Y;
			z=xyz.Z;
			//output
			{
				sprintf(output,"ATOM  %5d %4s %3s %c%4s    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					numb,atomname,dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
			}
			//real_output//__110430__//
			{
				output[77]=output[13];
			}
			fprintf(fp,"%s",output);
		}

		//--- OutHydr ---//start
		if(Hydro==1)
		{
			vector <XYZ> hydro_out=PDB.hydro_rec;
			int hydr_num=(int)hydro_out.size()<999?(int)hydro_out.size():999;
			for(k=0;k<hydr_num;k++)
			{
				//get hydrogen name
				char hname[5];
				sprintf(hname,"H%-3d",k+1);
				hname[4]='\0';
				string hname_str=hname;
				//get coordinate
				xyz=hydro_out[k];
				x=xyz.X;
				y=xyz.Y;
				z=xyz.Z;
				//get others
				numb=0;
				rfactor=1;
				temperature=20;
				//output
				{
					sprintf(output,"ATOM  %5d %4s %3s %c%4s    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
						numb,hname_str.c_str(),dummyaa,Real_Chain,pdbind.c_str(),x,y,z,rfactor,temperature,buf.c_str());
				}
				//real_output//__110430__//
				{
					output[77]=output[12];
				}
				fprintf(fp,"%s",output);
			}
		}
		//--- OutHydr ---//end

		cur_pos++;
		//terminal_process
		if(Chain_ID=='!'||Chain_ID==-1)
		{
			if(i<len-1)
			{
				PDB_Residue pdb;
				string pdb_ind;
				pdb=mol[i+1];
				pdb.get_PDB_residue_number(pdb_ind);
				nxt_chain=pdb_ind[0];
				if(cur_chain!=nxt_chain)
				{
					fprintf(fp,"%s\n",TER.c_str());
					cur_pos=0;
				}
			}
		}
	}
	//terminal
	return 1;
}



//================== All_Chain_Process =============//__110530__//
//-> calculate Sidechain Center (SC)
void Calculate_Sidechain_Center(PDB_Residue &PDB, XYZ &sc_center)
{
	int k;
	int number;
	XYZ xyz;
	char amino=PDB.get_AA();
	//check 'G'
	if(amino=='G')
	{
		PDB.get_backbone_atom(1,xyz);
		sc_center=xyz;
		return;
	}

	//sidechain_out
	int count=0;
	sc_center=0;
	number=PDB.get_sidechain_totnum();
	for(k=0;k<number;k++)
	{
		//get_sidechain
		if(PDB.get_sidechain_part_index(k)==0)continue;
		PDB.get_sidechain_atom(k,xyz);
		sc_center+=xyz;
		count++;
	}
	if(count==0)
	{
		PDB.get_backbone_atom(1,xyz);
		sc_center=xyz;
	}
	else sc_center/=count;
}

//-> Compare chain1 and chain2
int Compare_Chain1_and_Chain2_XYZ(vector <XYZ> & chain1, vector <XYZ> & chain2,
	double r_cut,vector < pair <int,int> > & result)
{
	int i,j;
	int ll=(int)chain1.size();
	int pl=(int)chain2.size();
	//distance check
	double dist2;
	double thres2=r_cut*r_cut;
	result.clear();
	int count=0;
	for(i=0;i<ll;i++)
	{
		//check miss
		if(chain1[i].X==-99999.0)continue;
		for(j=0;j<pl;j++)
		{
			//check miss
			if(chain2[j].X==-99999.0)continue;
			//calc distance
			dist2=chain1[i].distance_square(chain2[j]);
			if(dist2<thres2)
			{
				result.push_back(pair<int,int>(i,j));
				count++;
			}
		}
	}
	return count;
}
int Compare_Chain1_and_Chain2_PDB(vector <PDB_Residue> & chain1, vector <PDB_Residue> & chain2,
	double r_cut,vector < pair <int,int> > & result)
{
	int i,j;
	int ll=(int)chain1.size();
	int pl=(int)chain2.size();
	//distance check
	XYZ xyz;
	double dist2;
	double thres2=r_cut*r_cut;
	result.clear();
	int count=0;
	for(i=0;i<ll;i++)
	{
		//check miss
		chain1[i].get_backbone_atom(1,xyz);
		if(xyz.X==-99999.0)continue;
		for(j=0;j<pl;j++)
		{
			//check miss
			chain2[j].get_backbone_atom(1,xyz);
			if(xyz.X==-99999.0)continue;
			//calc distance
			dist2=PDB_Distance_Square(chain1[i],chain2[j]);
			if(dist2<thres2)
			{
				result.push_back(pair<int,int>(i,j));
				count++;
			}
		}
	}
	return count;
}

//-> Compare inner chain
int Compare_Inner_Chain_XYZ(vector <XYZ> & chain, double r_cut, int resi_gap,
	vector < pair <int,int> > & result)
{
	int i,j;
	int ll=(int)chain.size();
	int pl=(int)chain.size();
	//distance check
	double dist2;
	double thres2=r_cut*r_cut;
	result.clear();
	int count=0;
	for(i=0;i<ll;i++)
	{
		//check miss
		if(chain[i].X==-99999.0)continue;
		for(j=0;j<pl;j++)
		{
			//check miss
			if(chain[j].X==-99999.0)continue;
			//check separation
			if(abs(i-j)<resi_gap)continue;
			//calc distance
			dist2=chain[i].distance_square(chain[j]);
			if(dist2<thres2)
			{
				result.push_back(pair<int,int>(i,j));
				count++;
			}
		}
	}
	return count;
}
int Compare_Inner_Chain_PDB(vector <PDB_Residue> & chain, double r_cut, int resi_gap,
	vector < pair <int,int> > & result)
{
	int i,j;
	int ll=(int)chain.size();
	int pl=(int)chain.size();
	//distance check
	XYZ xyz;
	double dist2;
	double thres2=r_cut*r_cut;
	result.clear();
	int count=0;
	for(i=0;i<ll;i++)
	{
		//check miss
		chain[i].get_backbone_atom(1,xyz);
		if(xyz.X==-99999.0)continue;
		for(j=0;j<pl;j++)
		{
			//check miss
			chain[j].get_backbone_atom(1,xyz);
			if(xyz.X==-99999.0)continue;
			//check separation
			if(abs(i-j)<resi_gap)continue;
			//calc distance
			dist2=PDB_Distance_Square(chain[i],chain[j]);
			if(dist2<thres2)
			{
				result.push_back(pair<int,int>(i,j));
				count++;
			}
		}
	}
	return count;
}


//=================== Get Missing SEQRES !!!!! ===================//__2014_09_15__//
//[note]: from now on, we'll consider missing part in the PDB structure, using the following format !!
/*
MISS             MET     1       0.000   0.000   0.000  0.00  0.00
MISS             ASP     2       0.000   0.000   0.000  0.00  0.00
MISS             ASN     3       0.000   0.000   0.000  0.00  0.00
MISS             ASN     4       0.000   0.000   0.000  0.00  0.00
MISS             PRO     5       0.000   0.000   0.000  0.00  0.00
MISS             ASN     6       0.000   0.000   0.000  0.00  0.00
MISS             ILE     7       0.000   0.000   0.000  0.00  0.00
MISS             ASN     8       0.000   0.000   0.000  0.00  0.00
MISS             GLU     9       0.000   0.000   0.000  0.00  0.00
MISS             CYS    10       0.000   0.000   0.000  0.00  0.00
MISS             ILE    11       0.000   0.000   0.000  0.00  0.00
MISS             PRO    12       0.000   0.000   0.000  0.00  0.00
MISS             TYR    13       0.000   0.000   0.000  0.00  0.00
MISS             ASN    14       0.000   0.000   0.000  0.00  0.00
MISS             CYS    15       0.000   0.000   0.000  0.00  0.00
ATOM      1  N   LEU    16      24.358  -6.198  13.806  1.00 83.08           N
ATOM      2  CA  LEU    16      24.356  -7.651  13.730  1.00 83.08           C
ATOM      3  CB  LEU    16      24.725  -8.270  15.107  1.00 83.08           C
ATOM      4  CG  LEU    16      26.212  -7.994  15.562  1.00 83.08           C
ATOM      5  CD1 LEU    16      26.516  -6.517  15.765  1.00 83.08           C
ATOM      6  CD2 LEU    16      26.614  -8.822  16.776  1.00 83.08           C
*/

//---- part 0. data structure -----//
const int BufLen = 100000;
string AA3Coding[26]={"ALA","XXX","CYS","ASP","GLU","PHE","GLY","HIS","ILE","XXX","LYS","LEU","MET","ASN","XXX","PRO","GLN","ARG","SER","THR","XXX","VAL","TRP","XXX","TYR","XXX"};

//---- part 1. read SEQRES from original PDB file -------//
void Read_SEQRES_From_Orig_PDB(string &pdbfile,char chainID,string &out_str)
{
	ifstream pdbfile_in(pdbfile.c_str());
	char buf[BufLen];
	int start=0;
	int seqres=0;
	string sequence="";

	//---- read original PDB ----//
	while(pdbfile_in.getline(buf,BufLen))
	{
		// break condition
		if(start == 1 &&strncmp(buf,"TER",3) == 0)break;
		if(start == 1 &&strncmp(buf,"END",3) == 0)break;
		
		// extract sequence
		if(strncmp(buf,"SEQRES",6) == 0 && buf[11] == chainID)
		{
			seqres = 1;
			istringstream seq_in (buf+19);
			string aa3;
			while(seq_in >> aa3)
			{
				int aa1 = -1;
				for(int i=0;i<26;i++)
				if(aa3 == AA3Coding[i])aa1 = i;
				if(aa1<0) aa1=12;
				if(aa3=="XXX") aa1=12;
				sequence = sequence + (char)('A' + aa1);
			}
		}
	}
	pdbfile_in.close();

	//----- final process ----//
	out_str=sequence;
}


//============================ main process ======================//
void Main_Process(string &complex_file,string &contact_out, 
	int CAorCB, double radius, int resi_thres,int mode, int Hydro,int OutFifi)
{
	//get name
	string pdb_nam;
	getBaseName(complex_file,pdb_nam,'/','.');
	//class
	Mol_File mol_input;
	mol_input.MODRES=1;
	mol_input.HYDR_MODE=Hydro;   // [1: output hydrogen][0: don't out]

	//chain load
	vector <PDB_Chain> chains;
	int retv=mol_input.PDB_read_pdb_file(complex_file,chains);
	if(retv<0)
	{
		fprintf(stderr,"complex_file %s failed to open or format bad. \n",complex_file.c_str());
		return; //failed
	}

	//seqres load
	int i,j,k;
	int chain_size=(int)chains.size();
	PDB_Chain pdb_chain;
	int moln;
	char chain;
	PDB_Residue residue;
	//-> seqres data structure
	vector <string> seqres_rec;
	seqres_rec.resize(chain_size);
	for(i=0;i<chain_size;i++)
	{
		//-> pick chain
		pdb_chain=chains[i];
		moln=pdb_chain.get_length();
		chain=pdb_chain.get_chain_id();
		//-> load seqres
		string out_str;
		Read_SEQRES_From_Orig_PDB(complex_file,chain,out_str);
		if(out_str=="")
		{
			string amino_sequence;
			amino_sequence.resize(moln);
			for(k=0;k<moln;k++)
			{
				pdb_chain.get_residue(k,residue);
				amino_sequence[k]=residue.get_AA();
			}
			out_str=amino_sequence;
		}
		//-> record seqres
		seqres_rec[i]=out_str;
	}

	//do mapping
	vector<vector <XYZ> > xyz_rec;
	vector<vector <PDB_Residue> > pdb_rec;
	xyz_rec.resize(chain_size);
	pdb_rec.resize(chain_size);
	for(i=0;i<chain_size;i++)
	{
		//-> pick chain
		pdb_chain=chains[i];
		moln=pdb_chain.get_length();
		chain=pdb_chain.get_chain_id();
		//-> generate seqres data
		int seqres_len=seqres_rec[i].length();
		char *seqres=new char[seqres_len+1];
		for(k=0;k<seqres_len;k++)seqres[k]=seqres_rec[i][k];
		seqres[seqres_len]='\0';
		//-> generate atom data
		char *ami_=new char[moln+1];
		int *int_=new int[moln];
		char *ins_=new char[moln+1];
		char *tag_=new char[moln+1];
		vector <XYZ> in_xyz;
		vector <PDB_Residue> in_pdb;
		in_xyz.resize(moln);
		in_pdb.resize(moln);
		for(k=0;k<moln;k++)
		{
			pdb_chain.get_residue(k,residue);
			string pdbind_;
			residue.get_PDB_residue_number(pdbind_);
			string int_str=pdbind_.substr(1,4);
			//-> get sequence data
			ami_[k]=residue.get_AA();
			int_[k]=atoi(int_str.c_str());
			ins_[k]=pdbind_[5];
			tag_[k]=' ';
			//-> get xyz
			if(CAorCB==1)       //-> get CA
			{
				residue.get_backbone_atom(1,in_xyz[k]);
			}
			else if(CAorCB==0)  //-> get CB
			{
				if(residue.get_sidechain_part_index(0)==0)
				{
					residue.get_backbone_atom(1,in_xyz[k]);
				}
				else
				{
					residue.get_sidechain_atom(0,in_xyz[k]);
				}
			}
			else if(CAorCB==-1) //-> Sidechain Center
			{
				Calculate_Sidechain_Center(residue, in_xyz[k]);
			}
			else            //-> dummy record
			{
				in_xyz[k]=0;
			}
			//-> get pdb_residue
			in_pdb[k]=residue;
		}
		ami_[moln]='\0';
		ins_[moln]='\0';
		tag_[moln]='\0';
		//-> generate mapping
		vector <int> mapping;
		string out1,out2;
		process_oriami_record(seqres,ami_,int_,ins_,tag_,mapping,out1,out2);
		//-> calclate result
		vector <XYZ> out_xyz;
		Extract_Mapping_XYZ(mapping,in_xyz,out_xyz);
		vector <PDB_Residue> out_pdb;
		Extract_Mapping_PDB(seqres_rec[i],mapping,in_pdb,out_pdb,chain);
		//-> record
		xyz_rec[i]=out_xyz;
		pdb_rec[i]=out_pdb;
		//-> delete
		delete [] seqres;
		delete [] ami_;
		delete [] int_;
		delete [] ins_;
		delete [] tag_;
	}


	//---------- output -----------//
	FILE *fp=fopen(contact_out.c_str(),"wb");
	//calculate inner contact
	for(i=0;i<chain_size;i++)
	{
		//-> pick chain
		pdb_chain=chains[i];
		moln=pdb_chain.get_length();
		int relmoln=(int)xyz_rec[i].size();
		chain=pdb_chain.get_chain_id();
		//-> get data
		vector < pair <int,int> >  contact_inner;
		if(CAorCB>=-1)Compare_Inner_Chain_XYZ(xyz_rec[i], radius, resi_thres, contact_inner);
		else Compare_Inner_Chain_PDB(pdb_rec[i], radius, resi_thres, contact_inner);
		//-> output data
		fprintf(fp,"# %s%c %4d %4d : intra chain contact\n",pdb_nam.c_str(),chain,moln,relmoln);
		int size=(int)contact_inner.size();
		int ii,jj;
		for(k=0;k<size;k++)
		{
			ii=contact_inner[k].first;
			jj=contact_inner[k].second;
			fprintf(fp,"%4d %4d\n",ii+1,jj+1);
		}
		fprintf(fp,"\n");
		//-> output additional contact data
		if(mode==1 || mode==2)
		{
			//-> output contact/distance map in matrix format
			string out_file=pdb_nam+chain+".con";
			FILE *ff=fopen(out_file.c_str(),"wb");
			if(CAorCB>=-1)Output_XYZ_Miss_contact(ff,xyz_rec[i],xyz_rec[i]);
			else Output_PDB_Miss_contact(ff,pdb_rec[i],pdb_rec[i]);
			fclose(ff);
		}
	}
	//calculate inter contact
	for(i=0;i<chain_size;i++)
	{
		//-> pick chain 1
		pdb_chain=chains[i];
		int moln1=pdb_chain.get_length();
		int relmoln1=(int)xyz_rec[i].size();
		char chain1=pdb_chain.get_chain_id();
		for(j=i+1;j<chain_size;j++)
		{
			//-> pick chain 2
			pdb_chain=chains[j];
			int moln2=pdb_chain.get_length();
			int relmoln2=(int)xyz_rec[j].size();
			char chain2=pdb_chain.get_chain_id();
			//-> get data
			vector < pair <int,int> >  contact_inter;
			if(CAorCB>=-1)Compare_Chain1_and_Chain2_XYZ(xyz_rec[i], xyz_rec[j],radius,contact_inter);
			else Compare_Chain1_and_Chain2_PDB(pdb_rec[i], pdb_rec[j],radius,contact_inter);
			//-> output data
			fprintf(fp,"# %s%c %4d %4d <-> %s%c %4d %4d : inter chain contact\n",
				pdb_nam.c_str(),chain1,moln1,relmoln1,pdb_nam.c_str(),chain2,moln2,relmoln2);
			int size=(int)contact_inter.size();
			int ii,jj;
			for(k=0;k<size;k++)
			{
				ii=contact_inter[k].first;
				jj=contact_inter[k].second;
				fprintf(fp,"%4d %4d\n",ii+1,jj+1);
			}
			fprintf(fp,"\n");
			//-> output additional contact data
			if(mode==1 || mode==2)
			{
				//-> output contact/distance map in matrix format
				{
					string out_file=pdb_nam+"_"+chain1+"_"+chain2+".con";
					FILE *ff=fopen(out_file.c_str(),"wb");
					if(CAorCB>=-1)Output_XYZ_Miss_contact(ff,xyz_rec[i],xyz_rec[j]);
					else Output_PDB_Miss_contact(ff,pdb_rec[i],pdb_rec[j]);
					fclose(ff);
				}
				//-> output contact/distance map in matrix format
				{
					string out_file=pdb_nam+"_"+chain2+"_"+chain1+".con";
					FILE *ff=fopen(out_file.c_str(),"wb");
					if(CAorCB>=-1)Output_XYZ_Miss_contact(ff,xyz_rec[j],xyz_rec[i]);
					else Output_PDB_Miss_contact(ff,pdb_rec[j],pdb_rec[i]);
					fclose(ff);
				}
			}
		}
	}
	fclose(fp);

	//============ output each structure in 'MISS' foramt ===========//
	if(mode==0 || mode==2)
	{
		for(i=0;i<chain_size;i++)
		{
			//-> pick chain
			pdb_chain=chains[i];
			moln=pdb_chain.get_length();
			chain=pdb_chain.get_chain_id();
			//-> output PDB in 'MISS' format
			string out_file=pdb_nam+chain+".pdb";
			fp=fopen(out_file.c_str(),"wb");
			Output_PDB_Miss(fp,pdb_rec[i],chain,Hydro);
			fclose(fp);
			//-> output other features
			string feat_file=pdb_nam+chain;
			Feat_Main_Process(pdb_rec[i],feat_file,OutFifi);
		}
	}
}

//---------- usage ---------//
void Usage() 
{
	fprintf(stderr,"Version: 1.05 [2017-12-09] \n");
	fprintf(stderr,"Complex_Tool -i complex_file -o contact_out [-m mode] \n");
	fprintf(stderr,"             [-c CAorCB] [-r radius] [-R resi_thres] [-h Hydro] [-F feature] \n");
	fprintf(stderr,"Usage : \n\n");
	fprintf(stderr,"-i complex_file :      Input complex official PDB file in 4-char. \n");
	fprintf(stderr,"-o contact_out  :      Output file for intra/inter molecular contact. \n");
	fprintf(stderr,"-m mode         :     [0] output '.pdb' file for each chain only. \n");
	fprintf(stderr,"                       1  output '.con' matrix for intra/inter chain \n");
	fprintf(stderr,"                      -1  not output any additional file; 2 output all \n");
	fprintf(stderr,"-c CAorCB       :      use CA[1],CB[0],SC[-1],All[-2] to calculate contact. \n");
	fprintf(stderr,"                       by default, we use 0 to use CB-CB to define contact. \n");
	fprintf(stderr,"-r radius       :      within the radius for contact. (set to 8.0 by default)\n");
	fprintf(stderr,"-R resi_thres   :      inner contact residue separation (set to 6 by default)\n");
	fprintf(stderr,"-h Hydro        :     [0] DO NOT consider hydrogen atoms; 1 for considering. \n");
	fprintf(stderr,"-F feature      :      1 for AMI,CLE,SSE; 2 for ACC; 4 for FEAT. \n");
	fprintf(stderr,"                       these features could be combined (set to 0 by default)\n");
	fprintf(stderr,"-----------------------------------------------------------------------------\n");
	fprintf(stderr,"[note]: in this version, we consider MISS residue with repect to SEQRES \n");
	fprintf(stderr,"        we also renumber the residue number sequentially according to SEQRES \n");
}

//============== main ===============//
int main(int argc, char** argv)
{
	//---- Complex_To_ContactMap ---//process the homomer's contact for all other chains with the first one
	{
		if(argc<3)
		{
			Usage();
			exit(-1);
		}
		//required input
		string complex_file="";
		string contact_out="";
		//optional input
		int mode=0;
		int CAorCB=0;
		double radius=8.0;
		int resi_thres=6;
		int Hydro=0;
		int OutFifi=0;

		//command-line arguments process
		extern char* optarg;
		char c = 0;
		while ((c = getopt(argc, argv, "i:o:m:c:r:R:h:F:")) != EOF) 
		{
			switch (c) 
			{
				case 'i':
					complex_file = optarg;
					break;
				case 'o':
					contact_out = optarg;
					break;
				case 'm':
					mode = atoi(optarg);
					break;
				case 'c':
					CAorCB = atoi(optarg);
					break;
				case 'r':
					radius = atof(optarg);
					break;
				case 'R':
					resi_thres = atoi(optarg);
					break;
				case 'h':
					Hydro = atoi(optarg);
					break;
				case 'F':
					OutFifi = atoi(optarg);
					break;
				default:
					Usage();
					exit(-1);
			}
		}

		//---- check required input ----//
		if(complex_file=="")
		{
			fprintf(stderr,"complex_file is null !!\n");
			exit(-1);
		}
		if(contact_out=="")
		{
			fprintf(stderr,"contact_out is null !!\n");
			exit(-1);
		}

		//calculate
		Main_Process(complex_file,contact_out,CAorCB,radius,resi_thres,mode,Hydro,OutFifi);
		//exit
		exit(0);
	}
}
