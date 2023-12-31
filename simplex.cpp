/*
 * main.cpp
 *
 *  Created on: 21 Jul 2021
 *      Author: prakharraj1302
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
using namespace std;

class fraction
{
public:
	int numerator;
	int denominator;

	fraction();			// default constractor....check
	fraction(int, int); // default constractor with parameters ....check
	fraction AddedTo(fraction value) const;
	fraction MultipliedBy(fraction value) const; // binary observer type of operation
	fraction Subtract(fraction value) const;
	fraction DividedBy(fraction value) const;
	fraction isGreaterThan();
	fraction isEqualTo();
	fraction print() const;
};

fraction::fraction()
{
	numerator = 0;
	denominator = 1;
}
//**********************************
fraction::fraction(int newNumerator, int newDenominator)
{
	numerator = newNumerator;
	denominator = newDenominator;
}
//***************************
fraction fraction::AddedTo(fraction value) const
{
	fraction result;
	result.numerator = (numerator * value.denominator) + (value.numerator * denominator);
	result.denominator = denominator * value.denominator;
	return result;
}
//********************************
fraction fraction::Subtract(fraction value) const
{
	fraction result;
	result.numerator = (numerator * value.denominator) - (value.numerator * denominator);
	result.denominator = denominator * value.denominator;
	return result;
}
//**************************************
fraction fraction::MultipliedBy(fraction value) const
{
	fraction result;
	result.numerator = numerator * value.numerator;
	result.denominator = denominator * value.denominator;
	return result;
}
//************************************
fraction fraction::DividedBy(fraction value) const
{
	fraction result;
	result.numerator = numerator * value.denominator;
	result.denominator = denominator * value.numerator;
	return result;
}

// fraction fraction::print() const
// {
// int num=numerator;
// int den=denominator;
//  if(num>den)
//  {
//  for(int counter=2;counter<den;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//  else if(den>num)
//  {
//  for(int counter=2;counter<num;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//
//  cout<<num<<"/"<<den;
//  }

class frac
{

public:
	double num, den;
	frac()
	{
		cout << "frac cons called!!!!!!!!!!!!" << endl;

		num = 0;

		den = 1;
		cout << num << " " << den << endl;
	}

	friend class matrix;
};

class matrix : public frac
{

public:
	int rows, columns;

	float **p;
	frac **pf;

	void table(matrix &, matrix &, matrix &, int);

	matrix mx(matrix &, matrix &);
	matrix sm(matrix &, matrix &);
	matrix sb(matrix &, matrix &);

	void in();
	void inset(int, int);

	void disp();
	void dispstat();

	matrix trn();
	matrix inverse();
	void imat();
	void cost();
	;

	void initializer();
	void iterator(matrix &, matrix &, matrix &);
	int comp(matrix &);
	int minarray(matrix &);

	matrix(int = 1, int = 1);
	matrix(matrix &);
	friend class frac;
};

//----------------------------------------------------------------------------------------

matrix matrix::trn()
{
	cout << "trn fx called!!!" << endl;

	matrix tempmat(columns, rows);
	cout << columns << rows << endl;

	tempmat.dispstat();

	int n, m, i, j;
	n = rows;
	m = columns;
	cout << "!!!!!1";

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			cout << i << " " << j << " " << p[i][j] << " " << tempmat.p[i][j];

			tempmat.p[j][i] = p[i][j];
			tempmat.disp();
		}
	}
	return tempmat;
}

matrix matrix ::inverse()
{
	cout << "invrese fx called !!!" << endl;

	matrix B(rows, columns), C(columns, columns * 2);
	int row, column, step, i;
	double mult;

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			C.p[row][column] = p[row][column];

	for (i = 0; i < rows; i++)
		C.p[i][rows + i] = 1;

	for (step = 0; step < rows - 1; step++)
	{
		for (row = step + 1; row < rows; row++)
		{
			mult = C.p[row][step] / C.p[step][step];
			for (column = step; column < 2 * rows; column++)
				C.p[row][column] -= mult * C.p[step][column];
		}
	}

	for (step = 1; step <= rows - 1; step++)
	{
		for (row = rows - step - 1; row >= 0; row--)
		{
			mult = C.p[row][rows - step] / C.p[rows - step][rows - step];
			for (column = rows; column < 2 * rows; column++)
				C.p[row][column] -=
					mult * C.p[rows - step][column];
		}
	}

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			B.p[row][column] = C.p[row][rows + column] / C.p[row][row];
	return B;
}

matrix::matrix(int r, int c)
{

	// matrix  constructor

	cout << "mat constructor called !!!!!" << r << " " << c << endl;

	int i, j;
	rows = r;
	columns = c;

	p = new float *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		p[i] = new float[columns];
		// assert(p[i] != NULL);
	}

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			p[i][j] = 0;
	pf = NULL;
}

matrix::matrix(matrix &mf) : frac()
{

	// matrix  constructor

	cout << " ############  mat FRACTION  constructor called !!!!!" << endl;

	int i, j;
	rows = mf.rows;
	columns = mf.columns;

	pf = new frac *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		pf[i] = new frac[columns];
		// assert(p[i] != NULL);
	}

	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < columns; j++)
		{
			pf[i][j].num = mf.p[i][j];
		}
	}

	p = NULL;
}

void matrix ::imat()
{

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			if (i == j)
			{
				p[i][j] = 1;
			}
			else
			{
				p[i][j] = 0;
			}
		}
	}
}

void matrix ::cost()
{

	cout << "cost called !!" << endl;

	for (int i = 0; i < rows; i++)
	{
		p[i][1] = 0;
	}
}

matrix matrix ::mx(matrix &A, matrix &B)
{
	cout << "mx fx called!!!!\n";

	float sum;
	int i, j, k, n, m, l;
	n = A.rows;
	m = A.columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + A.p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix operator*(matrix &A, matrix &B)
{
	float sum;
	int i, j, k, n, m, l;
	n = A.rows;
	m = A.columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + A.p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix matrix ::sm(matrix &A, matrix &B)
{

	cout << "sm called!!!!!" << endl;

	int i, j, n, m;

	n = A.rows;
	m = A.columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = B.p[i][j] + A.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

matrix matrix ::sb(matrix &A, matrix &B)
{

	cout << "sb called!!!!!" << endl;

	int i, j, n, m;

	n = A.rows;
	m = A.columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = B.p[i][j] - A.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

void matrix ::in()
{
	cout << "matrix fx called!!!!\n";

	delete[] p;

	cout << "Enter the number of Rows and Columns: ";
	cin >> rows >> columns;
	p = new float *[rows];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new float[columns];
	}

	cout << "Enter the " << rows << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << columns << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::inset(int rset, int cset)
{
	cout << "matrix inset called!!!!\n";

	delete[] p;

	p = new float *[rset];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new float[cset];
	}

	cout << "Enter the " << rset << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << cset << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::disp()
{
	cout << "disp called!!!!!\n";

	for (int i = 0; i < rows; i++)
	{

		for (int j = 0; j < columns; j++)
		{
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
}
void matrix ::dispstat()
{
	cout << "dispstat called!!!!!\n";
	cout << rows << " " << columns << endl;
}

void matrix ::initializer()
{
	cout << "IN ITER !!!!" << endl;

	int nz, nv, ndesc;

	cout << "enter max mat" << endl;
	cout << "enter no of var" << endl;
	cin >> nz;
	matrix matz(1, nz);
	matz.inset(1, nz);
	matz.disp();

	cout << "no of desc" << endl;
	cin >> ndesc;
	matrix matv(ndesc, (nz + 1));
	matv.inset(ndesc, (nz + 1));
	matv.disp();

	matrix matB(ndesc, ndesc);
	matB.imat();
	matB.disp();

	iterator(matz, matv, matB);
}

void matrix ::iterator(matrix &z, matrix &v, matrix &B)
{
	cout << "IN ITERRATOR !!!!" << endl;

	matrix Xb;
	matrix DEL;
	matrix COST;
	matrix prc_1;
	matrix prc_2;
	matrix Aj;
	matrix Cj;
	matrix COST_trn;
	matrix temp_1;

	matrix B_inv;
	B_inv.inverse();
	B_inv.disp();

	matrix b(v.rows, 1);
	for (int i = 0; i < v.rows; i++)
	{
		b.p[i][0] = v.p[i][(v.columns - 1)];
	}
	b.disp();

	matrix cB(v.rows, 1);
	cB.cost();
	cB.disp();
	cB.trn();
	cB.disp();

	int opti_flag = 0;

	for (int j = 0; j < 1; j++)
	{

		cout << "ITERATION" << endl;

		Xb.mx(B_inv, b);
		prc_1.mx(cB, B_inv);
		prc_1.disp();
		prc_2.mx(prc_1, Aj);
		prc_2.disp();
		DEL.sb(prc_2, Cj);
		DEL.disp();
		opti_flag = temp_1.comp(DEL);
		cout << opti_flag << endl;
		DEL.disp();
	}
}

int matrix ::comp(matrix &delj)
{
	cout << "in comp !!!!" << endl;
	int minpos = 0;
	float minelem;
	int nflag = 0;

	for (int i = 0; i < delj.rows; i++)
	{
		if (delj.p[i][0] < 0)
		{
			nflag = 1;
			break;
		}
	}
	if (nflag == 1)
	{
		minpos = minarray(delj);
		cout << "index" << minpos << endl;
	}
	else
	{
		cout << "all +positive" << endl;
	}
	return minpos;
}

int matrix ::minarray(matrix &min)
{
	cout << "in minarray !!!!" << endl;

	int index = 0;
	int size = sizeof(min) / sizeof(min.p, min.p[0][0]);

	for (int i = 1; i < size; i++)
	{
		if (min.p[i] < min.p[index])
			index = i;
	}

	return index;
}

//----------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------

int main()
{
	cout << "start\n";

	matrix m1;
	m1.in();
	matrix m2(m1);

	m1.disp();

	cout << "main end\n";

	return 0;
}
//---------------------------------------------------------------------------------------
//_____________________________________V2____________________________________
/*
 * main.cpp
 *
 *  Created on: 21 Jul 2021
 *      Author: prakharraj1302
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
using namespace std;

// class fraction
//{
// public:
//
//	int numerator;
//	int denominator;
//
//
//
//
//
//	fraction (); //default constractor....check
//	fraction (int,int); //default constractor with parameters ....check
//	fraction AddedTo (fraction value) const;
//	fraction MultipliedBy (fraction value) const; //binary observer type of operation
//	fraction Subtract (fraction value) const;
//	fraction DividedBy (fraction value) const;
//	fraction  isGreaterThan ();
//	fraction isEqualTo ();
//	fraction print() const;
//
//
// };

// frac::fraction()
//{
//	numerator = 0;
//	denominator = 1;
// }
//**********************************
// fraction::fraction(int newNumerator, int newDenominator)
//{
//	numerator=newNumerator;
//	denominator=newDenominator;
// }
//***************************

////********************************

// fraction fraction::print() const
// {
// int num=numerator;
// int den=denominator;
//  if(num>den)
//  {
//  for(int counter=2;counter<den;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//  else if(den>num)
//  {
//  for(int counter=2;counter<num;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//
//  cout<<num<<"/"<<den;
//  }
//**************************************#####################################

class frac
{

public:
	double num, den;
	frac()
	{
		cout << "frac cons called!!!!!!!!!!!!" << endl;

		num = 0;

		den = 1;
		cout << num << " " << den << endl;
	}

	frac fadd(frac &);
	frac fsub(frac &);
	frac fmul(frac &);
	frac fdiv(frac &);
	frac fneg();
	void fdisp(frac &);
};

void frac::fdisp(frac &fval)
{
	cout << "[" << fval.num << " / " << fval.den << "]" << endl;
}
frac frac::fadd(frac &value)
{
	frac result;
	result.num = (num * value.den) + (value.num * den);
	result.den = den * value.den;
	return result;
}
frac frac::fsub(frac &value)
{
	frac result;
	result.num = (num * value.den) - (value.num * den);
	result.den = den * value.den;
	return result;
}
//**************************************
frac frac::fmul(frac &value)
{
	frac result;
	result.num = num * value.num;
	result.den = den * value.den;
	return result;
}
//************************************
frac frac::fdiv(frac &value)
{
	frac result;
	result.num = num * value.den;
	result.den = den * value.num;
	return result;
}
//************************************
frac frac ::fneg()
{
	num = -num;

	return (*this);
}

//**************************************#####################################

class matrix : public frac
{

public:
	int rows, columns;

	double **p;
	frac **pf;

	void table(matrix &, matrix &, matrix &, int);

	matrix mx(matrix &);
	matrix sm(matrix &, matrix &);
	matrix sb(matrix &);

	void in();
	void inset(int, int);

	matrix merger(matrix &, matrix &, matrix &, int);
	matrix A_indexer(int);
	int ratio(matrix &);

	void disp();
	void dispstat();

	void frac_disp();

	matrix trn();
	matrix inverse();
	void imat();
	void cost();
	;

	void initializer();
	void iterator(matrix &, matrix &, matrix &);
	int comp(matrix &);
	int minarray(matrix &);

	matrix extractor(matrix &);
	matrix matcopy(matrix &);
	matrix ftab_solver(matrix &, int);
	matrix ftab_updater(matrix &, int, int);
	matrix iter_result(matrix &);
	void del_j(matrix &);
	void cost_updt(matrix &);
	void Aj_updt(matrix &);
	void Cj_updt(matrix &);

	matrix(int = 1, int = 1);
	matrix(matrix &);
};
matrix global_V;
matrix global_B;
matrix global_Z;
matrix global_full_tab;
matrix global_A;

int global_flag;
int global_ratio;

void matrix ::frac_disp()
{
	cout << "frac_disp called!!!!!\n";

	for (int i = 0; i < rows; i++)
	{

		for (int j = 0; j < columns; j++)
		{
			// cout<<i<< " %  "<<j<<endl;
			cout << pf[i][j].num << "/" << pf[i][j].den << " ";
		}
		cout << endl;
	}
}

matrix matrix ::ftab_solver(matrix &ftab, int ent_row)
{
	cout << "ftab solver called !!!!!!111" << endl;

	//	matrix temp( (*this)  , 1);
	//	matrix mat( *this->rows , *this->columns);
	matrix mat(ftab);

	cout << "frac mat for table made $$$$$$$$$" << endl;

	for (int j = 0; j < ftab.rows; j++)
	{
		if (j == ent_row)
		{
			continue;
		}
		cout << j << endl;

		// double a ,b;
		// a = (pf[ent_row][columns - 1]).num / (pf[ent_row][columns - 1]).den;
		// b = (pf[j][columns - 1]).num / (pf[j][columns - 1]).den;

		// float x = temp.eq_solver(a , b);

		// cout<< a <<" ? "<< b <<" ? "  <<endl;

		mat = ftab_updater(ftab, ent_row, j);
	}

	return mat;
}
matrix matrix ::ftab_updater(matrix &tab, int p, int q)
{
	cout << "ftab updter called !!!!!!111" << endl;

	// double a = tab.p[p][columns - 1] ;
	// double b = tab.p[q][columns - 1];

	cout << "tab updater called !!!!!!!1!!!" << endl;
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	cout << " / " << p << " / " << q << endl;

	// cout<< "a= "<< a.num <<" "<<a.den<<" b = " <<b.num<<" "<<b.den <<endl;

	for (int i = 0; i < tab.columns; i++)
	{
		cout << " @@@@@@@@@@@@@@@@ TAB SOLVING @@@@@@@@@@@@@@@" << endl;

		cout << tab.pf[p][i].num << "/" << tab.pf[p][i].den << " $ " << tab.pf[q][i].num << "/" << tab.pf[q][i].den << endl;

		frac a, b, t1, t2;
		a = tab.pf[p][tab.columns - 1];
		b = tab.pf[q][tab.columns - 1];
		cout << "a= " << a.num << " " << a.den << " b = " << b.num << " " << b.den << endl
			 << "##" << endl;
		b = b.fneg();
		fdisp(b);
		t1 = b.fdiv(a);
		fdisp(t1);
		t2 = t1.fmul(tab.pf[p][i]);
		fdisp(t2);
		tab.pf[q][i] = (tab.pf[q][i]).fadd(t2);
		fdisp(tab.pf[q][i]);

		// tab.p[q][i]= tab.p[q][i] + ( (-b/a) * (tab.p[p][i])) ;			//dfrac

		cout << tab.pf[p][i].num << "/" << tab.pf[p][i].den << " $ " << tab.pf[q][i].num << "/" << tab.pf[q][i].den << endl;
	}
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

	for (int j = 0; j < tab.columns; j++)
	{
		frac a;
		a = tab.pf[p][tab.columns - 1];
		tab.pf[p][j] = (tab.pf[p][j]).fdiv(a);

		// tab.p[p][i] = tab.p[p][i] / tab.p[p][columns -1];				//dfrac
	}

	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

	cout << " ############SOLVED ###########" << endl;

	//	//extractor(tab);
	//	matrix dmat;
	//	dmat=extractor(tab);
	//	dmat.disp();
	iter_result(tab);
	return tab;
}
//----------------------------------------------------------------------------------------

matrix matrix ::extractor(matrix &fmat)
{

	cout << "in extractor ######################" << endl;
	matrix temp(fmat.rows, fmat.columns);
	int i, j;
	cout << fmat.rows << " " << fmat.columns << endl;

	for (i = 0; i < fmat.rows; i++)
	{

		for (j = 0; j < fmat.columns; j++)
		{

			temp.p[i][j] = ((fmat.pf[i][j].num) / (fmat.pf[i][j].den));
			cout << temp.p[i][j] << " ";
		}
		cout << endl;
	}
	return temp;
}
matrix matrix ::iter_result(matrix &fmat)
{
	cout << "in iterresult ######################" << endl;

	matrix Xb;
	cout << "NEW FEASIBLE SOLN Xb " << endl;
	for (int i = 0; i < fmat.rows; i++)
	{
		cout << fmat.pf[i][0].num << "/" << fmat.pf[i][0].den << endl;
	}

	matrix t1(fmat.rows, fmat.columns);
	matrix binv(global_V.rows, global_V.rows);
	t1 = extractor(fmat);
	t1.disp();

	cout << "BINVERSE EXTRACTION ###################################" << endl;
	for (int i = 0; i < fmat.rows; i++)
	{
		for (int j = 0; j < global_V.rows; j++)
		{
			cout << i << " & " << j << endl;
			(binv.p[i][j]) = (t1.p[i][j + 1]);
		}
	}

	binv.disp();
	matrix::del_j(binv);
}

void matrix ::del_j(matrix &Binv)
{
	cout << "in delj ######################" << endl;

	matrix prc_1;
	matrix prc_2;
	matrix DEL;
	matrix Aj(global_V.rows, global_V.columns - 1);
	matrix cost_trn_mat(global_V.rows, 1);
	matrix Cj(global_Z.rows, global_Z.columns);
	int op_flag = 0;

	cost_updt(cost_trn_mat);
	Aj_updt(Aj);
	Cj_updt(Cj);

	prc_1 = cost_trn_mat.mx(Binv);
	prc_1.disp();

	cout << "######" << endl;

	prc_2 = prc_1.mx(Aj);
	prc_2.disp();

	cout << "######" << endl;

	DEL = prc_2.sb(Cj);
	DEL.disp();
	op_flag = matrix::comp(DEL);
}

void matrix ::Aj_updt(matrix &Aj)
{
	cout << "in AJ updt ################" << endl;

	int fl = 0;
	for (int i = 0; i < global_A.rows; i++)
	{
		int c = 0;
		for (int j = 0; j < global_A.columns - 1; j++)
		{
			if (j == global_flag)
			{
				continue;
			}
			Aj.p[i][j] = global_A.p[i][j];
			c++;
			if (c == 3)
			{
				fl = 1;
				break;
			}
		}
		if (fl == 1)
		{
			break;
		}
	}
	Aj.disp();
}

void matrix ::Cj_updt(matrix &Cj)
{
	cout << "in z updt ################" << endl;

	cout << " in min config Z with added vars --- " << endl;

	matrix Z(global_Z.rows, global_Z.columns + global_V.rows);
	for (int i = 0; i < global_Z.columns; i++)
	{
		if (i == 3)
		{
			Z.p[0][i] = 0;
			continue;
		}
		Z.p[0][i] = global_Z.p[0][i];
	}
	Z.disp();

	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "Cj updation #########" << endl;

	int c = 0;
	for (int j = 0; j < Z.columns; j++)
	{
		if (j == global_flag)
		{
			continue;
		}
		Cj.p[0][j] = Z.p[0][j];
		c++;
		if (c == global_V.columns - 1)
		{
			break;
		}
	}
	cout << "##################" << endl;
}

void matrix ::cost_updt(matrix &costmat)
{

	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	cout << "cost_updt called ######" << endl;
	cout << " globalflag- " << global_flag << " globarat- " << global_ratio << endl;

	for (int i = 0; i < global_V.rows; i++)
	{
		if (i == global_ratio)
		{
			costmat.p[i][0] = global_Z.p[1][global_flag];
			continue;
		}
		costmat.p[i][0] = global_V.p[global_V.columns][i];
	}
	costmat.disp();
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
}

matrix matrix ::matcopy(matrix &mat)
{
	cout << "matcop[y called ######" << endl;

	matrix temp(mat);

	int i, j;
	int r, c;

	r = mat.rows;
	c = mat.columns;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			pf[i][j].num = mat.p[i][j];
		}
	}
	return temp;
}

matrix matrix::trn()
{
	cout << "trn fx called!!!" << endl;

	matrix tempmat(columns, rows);
	//	cout<<columns<<rows<<endl;

	// tempmat.dispstat();
	int n, m, i, j;
	n = rows;
	m = columns;
	//	cout<<"!!!!!1";

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			//	cout<<i<<" "<<j<<" "<<p[i][j]<<" "<<tempmat.p[i][j];

			tempmat.p[j][i] = p[i][j];
			tempmat.disp();
		}
	}
	return tempmat;
}

matrix matrix ::inverse()
{
	cout << "invrese fx called !!!" << endl;

	matrix B(rows, columns), C(columns, columns * 2);
	int row, column, step, i;
	double mult;

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			C.p[row][column] = p[row][column];

	for (i = 0; i < rows; i++)
		C.p[i][rows + i] = 1;

	for (step = 0; step < rows - 1; step++)
	{
		for (row = step + 1; row < rows; row++)
		{
			mult = C.p[row][step] / C.p[step][step];
			for (column = step; column < 2 * rows; column++)
				C.p[row][column] -= mult * C.p[step][column];
		}
	}

	for (step = 1; step <= rows - 1; step++)
	{
		for (row = rows - step - 1; row >= 0; row--)
		{
			mult = C.p[row][rows - step] / C.p[rows - step][rows - step];
			for (column = rows; column < 2 * rows; column++)
				C.p[row][column] -=
					mult * C.p[rows - step][column];
		}
	}

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			B.p[row][column] = C.p[row][rows + column] / C.p[row][row];
	return B;
}

matrix::matrix(int r, int c)
{

	// matrix  constructor

	cout << "mat constructor called !!!!!" << r << " " << c << endl;

	int i, j;
	rows = r;
	columns = c;

	p = new double *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		p[i] = new double[columns];
		// assert(p[i] != NULL);
	}

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			p[i][j] = 0;
	pf = NULL;
}

matrix::matrix(matrix &mf) : frac()
{

	// matrix  constructor

	cout << " ############  mat FRACTION  constructor called !!!!!" << endl;

	int i, j;
	rows = mf.rows;
	columns = mf.columns;

	pf = new frac *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		pf[i] = new frac[columns];
		// assert(p[i] != NULL);
	}
	cout << "frac alloc sucessful" << endl;

	if (mf.p == NULL)
	{

		for (i = 0; i < mf.rows; i++)
		{
			for (j = 0; j < mf.columns; j++)
			{
				cout << i << " % " << j << endl;
				pf[i][j] = mf.pf[i][j];
			}
		}
	}
	else
	{
		for (i = 0; i < mf.rows; i++)
		{
			for (j = 0; j < mf.columns; j++)
			{
				cout << i << " % " << j << endl;
				pf[i][j].num = mf.p[i][j];
			}
		}
	}
	p = NULL;
	cout << "mat FRACTION  constructor LEFT  " << endl;
}

void matrix ::imat()
{

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			if (i == j)
			{
				p[i][j] = 1;
			}
			else
			{
				p[i][j] = 0;
			}
		}
	}
}

void matrix ::cost()
{

	cout << "cost called !!" << endl;

	for (int i = 0; i < rows; i++)
	{
		p[i][1] = 0;
	}
}

matrix matrix ::mx(matrix &B)
{
	cout << "mx fx called!!!!\n";

	float sum;
	int i, j, k, n, m, l;
	n = rows;
	m = columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix operator*(matrix &A, matrix &B)
{
	double sum;
	int i, j, k, n, m, l;
	n = A.rows;
	m = A.columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + A.p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix matrix ::sm(matrix &A, matrix &B)
{

	cout << "sm called!!!!!" << endl;

	int i, j, n, m;

	n = A.rows;
	m = A.columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = B.p[i][j] + A.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

matrix matrix ::sb(matrix &B)
{

	cout << "sb called!!!!!" << endl;

	int i, j, n, m;

	n = rows;
	m = columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = p[i][j] - B.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

void matrix ::in()
{
	cout << "matrix fx called!!!!\n";

	delete[] p;

	cout << "Enter the number of Rows and Columns: ";
	cin >> rows >> columns;
	p = new double *[rows];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new double[columns];
	}

	cout << "Enter the " << rows << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << columns << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::inset(int rset, int cset)
{
	cout << "matrix inset called!!!!\n";

	delete[] p;

	p = new double *[rset];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new double[cset];
	}

	cout << "Enter the " << rset << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << cset << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::disp()
{
	cout << "disp called!!!!!\n";

	for (int i = 0; i < rows; i++)
	{

		for (int j = 0; j < columns; j++)
		{
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
}

void matrix ::dispstat()
{
	cout << "dispstat called!!!!!\n";
	cout << rows << " " << columns << endl;
}

void matrix ::initializer()
{
	cout << "IN ITER !!!!" << endl;

	int nz, nv, ndesc;

	cout << "enter max mat" << endl;
	cout << "enter no of var" << endl;
	cin >> nz;
	matrix matz(1, nz);
	matz.inset(1, nz);
	matz.disp();
	global_Z = matz;

	cout << "no of desc" << endl;
	cin >> ndesc;
	matrix matv(ndesc, (nz + 1));
	matv.inset(ndesc, (nz + 1));
	matv.disp();
	global_V = matv;

	matrix matB(ndesc, ndesc);
	matB.imat();
	matB.disp();
	global_B = matB;

	iterator(matz, matv, matB);
}

void matrix ::iterator(matrix &z, matrix &v, matrix &B)
{
	cout << "IN ITERRATOR !!!!" << endl;

	matrix Xb;
	matrix DEL;
	matrix COST;
	matrix prc_1;
	matrix prc_2;
	matrix prc_3;

	matrix Cj;
	matrix COST_trn;
	matrix temp_1;

	matrix Aj(v.rows, v.columns - 1);

	for (int i = 0; i < v.rows; i++)
	{
		for (int j = 0; j < v.columns - 1; j++)
		{
			cout << i << " " << j << endl;

			Aj.p[i][j] = v.p[i][j];
		}
	}
	Aj.disp();

	matrix B_inv;
	B_inv = B.inverse();
	B_inv.disp();

	matrix b(v.rows, 1);
	for (int i = 0; i < v.rows; i++)
	{
		b.p[i][0] = v.p[i][(v.columns - 1)];
	}
	b.disp();

	matrix cB(v.rows, 1);
	cB.cost();
	cB.disp();
	COST_trn = cB.trn();
	COST_trn.disp();

	Cj = z;
	Cj.disp();

	int opti_flag = 0;

	for (int j = 0; j < 1; j++)
	{

		cout << "ITERATION" << endl
			 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

		B_inv.disp();
		b.disp();
		COST_trn.disp();

		Xb = B_inv.mx(b);
		Xb.disp();
		cout << "######" << endl;

		prc_1 = COST_trn.mx(B_inv);
		prc_1.disp();

		cout << "######" << endl;

		prc_2 = prc_1.mx(Aj);
		prc_2.disp();

		cout << "######" << endl;
		DEL = prc_2.sb(Cj);
		DEL.disp();

		cout << "######" << endl;
		opti_flag = temp_1.comp(DEL);
		global_flag = opti_flag;

		cout << "opti_flag = " << opti_flag << endl;
		DEL.disp();

		cout << "######" << endl;
		matrix table(v.rows, B.columns + 2);

		table = prc_3.merger(Xb, B_inv, v, opti_flag);
		table.disp();

		matrix Xn(v.rows, 1);
		for (int k = 0; k < v.rows; k++)
		{
			Xn.p[k][0] = v.p[k][opti_flag];
		}
		matrix Atab(v.rows, v.columns - 1 + v.rows);

		Atab = v.A_indexer(1);

		cout << "######" << endl;

		int while_flag = 1;

		while (while_flag == 1)
		{

			cout << "in WHILE @@@@@@@@@@@@@@@@@@" << endl;

			int rat = Xb.ratio(Xn);
			global_ratio = rat;
			cout << "ratio is " << rat << endl;
			matrix ftab(table);
			ftab.frac_disp();

			ftab_solver(ftab, rat);

			cout << "######### END OF WHILE############" << endl;

			while_flag = 0;
		}
	}
}
int matrix ::ratio(matrix &X)
{

	cout << "in ratio !!!!!!!!!!1" << endl;

	matrix rat(X.rows, 1);

	for (int i = 0; i < X.rows; i++)
	{
		rat.p[i][0] = (p[i][0]) / (X.p[i][0]);
	}

	int index = 0;
	int size = sizeof(rat.p) / sizeof(rat.p[0]);

	for (int i = 1; i < size; i++)
	{
		if (rat.p[i] < rat.p[index])
			index = i;
	}
	return index;
}

matrix matrix ::merger(matrix &Xb, matrix &Bin, matrix &Xn, int n)
{
	cout << "merger called !!!!!!!!!!!!1" << endl;

	matrix merge_mat(Xn.rows, Bin.columns + 2);
	matrix full_tab(Xn.rows, Bin.columns + Xn.columns);

	for (int i = 0; i < Xb.rows; i++)
	{
		merge_mat.p[i][0] = Xb.p[i][0];
		for (int j = 0; j < Bin.columns; j++)
		{
			merge_mat.p[i][1 + j] = Bin.p[i][j];
		}
		for (int j = 0; j < Xn.columns; j++)
		{

			merge_mat.p[i][Bin.columns + 1] = Xn.p[i][j];
		}
	}
	global_full_tab = full_tab;

	for (int i = 0; i < Xb.rows; i++)
	{
		merge_mat.p[i][0] = Xb.p[i][0];
		for (int j = 0; j < Bin.columns; j++)
		{
			merge_mat.p[i][1 + j] = Bin.p[i][j];
		}
		merge_mat.p[i][Bin.columns + 1] = Xn.p[i][n];
	}
	return merge_mat;
}

matrix matrix ::A_indexer(int select)
{
	cout << "in INDEXER @@@@@@@@@" << endl;

	if (select == 1)
	{

		cout << "MIN CONFIG" << endl;

		matrix A(rows, columns - 1 + rows);
		A.disp();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < columns - 1; j++)
			{
				A.p[i][j] = p[i][j];
			}
			for (int k = 0; k < rows; k++)
			{

				if (i == k)
				{
					A.p[i][columns - 1 + k] = 1;
				}
				else
				{
					A.p[i][columns - 1 + k] = 0;
				}
			}
		}
		A.disp();
		global_A = A;

		return A;
	}
}

int matrix ::comp(matrix &delj)
{
	cout << "in comp !!!!" << endl;
	int minpos = 0;
	// float minelem;
	int nflag = 0;

	for (int i = 0; i < delj.columns; i++)
	{
		if (delj.p[0][i] < 0)
		{
			nflag = 1;
			break;
		}
	}
	if (nflag == 1)
	{
		minpos = minarray(delj);
		cout << "index = " << minpos << endl;
	}
	else
	{
		cout << "all +positive" << endl;
	}
	return minpos;
}

int matrix ::minarray(matrix &min)
{
	cout << "in minarray !!!!" << endl;

	int index = 0;
	int size = sizeof(min) / sizeof(min.p, min.p[0][0]);

	for (int i = 0; i < min.columns; i++)
	{

		if (min.p[0][i] < min.p[0][index])
		{
			index = i;
		}
	}
	return index;
}

//----------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------

int main()
{
	cout << "start\n";

	matrix MAT;
	MAT.initializer();

	cout << "main end\n";

	return 0;
}

//-------------------------------------------- V3    -------------------------------------------

/*
 * main.cpp
 *
 *  Created on: 21 Jul 2021
 *      Author: prakharraj1302
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
using namespace std;

// class fraction
//{
// public:
//
//	int numerator;
//	int denominator;
//
//
//
//
//
//	fraction (); //default constractor....check
//	fraction (int,int); //default constractor with parameters ....check
//	fraction AddedTo (fraction value) const;
//	fraction MultipliedBy (fraction value) const; //binary observer type of operation
//	fraction Subtract (fraction value) const;
//	fraction DividedBy (fraction value) const;
//	fraction  isGreaterThan ();
//	fraction isEqualTo ();
//	fraction print() const;
//
//
// };

// frac::fraction()
//{
//	numerator = 0;
//	denominator = 1;
// }
//**********************************
// fraction::fraction(int newNumerator, int newDenominator)
//{
//	numerator=newNumerator;
//	denominator=newDenominator;
// }
//***************************

////********************************

// fraction fraction::print() const
// {
// int num=numerator;
// int den=denominator;
//  if(num>den)
//  {
//  for(int counter=2;counter<den;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//  else if(den>num)
//  {
//  for(int counter=2;counter<num;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//
//  cout<<num<<"/"<<den;
//  }
//**************************************#####################################

class frac
{

public:
	double num, den;
	frac()
	{
		cout << "frac cons called!!!!!!!!!!!!" << endl;

		num = 0;

		den = 1;
		cout << num << " " << den << endl;
	}

	frac fadd(frac &);
	frac fsub(frac &);
	frac fmul(frac &);
	frac fdiv(frac &);
	frac fneg();
	void fdisp(frac &);
};

void frac::fdisp(frac &fval)
{
	cout << "[" << fval.num << " / " << fval.den << "]" << endl;
}
frac frac::fadd(frac &value)
{
	frac result;
	result.num = (num * value.den) + (value.num * den);
	result.den = den * value.den;
	return result;
}
frac frac::fsub(frac &value)
{
	frac result;
	result.num = (num * value.den) - (value.num * den);
	result.den = den * value.den;
	return result;
}
//**************************************
frac frac::fmul(frac &value)
{
	frac result;
	result.num = num * value.num;
	result.den = den * value.den;
	return result;
}
//************************************
frac frac::fdiv(frac &value)
{
	frac result;
	result.num = num * value.den;
	result.den = den * value.num;
	return result;
}
//************************************
frac frac ::fneg()
{
	num = -num;

	return (*this);
}

//**************************************#####################################

class matrix : public frac
{

public:
	int rows, columns;

	double **p;
	frac **pf;

	void table(matrix &, matrix &, matrix &, int);

	matrix mx(matrix &);
	matrix sm(matrix &, matrix &);
	matrix sb(matrix &);

	void in();
	void inset(int, int);

	matrix merger(matrix &, matrix &, matrix &);
	matrix A_indexer(int);
	int ratio(matrix &);

	void disp();
	void dispstat();

	void frac_disp();

	matrix trn();
	matrix inverse();
	void imat();
	void cost();
	;

	void initializer();
	void iterator(matrix &, matrix &, matrix &);
	int comp(matrix &);
	int minarray(matrix &);

	matrix extractor(matrix &);
	matrix matcopy(matrix &);
	matrix ftab_solver(matrix &, int);
	matrix ftab_updater(matrix &, int, int);
	int iter_result(matrix &);
	int del_j(matrix &);
	int ITER(matrix &);
	void cost_updt(matrix &);
	void Aj_updt(matrix &);
	void Cj_updt(matrix &);
	int rat_comp(int[]);

	matrix(int = 1, int = 1);
	matrix(matrix &);
};
matrix global_V;
matrix global_B;
matrix global_Binv;
matrix global_Z;
matrix global_full_tab;
matrix global_A;
matrix global_Cbt;
matrix global_Aj;
matrix global_b;

int global_flag;
int global_ratio;
int global_counter;

void matrix ::frac_disp()
{
	cout << "frac_disp called!!!!!\n";

	for (int i = 0; i < rows; i++)
	{

		for (int j = 0; j < columns; j++)
		{
			// cout<<i<< " %  "<<j<<endl;
			cout << pf[i][j].num << "/" << pf[i][j].den << " ";
		}
		cout << endl;
	}
}

matrix matrix ::ftab_solver(matrix &ftab, int ent_row)
{
	cout << "ftab solver called !!!!!!111" << endl;

	//	matrix temp( (*this)  , 1);
	//	matrix mat( *this->rows , *this->columns);
	matrix mat(ftab);

	cout << "frac mat for table made $$$$$$$$$" << endl;

	for (int j = 0; j < ftab.rows; j++)
	{
		if (j == ent_row)
		{
			continue;
		}
		cout << j << endl;

		// double a ,b;
		// a = (pf[ent_row][columns - 1]).num / (pf[ent_row][columns - 1]).den;
		// b = (pf[j][columns - 1]).num / (pf[j][columns - 1]).den;

		// float x = temp.eq_solver(a , b);

		// cout<< a <<" ? "<< b <<" ? "  <<endl;

		mat = ftab_updater(ftab, ent_row, j);
	}

	return mat;
}
matrix matrix ::ftab_updater(matrix &tab, int p, int q)
{
	cout << "ftab updter called !!!!!!111" << endl;

	// double a = tab.p[p][columns - 1] ;
	// double b = tab.p[q][columns - 1];

	cout << "tab updater called !!!!!!!1!!!" << endl;
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	cout << " / " << p << " / " << q << endl;

	// cout<< "a= "<< a.num <<" "<<a.den<<" b = " <<b.num<<" "<<b.den <<endl;

	for (int i = 0; i < tab.columns; i++)
	{
		cout << " @@@@@@@@@@@@@@@@ TAB SOLVING @@@@@@@@@@@@@@@" << endl;

		cout << tab.pf[p][i].num << "/" << tab.pf[p][i].den << " $ " << tab.pf[q][i].num << "/" << tab.pf[q][i].den << endl;

		frac a, b, t1, t2;
		a = tab.pf[p][tab.columns - 1];
		b = tab.pf[q][tab.columns - 1];
		cout << "a= " << a.num << " " << a.den << " b = " << b.num << " " << b.den << endl
			 << "##" << endl;
		b = b.fneg();
		fdisp(b);
		t1 = b.fdiv(a);
		fdisp(t1);
		t2 = t1.fmul(tab.pf[p][i]);
		fdisp(t2);
		tab.pf[q][i] = (tab.pf[q][i]).fadd(t2);
		fdisp(tab.pf[q][i]);

		// tab.p[q][i]= tab.p[q][i] + ( (-b/a) * (tab.p[p][i])) ;			//dfrac

		cout << tab.pf[p][i].num << "/" << tab.pf[p][i].den << " $ " << tab.pf[q][i].num << "/" << tab.pf[q][i].den << endl;
	}
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

	for (int j = 0; j < tab.columns; j++)
	{
		frac a;
		a = tab.pf[p][tab.columns - 1];
		tab.pf[p][j] = (tab.pf[p][j]).fdiv(a);

		// tab.p[p][i] = tab.p[p][i] / tab.p[p][columns -1];				//dfrac
	}

	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

	cout << " ############SOLVED ###########" << endl;

	//	//extractor(tab);
	//	matrix dmat;
	//	dmat=extractor(tab);
	//	dmat.disp();
	iter_result(tab);
	return tab;
}
//----------------------------------------------------------------------------------------

matrix matrix ::extractor(matrix &fmat)
{

	cout << "in extractor ######################" << endl;
	matrix temp(fmat.rows, fmat.columns);
	int i, j;
	cout << fmat.rows << " " << fmat.columns << endl;

	for (i = 0; i < fmat.rows; i++)
	{

		for (j = 0; j < fmat.columns; j++)
		{

			temp.p[i][j] = ((fmat.pf[i][j].num) / (fmat.pf[i][j].den));
			cout << temp.p[i][j] << " ";
		}
		cout << endl;
	}
	return temp;
}
int matrix ::iter_result(matrix &fmat)
{

	cout << "in iterresult ######################" << endl;
	int flag;

	cout << "NEW FEASIBLE SOLN Xb " << endl;
	for (int i = 0; i < fmat.rows; i++)
	{
		cout << fmat.pf[i][0].num << "/" << fmat.pf[i][0].den << endl;
	}

	matrix t1(fmat.rows, fmat.columns);
	matrix binv(global_V.rows, global_V.rows);
	t1 = extractor(fmat);
	t1.disp();

	cout << "BINVERSE EXTRACTION ###################################" << endl;
	for (int i = 0; i < fmat.rows; i++)
	{
		for (int j = 0; j < global_V.rows; j++)
		{
			cout << i << " & " << j << endl;
			(binv.p[i][j]) = (t1.p[i][j + 1]);
		}
	}

	binv.disp();
	global_Binv = binv;
	flag = del_j(binv);

	//	int whilef =0;
	//	while( whilef == 0){
	//
	//		cout<<" @@@@@@@@@@@@@   ITERATION 2.... ##############"<<endl;
	//
	//		matrix TAB(global_V.rows,global_B.columns+2);
	//		matrix V;
	//		V = global_V;
	//
	//		V.disp();
	//		Xb.disp();
	//		binv.disp();
	//
	//		TAB = merger( Xb, binv , V , flag );
	//		TAB.disp();
	//
	//
	//		matrix Xn(global_V.rows,1);
	//		for(int k = 0; k< global_V.rows ; k++){
	//			Xn.p[k][0] = global_V.p[k][flag];
	//
	//		}
	//
	//		int ratio = Xb.ratio(Xn);
	//		global_ratio = ratio;
	//		cout<<"ratio is "<<ratio<<endl;
	//		matrix ftab (TAB);
	//		ftab.frac_disp();
	//
	//		matrix solved_mat;
	//		solved_mat = ftab_solver(ftab, ratio);
	//		solved_mat.disp();
	//
	//		whilef = 1 ;
	//
	//	}
	//	cout<<" out of iteration "<<endl;
	//	return 0;
}
int matrix::ITER(matrix &Binv)
{

	matrix Xb;
	matrix V;
	matrix TAB(global_V.rows, global_B.columns + 2);
	matrix xn(global_V.rows, 1);
	matrix solved_mat;
	matrix Xn;

	Xb = Binv.mx(global_b);

	int whilef = 0;
	while (whilef == 0)
	{

		cout << " @@@@@@@@@@@@@   ITERATION 2.... ##############" << endl;

		V = global_V;

		V.disp();
		Xb.disp();
		Binv.disp();

		for (int k = 0; k < global_Aj.rows; k++)
		{
			xn.p[k][0] = global_Aj.p[k][global_flag];
		}
		Xn = Binv.mx(xn);

		TAB = merger(Xb, Binv, Xn);
		TAB.disp();

		int ratio = Xb.ratio(Xn);
		global_ratio = ratio;
		cout << "ratio is " << ratio << endl;

		global_counter++;
		if (global_counter == 2)
		{
			exit(EXIT_FAILURE);
		}

		matrix ftab(TAB);
		ftab.frac_disp();

		solved_mat = ftab_solver(ftab, ratio);
		solved_mat.disp();

		whilef = 1;
	}
	cout << " out of iteration " << endl;
	return 0;
}

int matrix ::del_j(matrix &Binv)
{
	cout << "in delj ######################" << endl;

	matrix prc_1;
	matrix prc_2;
	matrix DEL;
	matrix Aj(global_V.rows, global_V.columns - 1);
	matrix cost_mat(global_V.rows, 1);
	matrix Cj(global_Z.rows, global_Z.columns);
	matrix cost_mat_trn;
	int op_flag = 0;

	cost_updt(cost_mat);
	Aj_updt(Aj);
	Cj_updt(Cj);
	cost_mat_trn = cost_mat.trn();

	global_Aj = Aj;

	cout << " @@@@@@ ALL VALUES UPDATED @@@@@@2" << endl;

	prc_1 = cost_mat_trn.mx(Binv);
	prc_1.disp();

	cout << "######" << endl;

	prc_2 = prc_1.mx(Aj);
	prc_2.disp();

	cout << "######" << endl;

	DEL = prc_2.sb(Cj);
	DEL.disp();
	op_flag = matrix::comp(DEL);
	cout << op_flag << endl;
	global_flag = op_flag;

	int x;
	x = ITER(Binv);
	return op_flag;
}

void matrix ::Aj_updt(matrix &Aj)
{
	cout << "in AJ updt ################" << endl;
	global_A.disp();

	for (int j = 0; j < global_V.rows; j++)
	{
		int i = 0;
		int k = 0;

		while ((i < global_V.columns - 1) && (k < global_A.columns))
		{
			cout << "UPDATING Aj" << endl
				 << "i " << i << "j " << j << " k " << k << endl;

			if (k != global_flag)
			{
				Aj.p[j][i] = global_A.p[j][k];
				i++;
				k++;
				continue;
			}
			if (k == global_flag)
			{
				k++;
				continue;
			}
			else
			{
				cout << "erroe !!" << endl;
			}
		}
	}
	Aj.disp();

	//	int fl =0;
	//	for(int i=0; i < global_A.rows ; i++){
	//		int c=0;
	//		for(int j=0; j < global_A.columns - 1 ; j++){
	//			if (j == global_flag){
	//				continue;
	//			}
	//			Aj.p[i][j] = global_A.p[i][j];
	//			c++;
	//			if(c==3){
	//				fl=1;
	//				break;
	//			}
	//
	//
	//		}
	//		if (fl==1){
	//			break;
	//		}
	//
	//	}
}

void matrix ::Cj_updt(matrix &Cj)
{
	cout << "in z updt ################" << endl;

	cout << " in min config Z with added vars --- " << endl;

	matrix Zfull(global_Z.rows, global_Z.columns + global_V.rows);
	for (int i = 0; i < global_Z.columns; i++)
	{
		if (i == 3)
		{
			Zfull.p[0][i] = 0;
			continue;
		}
		Zfull.p[0][i] = global_Z.p[0][i];
	}
	Zfull.disp();

	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "Cj updation #########" << endl;

	int i = 0;
	int k = 0;
	while ((i < global_V.columns - 1) && (k < Zfull.columns))
	{
		if (k != global_flag)
		{
			Cj.p[0][i] = Zfull.p[0][k];
			i++;
			k++;
			continue;
		}
		if (k == global_flag)
		{

			k++;
			continue;
		}
		else
		{
			cout << "error!!" << endl;
		}
	}
	Cj.disp();

	//	int c =0 ;
	//	for (int j =0 ; j< Z.columns; j++){
	//		if( j == global_flag ){
	//			continue;
	//		}
	//		Cj.p[0][j] = Z.p[0][j];
	//		c++;
	//		if( c == global_V.columns-1){
	//			break;
	//		}
	//	}
	cout << "##################" << endl;
}

void matrix ::cost_updt(matrix &costmat)
{

	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	cout << "cost_updt called ######" << endl;
	cout << " globalflag- " << global_flag << " globarat- " << global_ratio << endl;
	cout << " global _v" << endl;
	global_V.disp();

	int i = 0;
	while (i < global_V.rows)
	{
		cout << "UPDATING" << endl;
		if (i != global_ratio)
		{
			costmat.p[0][i] = global_Cbt.p[0][i];
			i++;
			continue;
		}
		if (i == global_ratio)
		{
			costmat.p[0][i] = global_Z.p[0][global_flag];
			i++;
			continue;
		}
		else
		{
			cout << "error" << endl;
		}
	}

	//	for(int i=0; i < global_V.rows ; i++){
	//		if(i == global_ratio){
	//			costmat.p[i][0] = global_Z.p[1][global_flag];
	//			continue;
	//		}
	//		costmat.p[i][0] = global_V.p[global_V.columns][i];
	//
	//	}
	costmat.disp();
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
}

matrix matrix ::matcopy(matrix &mat)
{
	cout << "matcop[y called ######" << endl;

	matrix temp(mat);

	int i, j;
	int r, c;

	r = mat.rows;
	c = mat.columns;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			pf[i][j].num = mat.p[i][j];
		}
	}
	return temp;
}

matrix matrix::trn()
{
	cout << "trn fx called!!!" << endl;

	matrix tempmat(columns, rows);
	//	cout<<columns<<rows<<endl;

	// tempmat.dispstat();
	int n, m, i, j;
	n = rows;
	m = columns;
	//	cout<<"!!!!!1";

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			//	cout<<i<<" "<<j<<" "<<p[i][j]<<" "<<tempmat.p[i][j];

			tempmat.p[j][i] = p[i][j];
			tempmat.disp();
		}
	}
	return tempmat;
}

matrix matrix ::inverse()
{
	cout << "invrese fx called !!!" << endl;

	matrix B(rows, columns), C(columns, columns * 2);
	int row, column, step, i;
	double mult;

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			C.p[row][column] = p[row][column];

	for (i = 0; i < rows; i++)
		C.p[i][rows + i] = 1;

	for (step = 0; step < rows - 1; step++)
	{
		for (row = step + 1; row < rows; row++)
		{
			mult = C.p[row][step] / C.p[step][step];
			for (column = step; column < 2 * rows; column++)
				C.p[row][column] -= mult * C.p[step][column];
		}
	}

	for (step = 1; step <= rows - 1; step++)
	{
		for (row = rows - step - 1; row >= 0; row--)
		{
			mult = C.p[row][rows - step] / C.p[rows - step][rows - step];
			for (column = rows; column < 2 * rows; column++)
				C.p[row][column] -=
					mult * C.p[rows - step][column];
		}
	}

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			B.p[row][column] = C.p[row][rows + column] / C.p[row][row];
	return B;
}

matrix::matrix(int r, int c)
{

	// matrix  constructor

	cout << "mat constructor called !!!!!" << r << " " << c << endl;

	int i, j;
	rows = r;
	columns = c;

	p = new double *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		p[i] = new double[columns];
		// assert(p[i] != NULL);
	}

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			p[i][j] = 0;
	pf = NULL;
}

matrix::matrix(matrix &mf) : frac()
{

	// matrix  constructor

	cout << " ############  mat FRACTION  constructor called !!!!!" << endl;

	int i, j;
	rows = mf.rows;
	columns = mf.columns;

	pf = new frac *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		pf[i] = new frac[columns];
		// assert(p[i] != NULL);
	}
	cout << "frac alloc sucessful" << endl;

	if (mf.p == NULL)
	{

		for (i = 0; i < mf.rows; i++)
		{
			for (j = 0; j < mf.columns; j++)
			{
				cout << i << " % " << j << endl;
				pf[i][j] = mf.pf[i][j];
			}
		}
	}
	else
	{
		for (i = 0; i < mf.rows; i++)
		{
			for (j = 0; j < mf.columns; j++)
			{
				cout << i << " % " << j << endl;
				pf[i][j].num = mf.p[i][j];
			}
		}
	}
	p = NULL;
	cout << "mat FRACTION  constructor LEFT  " << endl;
}

void matrix ::imat()
{

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			if (i == j)
			{
				p[i][j] = 1;
			}
			else
			{
				p[i][j] = 0;
			}
		}
	}
}

void matrix ::cost()
{

	cout << "cost called !!" << endl;

	for (int i = 0; i < rows; i++)
	{
		p[i][1] = 0;
	}
}

matrix matrix ::mx(matrix &B)
{
	cout << "mx fx called!!!!\n";

	float sum;
	int i, j, k, n, m, l;
	n = rows;
	m = columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix operator*(matrix &A, matrix &B)
{
	double sum;
	int i, j, k, n, m, l;
	n = A.rows;
	m = A.columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + A.p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix matrix ::sm(matrix &A, matrix &B)
{

	cout << "sm called!!!!!" << endl;

	int i, j, n, m;

	n = A.rows;
	m = A.columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = B.p[i][j] + A.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

matrix matrix ::sb(matrix &B)
{

	cout << "sb called!!!!!" << endl;

	int i, j, n, m;

	n = rows;
	m = columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = p[i][j] - B.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

void matrix ::in()
{
	cout << "matrix fx called!!!!\n";

	delete[] p;

	cout << "Enter the number of Rows and Columns: ";
	cin >> rows >> columns;
	p = new double *[rows];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new double[columns];
	}

	cout << "Enter the " << rows << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << columns << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::inset(int rset, int cset)
{
	cout << "matrix inset called!!!!\n";

	delete[] p;

	p = new double *[rset];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new double[cset];
	}

	cout << "Enter the " << rset << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << cset << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::disp()
{
	cout << "disp called!!!!!\n";

	for (int i = 0; i < rows; i++)
	{

		for (int j = 0; j < columns; j++)
		{
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
}

void matrix ::dispstat()
{
	cout << "dispstat called!!!!!\n";
	cout << rows << " " << columns << endl;
}

void matrix ::initializer()
{
	cout << "IN ITER !!!!" << endl;

	int nz, nv, ndesc;

	cout << "enter max mat" << endl;
	cout << "enter no of var" << endl;
	cin >> nz;
	matrix matz(1, nz);
	matz.inset(1, nz);
	matz.disp();
	global_Z = matz;

	cout << "no of desc" << endl;
	cin >> ndesc;
	matrix matv(ndesc, (nz + 1));
	matv.inset(ndesc, (nz + 1));
	matv.disp();
	global_V = matv;

	matrix matB(ndesc, ndesc);
	matB.imat();
	matB.disp();
	global_B = matB;

	iterator(matz, matv, matB);
}

void matrix ::iterator(matrix &z, matrix &v, matrix &B)
{
	cout << "IN ITERRATOR !!!!" << endl;

	matrix Xb;
	matrix DEL;
	matrix COST;
	matrix prc_1;
	matrix prc_2;
	matrix prc_3;

	matrix Cj;
	matrix COST_trn;
	matrix temp_1;

	matrix Aj(v.rows, v.columns - 1);

	for (int i = 0; i < v.rows; i++)
	{
		for (int j = 0; j < v.columns - 1; j++)
		{
			cout << i << " " << j << endl;

			Aj.p[i][j] = v.p[i][j];
		}
	}
	Aj.disp();

	matrix Binv;
	Binv = B.inverse();
	Binv.disp();

	matrix b(v.rows, 1);
	for (int i = 0; i < v.rows; i++)
	{
		b.p[i][0] = v.p[i][(v.columns - 1)];
	}
	b.disp();

	global_b = b;

	matrix cB(v.rows, 1);
	cB.cost();
	cB.disp();
	COST_trn = cB.trn();
	COST_trn.disp();
	global_Cbt = COST_trn;

	Cj = z;
	Cj.disp();

	int opti_flag = 0;

	for (int j = 0; j < 1; j++)
	{

		cout << "ITERATION" << endl
			 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

		Binv.disp();
		b.disp();
		COST_trn.disp();

		Xb = Binv.mx(b);
		Xb.disp();
		cout << "######" << endl;

		prc_1 = COST_trn.mx(Binv);
		prc_1.disp();

		cout << "######" << endl;

		prc_2 = prc_1.mx(Aj);
		prc_2.disp();

		cout << "######" << endl;
		DEL = prc_2.sb(Cj);
		DEL.disp();

		cout << "######" << endl;
		opti_flag = temp_1.comp(DEL);
		global_flag = opti_flag;

		cout << "opti_flag = " << opti_flag << endl;
		DEL.disp();

		cout << "######" << endl;
		matrix table(v.rows, B.columns + 2);

		matrix xn(v.rows, 1);
		for (int k = 0; k < v.rows; k++)
		{
			xn.p[k][0] = v.p[k][opti_flag];
		}
		matrix Xn;
		Xn = Binv.mx(xn);

		table = prc_3.merger(Xb, Binv, Xn);
		table.disp();

		matrix Atab(v.rows, v.columns - 1 + v.rows);

		Atab = v.A_indexer(1);

		cout << "######" << endl;

		int while_flag = 1;

		while (while_flag == 1)
		{

			cout << "in WHILE @@@@@@@@@@@@@@@@@@" << endl;

			int rat = Xb.ratio(Xn);
			global_ratio = rat;
			cout << "ratio is " << rat << endl;
			matrix ftab(table);
			ftab.frac_disp();

			ftab_solver(ftab, rat);

			cout << "######### END OF WHILE  ############" << endl;

			while_flag = 0;
		}
	}
}
int matrix ::ratio(matrix &X)
{

	cout << "in ratio !!!!!!!!!!1" << endl;

	matrix rat(X.rows, 1);
	int rat_arr[X.rows];

	for (int i = 0; i < X.rows; i++)
	{
		rat.p[i][0] = (p[i][0]) / (X.p[i][0]);
		rat_arr[i] = rat.p[i][0];
		cout << "RATIO ======" << rat.p[i][0] << rat_arr[i] << endl;
	}

	int pos = 0;
	//	//int size = sizeof(rat.p)/sizeof(rat.p[0][0]);
	//
	//
	//	for(int j = 0; j < rat.rows ; j++)
	//	{
	//		if(rat.p[j][0] < rat.p[0][0] && rat.p[j][0] > 0 )
	//			index = j;
	//	}

	// pos = rat_comp(rat_arr);
	cout << pos;

	return pos;
}
int matrix ::rat_comp(int arr[])
{

	int n = sizeof(arr) / sizeof(arr[0]);
	int x = 0;
	int pos = 0;
	while (x < n)
	{
		if (arr[x] < 0)
		{
			x++;
		}
	}

	for (int j = x; j < n; j++)
	{
		cout << x << " $ " << j << endl;
		if (arr[j] < arr[x] && arr[j] > 0)
		{
			x = j;
			pos = j;
		}
	}
	return pos;
}

matrix matrix ::merger(matrix &Xb, matrix &Bin, matrix &Xn)
{
	cout << "merger called !!!!!!!!!!!!1" << endl;

	matrix merge_mat(Xn.rows, Bin.columns + 2);
	//	matrix full_tab(Xn.rows,Bin.columns+ Xn.columns);
	//
	//	for (int i =0 ; i < Xb.rows ; i++){
	//		merge_mat.p[i][0] = Xb.p[i][0];
	//		for (int j = 0; j < Bin.columns ; j++ ){
	//			merge_mat.p[i][1+j] = Bin.p[i][j] ;
	//		}
	//		for (int j = 0; j < Xn.columns ; j++ ){
	//
	//			merge_mat.p[i][Bin.columns+1] = Xn.p[i][j];
	//
	//		}
	//	}
	//	global_full_tab = full_tab;

	for (int i = 0; i < Xb.rows; i++)
	{
		merge_mat.p[i][0] = Xb.p[i][0];
		for (int j = 0; j < Bin.columns; j++)
		{
			merge_mat.p[i][1 + j] = Bin.p[i][j];
		}

		merge_mat.p[i][Bin.columns + 1] = Xn.p[i][0];
	}
	return merge_mat;
}

matrix matrix ::A_indexer(int select)
{
	cout << "in INDEXER @@@@@@@@@" << endl;

	if (select == 1)
	{

		cout << "MIN CONFIG" << endl;

		matrix A(rows, columns - 1 + rows);
		A.disp();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < columns - 1; j++)
			{
				A.p[i][j] = p[i][j];
			}
			for (int k = 0; k < rows; k++)
			{

				if (i == k)
				{
					A.p[i][columns - 1 + k] = 1;
				}
				else
				{
					A.p[i][columns - 1 + k] = 0;
				}
			}
		}
		A.disp();
		global_A = A;

		return A;
	}
}

int matrix ::comp(matrix &delj)
{
	cout << "in comp !!!!" << endl;
	int minpos = 0;
	// float minelem;
	int nflag = 0;

	for (int i = 0; i < delj.columns; i++)
	{
		if (delj.p[0][i] < 0)
		{
			nflag = 1;
			break;
		}
	}
	if (nflag == 1)
	{
		minpos = minarray(delj);
		cout << "index = " << minpos << endl;
	}
	else
	{
		cout << "all +positive" << endl;
	}
	return minpos;
}

int matrix ::minarray(matrix &min)
{
	cout << "in minarray !!!!" << endl;

	int index = 0;
	int size = sizeof(min) / sizeof(min.p, min.p[0][0]);

	for (int i = 0; i < min.columns; i++)
	{

		if (min.p[0][i] < min.p[0][index])
		{
			index = i;
		}
	}
	return index;
}

//----------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------

int main()
{
	cout << "start\n";

	matrix MAT;
	MAT.initializer();

	cout << "main end\n";

	return 0;
}

//-------------------------------------- V 3 -------------------------------------------------

/*
 * main.cpp
 *
 *  Created on: 21 Jul 2021
 *      Author: prakharraj1302
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
using namespace std;

// class fraction
//{
// public:
//
//	int numerator;
//	int denominator;
//
//
//
//
//
//	fraction (); //default constractor....check
//	fraction (int,int); //default constractor with parameters ....check
//	fraction AddedTo (fraction value) const;
//	fraction MultipliedBy (fraction value) const; //binary observer type of operation
//	fraction Subtract (fraction value) const;
//	fraction DividedBy (fraction value) const;
//	fraction  isGreaterThan ();
//	fraction isEqualTo ();
//	fraction print() const;
//
//
// };

// frac::fraction()
//{
//	numerator = 0;
//	denominator = 1;
// }
//**********************************
// fraction::fraction(int newNumerator, int newDenominator)
//{
//	numerator=newNumerator;
//	denominator=newDenominator;
// }
//***************************

////********************************

// fraction fraction::print() const
// {
// int num=numerator;
// int den=denominator;
//  if(num>den)
//  {
//  for(int counter=2;counter<den;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//  else if(den>num)
//  {
//  for(int counter=2;counter<num;counter++)
//          {
//          while(num%counter==0 & den%counter==0)
//          {
//          num=(num/counter);
//          den=(den/counter);
//          }
//          }
//  }
//
//  cout<<num<<"/"<<den;
//  }
//**************************************#####################################

class frac
{

public:
	double num, den;
	frac()
	{
		// cout<<"frac cons called!!!!!!!!!!!!"<<endl;

		num = 0;

		den = 1;
		// cout<<num<<" "<<den<<endl;
	}

	frac fadd(frac &);
	frac fsub(frac &);
	frac fmul(frac &);
	frac fdiv(frac &);
	frac fneg();
	void fdisp(frac &);
};

void frac::fdisp(frac &fval)
{
	cout << "[" << fval.num << " / " << fval.den << "]" << endl;
}
frac frac::fadd(frac &value)
{
	frac result;
	result.num = (num * value.den) + (value.num * den);
	result.den = den * value.den;
	return result;
}
frac frac::fsub(frac &value)
{
	frac result;
	result.num = (num * value.den) - (value.num * den);
	result.den = den * value.den;
	return result;
}
//**************************************
frac frac::fmul(frac &value)
{
	frac result;
	result.num = num * value.num;
	result.den = den * value.den;
	return result;
}
//************************************
frac frac::fdiv(frac &value)
{
	frac result;
	result.num = num * value.den;
	result.den = den * value.num;
	return result;
}
//************************************
frac frac ::fneg()
{
	num = -num;

	return (*this);
}

//**************************************#####################################

class matrix : public frac
{

public:
	int rows, columns;

	double **p;
	frac **pf;

	void table(matrix &, matrix &, matrix &, int);

	matrix mx(matrix &);
	matrix sm(matrix &, matrix &);
	matrix sb(matrix &);

	void in();
	void inset(int, int);

	matrix merger(matrix &, matrix &, matrix &);
	matrix A_indexer(int);
	int ratio(matrix &);

	void disp();
	void dispstat();

	void frac_disp();

	matrix trn();
	matrix inverse();
	void imat();
	void cost();
	;

	int initializer();
	int iterator(matrix &, matrix &, matrix &);
	int comp(matrix &);
	int minarray(matrix &);

	matrix extractor(matrix &);
	matrix matcopy(matrix &);
	int ftab_solver(matrix &, int);
	int ftab_updater(matrix &, int, int);
	int iter_result(matrix &);
	int del_j(matrix &);
	int ITER(matrix &);
	void cost_updt(matrix &);
	void Aj_updt(matrix &);
	void Cj_updt(matrix &);
	int rat_comp(int[]);
	void head();
	void head_ch(int, int);
	int result();

	matrix(int = 1, int = 1);
	matrix(matrix &);
};
matrix global_V;
matrix global_B;
matrix global_Binv;
matrix global_Z;
matrix global_full_tab;
matrix global_A;
matrix global_Cbt;
matrix global_Aj;
matrix global_b;
matrix global_Cj;
matrix global_head;
matrix global_Xb;

int global_flag;
int global_ratio;
int global_counter;
int run = 0;

void matrix ::head()
{

	cout << "************** in head  **********8**" << endl;
	matrix head(1, global_Z.columns + global_V.rows);

	for (int i = 0; i < head.columns; i++)
	{
		if (i < global_Z.columns)
		{
			head.p[0][i] = 1;
		}
		if (i >= head.columns)
		{
			head.p[0][i] = 0;
		}
	}
	cout << "init head  -------";
	head.disp();
	global_head = head;
}

void matrix ::head_ch(int f, int r)
{
	cout << "************** in head_ch  **********8**" << endl;
	cout << f << " ***** " << r;

	global_head.disp();

	int fl = -1;
	int rt = -1;
	int h0 = 0;
	int h1 = 0;

	for (int i = 0; i < global_head.columns; i++)
	{

		if (global_head.p[0][i] == 1)
		{
			fl++;
			cout << fl << i << endl;
			if (fl == f)
			{
				h0 = i; // global_head.p[0][i] = 0;
			}
		}
	}
	for (int j = 0; j < global_head.columns; j++)
	{
		if (global_head.p[0][j] == 0)
		{
			rt++;
			cout << rt << j << endl;
			if (rt == r)
			{
				h1 = j;
			}
		}
	}

	global_head.p[0][h0] = 0;

	global_head.p[0][h1] = 1;

	global_head.disp();
}

void matrix ::frac_disp()
{
	cout << "frac_disp called!!!!!\n";

	for (int i = 0; i < rows; i++)
	{

		for (int j = 0; j < columns; j++)
		{
			// cout<<i<< " %  "<<j<<endl;
			cout << pf[i][j].num << "/" << pf[i][j].den << " ";
		}
		cout << endl;
	}
}

int matrix ::ftab_solver(matrix &ftab, int ent_row)
{
	cout << "ftab solver called !!!!!!111" << endl;

	//	matrix temp( (*this)  , 1);
	//	matrix mat( *this->rows , *this->columns);
	matrix mat(ftab);
	int check;

	cout << "frac mat for table made $$$$$$$$$" << endl;

	for (int j = 0; j < ftab.rows; j++)
	{
		if (j == ent_row)
		{
			continue;
		}
		cout << j << endl;

		// double a ,b;
		// a = (pf[ent_row][columns - 1]).num / (pf[ent_row][columns - 1]).den;
		// b = (pf[j][columns - 1]).num / (pf[j][columns - 1]).den;

		// float x = temp.eq_solver(a , b);

		// cout<< a <<" ? "<< b <<" ? "  <<endl;

		check = ftab_updater(ftab, ent_row, j);
		if (check == -1)
		{
			return -1;
		}
	}

	return 0;
}
int matrix ::ftab_updater(matrix &tab, int p, int q)
{
	cout << "ftab updter called !!!!!!111" << endl;

	// double a = tab.p[p][columns - 1] ;
	// double b = tab.p[q][columns - 1];

	cout << "tab updater called !!!!!!!1!!!" << endl;
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	cout << " / " << p << " / " << q << endl;

	// cout<< "a= "<< a.num <<" "<<a.den<<" b = " <<b.num<<" "<<b.den <<endl;

	for (int i = 0; i < tab.columns; i++)
	{
		cout << " @@@@@@@@@@@@@@@@ TAB SOLVING @@@@@@@@@@@@@@@" << endl;

		cout << tab.pf[p][i].num << "/" << tab.pf[p][i].den << " $ " << tab.pf[q][i].num << "/" << tab.pf[q][i].den << endl;

		frac a, b, t1, t2;
		a = tab.pf[p][tab.columns - 1];
		b = tab.pf[q][tab.columns - 1];
		cout << "a= " << a.num << " " << a.den << " b = " << b.num << " " << b.den << endl
			 << "##" << endl;
		b = b.fneg();
		fdisp(b);
		t1 = b.fdiv(a);
		fdisp(t1);
		t2 = t1.fmul(tab.pf[p][i]);
		fdisp(t2);
		tab.pf[q][i] = (tab.pf[q][i]).fadd(t2);
		fdisp(tab.pf[q][i]);

		// tab.p[q][i]= tab.p[q][i] + ( (-b/a) * (tab.p[p][i])) ;			//dfrac

		cout << tab.pf[p][i].num << "/" << tab.pf[p][i].den << " $ " << tab.pf[q][i].num << "/" << tab.pf[q][i].den << endl;
	}
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

	for (int j = 0; j < tab.columns; j++)
	{
		frac a;
		a = tab.pf[p][tab.columns - 1];
		tab.pf[p][j] = (tab.pf[p][j]).fdiv(a);

		// tab.p[p][i] = tab.p[p][i] / tab.p[p][columns -1];				//dfrac
	}

	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	tab.frac_disp();
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

	cout << " ############SOLVED ###########" << endl;

	//	//extractor(tab);
	//	matrix dmat;
	//	dmat=extractor(tab);

	//	dmat.disp();

	int check;
	check = iter_result(tab);
	if (check == -1)
	{
		return 0;
	}
	return 0;
}
//----------------------------------------------------------------------------------------

matrix matrix ::extractor(matrix &fmat)
{

	cout << "in extractor ######################" << endl;
	matrix temp(fmat.rows, fmat.columns);
	int i, j;
	cout << fmat.rows << " " << fmat.columns << endl;

	for (i = 0; i < fmat.rows; i++)
	{

		for (j = 0; j < fmat.columns; j++)
		{

			temp.p[i][j] = ((fmat.pf[i][j].num) / (fmat.pf[i][j].den));
			cout << temp.p[i][j] << " ";
		}
		cout << endl;
	}
	return temp;
}
int matrix ::iter_result(matrix &fmat)
{

	cout << "in iterresult ######################" << endl;
	int flag;

	cout << "NEW FEASIBLE SOLN Xb " << endl;
	for (int i = 0; i < fmat.rows; i++)
	{
		cout << fmat.pf[i][0].num << "/" << fmat.pf[i][0].den << endl;
	}

	matrix t1(fmat.rows, fmat.columns);
	matrix t2(fmat.rows, 1);

	matrix binv(global_V.rows, global_V.rows);
	t1 = extractor(fmat);
	t1.disp();

	for (int i = 0; i < fmat.rows; i++)
	{
		t2.p[i][0] = t1.p[i][0];
	}
	global_Xb = t2;

	cout << "BINVERSE EXTRACTION ###################################" << endl;
	for (int i = 0; i < fmat.rows; i++)
	{
		for (int j = 0; j < global_V.rows; j++)
		{
			cout << i << " & " << j << endl;
			(binv.p[i][j]) = (t1.p[i][j + 1]);
		}
	}

	binv.disp();
	global_Binv = binv;
	flag = del_j(binv);
	if (flag == -1)
	{
		return -1;
	}

	//	int whilef =0;
	//	while( whilef == 0){
	//
	//		cout<<" @@@@@@@@@@@@@   ITERATION 2.... ##############"<<endl;
	//
	//		matrix TAB(global_V.rows,global_B.columns+2);
	//		matrix V;
	//		V = global_V;
	//
	//		V.disp();
	//		Xb.disp();
	//		binv.disp();
	//
	//		TAB = merger( Xb, binv , V , flag );
	//		TAB.disp();
	//
	//
	//		matrix Xn(global_V.rows,1);
	//		for(int k = 0; k< global_V.rows ; k++){
	//			Xn.p[k][0] = global_V.p[k][flag];
	//
	//		}
	//
	//		int ratio = Xb.ratio(Xn);
	//		global_ratio = ratio;
	//		cout<<"ratio is "<<ratio<<endl;
	//		matrix ftab (TAB);
	//		ftab.frac_disp();
	//
	//		matrix solved_mat;
	//		solved_mat = ftab_solver(ftab, ratio);
	//		solved_mat.disp();
	//
	//		whilef = 1 ;
	//
	//	}
	//	cout<<" out of iteration "<<endl;
	//	return 0;
}
int matrix::ITER(matrix &Binv)
{

	matrix Xb;
	matrix V;
	matrix TAB(global_V.rows, global_B.columns + 2);
	matrix xn(global_V.rows, 1);
	matrix solved_mat;
	matrix Xn;

	Xb = Binv.mx(global_b);

	int whilef = 0;
	while (whilef == 0)
	{

		cout << " @@@@@@@@@@@@@   ITERATION 2.... ##############" << endl;

		V = global_V;

		V.disp();
		Xb.disp();
		Binv.disp();

		for (int k = 0; k < global_Aj.rows; k++)
		{
			xn.p[k][0] = global_Aj.p[k][global_flag];
		}
		Xn = Binv.mx(xn);

		TAB = merger(Xb, Binv, Xn);
		TAB.disp();

		int ratio = Xb.ratio(Xn);
		global_ratio = ratio;
		cout << "ratio is " << ratio << endl;

		head_ch(global_flag, ratio);

		global_counter++;
		if (global_counter == 2)
		{
			exit(EXIT_FAILURE);
		}

		matrix ftab(TAB);
		ftab.frac_disp();

		solved_mat = ftab_solver(ftab, ratio);
		solved_mat.disp();

		whilef = 1;
	}
	cout << " out of iteration " << endl;
	return 0;
}

int matrix ::del_j(matrix &Binv)
{
	cout << "in delj ######################" << endl;

	matrix prc_1;
	matrix prc_2;
	matrix DEL;
	matrix Aj(global_V.rows, global_V.columns - 1);
	matrix cost_mat(global_V.rows, 1);
	matrix Cj(global_Z.rows, global_Z.columns);
	matrix cost_mat_trn;
	int op_flag = 0;

	cost_updt(cost_mat);
	Cj_updt(Cj);
	Aj_updt(Aj);

	cost_mat_trn = cost_mat.trn();

	global_Aj = Aj;

	cout << " @@@@@@ ALL VALUES UPDATED @@@@@@2" << endl;

	cost_mat_trn.disp();
	Binv.disp();
	Aj.disp();
	Cj.disp();

	prc_1 = cost_mat_trn.mx(Binv);
	prc_1.disp();

	cout << "######" << endl;

	prc_2 = prc_1.mx(Aj);
	prc_2.disp();

	cout << "######" << endl;

	DEL = prc_2.sb(Cj);
	DEL.disp();
	op_flag = matrix::comp(DEL);
	if (run == 0)
	{
		return -1;
	}
	cout << op_flag << endl;
	global_flag = op_flag;

	int x;
	x = ITER(Binv);
	return op_flag;
}

void matrix ::Aj_updt(matrix &Aj)
{
	cout << "in AJ updt ################" << endl;
	global_A.disp();
	global_Aj.disp();
	global_head.disp();

	for (int i = 0; i < global_V.rows; i++)
	{
		int j = 0;
		int k = 0;
		cout << i << j << k << endl;
		while (j < Aj.columns && k < global_A.columns)
		{

			cout << i << j << k << " ?" << endl;
			cout << global_head.p[0][j] << " ^ " << global_A.p[i][k] << endl;

			if (global_head.p[0][k] == 1)
			{
				Aj.p[i][j] = global_A.p[i][k];
				j++;
				k++;
				continue;
			}
			if (global_head.p[0][k] == 0)
			{
				k++;
				continue;
			}
		}
	}

	//			while ( (i < global_V.columns-1) && (k < global_A.columns) ){
	//				cout<<"UPDATING Aj"<<endl<<"i "<<i<<"j "<<j<<" k " <<k<<endl;
	//
	//
	//				if( k != global_flag){
	//					if( k <= global_Aj.columns-1){
	//						Aj.p[j][i] = global_Aj.p[j][k];
	//						i++;
	//						k++;
	//						continue;
	//					}
	//					if (k > global_Aj.columns-1){
	//						Aj.p[j][i] = global_A.p[j][k];
	//						i++;
	//						k++;
	//						continue;
	//					}
	//				}
	//				if( k == global_flag ){
	//					k++;
	//					continue;
	//				}
	//				else{
	//					cout<<"erroe !!"<<endl;
	//				}
	//			}

	Aj.disp();

	//	int fl =0;
	//	for(int i=0; i < global_A.rows ; i++){
	//		int c=0;
	//		for(int j=0; j < global_A.columns - 1 ; j++){
	//			if (j == global_flag){
	//				continue;
	//			}
	//			Aj.p[i][j] = global_A.p[i][j];
	//			c++;
	//			if(c==3){
	//				fl=1;
	//				break;
	//			}
	//
	//
	//		}
	//		if (fl==1){
	//			break;
	//		}
	//
	//	}
}

void matrix ::Cj_updt(matrix &Cj)
{
	cout << "in z updt ################" << endl;

	cout << " in min config Z with added vars --- " << endl;

	matrix Zfull(global_Z.rows, global_Z.columns + global_V.rows);
	for (int i = 0; i < global_Z.columns; i++)
	{
		if (i == 3)
		{
			Zfull.p[0][i] = 0;
			continue;
		}
		Zfull.p[0][i] = global_Z.p[0][i];
	}
	Zfull.disp();

	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "Cj updation #########" << endl;
	global_Cj.disp();

	int i = 0;
	int k = 0;

	int j = 0;
	cout << i << j << k << endl;
	while (j < Cj.columns && k < Zfull.columns)
	{

		cout << i << j << k << " ?" << endl;
		cout << global_head.p[0][j] << " ^ " << Zfull.p[0][k] << endl;

		if (global_head.p[0][k] == 1)
		{
			Cj.p[0][j] = Zfull.p[0][k];
			j++;
			k++;
			continue;
		}
		if (global_head.p[0][k] == 0)
		{
			k++;
			continue;
		}
	}

	//	while((i < global_V.columns-1) && (k < Zfull.columns )){
	//		if( k != global_flag)
	//		{
	//			if (k > global_Cj.columns-1){
	//				Cj.p[0][i] = Zfull.p[0][k];
	//				i++;
	//				k++;
	//				continue;
	//
	//			}
	//			if(k<=global_Cj.columns-1){
	//				Cj.p[0][i] = global_Cj.p[0][k];
	//
	//				i++;
	//				k++;
	//				continue;
	//			}
	//		}
	//		if(k == global_flag)
	//		{
	//
	//			k++;
	//			continue;
	//
	//		}
	//		else{
	//			cout<<"error!!"<<endl;
	//
	//		}
	//	}
	global_Cj = Cj;
	Cj.disp();

	//	int c =0 ;
	//	for (int j =0 ; j< Z.columns; j++){
	//		if( j == global_flag ){
	//			continue;
	//		}
	//		Cj.p[0][j] = Z.p[0][j];
	//		c++;
	//		if( c == global_V.columns-1){
	//			break;
	//		}
	//	}
	cout << "##################" << endl;
}

void matrix ::cost_updt(matrix &costmat)
{

	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	cout << "cost_updt called ######" << endl;
	cout << " globalflag- " << global_flag << " globarat- " << global_ratio << endl;
	cout << " global _v" << endl;
	global_V.disp();
	global_Cbt.disp();
	global_Z.disp();
	global_Cj.disp();

	int i = 0;
	while (i < global_V.rows)
	{
		cout << "UPDATING" << endl;
		if (i != global_ratio)
		{
			costmat.p[i][0] = global_Cbt.p[0][i];
			i++;
			continue;
		}
		if (i == global_ratio)
		{
			costmat.p[i][0] = global_Cj.p[0][global_flag];
			i++;
			continue;
		}
		else
		{
			cout << "error" << endl;
		}
	}

	//	for(int i=0; i < global_V.rows ; i++){
	//		if(i == global_ratio){
	//			costmat.p[i][0] = global_Z.p[1][global_flag];
	//			continue;
	//		}
	//		costmat.p[i][0] = global_V.p[global_V.columns][i];
	//
	//	}
	global_Cbt = costmat.trn();
	costmat.disp();
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
}

matrix matrix ::matcopy(matrix &mat)
{
	cout << "matcop[y called ######" << endl;

	matrix temp(mat);

	int i, j;
	int r, c;

	r = mat.rows;
	c = mat.columns;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			pf[i][j].num = mat.p[i][j];
		}
	}
	return temp;
}

matrix matrix::trn()
{
	cout << "trn fx called!!!" << endl;

	matrix tempmat(columns, rows);
	//	cout<<columns<<rows<<endl;

	// tempmat.dispstat();
	int n, m, i, j;
	n = rows;
	m = columns;
	//	cout<<"!!!!!1";

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			//	cout<<i<<" "<<j<<" "<<p[i][j]<<" "<<tempmat.p[i][j];

			tempmat.p[j][i] = p[i][j];
			tempmat.disp();
		}
	}
	return tempmat;
}

matrix matrix ::inverse()
{
	cout << "invrese fx called !!!" << endl;

	matrix B(rows, columns), C(columns, columns * 2);
	int row, column, step, i;
	double mult;

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			C.p[row][column] = p[row][column];

	for (i = 0; i < rows; i++)
		C.p[i][rows + i] = 1;

	for (step = 0; step < rows - 1; step++)
	{
		for (row = step + 1; row < rows; row++)
		{
			mult = C.p[row][step] / C.p[step][step];
			for (column = step; column < 2 * rows; column++)
				C.p[row][column] -= mult * C.p[step][column];
		}
	}

	for (step = 1; step <= rows - 1; step++)
	{
		for (row = rows - step - 1; row >= 0; row--)
		{
			mult = C.p[row][rows - step] / C.p[rows - step][rows - step];
			for (column = rows; column < 2 * rows; column++)
				C.p[row][column] -=
					mult * C.p[rows - step][column];
		}
	}

	for (row = 0; row < rows; row++)
		for (column = 0; column < rows; column++)
			B.p[row][column] = C.p[row][rows + column] / C.p[row][row];
	return B;
}

matrix::matrix(int r, int c)
{

	// matrix  constructor

	// cout<<"mat constructor called !!!!!"<<r<<" "<<c<<endl;

	int i, j;
	rows = r;
	columns = c;

	p = new double *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		p[i] = new double[columns];
		// assert(p[i] != NULL);
	}

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			p[i][j] = 0;
	pf = NULL;
}

matrix::matrix(matrix &mf) : frac()
{

	// matrix  constructor

	cout << " ############  mat FRACTION  constructor called !!!!!" << endl;

	int i, j;
	rows = mf.rows;
	columns = mf.columns;

	pf = new frac *[rows];
	// assert(p != NULL);								//analyse
	for (i = 0; i < rows; i++)
	{
		pf[i] = new frac[columns];
		// assert(p[i] != NULL);
	}
	cout << "frac alloc sucessful" << endl;

	if (mf.p == NULL)
	{

		for (i = 0; i < mf.rows; i++)
		{
			for (j = 0; j < mf.columns; j++)
			{
				cout << i << " % " << j << endl;
				pf[i][j] = mf.pf[i][j];
			}
		}
	}
	else
	{
		for (i = 0; i < mf.rows; i++)
		{
			for (j = 0; j < mf.columns; j++)
			{
				cout << i << " % " << j << endl;
				pf[i][j].num = mf.p[i][j];
			}
		}
	}
	p = NULL;
	cout << "mat FRACTION  constructor LEFT  " << endl;
}

void matrix ::imat()
{

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			if (i == j)
			{
				p[i][j] = 1;
			}
			else
			{
				p[i][j] = 0;
			}
		}
	}
}

void matrix ::cost()
{

	cout << "cost called !!" << endl;

	for (int i = 0; i < rows; i++)
	{
		p[i][1] = 0;
	}
}

matrix matrix ::mx(matrix &B)
{
	cout << "mx fx called!!!!\n";

	float sum;
	int i, j, k, n, m, l;
	n = rows;
	m = columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix operator*(matrix &A, matrix &B)
{
	double sum;
	int i, j, k, n, m, l;
	n = A.rows;
	m = A.columns;
	l = B.columns;
	matrix tempMatrix(n, l);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < l; j++)
		{
			sum = 0;
			for (k = 0; k < m; k++)
				sum = sum + A.p[i][k] * B.p[k][j];
			tempMatrix.p[i][j] = sum;
		}
	}

	return tempMatrix;
}

matrix matrix ::sm(matrix &A, matrix &B)
{

	cout << "sm called!!!!!" << endl;

	int i, j, n, m;

	n = A.rows;
	m = A.columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = B.p[i][j] + A.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

matrix matrix ::sb(matrix &B)
{

	cout << "sb called!!!!!" << endl;

	int i, j, n, m;

	n = rows;
	m = columns;
	matrix tempMatrix(n, m);
	tempMatrix.dispstat();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{

			tempMatrix.p[i][j] = p[i][j] - B.p[i][j];
			cout << i << " " << j << endl;
		}
	}

	return tempMatrix;
};

void matrix ::in()
{
	cout << "matrix fx called!!!!\n";

	delete[] p;

	cout << "Enter the number of Rows and Columns: ";
	cin >> rows >> columns;
	p = new double *[rows];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new double[columns];
	}

	cout << "Enter the " << rows << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << columns << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::inset(int rset, int cset)
{
	cout << "matrix inset called!!!!\n";

	delete[] p;

	p = new double *[rset];

	for (int i = 0; i < rows; i++)
	{
		p[i] = new double[cset];
	}

	cout << "Enter the " << rset << " rows of the matrix:" << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << "Enter the " << cset << " elements of row number " << i + 1 << ": ";
		for (int j = 0; j < columns; j++)
		{

			cin >> p[i][j];
			cout << p[i][j] << endl;
		}
	}
}

void matrix ::disp()
{
	cout << "disp called!!!!!\n";

	for (int i = 0; i < rows; i++)
	{

		for (int j = 0; j < columns; j++)
		{
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
}

void matrix ::dispstat()
{
	cout << "dispstat called!!!!!\n";
	cout << rows << " " << columns << endl;
}

int matrix ::initializer()
{
	cout << "IN ITER !!!!" << endl;

	int nz, ndesc;

	cout << "enter max mat" << endl;
	cout << "enter no of var" << endl;
	cin >> nz;
	matrix matz(1, nz);
	matz.inset(1, nz);
	matz.disp();
	global_Z = matz;

	cout << "no of desc" << endl;
	cin >> ndesc;
	matrix matv(ndesc, (nz + 1));
	matv.inset(ndesc, (nz + 1));
	matv.disp();
	global_V = matv;

	matrix matB(ndesc, ndesc);
	matB.imat();
	matB.disp();
	global_B = matB;

	head();
	int check = iterator(matz, matv, matB);
	if (check == -1)
	{
		return 0;
	}
}

int matrix ::iterator(matrix &z, matrix &v, matrix &B)
{
	cout << "IN ITERRATOR !!!!" << endl;

	matrix Xb;
	matrix DEL;
	matrix COST;
	matrix prc_1;
	matrix prc_2;
	matrix prc_3;

	matrix Cj;
	matrix COST_trn;
	matrix temp_1;

	matrix Aj(v.rows, v.columns - 1);

	for (int i = 0; i < v.rows; i++)
	{
		for (int j = 0; j < v.columns - 1; j++)
		{
			cout << i << " " << j << endl;

			Aj.p[i][j] = v.p[i][j];
		}
	}
	Aj.disp();

	matrix Binv;
	Binv = B.inverse();
	Binv.disp();

	matrix b(v.rows, 1);
	for (int i = 0; i < v.rows; i++)
	{
		b.p[i][0] = v.p[i][(v.columns - 1)];
	}
	b.disp();

	global_b = b;

	matrix cB(v.rows, 1);
	cB.cost();
	cB.disp();
	COST_trn = cB.trn();
	COST_trn.disp();
	global_Cbt = COST_trn;

	Cj = z;
	global_Cj = Cj;
	Cj.disp();

	int opti_flag = 0;

	for (int j = 0; j < 1; j++)
	{

		cout << "ITERATION" << endl
			 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

		Binv.disp();
		b.disp();
		COST_trn.disp();

		Xb = Binv.mx(b);
		Xb.disp();
		cout << "######" << endl;

		prc_1 = COST_trn.mx(Binv);
		prc_1.disp();

		cout << "######" << endl;

		prc_2 = prc_1.mx(Aj);
		prc_2.disp();

		cout << "######" << endl;
		DEL = prc_2.sb(Cj);
		DEL.disp();

		cout << "######" << endl;
		opti_flag = temp_1.comp(DEL);
		global_flag = opti_flag;

		cout << "opti_flag = " << opti_flag << endl;
		DEL.disp();

		cout << "######" << endl;
		matrix table(v.rows, B.columns + 2);

		matrix xn(v.rows, 1);
		for (int k = 0; k < v.rows; k++)
		{
			xn.p[k][0] = v.p[k][opti_flag];
		}
		matrix Xn;
		Xn = Binv.mx(xn);

		table = prc_3.merger(Xb, Binv, Xn);
		table.disp();

		matrix Atab(v.rows, v.columns - 1 + v.rows);

		Atab = v.A_indexer(1);

		cout << "######" << endl;

		int while_flag = 1;

		while (while_flag == 1)
		{

			cout << "in WHILE @@@@@@@@@@@@@@@@@@" << endl;

			int rat = Xb.ratio(Xn);
			global_ratio = rat;
			cout << "ratio is " << rat << endl;
			matrix ftab(table);
			ftab.frac_disp();

			head_ch(opti_flag, rat);
			int check = ftab_solver(ftab, rat);
			if (check == -1)
			{
				return -1;
			}

			cout << "######### END OF WHILE  ############" << endl;

			while_flag = 0;
		}
	}
}
int matrix ::ratio(matrix &X)
{

	cout << "in ratio !!!!!!!!!!1" << endl;

	matrix rat(X.rows, 1);
	int rat_arr[X.rows];

	for (int i = 0; i < X.rows; i++)
	{
		rat.p[i][0] = (p[i][0]) / (X.p[i][0]);
		rat_arr[i] = rat.p[i][0];
		cout << "RATIO ======" << rat.p[i][0] << rat_arr[i] << endl;
	}

	//	//int size = sizeof(rat.p)/sizeof(rat.p[0][0]);
	//
	//
	//	for(int j = 0; j < rat.rows ; j++)
	//	{
	//		if(rat.p[j][0] < rat.p[0][0] && rat.p[j][0] > 0 )
	//			index = j;
	//	}

	// pos = rat_comp(rat_arr);
	int x = 0, c = 0;
	int pos = 0;
	cout << rat.rows << endl;
	for (x = 0; x < rat.rows; x++)
	{
		if (rat.p[x][0] < 0)
		{
			c++;
			cout << c << x << endl;
		}
	}
	cout << "here" << endl;

	for (int j = c; j < rat.rows; j++)
	{
		cout << c << " $ " << j << endl;
		if (rat.p[j][0] <= rat.p[c][0] && rat.p[j][0] > 0)
		{
			c = j;
			pos = j;
		}
	}
	cout << pos;

	return pos;
}

matrix matrix ::merger(matrix &Xb, matrix &Bin, matrix &Xn)
{
	cout << "merger called !!!!!!!!!!!!1" << endl;

	matrix merge_mat(Xn.rows, Bin.columns + 2);
	//	matrix full_tab(Xn.rows,Bin.columns+ Xn.columns);
	//
	//	for (int i =0 ; i < Xb.rows ; i++){
	//		merge_mat.p[i][0] = Xb.p[i][0];
	//		for (int j = 0; j < Bin.columns ; j++ ){
	//			merge_mat.p[i][1+j] = Bin.p[i][j] ;
	//		}
	//		for (int j = 0; j < Xn.columns ; j++ ){
	//
	//			merge_mat.p[i][Bin.columns+1] = Xn.p[i][j];
	//
	//		}
	//	}
	//	global_full_tab = full_tab;

	for (int i = 0; i < Xb.rows; i++)
	{
		merge_mat.p[i][0] = Xb.p[i][0];
		for (int j = 0; j < Bin.columns; j++)
		{
			merge_mat.p[i][1 + j] = Bin.p[i][j];
		}

		merge_mat.p[i][Bin.columns + 1] = Xn.p[i][0];
	}
	return merge_mat;
}

matrix matrix ::A_indexer(int select)
{
	cout << "in INDEXER @@@@@@@@@" << endl;

	if (select == 1)
	{

		cout << "MIN CONFIG" << endl;

		matrix A(rows, columns - 1 + rows);
		A.disp();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < columns - 1; j++)
			{
				A.p[i][j] = p[i][j];
			}
			for (int k = 0; k < rows; k++)
			{

				if (i == k)
				{
					A.p[i][columns - 1 + k] = 1;
				}
				else
				{
					A.p[i][columns - 1 + k] = 0;
				}
			}
		}
		A.disp();
		global_A = A;

		return A;
	}
}

int matrix ::comp(matrix &delj)
{
	cout << "in comp !!!!" << endl;
	int minpos = 0;
	// float minelem;
	int nflag = 0;

	for (int i = 0; i < delj.columns; i++)
	{
		if (delj.p[0][i] < 0)
		{
			nflag = 1;
			cout << delj.p[0][i] << endl;

			break;
		}
	}
	if (nflag == 1)
	{
		minpos = minarray(delj);
		cout << "index = " << minpos << endl;
	}
	else
	{

		cout << " ALL DELj ARE POSITIVE " << endl;
		result();
	}
	return minpos;
}
int matrix ::result()
{
	cout << " optimality reached " << endl;

	double Z = 0;

	global_Xb.disp();

	for (int i = 0; i < global_Xb.rows; i++)
	{

		cout << "X" << i + 1 << " -> " << global_Xb.p[i][0] << endl;
	}

	for (int i = 0; i < global_Z.columns; i++)
	{
		if (global_head.p[0][i] == 0)
			Z = Z + ((global_Xb.p[0][i]) * global_Z.p[0][i]);
	}
	cout << " MAX Z -> " << Z << endl;
	run = 0;

	return 0;
}

int matrix ::minarray(matrix &min)
{
	cout << "in minarray !!!!" << endl;

	int index = 0;
	//	int size = sizeof(min)/sizeof(min.p , min.p[0][0]);

	for (int i = 0; i < min.columns; i++)
	{

		if (min.p[0][i] < min.p[0][index])
		{
			index = i;
		}
	}
	return index;
}

//----------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------

int main()
{
	cout << "start\n";

	matrix MAT;
	run = 1;

	MAT.initializer();

	cout << "main end\n";

	return 0;
}
