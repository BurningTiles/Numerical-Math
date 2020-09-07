/**
 * Author  : BurningTiles
 * Created : 2020-09-07 17:50:58
 * Link    : GitHub.com/BurningTiles
 * Header  : Numerical Math
**/

#ifndef NMATH_H
#define NMATH_H

#include "function.h"
#include <iostream>
#include <iomanip>

#ifndef accuracy
#define accuracy 0.00001
#endif

#ifndef error
#define error 0.000000001
#endif

using namespace std;




//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Bisection{
	 private:
	 	double a,b;
	 public:
	 	bool input_range(double,double);
	 	void print();
	 	bool main();
};

bool Bisection::input_range(double x, double y){
	double fx=f(x) , fy=f(y);
	if ( fx == 0 ) {
		cout << "Root is " << x << endl;
		return 1;
	}
	else if( fy == 0 ){
		cout << "Root is " << y << endl;
		return 1;
	}
	else if( ( fx < 0 && fy < 0 ) || ( fx > 0 && fy > 0 ) ){
		cerr << "There is no root in given range." << endl;
		return 1;
	}
	else{
		if( fx < 0 ){
			a=x; b=y;
		}
		else {
			a=y; b=x;
		}
		return 0;
	}
}

void Bisection::print(){
	double average=0,prev_average,value;
	int i = 1;
	cout << "Iteration   Approx_Root   Value_of__f(x)" << endl;
	do {
		prev_average = average;
		average = (a+b)/2;
		value = f(average);
		if ( value < 0 ) 
			a = average;
		else if ( value > 0 )
			b = average;
		else{
			cout << endl << endl << "Root of given function is " << average << endl;
			return;
		}
		cout << setw(9) << left << i;
		cout.setf(ios::floatfield,ios::fixed);
		cout << setprecision(8)  << setw(14) << right << average;
		cout << setprecision(10) << setw(17) << right << value << endl;
		i++;
	}while( mod(average,prev_average) > accuracy);
	cout << endl << endl << "Root of given function is " << setprecision(8) << average << endl;
}

bool Bisection::main(){
	double temp1,temp2;
	cout << endl << "Enter range in which root lies : ";
	if( !(cin >> temp1 >> temp2)){
		cerr << "Enter valid input." << endl;
		return 1;
	}

	if(input_range(temp1,temp2))
		return 0;
	print();
	return 0;
}







//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Correlation
{
private:
	double *x, *y;
	int length;

public:
	bool input_data(double[], double[], int);
	void print();
	bool main();
};

bool Correlation::input_data(double tempx[], double tempy[], int templ)
{
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i = 0; i < length; i++)
	{
		*(x + i) = tempx[i];
		*(y + i) = tempy[i];
	}
	return 0;
}

void Correlation::print()
{
	double tx = 0, ty = 0, txy = 0, txx = 0, tyy = 0, r;

	for (int i = 0; i < length; i++)
	{
		tx += *(x + i);
		ty += *(y + i);
		txy += *(x + i) * *(y + i);
		txx += *(x + i) * *(x + i);
		tyy += *(y + i) * *(y + i);
	}
	r = (txy - tx * ty / length) / (sqrt(txx - tx * tx / length) * sqrt(tyy - ty * ty / length));
	cout.setf(ios::floatfield, ios::fixed);
	cout << endl
		 << setprecision(8);
	cout << endl
		 << "Total of x  : " << setw(12) << tx << endl
		 << "Total of y  : " << setw(12) << ty << endl
		 << "Total of xy : " << setw(12) << txy << endl
		 << "Total of xx : " << setw(12) << txx << endl
		 << "Total of yy : " << setw(12) << tyy << endl
		 << endl
		 << "Correlation coefficient (r) = " << r << endl;
}

bool Correlation::main()
{
	unsigned int templ;
	cout << endl
		 << "Enter length of data : ";
	if (!(cin >> templ) || templ < 2)
	{
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ], tempy[templ];
	cout << endl
		 << "Enter values of x : ";
	for (int i = 0; i < templ; i++)
	{
		if (!(cin >> tempx[i]))
		{
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}

	cout << endl
		 << "Enter values of y corresponding to x : ";
	for (int i = 0; i < templ; i++)
	{
		if (!(cin >> tempy[i]))
		{
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}

	if (input_data(tempx, tempy, templ))
		return 0;
	print();
	return 0;
}







//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Euler{
	 private:
	 	double *x,*y,h,in_x;
	 	int n;
	 public:
	 	bool input_data(double,double,double);
	 	void result();
	 	bool main();
};

bool Euler::input_data(double tempx0, double tempy0, double tempin){
	x = new double[n+1];
	y = new double[n+1];
	*x = tempx0;
	*y = tempy0;
	in_x = tempin;
	for ( int i=0 ; i<=n ; i++)
		*(x+i) = *x+i*h;
	return 0;
}

void Euler::result(){
	cout.setf( ios::floatfield , ios::fixed );
	cout << endl << setprecision(8) << "values_of_x     value_of_y/f(x)" << endl
	             << setw(11) << *x << setw(20) << *y  << endl;
	for ( int i=1 ; i<=n ; i++ ){
		*(y+i) = *(y+i-1) + h * fxy( *(x+i-1) , *(y+i-1) );
		cout << setw(11) << *(x+i) << setw(20) << *(y+i) << endl;
	}
}

bool Euler::main(){
	double tempx0,tempy0,tempin;
	bool choice;
	cout << endl << "Enter x0 and y0 : ";
	if( !(cin >> tempx0 >> tempy0)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Enter value of x : ";
	if( !(cin >> tempin) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Chose steps h(0) or n(1) ? : ";
	if ( !(cin >> choice) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if ( choice ) {
		cout << endl << "Enter value of n (step number) : ";
		if( !(cin >> n) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		h = ( tempin - tempx0 ) / n;
	}
	else {
		cout << endl << "Enter value of h (step size) : ";
		if( !(cin >> h) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		n = round(( tempin - tempx0 ) / h);
	}
	
	if ( mod( (tempx0+n*h) , tempin ) > error ){
		cerr << "y(x) can\'t be found using given data.";
		return 1;
	}
	
	if( input_data(tempx0,tempy0,tempin) )
		return 0;
	result();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Falseposition{
	 private:
	 	double a,b,fa,fb;
	 public:
	 	bool input_range(double,double);
	 	void print();
	 	bool main();
};

bool Falseposition::input_range(double x, double y){
	double fx=f(x) , fy=f(y);
	if ( fx == 0 ) {
		cout << "Root is " << x << endl;
		return 1;
	}
	else if( fy == 0 ){
		cout << "Root is " << y << endl;
		return 1;
	}
	else if( ( fx < 0 && fy < 0 ) || ( fx > 0 && fy > 0 ) ){
		cerr << "There is no root in given range." << endl;
		return 1;
	}
	else{
		if( fx < 0 ){
			a=x; b=y; fa=fx; fb=fy;
		}
		else {
			a=y; b=x; fa=fy; fb=fx;
		}
		return 0;
	}
}

void Falseposition::print(){
	double position=0,prev_position,value;
	int i = 1;
	cout << "Iteration   Approx_Root   Value_of__f(x)" << endl;
	do {
		prev_position = position;
		position = ((a*fb)-(b*fa))/(fb-fa);
		value = f(position);
		if ( value < 0 ) {
			a = position;
			fa = value;
		}
		else if ( value > 0 ) {
			b = position;
			fb = value;
		}
		else{
			cout << endl << endl << "Root of given function is " << position << endl;
			return;
		}
		cout << setw(9) << left << i;
		cout.setf(ios::floatfield,ios::fixed);
		cout << setprecision(8)  << setw(14) << right << position;
		cout << setprecision(10) << setw(17) << right << value << endl;
		i++;
	}while( mod(position,prev_position) > accuracy);
	cout << endl << endl << "Root of given function is " << setprecision(8) << position << endl;
}

bool Falseposition::main(){
	double temp1,temp2;
	cout << endl << "Enter range in which root lies : ";
	if( !(cin >> temp1 >> temp2)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if(input_range(temp1,temp2))
		return 0;
	print();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Gaussjacobi{
	 private:
	 	double data[3][4];
	 public:
	 	bool input_data(double[][4]);
	 	void print();
	 	bool main();
};

bool Gaussjacobi::input_data(double temp[][4]){
	bool flag[3]={0,0,0};
	for (int i=0 ; i<3 ; i++){
		int max = (temp[i][0] > temp[i][1])? (temp[i][0] > temp[i][2])? 0 : 2 : (temp[i][1] > temp[i][2])? 1 : 2 ;
		for (int j=0 ; j<4 ; j++)
			data[max][j] = temp[i][j];
		flag[max] = 1;
	}
	if ( (flag[0]+flag[1]+flag[2]) != 3 ) {
		cerr<<"Given system is not diagonally dominant."<<endl;
		return 1;
	}
	return 0;
}

void Gaussjacobi::print(){
	double ans[3]={0,0,0} ,prev_ans[3];
	int i = 1;
	
	cout << endl << "Iteration   Approx_Root_x   Approx_Root_y   Approx_Root_z" << endl
	             << "0              0.00000000      0.00000000      0.00000000" << endl;
	
	do {
		for (int i=0 ; i<3 ; i++)
			prev_ans[i]=ans[i];
		
		ans[0] = (data[0][3] - (data[0][1]*prev_ans[1]) - (data[0][2]*prev_ans[2]) ) / data[0][0];
		ans[1] = (data[1][3] - (data[1][0]*prev_ans[0]) - (data[1][2]*prev_ans[2]) ) / data[1][1];
		ans[2] = (data[2][3] - (data[2][1]*prev_ans[1]) - (data[2][0]*prev_ans[0]) ) / data[2][2];
		
		cout << setw(9) << left << i;
		cout.setf(ios::floatfield,ios::fixed);
		cout << setprecision(8)  << setw(16) << right << ans[0] << setw(16) << ans[1] << setw(16) << ans[2] << endl;
		i++;
	}while( mod(ans[0],prev_ans[0]) > accuracy || mod(ans[0],prev_ans[0]) > accuracy || mod(ans[0],prev_ans[0]) > accuracy );
	
	cout << endl << endl << "Roots of given equations are : " << endl << setprecision(8) << "x = " << ans[0] << "   y = " << ans[1] << "   z = " << ans[2] << endl;
}

bool Gaussjacobi::main(){
	double t[3][4];
	cout << endl << "Input format : " << endl << "If equation is ax + by + cz = d " << endl << "Then input will be a >b >c >d." <<endl;
	for(int i=0;i<3;i++){
		cout << endl << "Enter coefficients of equation" << i+1 << " : ";
		if( !(cin >> t[i][0] >> t[i][1] >> t[i][2] >> t[i][3] ) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
	}
	
	if( input_data(t) )
		return 0;
	print();
	return 0;
}







//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Gaussseidel{
	 private:
	 	double data[3][4];
	 public:
	 	bool input_data(double[][4]);
	 	void print();
	 	bool main();
};

bool Gaussseidel::input_data(double temp[][4]){
	bool flag[3]={0,0,0};
	for (int i=0 ; i<3 ; i++){
		int max = (temp[i][0] > temp[i][1])? (temp[i][0] > temp[i][2])? 0 : 2 : (temp[i][1] > temp[i][2])? 1 : 2 ;
		for (int j=0 ; j<4 ; j++)
			data[max][j] = temp[i][j];
		flag[max] = 1;
	}
	if ( (flag[0]+flag[1]+flag[2]) != 3 ) {
		cerr<<"Given system is not diagonally dominant."<<endl;
		return 1;
	}
	return 0;
}

void Gaussseidel::print(){
	double ans[3]={0,0,0} ,prev_ans[3];
	int i = 1;
	
	cout << endl << "Iteration   Approx_Root_x   Approx_Root_y   Approx_Root_z" << endl
	             << "0              0.00000000      0.00000000      0.00000000" << endl;
	
	do {
		for (int i=0 ; i<3 ; i++)
			prev_ans[i]=ans[i];
		
		ans[0] = (data[0][3] - (data[0][1]*prev_ans[1]) - (data[0][2]*prev_ans[2]) ) / data[0][0];
		ans[1] = (data[1][3] - (data[1][0]*     ans[0]) - (data[1][2]*prev_ans[2]) ) / data[1][1];
		ans[2] = (data[2][3] - (data[2][1]*     ans[1]) - (data[2][0]*     ans[0]) ) / data[2][2];
		
		cout << setw(9) << left << i;
		cout.setf(ios::floatfield,ios::fixed);
		cout << setprecision(8)  << setw(16) << right << ans[0] << setw(16) << ans[1] << setw(16) << ans[2] << endl;
		i++;
	}while( mod(ans[0],prev_ans[0]) > accuracy || mod(ans[0],prev_ans[0]) > accuracy || mod(ans[0],prev_ans[0]) > accuracy );
	
	cout << endl << endl << "Roots of given equations are : " << endl << setprecision(8) << "x = " << ans[0] << "   y = " << ans[1] << "   z = " << ans[2] << endl;
}

bool Gaussseidel::main(){
	double t[3][4];
	cout << endl << "Input format : " << endl << "If equation is ax + by + cz = d " << endl << "Then input will be a >b >c >d." <<endl;
	for(int i=0;i<3;i++){
		cout << endl << "Enter coefficients of equation" << i+1 << " : ";
		if( !(cin >> t[i][0] >> t[i][1] >> t[i][2] >> t[i][3] ) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
	}
	
	if( input_data(t) )
		return 0;
	print();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Lagrange{
	 private:
	 	double *x,*y;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void get();
	 	bool main();
};

bool Lagrange::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Lagrange::get(){
	double ans , tempu , tempd , in_x;
	
	cout << endl << "For exit insert eof." <<endl;
	while ( cout << endl << "Enter value of x : " && cin >> in_x ) {
		ans = 0;
		for ( int i=0 ; i<length ; i++ ){
			tempu = *(y+i);
			for ( int j=0 ; j<length ; j++ ){
				if (i==j)
					continue;
				tempu *= in_x - *(x+j);
			}
			tempd = 1;
			for ( int j=0 ; j<length ; j++ ){
				if (i==j)
					continue;
				tempd *= *(x+i) - *(x+j);
			}
			ans += tempu / tempd;
		}
		cout << "Value of y at x : " << ans << endl;
	}
}

bool Lagrange::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ < 2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ] , tempy[templ];
	
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	get();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Linefitting{
	 private:
	 	double *x,*y;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void print();
	 	bool main();
};

bool Linefitting::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Linefitting::print(){
	double tx=0,ty=0,txy=0,txx=0,a,b;
	
	for (int i=0 ; i<length ; i++) {
		tx  += *(x+i);
		ty  += *(y+i);
		txy += *(x+i) * *(y+i);
		txx += *(x+i) * *(x+i);
	}
	b = ( length * txy - tx * ty ) / ( length * txx - tx * tx );
	a = ty / length - b * tx / length;
	
	cout.setf(ios::floatfield,ios::fixed);
	cout << setprecision(6);
	cout << endl << "Total of x  : " << setw(12) << tx  << endl
	             << "Total of y  : " << setw(12) << ty  << endl
	             << "Total of xy : " << setw(12) << txy << endl
	             << "Total of xx : " << setw(12) << txx << endl
	     << endl << "Normal equations are : " << endl
	             << setw(12) << ty  << " = " << setw(12) << length << "a + " << setw(12) << tx  << endl
	             << setw(12) << txy << " = " << setw(12) << tx     << "a + " << setw(12) << txx << endl
	     << endl << "Required straight line is : " << endl
	             << "y = " << a << " + " << b << "x" << endl;
}

bool Linefitting::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ<2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ],tempy[templ];
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	print();
	return 0;
}





//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Modifiedeuler{
	 private:
	 	double *x,*y,h,in_x;
	 	int n;
	 public:
	 	bool input_data(double,double,double);
	 	void result();
	 	bool main();
};

bool Modifiedeuler::input_data(double tempx0, double tempy0, double tempin){
	x = new double[n+1];
	y = new double[n+1];
	*x = tempx0;
	*y = tempy0;
	in_x = tempin;
	for ( int i=0 ; i<=n ; i++)
		*(x+i) = *x+i*h;
	return 0;
}

void Modifiedeuler::result(){
	double tempvalue;
	cout.setf( ios::floatfield , ios::fixed );
	cout << endl << setprecision(8) << "values_of_x     value_of_y/f(x)" << endl
	             << setw(11) << *x << setw(20) << *y  << endl;
	for ( int i=1 ; i<=n ; i++ ){
		tempvalue = *(y+i-1) + h * fxy( *(x+i-1) , *(y+i-1) );
		*(y+i)    = *(y+i-1) + h * ( fxy(*(x+i-1),*(y+i-1)) + fxy(*(x+i),tempvalue) ) / 2;
		cout << setw(11) << *(x+i) << setw(20) << *(y+i) << endl;
	}
}

bool Modifiedeuler::main(){
	double tempx0,tempy0,tempin;
	bool choice;
	cout << endl << "Enter x0 and y0 : ";
	if( !(cin >> tempx0 >> tempy0)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Enter value of x : ";
	if( !(cin >> tempin) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Chose steps h(0) or n(1) ? : ";
	if ( !(cin >> choice) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if ( choice ) {
		cout << endl << "Enter value of n (step number) : ";
		if( !(cin >> n) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		h = ( tempin - tempx0 ) / n;
	}
	else {
		cout << endl << "Enter value of h (step size) : ";
		if( !(cin >> h) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		n = round(( tempin - tempx0 ) / h);
	}
	
	if ( mod( (tempx0+n*h) , tempin ) > error ){
		cerr << "y(x) can\'t be found using given data.";
		return 1;
	}
	
	if( input_data(tempx0,tempy0,tempin) )
		return 0;
	result();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Newtonbackward{
	 private:
	 	double *x,*y,h;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void get();
	 	bool main();
	 	inline double countp(double,int);
};

bool Newtonbackward::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Newtonbackward::get(){
	double diff_table[length][length] , in_x , ans , p;
	
	for ( int i=0 ; i<length ; i++ )
		diff_table[i][0] = *(y+i);
		
	for ( int j=1 ; j<length ; j++ ) 
		for ( int i=0 ; i<length-j ; i++ )
			diff_table[i][j] = diff_table[i+1][j-1] - diff_table[i][j-1];
	
	cout << endl << "For exit insert eof." <<endl;
	while ( cout << endl << "Enter value of x : " && cin >> in_x ) {
		ans = *(y+length-1);
		p = (in_x - *(x+length-1)) / h;
		for ( int i=1 ; i<length ; i++ )
			ans += countp( p , i-1 ) * diff_table[length-i-1][i] / double(fact(i));
		cout << "Value of y at x : " << ans << endl;
	}
}

inline double Newtonbackward::countp(double p,int n) {
	return (n<1)? p : double( (p+n) * countp (p,n-1) );
}

bool Newtonbackward::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ < 2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ] , tempy[templ];
	
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	h = tempx[1]-tempx[0];
	for(int i=templ-1 ; i>0 ; i--) {
		if( mod(tempx[i]-tempx[i-1],h) > error ) {
			cerr << endl << "Sorry given data is not sequential." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	get();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Newtondivided{
	 private:
	 	double *x,*y;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void get();
	 	bool main();
	 	inline double countx(double,int);
};

bool Newtondivided::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Newtondivided::get(){
	double diff_table[length][length] , in_x , ans , p;
	
	for ( int i=0 ; i<length ; i++ )
		diff_table[i][0] = *(y+i);
		
	for ( int j=1 ; j<length ; j++ ) 
		for ( int i=0 ; i<length-j ; i++ )
			diff_table[i][j] = ( diff_table[i+1][j-1] - diff_table[i][j-1] ) / ( *(x+i+j) - *(x+i) );
	
	cout << endl << "For exit insert eof." <<endl;
	while ( cout << endl << "Enter value of x : " && cin >> in_x ) {
		ans = *y;
		for ( int i=1 ; i<length ; i++ )
			ans += countx( in_x , i-1 ) * diff_table[0][i];
		cout << "Value of y at x : " << ans << endl;
	}
}

inline double Newtondivided::countx(double dx,int n) {
	return (n<0)? 1 : double( (dx-*(x+n)) * countx (dx,n-1) );
}

bool Newtondivided::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ < 2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ] , tempy[templ];
	
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	get();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Newtonforward{
	 private:
	 	double *x,*y,h;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void get();
	 	bool main();
	 	inline double countp(double,int);
};

bool Newtonforward::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Newtonforward::get(){
	double diff_table[length][length] , in_x , ans , p;
	
	for ( int i=0 ; i<length ; i++ )
		diff_table[i][0] = *(y+i);
		
	for ( int j=1 ; j<length ; j++ ) 
		for ( int i=0 ; i<length-j ; i++ )
			diff_table[i][j] = diff_table[i+1][j-1] - diff_table[i][j-1];
	
	cout << endl << "For exit insert eof." <<endl;
	while ( cout << endl << "Enter value of x : " && cin >> in_x ) {
		ans = *y;
		p = (in_x - *x) / h;
		for ( int i=1 ; i<length ; i++ )
			ans += countp( p , i-1 ) * diff_table[0][i] / double(fact(i));
		cout << "Value of y at x : " << ans << endl;
	}
}

inline double Newtonforward::countp(double p,int n) {
	return (n<1)? p : double( (p-n) * countp (p,n-1) );
}

bool Newtonforward::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ < 2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ] , tempy[templ];
	
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	h = tempx[1]-tempx[0];
	for(int i=templ-1 ; i>0 ; i--) {
		if( mod(tempx[i]-tempx[i-1],h) < error ) {
			cerr << endl << "Sorry given data is not sequential." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	get();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Newtonraphson{
	 private:
	 	double a,fa,f1a;
	 public:
	 	bool input_approx(double);
	 	void print();
	 	bool main();
};

bool Newtonraphson::input_approx(double x){
	a=x; fa=f(x); f1a=f1(x);
	if ( fa == 0 ) {
		cout << "Root is " << x << endl;
		return 1;
	}
	return 0;
}

void Newtonraphson::print(){
	double root=0,value,prev_root;
	int i = 1;
	cout << "Iteration   Approx_Root   Value_of__f(x)" << endl;
	do {
		prev_root = root;
		root = (a - fa/f1a);
		value = f(root);
		if ( value == 0 ) {
			cout << endl << endl << "Root of given function is " << root << endl;
			return;
		}
		
		a=root; fa=value; f1a=f1(root);
		
		cout << setw(9) << left << i;
		cout.setf(ios::floatfield,ios::fixed);
		cout << setprecision(8)  << setw(14) << right << root;
		cout << setprecision(10) << setw(17) << right << value << endl;
		i++;
	}while( mod(root,prev_root) > accuracy);
	cout << endl << endl << "Root of given function is " << setprecision(8) << root << endl;
}

bool Newtonraphson::main(){
	double temp1;
	cout << endl << "Enter approximated root : ";
	if( !(cin >> temp1)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if(input_approx(temp1))
		return 0;
	print();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Regressionlines{
	 private:
	 	double *x,*y;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void print();
	 	bool main();
};

bool Regressionlines::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Regressionlines::print(){
	double tx=0,ty=0,txy=0,txx=0,tyy=0,r,bxy,byx;
	
	for (int i=0 ; i<length ; i++) {
		tx  += *(x+i);
		ty  += *(y+i);
		txy += *(x+i) * *(y+i);
		txx += *(x+i) * *(x+i);
		tyy += *(y+i) * *(y+i);
	}
	bxy = ( txy - tx * ty / length ) / ( tyy - ty * ty / length );
	byx = ( txy - tx * ty / length ) / ( txx - tx * tx / length );
	r   = sqrt( bxy * byx );
	
	cout.setf(ios::floatfield,ios::fixed);
	cout << endl << setprecision(8);
	cout << endl << "Total of x  : " << setw(12) << tx  << endl
	             << "Total of y  : " << setw(12) << ty  << endl
	             << "Total of xy : " << setw(12) << txy << endl
	             << "Total of xx : " << setw(12) << txx << endl
	             << "Total of yy : " << setw(12) << tyy << endl << endl
	             << "Regression coeffient of x on y (bxy) : " << setw(12) << bxy << endl
	             << "Regression coeffient of y on x (byx) : " << setw(12) << byx << endl << endl
	             << "Regression line of x on y :   x = " << setw(12) << bxy << "y + " << (tx / length - ty * bxy / length) <<endl
	             << "Regression line of y on x :   y = " << setw(12) << byx << "x + " << (ty / length - tx * byx / length) <<endl
	     << endl << "Correlation coefficient (r) = " << r << endl;
}

bool Regressionlines::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ<2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ],tempy[templ];
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	print();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Rkmethod2{
	 private:
	 	double *x,*y,h,in_x;
	 	int n;
	 public:
	 	bool input_data(double,double,double);
	 	void result();
	 	bool main();
};

bool Rkmethod2::input_data(double tempx0, double tempy0, double tempin){
	x = new double[n+1];
	y = new double[n+1];
	*x = tempx0;
	*y = tempy0;
	in_x = tempin;
	for ( int i=0 ; i<=n ; i++)
		*(x+i) = *x+i*h;
	return 0;
}

void Rkmethod2::result(){
	double k,k1,k2;
	cout.setf( ios::floatfield , ios::fixed );
	cout << endl << setprecision(8) << "values_of_x     value_of_y/f(x)" << endl
	             << setw(11) << *x << setw(20) << *y  << endl;
	for ( int i=1 ; i<=n ; i++ ){
		k1 = fxy( *(x+i-1)   , *(y+i-1)      );
		k2 = fxy( *(x+i-1)+h , *(y+i-1)+h*k1 );
		k  = (k1 + k2) / 2;
		*(y+i) = *(y+i-1) + h * k;
		cout << setw(11) << *(x+i) << setw(20) << *(y+i) << endl;
	}
}

bool Rkmethod2::main(){
	double tempx0,tempy0,tempin;
	bool choice;
	cout << endl << "Enter x0 and y0 : ";
	if( !(cin >> tempx0 >> tempy0)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Enter value of x : ";
	if( !(cin >> tempin) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Chose steps h(0) or n(1) ? : ";
	if ( !(cin >> choice) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if ( choice ) {
		cout << endl << "Enter value of n (step number) : ";
		if( !(cin >> n) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		h = ( tempin - tempx0 ) / n;
	}
	else {
		cout << endl << "Enter value of h (step size) : ";
		if( !(cin >> h) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		n = round(( tempin - tempx0 ) / h);
	}
	
	if ( mod( (tempx0+n*h) , tempin ) > error ){
		cerr << "y(x) can\'t be found using given data.";
		return 1;
	}
	
	if( input_data(tempx0,tempy0,tempin) )
		return 0;
	result();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Rkmethod4{
	 private:
	 	double *x,*y,h,in_x;
	 	int n;
	 public:
	 	bool input_data(double,double,double);
	 	void result();
	 	bool main();
};

bool Rkmethod4::input_data(double tempx0, double tempy0, double tempin){
	x = new double[n+1];
	y = new double[n+1];
	*x = tempx0;
	*y = tempy0;
	in_x = tempin;
	for ( int i=0 ; i<=n ; i++)
		*(x+i) = *x+i*h;
	return 0;
}

void Rkmethod4::result(){
	double k,k1,k2,k3,k4;
	cout.setf( ios::floatfield , ios::fixed );
	cout << endl << setprecision(8) << "values_of_x     value_of_y/f(x)" << endl
	             << setw(11) << *x << setw(20) << *y  << endl;
	for ( int i=1 ; i<=n ; i++ ){
		k1 = fxy( *(x+i-1)     , *(y+i-1)        );
		k2 = fxy( *(x+i-1)+h/2 , *(y+i-1)+h*k1/2 );
		k3 = fxy( *(x+i-1)+h/2 , *(y+i-1)+h*k2/2 );
		k4 = fxy( *(x+i-1)+h   , *(y+i-1)+h*k3   );
		k  = (k1 + 2*k2 + 2*k3 + k4) / 6;
		*(y+i) = *(y+i-1) + h * k;
		cout << setw(11) << *(x+i) << setw(20) << *(y+i) << endl;
	}
}

bool Rkmethod4::main(){
	double tempx0,tempy0,tempin;
	bool choice;
	cout << endl << "Enter x0 and y0 : ";
	if( !(cin >> tempx0 >> tempy0)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Enter value of x : ";
	if( !(cin >> tempin) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Chose steps h(0) or n(1) ? : ";
	if ( !(cin >> choice) ){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if ( choice ) {
		cout << endl << "Enter value of n (step number) : ";
		if( !(cin >> n) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		h = ( tempin - tempx0 ) / n;
	}
	else {
		cout << endl << "Enter value of h (step size) : ";
		if( !(cin >> h) ){
			cerr << "Enter valid input." << endl;
			return 1;
		}
		n = round(( tempin - tempx0 ) / h);
	}
	
	if ( mod( (tempx0+n*h) , tempin ) > error ){
		cerr << "y(x) can\'t be found using given data.";
		return 1;
	}
	
	if( input_data(tempx0,tempy0,tempin) )
		return 0;
	result();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Secant{
	 private:
	 	double a,b,fa,fb;
	 public:
	 	bool input_range(double,double);
	 	void print();
	 	bool main();
};

bool Secant::input_range(double x, double y){
	double fx=f(x) , fy=f(y);
	if ( fx == 0 ) {
		cout << "Root is " << x << endl;
		return 1;
	}
	else if( fy == 0 ){
		cout << "Root is " << y << endl;
		return 1;
	}
	else if( ( fx < 0 && fy < 0 ) || ( fx > 0 && fy > 0 ) ){
		cerr << "There is no root in given range." << endl;
		return 1;
	}
	else{
		if( fx < 0 ){
			a=x; b=y; fa=fx; fb=fy;
		}
		else {
			a=y; b=x; fa=fy; fb=fx;
		}
		return 0;
	}
}

void Secant::print(){
	double position=0,prev_position,value;
	int i = 1;
	cout << "Iteration   Approx_Root   Value_of__f(x)" << endl;
	do {
		prev_position = position;
		position = ((a*fb)-(b*fa))/(fb-fa);
		value = f(position);
		if ( value == 0 ) {
			cout << endl << endl << "Root of given function is " << position << endl;
			return;
		}
		
		a=b; fa=fb;
		b=position; fb=value;
		
		cout << setw(9) << left << i;
		cout.setf(ios::floatfield,ios::fixed);
		cout << setprecision(8)  << setw(14) << right << position;
		cout << setprecision(10) << setw(17) << right << value << endl;
		i++;
	}while( mod(position,prev_position) > accuracy);
	cout << endl << endl << "Root of given function is " << setprecision(8) << position << endl;
}

bool Secant::main(){
	double temp1,temp2;
	cout << endl << "Enter range in which root lies : ";
	if( !(cin >> temp1 >> temp2)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if(input_range(temp1,temp2))
		return 0;
	print();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Simpson1by3{
	 private:
	 	double *x,*y,h;
	 	int n;
	 public:
	 	bool input_data(double,double,int);
	 	void result();
	 	bool main();
};

bool Simpson1by3::input_data(double tempa, double tempb, int tempn){
	n = tempn;
	x = new double[n+1];
	y = new double[n+1];
	h = ( tempb - tempa ) / tempn;
	for (int i=0 ; i<=n ; i++){
		*(x+i) = tempa + i * h;
		*(y+i) = f( *(x+i) );
	}
	return 0;
}

void Simpson1by3::result(){
	double ans , odd=0 , even=0;
	for ( int i=1 ; i<n ; i++ ) {
		if (i%2)
			odd += *(y+i);
		else
			even += *(y+i);
	}
	ans = (*(y) + *(y+n) + 2 * even + 4 * odd ) * h / 3; 
	cout << endl << "I = " << ans << endl;
}

bool Simpson1by3::main(){
	double temp1,temp2;
	unsigned int tempn;
	cout << endl << "Enter limits : ";
	if( !(cin >> temp1 >> temp2)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Enter iterations (n) : ";
	if( !(cin >> tempn) || tempn < 2){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if( input_data(temp1,temp2,tempn) )
		return 0;
	result();
	return 0;
}






//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Simpson1by3nf{
	 private:
	 	double *x,*y,h;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void get();
	 	bool main();
};

bool Simpson1by3nf::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Simpson1by3nf::get(){
	double ans , odd=0 , even=0;
	for ( int i=1 ; i<length-1 ; i++ ) {
		if (i%2)
			odd += *(y+i);
		else
			even += *(y+i);
	}
	ans = (*(y) + *(y+length-1) + 2 * even + 4 * odd ) * h / 3; 
	cout << endl << "I = " << ans << endl;
}

bool Simpson1by3nf::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ < 2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ] , tempy[templ];
	
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	h = tempx[1]-tempx[0];
	for(int i=templ-1 ; i>0 ; i--) {
		if( mod(tempx[i]-tempx[i-1],h) > error) {
			cerr << endl << "Sorry given data is not sequential." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	get();
	return 0;
}







//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Simpson3by8{
	 private:
	 	double *x,*y,h;
	 	int n;
	 public:
	 	bool input_data(double,double,int);
	 	void result();
	 	bool main();
};

bool Simpson3by8::input_data(double tempa, double tempb, int tempn){
	n = tempn;
	x = new double[n+1];
	y = new double[n+1];
	h = ( tempb - tempa ) / tempn;
	for (int i=0 ; i<=n ; i++){
		*(x+i) = tempa + i * h;
		*(y+i) = f( *(x+i) );
	}
	return 0;
}

void Simpson3by8::result(){
	double ans , remain=0 , mul3=0;
	for ( int i=1 ; i<n ; i++ ) {
		if ( i%3 == 0 )
			mul3 += *(y+i);
		else
			remain += *(y+i);
	}
	ans = (*(y) + *(y+n) + 2 * mul3 + 3 * remain ) * h * 3 / 8; 
	cout << endl << "I = " << ans << endl;
}

bool Simpson3by8::main(){
	double temp1,temp2;
	unsigned int tempn;
	cout << endl << "Enter limits : ";
	if( !(cin >> temp1 >> temp2)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Enter iterations (n) : ";
	if( !(cin >> tempn) || tempn < 2){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if( input_data(temp1,temp2,tempn) )
		return 0;
	result();
	return 0;
}








//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Simpson3by8nf{
	 private:
	 	double *x,*y,h;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void get();
	 	bool main();
};

bool Simpson3by8nf::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Simpson3by8nf::get(){
	double ans , remain=0 , mul3=0;
	for ( int i=1 ; i<length-1 ; i++ ) {
		if ( (i+1)%3 == 0 )
			mul3 += *(y+i);
		else
			remain += *(y+i);
	}
	ans = (*(y) + *(y+length-1) + 2 * mul3 + 3 * remain ) * h * 3 / 8; 
	cout << endl << "I = " << ans << endl;
}

bool Simpson3by8nf::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ < 2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ] , tempy[templ];
	
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	h = tempx[1]-tempx[0];
	for(int i=templ-1 ; i>0 ; i--) {
		if( mod(tempx[i]-tempx[i-1],h) > error) {
			cerr << endl << "Sorry given data is not sequential." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	get();
	return 0;
}








//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Trapezoidal{
	 private:
	 	double *x,*y,h;
	 	int n;
	 public:
	 	bool input_data(double,double,int);
	 	void result();
	 	bool main();
};

bool Trapezoidal::input_data(double tempa, double tempb, int tempn){
	n = tempn;
	x = new double[n+1];
	y = new double[n+1];
	h = ( tempb - tempa ) / tempn;
	for (int i=0 ; i<=n ; i++){
		*(x+i) = tempa + i * h;
		*(y+i) = f( *(x+i) );
	}
	return 0;
}

void Trapezoidal::result(){
	double ans , temp=0;
	for ( int i=1 ; i<n ; i++ ) {
		temp += *(y+i);
	}
	ans = (*(y) + *(y+n) + 2 * temp ) * h / 2; 
	cout << endl << "I = " << ans << endl;
}

bool Trapezoidal::main(){
	double temp1,temp2;
	unsigned int tempn;
	cout << endl << "Enter limits : ";
	if( !(cin >> temp1 >> temp2)){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	cout << endl << "Enter iterations (n) : ";
	if( !(cin >> tempn) || tempn < 2){
		cerr << "Enter valid input." << endl;
		return 1;
	}
	
	if( input_data(temp1,temp2,tempn) )
		return 0;
	result();
	return 0;
}








//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class Trapezoidalnf{
	 private:
	 	double *x,*y,h;
	 	int length;
	 public:
	 	bool input_data(double[],double[],int);
	 	void get();
	 	bool main();
};

bool Trapezoidalnf::input_data(double tempx[], double tempy[], int templ){
	length = templ;
	x = new double[length];
	y = new double[length];
	for (int i=0 ; i<length ; i++){
		*(x+i) = tempx[i];
		*(y+i) = tempy[i];
	}
	return 0;
}

void Trapezoidalnf::get(){
	double ans , temp=0;
	for ( int i=1 ; i<length-1 ; i++ ) {
		temp += *(y+i);
	}
	ans = (*(y) + *(y+length-1) + 2 * temp ) * h / 2; 
	cout << endl << "I = " << ans << endl;
}

bool Trapezoidalnf::main(){
	unsigned int templ;
	cout << endl << "Enter length of data : ";
	if( !(cin >> templ) || templ < 2){
		cerr << "Enter valid length" << endl;
		return 1;
	}
	double tempx[templ] , tempy[templ];
	
	cout << endl << "Enter values of x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempx[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	h = tempx[1]-tempx[0];
	for(int i=templ-1 ; i>0 ; i--) {
		if( mod(tempx[i]-tempx[i-1],h) > error) {
			cerr << endl << "Sorry given data is not sequential." << endl;
			return 1;
		}
	}
	
	cout << endl << "Enter values of y corresponding to x : ";
	for (int i=0; i<templ ; i++){
		if( !(cin >> tempy[i]) ){
			cerr << "Enter valid data." << endl;
			return 1;
		}
	}
	
	if( input_data(tempx,tempy,templ) )
		return 0;
	get();
	return 0;
}






#endif