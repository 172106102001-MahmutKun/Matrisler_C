#include<iostream>
#include<iomanip>
#include "Eigen/Eigenvalues"
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

double mat[5][5];

void matrisCarp(double (&res)[5][5], double mat1[5][5], double mat2[5][5]) {
	for (int i=0; i<5; i++)
		for (int j=0;j<5;j++)
			for (int k=0; k<5; k++)
				res[i][j] += mat1[i][k]*mat2[k][j];
}

double determinant( double mat[5][5], int n) {
   double det = 0;
   double submat[5][5];
	if (n == 2)
		return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));

	for (int k = 0; k < n; k++) {
		int subi = 0;
		for (int i = 1; i < n; i++) {
			int subj = 0;
			for (int j = 0; j < n; j++) {
				if (j == k)
					continue;
				submat[subi][subj] = mat[i][j];
				subj++;
			}
			subi++;
		}
		det = det + (pow(-1, k) * mat[0][k] * determinant( submat, n - 1 ));
	}
   
	return det;
}

double ATAiz(double mat[5][5]) {
	double iz = 0;
	double matT[5][5];
	for (int i=0; i<5; i++)
		for (int j=0; j<5; j++)
			matT[i][j] = mat[j][i];

	double res[5][5];
	matrisCarp(res, mat, matT);

	for (int i=0; i<5; i++)
		iz += res[i][i];
	return iz;
}

void satirNormlari(double mat[5][5]) {
	for (int i=0; i<5; i++) {
		double curNorm = 0;
		for (int j=0; j<5; j++)
			curNorm += mat[i][j]*mat[i][j];
		printf("%d. satirin normu: %.4lf\n",i+1, sqrt(curNorm));
	}
}
void sutunNormlari(double mat[5][5]) {
	for (int j=0; j<5; j++) {
		double curNorm = 0;
		for (int i=0; i<5; i++)
			curNorm += mat[i][j]*mat[i][j];
		printf("%d. sutunun normu: %.4lf\n",j+1, sqrt(curNorm));
	}
}

double oklidNorm(double mat[5][5]) {
	double curNorm = 0;
	for (int i=0; i<5; i++)
		for (int j=0; j<5; j++)
			curNorm += mat[i][j]*mat[i][j];

	return sqrt(curNorm);
}

string oklidIzKontrol(double mat[5][5]) {
	if (sqrt(ATAiz(mat)) == oklidNorm(mat))
		return "ATA izinin karekoku ile oklid normlari birbirine esit.";
	else
		return "ATA izinin karekoku ile oklid normlari birbirine esit degil!";
}

void oklidNormlastir(double mat[5][5]) {
	double oklidN = oklidNorm(mat);
	cout << "Oklid Normaline gore normalastirilmis matris:" << endl;
	for (int i=0; i<5; i++){
		for (int j=0; j<5; j++)
			cout << mat[i][j]/oklidN << " ";
		cout << endl;
	}
}

void ozdeger(double mat[5][5]) {
	Matrix<double,5,5> mat2;
	for (int i=0; i<5; i++)
		for (int j=0; j<5; j++)
			mat2(i,j) = mat[i][j];

	ComplexEigenSolver<Matrix<double,5,5> > res;
	res.compute(mat2);
	cout << "Matrisin ozdegerleri:" << endl << res.eigenvalues() << endl;
}

void ters(double mat[5][5]) {
	Matrix<double,5,5> mat2;
	for (int i=0; i<5; i++)
		for (int j=0; j<5; j++)
			mat2(i,j) = mat[i][j];

	mat2 = mat2.inverse();

	cout << "Matrisin tersi:" << endl;
	for (int i=0; i<5; i++){
		for (int j=0; j<5; j++)
			cout << mat2(i,j) << " ";
		cout << endl;
	}
}

void xBilinmeyenler(double mat[5][5]) {
	MatrixXd A(5,4);
	VectorXd b(5);

	for (int i=0; i<5; i++)
		for (int j=0; j<4; j++)
			A(i,j) = mat[i][j];

	for (int i=0; i<5; i++)
		b(i) = mat[i][4];

	VectorXd x = A.colPivHouseholderQr().solve(b);
	cout << "X Bilinmeyenler vektoru:\n" << x << endl;
}

int main() {
	for (int i=0; i<5; i++)
		for(int j=0; j<5; j++)
			cin >> mat[i][j];
			// mat[i][j] = j;

    cout << fixed << showpoint;
    cout << setprecision(4);

    while(true) {
    	cout << "0) Cikis" << endl;
    	cout << "1) A Matsisinin Determinanti" << endl;
    	cout << "2) ATA Matrisinin Izi" << endl;
    	cout << "3) A Matrisinin Satir Normlari" << endl;
    	cout << "4) A Matrisinin Sutun Normlari" << endl;
    	cout << "5) A Matrisinin Oklid Normu (N(A))" << endl;
    	cout << "6) Oklid ve ATA Izi Kontrolu" << endl;
    	cout << "7) Oklid Normu Ile Normlastirma" << endl;
    	cout << "8) A Matrisinin Ozdegerleri" << endl;
    	cout << "11) A Matrisinin Tersi" << endl;
    	cout << "14) x Bilinmeyenler Vektoru" << endl;
    	cout << "Seciminz: ";
    	int x;
    	cin >> x;

    	if (x == 0)
    		break;
    	else if (x == 1)
			cout <<"A Matsisinin Determinanti: "<< determinant(mat, 5) << endl << endl;
    	else if (x == 2)
			cout << "ATA Matrisinin Izi: " << ATAiz(mat) << endl << endl;

    	else if (x == 3)
			satirNormlari(mat);

    	else if (x == 4)
			sutunNormlari(mat);

    	else if (x == 5)
			cout << "A Matrisinin Oklid Normu (N(A)): " << oklidNorm(mat) << endl << endl;

    	else if (x == 6)
			cout << oklidIzKontrol(mat) << endl << endl;	

    	else if (x == 7)
			oklidNormlastir(mat);

    	else if (x == 8)
			ozdeger(mat);

    	else if (x == 11)
			ters(mat);

    	else if (x == 14)
			xBilinmeyenler(mat);



    }

}

/*
1 2 3 4 5
5 4 3 2 1
1 1 3 45 4
4 2 2 2 2
5 3 12 1 1
*/










