#include"Variable.h"
#include<algorithm>
#include"eigen3/Eigen/Eigen"

using namespace Eigen;

double radicalSign_g(int y) {                      //�βΣ���ǰ���y����
	int n = y - ny / 2;                            //��10.05Ϊ�е㣬��Ϊ������Ϊ��
	double length = dy * n;
	double ans = 2.0 * M_PI * (10.05 - length);
	return ans;
}

//��������Է��̺�����ܶȸ��£�Vx��Vy
void VxVy_update() {
	//�������������ڲ���
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //�����ų����߽�㣬���ų��ڲ��߽��
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;
			
			velocity_x[i][j] = pitch * velocity_parallel_i[i][j] - 1.0 / density[i][j] * (density[i + 1][j] - density[i - 1][j]) / (2*dx);

			velocity_y[i][j] = -1.0 / density[i][j] * (density[i][j + 1] - density[i][j - 1]) / (2*dy);
		}
	}

	//core plasma (GG')    ���õ�����ߵ�ֵ
	for (int j = int(l_percentage_private_flux * nx + 1); j < (1 - l_percentage_private_flux) * nx - 1; j++) {
		velocity_x[0][j] = velocity_x[1][j];

		velocity_y[0][j] = velocity_y[1][j];
	}
	//target conditions	(AB and CD)    ����Ϊ�̶�ֵxΪpitch*v//,yΪ0
	for (int j = 1; j < ny - 1; j++) {
		velocity_x[0][j] = pitch * velocity_parallel_i[0][j];

		velocity_y[0][j] = 0;
	}
	for (int j= 1; j < ny - 1; j++) {
		velocity_x[nx-1][j] = pitch * velocity_parallel_i[nx-1][j];

		velocity_y[nx-1][j] = 0;
	}
	//private flux	(AE and E'C) and wall (BD)      ���õ������ڵ��ֵ
	for (int i = 0; i < l_percentage_private_flux * nx - 1; i++) {
		velocity_x[i][0] = velocity_x[i][1];

		velocity_y[i][0] = velocity_y[i][1];
	}
	for (int i = int((1 - l_percentage_private_flux) * nx + 1); i <= nx - 1; i++) {             
		velocity_x[i][0] = velocity_x[i][1];

		velocity_y[i][0] = velocity_y[i][1];
	}
	for (int i = 0; i <= nx - 1; i++) {
		velocity_x[i][ny-1] = velocity_x[i][ny-2];

		velocity_y[i][ny-1] = velocity_y[i][ny-2];
	}
	//�ڲ����µ��������߽��߷ֱ�����         ���õ������ڵ��ֵ
	for (int i = 0; i <= w_percentage_private_flux * ny - 1; i++) {
		velocity_x[int(l_percentage_private_flux * nx - 1)][i] = velocity_x[int(l_percentage_private_flux * nx - 2)][i];
		velocity_y[int(l_percentage_private_flux * nx - 1)][i] = velocity_y[int(l_percentage_private_flux * nx - 2)][i];
		
		velocity_x[int(l_percentage_private_flux * nx )][i] = velocity_x[int(l_percentage_private_flux * nx + 1)][i];
		velocity_y[int(l_percentage_private_flux * nx )][i] = velocity_y[int(l_percentage_private_flux * nx + 1)][i];

		velocity_x[int((1 - l_percentage_private_flux) * nx - 1)][i] = velocity_x[int((1 - l_percentage_private_flux) * nx - 2)][i];
		velocity_y[int((1 - l_percentage_private_flux) * nx - 1)][i] = velocity_y[int((1 - l_percentage_private_flux) * nx - 2)][i];

		velocity_x[int((1 - l_percentage_private_flux) * nx)][i] = velocity_x[int((1 - l_percentage_private_flux) * nx + 1)][i];
		velocity_y[int((1 - l_percentage_private_flux) * nx)][i] = velocity_y[int((1 - l_percentage_private_flux) * nx + 1)][i];
	}

}


//��������ion-ion��ײBraginskii��ײʱ��
double tau_i_Brag(int i, int j) {       //ֻ��Ҫ�����������¶�Ti/�ܶȵ�����
	double ans = 3.0 / 4.0 * (pow(m_i, 0.5) / pow(M_PI, 0.5)) *
		(pow(temperature_i[i][j], 1.5) / pow(Z_i, 4.0) / density[i][j] / 12.0) * 
		pow((4.0 * M_PI * permittivity / pow(electron, 2.0)), 2.0);
	return ans;
}
//�������̼���ճ��ϵ��
double eta_i_x(int i, int j) {
	double ans = pow(pitch, 2.0) * 0.96 * density[i][j] * temperature_i[i][j] * tau_i_Brag(i, j);
	return 0.9924*ans;
}
//�������̾���ճ��ϵ��
double eta_i_y(int i, int j) {
	double ans = 0.2 * m_i * density[i][j];
	return 0.9924*ans;
}

//����-��������ƽ����
double k(int i,int j) {
	double res = 4.8e-15 * 12.0 * pow(electron, 1.5) * pow(density[i][j], 2.0) / pow(temperature_e[i][j], 1.5);    //��mp=mi�������Ӻ���������
	return res;
}
//electrons-ions��ײʱ���e
double tau_e(int i,int j) {
	double ans = 3.0 / 4.0 * pow(m_e, 0.5) / pow(2.0 * M_PI, 0.5) *
		(pow(temperature_e[i][j], 1.5) / density[i][j] / 12.0) *
		pow((4.0 * M_PI * permittivity / pow(electron , 2.0)), 2.0);                  //electron ?????????????
	return ans;
}
//���Ӽ�����ϵ��
double kappa_e_x(int i,int j) {
	double ans = pow(pitch, 2.0) * 3.2 * density[i][j] * temperature_e[i][j] * tau_e(i, j) / m_e;
	return 1.3288*ans;
}
//���Ӿ�����ϵ��
double kappa_e_y(int i,int j) {
	return 1.3288*density[i][j];
}

//���Ӽ�����ϵ��
double kappa_i_x(int i, int j) {
	double ans = pow(pitch, 2.0) * 3.9 * density[i][j] * temperature_i[i][j] * tau_i_Brag(i,j) / m_i;
	return 1.0140*ans;
}
//���Ӿ�����ϵ��
double kappa_i_y(int i, int j) {
	return 1.0140*density[i][j];
}

void Eigen_use() {
	SparseMatrix<double> sparse(nx * ny, nx * ny);
	VectorXd xx = VectorXd::Zero(nx * ny), bb = VectorXd::Zero(nx * ny);
	std::vector<Eigen::Triplet<double>> tripletlist;

	for (int i = ny; i < nx * ny; i++) {
		if (a[i] != 0) {
			tripletlist.push_back(Triplet<double>(i, i - ny, a[i]));
		}
	}

	for (int i = 1; i < nx * ny; i++) {
		if (a[i] != 0) {
			tripletlist.push_back(Triplet<double>(i, i - 1, b[i]));
		}
	}
	for (int i = 0; i < nx * ny; i++) {
		tripletlist.push_back(Triplet<double>(i, i, c[i]));

	}
	for (int i = 0; i < nx * ny - 1; i++) {
		if (a[i] != 0) {
			tripletlist.push_back(Triplet<double>(i, i + 1, d[i]));
		}
	}
	for (int i = 0; i < nx * ny - ny; i++) {
		if (a[i] != 0) {
			tripletlist.push_back(Triplet<double>(i, i + ny, e[i]));
		}
	}
	sparse.setFromTriplets(tripletlist.begin(), tripletlist.end());

	sparse.makeCompressed();
	BiCGSTAB<SparseMatrix<double>> solver_sparse;

	solver_sparse.setTolerance(0.001);
	solver_sparse.compute(sparse);

	for (int i = 0; i < nx * ny; i++) {
		bb[i] = source[i];
	}
	xx = solver_sparse.solve(bb);

	for (int i = 0; i < nx * ny; i++) {
		x[i] = xx[i];
	}
}