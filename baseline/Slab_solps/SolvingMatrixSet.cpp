#include"Variable.h"
#include"SolvingMatrixSet.h"
#include"ohters.h"
#include<algorithm>
using namespace std;


double max_upstream(double D, double F) {
	return D + max(F, 0.0);
}

double max_downstream(double D, double F) {
	return D + max(-F, 0.0);
}


void ContinuityMatrixRenew() {
	//�ų����б߽�������
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //�����ų����߽�㣬���ų��ڲ��߽��
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//����һ��ӭ���ʽ
			//��������Խ��߶��壬�ֱ�������ϡ������С���������
			//S
			a[num] = -(
				max_upstream(radicalSign_g(j) * dy / dx,
					(radicalSign_g(j) * pitch * (velocity_parallel_i[i][j] + velocity_parallel_i[i - 1][j]) / 2.0 * dy))
				);
			//E
			b[num] = -((radicalSign_g(j) + radicalSign_g(j - 1)) / 2.0 * dx / dy);
			//W
			d[num] = -((radicalSign_g(j) + radicalSign_g(j + 1)) / 2.0 * dx / dy);
			//N
			e[num] = -(
				max_downstream(radicalSign_g(j) * dy / dx,
					(radicalSign_g(j) * pitch * (velocity_parallel_i[i][j] + velocity_parallel_i[i + 1][j]) / 2.0 * dy))
				);

			//P
			c[num] = -a[num] - b[num] - d[num] - e[num] +
				(radicalSign_g(j) * pitch * (velocity_parallel_i[i][j] + velocity_parallel_i[i + 1][j]) / 2.0 * dy) -               //Fn
				(radicalSign_g(j) * pitch * (velocity_parallel_i[i][j] + velocity_parallel_i[i - 1][j]) / 2.0 * dy) +               //-Fs
				(dx * dy / delta_t) * radicalSign_g(j);



			source[num] = (dx * dy / delta_t) * radicalSign_g(j) * density[i][j];      //������ap0*np0������һʱ���Ӱ��Դ��
		}
	}
}

void MomentumMatrixRenew() {
	//�ų����б߽�������
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //�����ų����߽�㣬���ų��ڲ��߽��
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//����һ��ӭ���ʽ
			//��������Խ��߶��壬�ֱ�������ϡ������С���������
			//S
			a[num] = -(
				max_upstream((radicalSign_g(j) * 4.0 / 3.0 * eta_i_x(i, j) * dy + radicalSign_g(j) * 4.0 / 3.0 * eta_i_x(i - 1, j) * dy) / 2.0 / dx,
					(radicalSign_g(j) * m_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * m_i * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0)
				);
			//E
			b[num] = -(
				max_upstream((radicalSign_g(j) * eta_i_y(i, j) * dx + radicalSign_g(j - 1) * eta_i_y(i, j - 1) * dx) / 2.0 / dy,
					(radicalSign_g(j) * m_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * m_i * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0)
				);
			//W
			d[num] = -(
				max_downstream((radicalSign_g(j) * eta_i_y(i, j) * dx + radicalSign_g(j + 1) * eta_i_y(i, j + 1) * dx) / 2.0 / dy,
					((radicalSign_g(j) * m_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * m_i * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0))
				);
			//N
			e[num] = -(
				max_downstream((radicalSign_g(j) * 4.0 / 3.0 * eta_i_x(i, j) * dy + radicalSign_g(j) * 4.0 / 3.0 * eta_i_x(i + 1, j) * dy) / 2.0 / dx,
					((radicalSign_g(j) * m_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * m_i * density[i + 1][j] * velocity_x[i + 1][j] * dy) / 2.0))
				);

			//P
			c[num] = -a[num] - b[num] - d[num] - e[num] +
				(radicalSign_g(j) * m_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * m_i * density[i + 1][j] * velocity_x[i + 1][j] * dy) / 2.0 +                 //deltaF,Fn+Fw-Fs-Fe
				(radicalSign_g(j) * m_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * m_i * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0 -
				(radicalSign_g(j) * m_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * m_i * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0 -
				(radicalSign_g(j) * m_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * m_i * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0 +
				m_i * density[i][j] * dx * dy * radicalSign_g(j) / delta_t;


			source[num] = (m_i * density[i][j] * dx * dy * radicalSign_g(j) / delta_t) * velocity_parallel_i[i][j] -          //������ap0*np0������һʱ���Ӱ��Դ��
				radicalSign_g(j) * pitch * (pressure[i + 1][j] - pressure[i - 1][j]) / (2.0 * dx) * (dx * dy);
		}
	}
}

void ElectronEnergyMatrixRenew() {
	//�ų����б߽�������
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //�����ų����߽�㣬���ų��ڲ��߽��
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//����һ��ӭ���ʽ
			//��������Խ��߶��壬�ֱ�������ϡ������С���������
			 
			D_s[num] = (radicalSign_g(j) * kappa_e_x(i, j) * dy / dx + radicalSign_g(j) * kappa_e_x(i - 1, j) * dy / dx) / 2.0;
			F_s[num] = max(0.0, (radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0);

			D_n[num] = (radicalSign_g(j) * kappa_e_x(i, j) * dy / dx + radicalSign_g(j) * kappa_e_x(i + 1, j) * dy / dx) / 2.0;
			F_n[num] = max(0.0, -((radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i + 1][j] * velocity_x[i + 1][j] * dy) / 2.0));

			D_e[num] = (radicalSign_g(j) * kappa_e_y(i, j) * dx / dy + radicalSign_g(j - 1) * kappa_e_y(i, j - 1) * dx / dy) / 2.0;
			F_e[num] = max(0.0, (radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * 5.0 / 2.0 * Z_i * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0);
			D_w[num] = (radicalSign_g(j) * kappa_e_y(i, j) * dx / dy + radicalSign_g(j + 1) * kappa_e_y(i, j + 1) * dx / dy) / 2.0;
			F_w[num] = max(0.0, -((radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * 5.0 / 2.0 * Z_i * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0));

			
			//S
			a[num] = -(
				max_upstream((radicalSign_g(j) * kappa_e_x(i, j) * dy / dx + radicalSign_g(j) * kappa_e_x(i - 1, j) * dy / dx) / 2.0,
					(radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0)
				);
			//E
			b[num] = -(
				max_upstream((radicalSign_g(j) * kappa_e_y(i, j) * dx / dy + radicalSign_g(j - 1) * kappa_e_y(i, j - 1) * dx / dy) / 2.0,
					(radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * 5.0 / 2.0 * Z_i * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0)
				);
			//W
			d[num] = -(
				max_downstream((radicalSign_g(j) * kappa_e_y(i, j) * dx / dy + radicalSign_g(j + 1) * kappa_e_y(i, j + 1) * dx / dy) / 2.0,
					((radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * 5.0 / 2.0 * Z_i * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0))
				);
			//N
			e[num] = -(
				max_downstream((radicalSign_g(j) * kappa_e_x(i, j) * dy / dx + radicalSign_g(j) * kappa_e_x(i + 1, j) * dy / dx) / 2.0,
					((radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i + 1][j] * velocity_x[i + 1][j] * dy) / 2.0))
				);

			//P
			c[num] = -a[num] - b[num] - d[num] - e[num] +
				(radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i + 1][j] * velocity_x[i + 1][j] * dy) / 2.0 +   //deltaF,Fn+Fw-Fs-Fe
				(radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * 5.0 / 2.0 * Z_i * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0 -
				(radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0 -
				(radicalSign_g(j) * 5.0 / 2.0 * Z_i * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * 5.0 / 2.0 * Z_i * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0 +
				radicalSign_g(j) * dx * dy * k(i, j) +
				3.0 / 2.0 * Z_i * density[i][j] * dx * dy * radicalSign_g(j) / delta_t;


			source[num] = temperature_e[i][j] * 3.0 / 2.0 * Z_i * density[i][j] * dx * dy * radicalSign_g(j) / delta_t +
				radicalSign_g(j) * dx * dy * k(i, j) * temperature_i[i][j] +                                                       //k(Te-Ti),Te�Ƶ�c��ȥ��
				radicalSign_g(j) * dx * dy * pitch * velocity_parallel_i[i][j] * (pressure_e[i + 1][j] - pressure_e[i - 1][j]) / (2.0 * dx);
		}
	}
}

void IonEnergyMatrixRenew() {
	//�ų����б߽�������
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //�����ų����߽�㣬���ų��ڲ��߽��
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//����һ��ӭ���ʽ
			//��������Խ��߶��壬�ֱ�������ϡ������С���������
			
			D_s[num] = (radicalSign_g(j) * kappa_i_x(i, j) * dy + radicalSign_g(j) * kappa_i_x(i - 1, j) * dy) / 2.0 / dx;
			F_s[num] = max(0.0, (radicalSign_g(j) * 5.0 / 2.0 * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0);

			D_n[num] = (radicalSign_g(j) * kappa_i_x(i, j) * dy + radicalSign_g(j) * kappa_i_x(i + 1, j) * dy) / 2.0 / dx;
			F_n[num] = max(0.0,-((radicalSign_g(j) * 5 / 2 * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * 5.0 / 2.0 * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0));
			
			D_e[num] = (radicalSign_g(j) * kappa_i_y(i, j) * dx + radicalSign_g(j - 1) * kappa_i_y(i, j - 1) * dx) / 2.0 / dy;
			F_e[num] = max(0.0, (radicalSign_g(j) * 5.0 / 2.0 * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * 5.0 / 2.0 * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0);
			D_w[num] = (radicalSign_g(j) * kappa_i_y(i, j) * dx + radicalSign_g(j + 1) * kappa_i_y(i, j + 1) * dx) / 2.0 / dy;
			F_w[num] = max(0.0,-((radicalSign_g(j) * 5 / 2 * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * 5.0 / 2.0 * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0));

			//S
			a[num] = -(
				max_upstream((radicalSign_g(j) * kappa_i_x(i, j) * dy + radicalSign_g(j) * kappa_i_x(i - 1, j) * dy) / 2.0 / dx,
					(radicalSign_g(j) * 5.0 / 2.0 * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0)
				);
			//E
			b[num] = -(
				max_upstream((radicalSign_g(j) * kappa_i_y(i, j) * dx + radicalSign_g(j - 1) * kappa_i_y(i, j - 1) * dx) / 2.0 / dy,
					(radicalSign_g(j) * 5.0 / 2.0 * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * 5.0 / 2.0 * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0)
				);
			//W
			d[num] = -(
				max_downstream((radicalSign_g(j) * kappa_i_y(i, j) * dx + radicalSign_g(j + 1) * kappa_i_y(i, j + 1) * dx) / 2.0 / dy,
					((radicalSign_g(j) * 5 / 2 * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * 5.0 / 2.0 * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0))
				);
			//N
			e[num] = -(
				max_downstream((radicalSign_g(j) * kappa_i_x(i, j) * dy + radicalSign_g(j) * kappa_i_x(i + 1, j) * dy) / 2.0 / dx,
					((radicalSign_g(j) * 5 / 2 * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * density[i + 1][j] * velocity_x[i + 1][j] * dy) / 2.0))
				);

			//P
			c[num] = -a[num] - b[num] - d[num] - e[num] +
				(radicalSign_g(j) * 5 / 2 * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * density[i + 1][j] * velocity_x[i + 1][j] * dy) / 2.0 +           //deltaF,Fn+Fw-Fs-Fe
				(radicalSign_g(j) * 5 / 2 * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j + 1) * 5.0 / 2.0 * density[i][j + 1] * velocity_y[i][j + 1] * dx) / 2.0 -
				(radicalSign_g(j) * 5.0 / 2.0 * density[i][j] * velocity_x[i][j] * dy + radicalSign_g(j) * 5.0 / 2.0 * density[i - 1][j] * velocity_x[i - 1][j] * dy) / 2.0 -
				(radicalSign_g(j) * 5.0 / 2.0 * density[i][j] * velocity_y[i][j] * dx + radicalSign_g(j - 1) * 5.0 / 2.0 * density[i][j - 1] * velocity_y[i][j - 1] * dx) / 2.0 +
				3.0 / 2.0 * density[i][j] * dx * dy * radicalSign_g(j) / delta_t;

			source[num] = temperature_i[i][j] * 3.0 / 2.0 * density[i][j] * dx * dy * radicalSign_g(j) / delta_t +     //ap0
				radicalSign_g(j) * dx * dy * k(i, j) * (temperature_e[i][j] - temperature_i[i][j]) -      //k(Te-Ti)
				radicalSign_g(j) * dx * dy * pitch * velocity_parallel_i[i][j] * (pressure_e[i + 1][j] - pressure_e[i - 1][j]) / (2.0 * dx)    //partial p/x
				//n-s
				- radicalSign_g(j) / 2.0 * m_i * (
					(density[i + 1][j] + density[i][j]) / 2.0 * (velocity_x[i + 1][j] + velocity_x[i][j]) / 2.0 * (pow((velocity_parallel_i[i + 1][j] + velocity_parallel_i[i][j]) / 2.0, 2.0)) * dy -
					(density[i - 1][j] + density[i][j]) / 2.0 * (velocity_x[i - 1][j] + velocity_x[i][j]) / 2.0 * (pow((velocity_parallel_i[i - 1][j] + velocity_parallel_i[i][j]) / 2.0, 2.0)) * dy
					) 
				//n-s
				+ radicalSign_g(j) * (
					(eta_i_x(i + 1, j) + eta_i_x(i, j)) / 2.0 * (pow(velocity_parallel_i[i + 1][j], 2.0) - pow(velocity_parallel_i[i][j], 2.0)) / dx * dy -
					(eta_i_x(i - 1, j) + eta_i_x(i, j)) / 2.0 * (pow(velocity_parallel_i[i][j], 2.0) - pow(velocity_parallel_i[i - 1][j], 2.0)) / dx * dy
					)
				//w-e
				- 1.0 / 2.0 * m_i * (
					(radicalSign_g(j + 1) + radicalSign_g(j)) / 2.0 * (density[i][j + 1] + density[i][j]) / 2.0 * (velocity_y[i][j + 1] + velocity_y[i][j]) / 2.0 * (pow((velocity_parallel_i[i][j + 1] + velocity_parallel_i[i][j]) / 2.0, 2.0)) * dx -
					(radicalSign_g(j - 1) + radicalSign_g(j)) / 2.0 * (density[i][j - 1] + density[i][j]) / 2.0 * (velocity_y[i][j - 1] + velocity_y[i][j]) / 2.0 * (pow((velocity_parallel_i[i][j - 1] + velocity_parallel_i[i][j]) / 2.0, 2.0)) * dx
					) 
				//w-e
				+ 1.0 / 2.0 * (
					(radicalSign_g(j + 1) + radicalSign_g(j)) / 2.0 * (eta_i_y(i, j + 1) + eta_i_y(i, j)) / 2.0 * (pow(velocity_parallel_i[i][j + 1], 2.0) - pow(velocity_parallel_i[i][j], 2.0)) / dy * dx -
					(radicalSign_g(j - 1) + radicalSign_g(j)) / 2.0 * (eta_i_y(i, j - 1) + eta_i_y(i, j)) / 2.0 * (pow(velocity_parallel_i[i][j], 2.0) - pow(velocity_parallel_i[i][j - 1], 2.0)) / dy * dx
					)

				- 0.5 * m_i * (density[i][j] * pow(velocity_parallel_i[i][j], 2.0) - density_old[i][j] * pow(velocity_parallel_i_old[i][j], 2.0)) * dx * dy * radicalSign_g(j) / delta_t;  //������
		}
	}
}