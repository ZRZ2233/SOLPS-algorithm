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
	//排除所有边界点的设置
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //首先排除外界边界点，再排除内部边界点
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//采用一阶迎风格式
			//求解矩阵五对角线定义，分别代表点的南、东、中、西、北侧
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



			source[num] = (dx * dy / delta_t) * radicalSign_g(j) * density[i][j];      //这里是ap0*np0代表上一时层的影响源项
		}
	}
}

void MomentumMatrixRenew() {
	//排除所有边界点的设置
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //首先排除外界边界点，再排除内部边界点
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//采用一阶迎风格式
			//求解矩阵五对角线定义，分别代表点的南、东、中、西、北侧
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


			source[num] = (m_i * density[i][j] * dx * dy * radicalSign_g(j) / delta_t) * velocity_parallel_i[i][j] -          //这里是ap0*np0代表上一时层的影响源项
				radicalSign_g(j) * pitch * (pressure[i + 1][j] - pressure[i - 1][j]) / (2.0 * dx) * (dx * dy);
		}
	}
}

void ElectronEnergyMatrixRenew() {
	//排除所有边界点的设置
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //首先排除外界边界点，再排除内部边界点
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//采用一阶迎风格式
			//求解矩阵五对角线定义，分别代表点的南、东、中、西、北侧
			 
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
				radicalSign_g(j) * dx * dy * k(i, j) * temperature_i[i][j] +                                                       //k(Te-Ti),Te移到c中去了
				radicalSign_g(j) * dx * dy * pitch * velocity_parallel_i[i][j] * (pressure_e[i + 1][j] - pressure_e[i - 1][j]) / (2.0 * dx);
		}
	}
}

void IonEnergyMatrixRenew() {
	//排除所有边界点的设置
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {           //首先排除外界边界点，再排除内部边界点
			if ((i == l_percentage_private_flux * nx - 1 || i == l_percentage_private_flux * nx || i == (1 - l_percentage_private_flux) * nx - 1 || i == (1 - l_percentage_private_flux) * nx) && (j - (w_percentage_private_flux * ny - 1) <= 1e-6)) continue;

			int num = i * ny + j;
			//采用一阶迎风格式
			//求解矩阵五对角线定义，分别代表点的南、东、中、西、北侧
			
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

				- 0.5 * m_i * (density[i][j] * pow(velocity_parallel_i[i][j], 2.0) - density_old[i][j] * pow(velocity_parallel_i_old[i][j], 2.0)) * dx * dy * radicalSign_g(j) / delta_t;  //动力项
		}
	}
}