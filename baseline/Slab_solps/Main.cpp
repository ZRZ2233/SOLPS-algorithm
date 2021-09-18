#include<iostream>
#include<fstream>
#include"Variable.h"
#include"EquationSolve.h"
#include"SolvingMatrixSet.h"
#include"BoundarySet.h"
#include"ohters.h"
#include<string>
#include"eigen3/Eigen/Eigen"

using namespace Eigen;
using namespace std;

//变量存储矩阵初始值设置
void initialValueSet() {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			density[i][j] = 3e18;
			density_old[i][j] = 3e18;

			velocity_x[i][j] = 0;
			velocity_y[i][j] = 0;
			velocity_parallel_i[i][j] = 1;
			velocity_parallel_i_old[i][j] = 1;

			temperature_e[i][j] = 40.0 * electron;
			temperature_i[i][j] = 40.0 * electron;

			pressure_e[i][j] = density[i][j] * temperature_e[i][j];
			pressure_i[i][j] = density[i][j] * temperature_i[i][j];
			pressure[i][j] = pressure_e[i][j] + pressure_i[i][j];

		
		}
	}
}
//把Ax=b中的所有变量都归零
void setZero() {
	for (int i = 0; i < nx * ny; i++) {
		a[i] = 0;
		b[i] = 0;
		c[i] = 0;
		d[i] = 0;
		e[i] = 0;

		x[i] = 0;

		source[i] = 0;
	}
}
//将计算出的值更新到相应变量存储矩阵中去
void resultUpdate(double* new_result, double(*old_result)[ny]) {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			old_result[i][j] = new_result[i * ny + j];
		}
	}
}
//输出最后结果到文件中
void output(string str, double(*matrix)[ny]) {
	ofstream out("C:\\Users\\Admin\\Desktop\\baseline_data\\" + str + ".dat", ios::app);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			out << matrix[i][j] << " ";
		}
		out << endl;
	}
	out.close();
}
void output1(string str, double(*matrix)[ny]) {                //用于温度，自动转换能量为温度
	ofstream out("C:\\Users\\Admin\\Desktop\\baseline_data\\" + str + ".dat", ios::app);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			out << matrix[i][j] / electron << " ";
		}
		out << endl;
	}
	out.close();
}
void output2(string str, double* matrix) {
	ofstream out("C:\\Users\\Admin\\Desktop\\baseline_data\\" + str + ".dat", ios::app);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			out << matrix[i * ny + j] << "||";
		}
		out << endl;
	}
	out.close();
}

//output1("a", a);
//output1("b", b);
//output1("c", c);
//output1("d", d);
//output1("e", e);
//output1("source", source);

//string s = "density";
//s = to_string(i) + s;
//output(s, density);

//output("velocityx", velocity_x);
//output("velocityy", velocity_y);
//output("density", density);
void main() {
	//前期工作，初始值设置
	initialValueSet();
	clock_t start, end;
	start = clock();
	double time = 0;
	int count = 1;

	//收敛判断
	double num_convergence = 1e-3;
	double num_int = 3;
	bool convergence = true;
	double A, B, C, D; //四个值保存均方根
	double A_old = 1, B_old = 1, C_old = 1, D_old = 1;

	//residual保存
	ofstream outnew("C:\\Users\\Admin\\Desktop\\baseline_data\\residual.dat", ios::app);

	do {

		cout << "baseline计算时间：" << time << "s" << endl;

		//连续性方程求解密度////////////////////////////////////////////////////////////////
		{
			setZero();                                          //把Ax=b中的所有变量都归零

			ContinuityMatrixRenew();                            //求解区矩阵设置

			ContinuityMatrixBoundary();                         //边界矩阵设置

			Eigen_use();                                        //矩阵求解 

			resultUpdate(x, density);                          //将计算出的值更新到相应变量存储矩阵中去
		}

		VxVy_update();                                     //利用刚刚计算出来的密度更新Vx、Vy


		//离子平行动量方程求解平行速度/////////////////////////////////////////////////////////
		{
			setZero();

			MomentumMatrixRenew();

			MomentumMatrixBoundary();

			Eigen_use();

			resultUpdate(x, velocity_parallel_i);
		}

		VxVy_update();


		//电子能量方程，求解电子温度/////////////////////////////////////////////////////////  
		{
			setZero();

			ElectronEnergyMatrixRenew();

			ElectronEnergyMatrixBoundary();

			Eigen_use();

			resultUpdate(x, temperature_e);
		}


		//离子能量方程，求解离子温度///////////////////////////////////////////////////////
		{
			setZero();

			IonEnergyMatrixRenew();

			IonEnergyMatrixBoundary();

			Eigen_use();

			resultUpdate(x, temperature_i);
		}
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				density_old[i][j] = density[i][j];
				velocity_parallel_i_old[i][j] = velocity_parallel_i_old[i][j];
			}
		}
		//利用刚刚算出来的温度更新压强////////////////////////////////////////////////////// 
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				pressure_e[i][j] = density[i][j] * temperature_e[i][j];
				pressure_i[i][j] = density[i][j] * temperature_i[i][j];
				pressure[i][j] = pressure_i[i][j] + pressure_e[i][j];
			}
		}

		//进入下一个时间步
		time += delta_t;

		//收敛判断
		A = 0, B = 0, C = 0, D = 0;
		convergence = true;

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				A += pow(density[i][j], 2.0);
				B += pow(velocity_parallel_i[i][j], 2.0);
				C += pow(temperature_e[i][j], 2.0);
				D += pow(temperature_i[i][j], 2.0);
			}
		}
		A = pow(A / (nx * ny), 0.5);
		B = pow(B / (nx * ny), 0.5);
		C = pow(C / (nx * ny), 0.5);
		D = pow(D / (nx * ny), 0.5);

		outnew << time << " " << fabs(A - A_old) / A_old << " " << fabs(B - B_old) / B_old << " " << fabs(C - C_old) / C_old << " " << fabs(D - D_old) / D_old << endl;

		if (fabs(A - A_old) / A_old > num_convergence) convergence = false;
		if (fabs(B - B_old) / B_old > num_convergence) convergence = false;
		if (fabs(C - C_old) / C_old > num_convergence) convergence = false;
		if (fabs(D - D_old) / D_old > num_convergence) convergence = false;

		if (convergence || time == 0.2691) {
			cout << "baseline收敛了" << endl;

			end = clock();
			double sum_time = (double)(end - start) / CLOCKS_PER_SEC;

			string folder = to_string(num_int) + "_" + to_string(time) + "time" + to_string(sum_time);
			string command = "mkdir -p C:\\Users\\Admin\\Desktop\\baseline_data\\" + folder;
			system(command.c_str());

			output(folder + "\\density", density);
			output(folder + "\\velocity", velocity_parallel_i);
			output1(folder + "\\e_temperature", temperature_e);
			output1(folder + "\\i_temperature", temperature_i);
			output(folder + "\\pressure_e", pressure_e);
			output(folder + "\\pressure_i", pressure_i);
			output(folder + "\\pressure", pressure);
			output(folder + "\\v_x", velocity_x);
			output(folder + "\\x_y", velocity_y);

			num_convergence = num_convergence / 10.0;
			num_int++;
		}

		A_old = A;
		B_old = B;
		C_old = C;
		D_old = D;

	} while (time < end_time);

	outnew.close();

	//将能量转换为温度
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			temperature_e[i][j] = temperature_e[i][j] / electron;
			temperature_i[i][j] = temperature_i[i][j] / electron;
		}
	}

	//输出最后结果到文件中
	output("density", density);
	output("velocity", velocity_parallel_i);
	output("e_temperature", temperature_e);
	output("i_temperature", temperature_i);
	output("pressure_e", pressure_e);
	output("pressure_i", pressure_i);
	output("pressure", pressure);

}

