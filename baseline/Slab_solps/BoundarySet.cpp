#include"Variable.h"
#include"BoundarySet.h"


void ContinuityMatrixBoundary() {
	//core plasma (GG')
	for (int i = int((l_percentage_private_flux * nx + 1) * ny); i < ((1 - l_percentage_private_flux) * nx - 1) * ny;) {                 //c[i]=1,source[i]=0;代表边界值固定为0
		c[i] = 1;
		source[i] = 3.0e19;
		i += ny;
	}
	//target conditions	(AB and CD)
	for (int i = 1; i < ny - 1; i++) {
		c[i] = 1;                                       //通过将边界值设置位相邻层的值，实现一阶导数为0
		source[i] = density[1][i];
	}
	for (int i = 1; i < ny - 1; i++) {
		c[(nx - 1) * ny + i] = 1;
		source[(nx - 1) * ny + i] = density[nx - 2][i];
	}
	//private flux	(AE and E'C) and wall (BD)
	for (int i = 0; i < (l_percentage_private_flux * nx - 1) * ny;) {
		c[i] = 1;
		source[i] = 3.0e18;
		i += ny;
	}
	for (int i = int(((1 - l_percentage_private_flux) * nx + 1) * ny); i <= (nx - 1) * ny;) {             //??????????????????
		c[i] = 1;
		source[i] = 3.0e18;
		i += ny;
	}
	for (int i = ny - 1; i <= nx * ny - 1; ) {
		c[i] = 1;
		source[i] = 3.0e18;
		i += ny;
	}
	//内部四条边界线分别设置
	for (int i = 0; i - (w_percentage_private_flux * ny - 1) <= 1e-7; i++) {
		c[int(l_percentage_private_flux * nx * ny + i)] = 1;                                              //core plasma (GG')两个边界，采用对称边界条件，等于相邻值
		source[int(l_percentage_private_flux * nx * ny + i)] = density[int(l_percentage_private_flux * nx + 1)][i];

		c[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = 1;
		source[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = density[int((1 - l_percentage_private_flux) * nx - 2)][i];

		c[int((l_percentage_private_flux * nx - 1) * ny + i)] = 1;                                              //private flux (FE F'E')两个边界，采用对称边界条件，等于相邻值
		source[int((l_percentage_private_flux * nx - 1) * ny + i)] = density[int(l_percentage_private_flux * nx - 2)][i];

		c[int((1 - l_percentage_private_flux) * nx * ny + i)] = 1;
		source[int((1 - l_percentage_private_flux) * nx * ny + i)] = density[int((1 - l_percentage_private_flux) * nx + 1)][i];

	}

}

void MomentumMatrixBoundary() {
	//core plasma (GG')
	for (int i = int(l_percentage_private_flux * nx + 1); i < (1 - l_percentage_private_flux) * nx - 1; i++) {                 
		c[i * ny] = 1.0e-5;
		source[i * ny] = velocity_parallel_i[i][1]*1.0e-5;
	}
	//target conditions	(AB and CD)
	for (int i = 1; i < ny - 1; i++) {
		c[i] = 1.0e-5;
		source[i] = -6.1899383e4* 1.0e-5;
	}
	for (int i = 1; i < ny - 1; i++) {
		c[(nx - 1) * ny + i] = 1.0e-5;
		source[(nx - 1) * ny + i] = 6.1899383e4* 1.0e-5;
	}
	//private flux	(AE and E'C) and wall (BD)
	for (int i = 0; i < l_percentage_private_flux * nx - 1;i++) {
		c[i*ny] = 1.0e-5;
		source[i * ny] = velocity_parallel_i[i][1]* 1.0e-5;
	}
	for (int i = int((1 - l_percentage_private_flux) * nx + 1); i <= nx - 1;i++) {
		c[i*ny] = 1.0e-5;
		source[i * ny] = velocity_parallel_i[i][1]* 1.0e-5;
	}
	for (int i = 0; i <= nx - 1;i++ ) {
		c[(i+1)*ny-1] = 1.0e-5;
		source[(i + 1) * ny - 1] = velocity_parallel_i[i][ny-2]* 1.0e-5;
	}
	//内部四条边界线分别设置
	for (int i = 0; i - (w_percentage_private_flux * ny - 1) <= 1e-7; i++) {
		c[int(l_percentage_private_flux * nx * ny + i)] = 1.0e-5;                                              //core plasma (GG')两个边界，采用对称边界条件，等于相邻值
		//source[int(l_percentage_private_flux * nx * ny + i)] = velocity_parallel_i[int(l_percentage_private_flux * nx + 1)][i]* 1.0e-5;
		source[int(l_percentage_private_flux * nx * ny + i)] = 0;

		c[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = 1.0e-5;
		//source[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = velocity_parallel_i[int((1 - l_percentage_private_flux) * nx - 2)][i]* 1.0e-5;
		source[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = 0;

		c[int((l_percentage_private_flux * nx - 1) * ny + i)] = 1.0e-5;                                              //private flux (FE F'E')两个边界，采用对称边界条件，等于相邻值
		//source[int((l_percentage_private_flux * nx - 1) * ny + i)] = velocity_parallel_i[int(l_percentage_private_flux * nx /*- 2*/)][i]* 1.0e-5;
		source[int((l_percentage_private_flux * nx - 1) * ny + i)] = 0;

		c[int((1 - l_percentage_private_flux) * nx * ny + i)] = 1.0e-5;
		//source[int((1 - l_percentage_private_flux) * nx * ny + i)] = velocity_parallel_i[int((1 - l_percentage_private_flux) * nx /*+ 1*/-1)][i]* 1.0e-5;
		source[int((1 - l_percentage_private_flux) * nx * ny + i)] = 0;
	}
}

void ElectronEnergyMatrixBoundary() {
	//core plasma (GG')
	for (int i = int((l_percentage_private_flux * nx + 1) * ny); i < ((1 - l_percentage_private_flux) * nx - 1) * ny;) {                 //c[i]=1,source[i]=0;代表边界值固定为0
		c[i] = 1.0/electron;
		source[i] = 400 ;
		i += ny;
	}
	//target conditions	(AB and CD)
	for (int i = 1; i < ny - 1; i++) {
		c[i] = 1.0/electron;                                       //通过将边界值设置位相邻层的值，实现一阶导数为0
		source[i] = 40 ;
	}
	for (int i = 1; i < ny - 1; i++) {
		c[(nx - 1) * ny + i] = 1.0/electron;
		source[(nx - 1) * ny + i] = 40 ;
	}
	//private flux	(AE and E'C) and wall (BD)
	for (int i = 0; i < (l_percentage_private_flux * nx - 1) * ny;) {
		c[i] = 1.0/electron;
		source[i] = 40 ;
		i += ny;
	}
	for (int i = int(((1 - l_percentage_private_flux) * nx + 1) * ny); i <= (nx - 1) * ny;) {
		c[i] = 1.0/electron;
		source[i] = 40 ;
		i += ny;
	}
	for (int i = ny - 1; i <= nx * ny - 1; ) {
		c[i] = 1.0/electron;
		source[i] = 40 ;
		i += ny;
	}
	//内部四条边界线分别设置
	for (int i = 0; i - (w_percentage_private_flux * ny - 1) <= 1e-7; i++) {
		c[int(l_percentage_private_flux * nx * ny + i)] = 1.0/electron;                                              //core plasma (GG')两个边界，采用对称边界条件，等于相邻值
		source[int(l_percentage_private_flux * nx * ny + i)] = temperature_e[int(l_percentage_private_flux * nx + 1)][i]/electron;

		c[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = 1.0 / electron;
		source[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = temperature_e[int((1 - l_percentage_private_flux) * nx - 2)][i] / electron;

		c[int((l_percentage_private_flux * nx - 1) * ny + i)] = 1.0 / electron;                                              //private flux (FE F'E')两个边界，采用对称边界条件，等于相邻值
		source[int((l_percentage_private_flux * nx - 1) * ny + i)] = temperature_e[int(l_percentage_private_flux * nx - 2)][i] / electron;

		c[int((1 - l_percentage_private_flux) * nx * ny + i)] = 1.0 / electron;
		source[int((1 - l_percentage_private_flux) * nx * ny + i)] = temperature_e[int((1 - l_percentage_private_flux) * nx + 1)][i] / electron;
	}
}

void IonEnergyMatrixBoundary() {
	//core plasma (GG')
	for (int i = int((l_percentage_private_flux * nx + 1) * ny); i < ((1 - l_percentage_private_flux) * nx - 1) * ny;) {                 //c[i]=1,source[i]=0;代表边界值固定为0
		c[i] = 1/electron;
		source[i] = 400;
		i += ny;
	}
	//target conditions	(AB and CD)
	for (int i = 1; i < ny - 1; i++) {
		c[i] = 1/electron;                                       //通过将边界值设置位相邻层的值，实现一阶导数为0
		source[i] = 40;
	}
	for (int i = 1; i < ny - 1; i++) {
		c[(nx - 1) * ny + i] = 1/electron;
		source[(nx - 1) * ny + i] = 40;
	}
	//private flux	(AE and E'C) and wall (BD)
	for (int i = 0; i < (l_percentage_private_flux * nx - 1) * ny;) {
		c[i] = 1/electron;
		source[i] = 40;
		i += ny;
	}
	for (int i = int(((1 - l_percentage_private_flux) * nx + 1) * ny); i <= (nx - 1) * ny;) {
		c[i] = 1/electron;
		source[i] = 40;
		i += ny;
	}
	for (int i = ny - 1; i <= nx * ny - 1; ) {
		c[i] = 1/electron;
		source[i] = 40;
		i += ny;
	}
	//内部四条边界线分别设置
	for (int i = 0; i - (w_percentage_private_flux * ny - 1) <= 1e-7; i++) {
		c[int(l_percentage_private_flux * nx * ny + i)] = 1/electron;                                              //core plasma (GG')两个边界，采用对称边界条件，等于相邻值
		source[int(l_percentage_private_flux * nx * ny + i)] = temperature_i[int(l_percentage_private_flux * nx + 1)][i]/electron;

		c[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = 1/electron;
		source[int(((1 - l_percentage_private_flux) * nx - 1) * ny + i)] = temperature_i[int((1 - l_percentage_private_flux) * nx - 2)][i]/electron;

		c[int((l_percentage_private_flux * nx - 1) * ny + i)] = 1/electron;                                              //private flux (FE F'E')两个边界，采用对称边界条件，等于相邻值
		source[int((l_percentage_private_flux * nx - 1) * ny + i)] = temperature_i[int(l_percentage_private_flux * nx - 2)][i]/electron;

		c[int((1 - l_percentage_private_flux) * nx * ny + i)] = 1/electron;
		source[int((1 - l_percentage_private_flux) * nx * ny + i)] = temperature_i[int((1 - l_percentage_private_flux) * nx + 1)][i]/electron;
	}

}