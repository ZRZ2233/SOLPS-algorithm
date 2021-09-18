#include"Variable.h"

double delta_t = 0.0001;              //时间间隔
double end_time = 10.0;               //首次总收敛时间

double pitch = 0.1;                        //bx=Bx/B
double m_i = 1.6726231e-27;                //离子质量,kg
double Z_i = 1;                            //离子电荷数
double m_e = 9.10938356e-31;               //电子质量
double permittivity = 8.8541878e-12;       //真空电容率
double M_PI = 3.14159265358979323846;      //pi
double electron = 1.602176634e-19;         //单位电荷，C
double K = 1.0;                        //1eV=11605K,尝试里面的温度单位，最后应该还是eV

double a[nx * ny], b[nx * ny], c[nx * ny], d[nx * ny], e[nx * ny];       //求解矩阵五对角线定义，分别代表点的南、东、西、北侧
double source[nx * ny];                                                  //矩阵求解Ax=b中的源项          a,b,c,d,e  * x = source
double x[nx * ny];

double D_n[nx * ny];
double D_s[nx * ny];
double F_n[nx * ny];
double F_s[nx * ny];
double D_w[nx * ny];
double D_e[nx * ny];
double F_w[nx * ny];
double F_e[nx * ny];

//定义变量存储矩阵，密度、离子平行速度、压力、温度
double density[nx][ny];
double density_old[nx][ny];

double velocity_x[nx][ny];
double velocity_y[nx][ny];
double velocity_parallel_i[nx][ny];
double velocity_parallel_i_old[nx][ny];

double pressure_e[nx][ny];
double pressure_i[nx][ny];
double pressure[nx][ny];

double temperature_e[nx][ny];
double temperature_i[nx][ny];