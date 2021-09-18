#pragma once

//求解模型大小，网格设置
const double lx = 20.0, ly = 0.1;
const int nx = 150, ny = 30;
const double dx = lx / (nx-1);
const double dy = ly / (ny-1);

const double width_private_flux = 0.04;                              //定义私有通量区宽与长
const double length_private_flux = 3.2;
const double w_percentage_private_flux = width_private_flux / ly;   //0.4 下面宽
const double l_percentage_private_flux = length_private_flux / lx;  //0.16 右边长

extern double delta_t;              //时间间隔
extern double end_time;               //首次总收敛时间

extern double pitch;                        //bx=Bx/B
extern double m_i;                //离子质量,kg
extern double Z_i;                            //离子电荷数
extern double m_e;               //电子质量
extern double permittivity;       //真空电容率
extern double M_PI;      //pi
extern double electron;         //单位电荷，C
extern double K;                 //eV转K

extern double a[], b[], c[], d[], e[];       //求解矩阵五对角线定义，分别代表点的南、东、西、北侧
extern double source[];                                                  //矩阵求解Ax=b中的源项          a,b,c,d,e  * x = source
extern double x[];

extern double D_n[];
extern double D_s[];
extern double F_n[];
extern double F_s[];
extern double D_w[];
extern double D_e[];
extern double F_w[];
extern double F_e[];

//定义变量存储矩阵，密度、离子平行速度、压力、温度
extern double density[][ny];
extern double density_old[][ny];

extern double velocity_x[][ny];
extern double velocity_y[][ny];
extern double velocity_parallel_i[][ny];
extern double velocity_parallel_i_old[][ny];

extern double pressure_e[][ny];
extern double pressure_i[][ny];
extern double pressure[][ny];

extern double temperature_e[][ny];
extern double temperature_i[][ny];

