#include"Variable.h"

double delta_t = 0.0001;              //ʱ����
double end_time = 10.0;               //�״�������ʱ��

double pitch = 0.1;                        //bx=Bx/B
double m_i = 1.6726231e-27;                //��������,kg
double Z_i = 1;                            //���ӵ����
double m_e = 9.10938356e-31;               //��������
double permittivity = 8.8541878e-12;       //��յ�����
double M_PI = 3.14159265358979323846;      //pi
double electron = 1.602176634e-19;         //��λ��ɣ�C
double K = 1.0;                        //1eV=11605K,����������¶ȵ�λ�����Ӧ�û���eV

double a[nx * ny], b[nx * ny], c[nx * ny], d[nx * ny], e[nx * ny];       //��������Խ��߶��壬�ֱ�������ϡ�������������
double source[nx * ny];                                                  //�������Ax=b�е�Դ��          a,b,c,d,e  * x = source
double x[nx * ny];

double D_n[nx * ny];
double D_s[nx * ny];
double F_n[nx * ny];
double F_s[nx * ny];
double D_w[nx * ny];
double D_e[nx * ny];
double F_w[nx * ny];
double F_e[nx * ny];

//��������洢�����ܶȡ�����ƽ���ٶȡ�ѹ�����¶�
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