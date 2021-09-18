#pragma once

//���ģ�ʹ�С����������
const double lx = 20.0, ly = 0.1;
const int nx = 150, ny = 30;
const double dx = lx / (nx-1);
const double dy = ly / (ny-1);

const double width_private_flux = 0.04;                              //����˽��ͨ�������볤
const double length_private_flux = 3.2;
const double w_percentage_private_flux = width_private_flux / ly;   //0.4 �����
const double l_percentage_private_flux = length_private_flux / lx;  //0.16 �ұ߳�

extern double delta_t;              //ʱ����
extern double end_time;               //�״�������ʱ��

extern double pitch;                        //bx=Bx/B
extern double m_i;                //��������,kg
extern double Z_i;                            //���ӵ����
extern double m_e;               //��������
extern double permittivity;       //��յ�����
extern double M_PI;      //pi
extern double electron;         //��λ��ɣ�C
extern double K;                 //eVתK

extern double a[], b[], c[], d[], e[];       //��������Խ��߶��壬�ֱ�������ϡ�������������
extern double source[];                                                  //�������Ax=b�е�Դ��          a,b,c,d,e  * x = source
extern double x[];

extern double D_n[];
extern double D_s[];
extern double F_n[];
extern double F_s[];
extern double D_w[];
extern double D_e[];
extern double F_w[];
extern double F_e[];

//��������洢�����ܶȡ�����ƽ���ٶȡ�ѹ�����¶�
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

