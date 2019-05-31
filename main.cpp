#include <iostream>
#include <math.h>
#include <iomanip> //控制输出精度
using namespace std;
double getpi(double lambda,double k);
double getbeta(double Ma2,double k);
double getalpha(double Ma,double k);
double lambda_0,lambda_1,lambda_2w;
double k2,k,k0 ;
//循环优化参数
float alpha_1 = 3.1415926/9,beta_2 = 3.1415926/9,beta_2pre=beta_2,alpha_1pre = alpha_1; // pi /(180/20)
float sigma_R = 1,sigma_Rpre = sigma_R,sigma_s = 1,sigma_spre = sigma_s;//循环优化参数

double p_1,p_2,w_2;
double T_0,p_0;
double psi;
double m ;

int main()
{
    //
    double T_0x = 1144.9,p_0x=521276;
    //
    double p_1x;
    double T_1x = T_0x;
    //
    lambda_1 = 0.5;
    //
    double q_lambda_2w,q_lambda1,q_lambda0,p_1;
    double K2,K,K0;
    double T_2pre,T_2 = 330;
    double T_1pre,T_1 = 1200  ;//初始静温 待使用参数
    double T_0pre,T_0 = 1200;


    T_1pre = 1105;
                double C = 1.011554*pow(1.8*T_1pre,7)*1e-25 -
                           1.452677*pow(1.8*T_1pre,6)*1e-21 +
                           7.6215767*1e-18*pow(1.8*T_1pre,5) -
                           1.5128259*1e-14*pow(1.8*T_1pre,4) -
                           6.67178376 * 1e-12 *pow(1.8*T_1pre,3) +
                           6.5519486 *1e-8*pow(1.8*T_1pre,2) -
                           5.1536879 * 1e-5 * 1.8*T_1pre + 0.25020051;

                k = C/(C-29.56738/437);
                cout<<k<<endl;



    while(lambda_1<1)
    {
        while(1) //静叶sigma 和 alpha 的循环过程
        {
            do
            {
                T_1pre = T_1;
                double C = 1.011554*pow(1.8*T_1pre,7)*1e-25 -
                           1.452677*pow(1.8*T_1pre,6)*1e-21 +
                           7.6215767*1e-18*pow(1.8*T_1pre,5) -
                           1.5128259*1e-14*pow(1.8*T_1pre,4) -
                           6.67178376 * 1e-12 *pow(1.8*T_1pre,3) +
                           6.5519486 *1e-8*pow(1.8*T_1pre,2) -
                           5.1536879 * 1e-5 * 1.8*T_1pre + 0.25020051;

                k = C/(C-29.56738/437);
                cout<<"----k = "<<k<<endl;
                //气动函数
                q_lambda1 = pow((k+1)/2,1/(k-1)) * lambda_1 *
                            pow( ( 1 - (k-1)/(k+1)*pow(lambda_1,2) ),1/(k-1));
                cout<<"----q_lambda1 = " << q_lambda1<<endl;
                //静叶出口静压
                p_1 = p_1x * pow( ( 1 - (k-1)/(k+1)*pow( lambda_1 , 2 ) ),1/(k-1));
                T_1 = T_1x *( 1 - (k-1) / (k+1)*lambda_1*lambda_1  );
                //
            }
            while(abs(T_1-T_1pre)>0.002);

            cout<<"T_1="<<T_1<<endl;
            K = sqrt( k/290.056 * pow(2/(k+1),(k+1)/(k-1)) );
            m = 0.0013046658 * K * p_1x* q_lambda1 *sin(alpha_1);

            double x_n1=0,x_n=0;
            do
            {
                T_0pre = T_0;
                double C0 = 1.011554*pow(1.8*T_0pre,7)*1e-25 -
                            1.452677*pow(1.8*T_0pre,6)*1e-21 +
                            7.6215767*1e-18*pow(1.8*T_0pre,5) -
                            1.5128259*1e-14*pow(1.8*T_0pre,4) -
                            6.67178376 * 1e-12 *pow(1.8*T_0pre,3) +
                            6.5519486 *1e-8*pow(1.8*T_0pre,2) -
                            5.1536879 * 1e-5 * 1.8*T_0pre + 0.25020051;

                k0 = C0/(C0-29.56738/437);
                cout<<"----k0 = "<<k0<<endl;
                K0 = sqrt( k0/290.056 * pow(2/(k0+1),(k0+1)/(k0-1)) );
                //
                q_lambda0 = 1.035813 * q_lambda1 * sin(alpha_1) * K / K0 *sigma_s;
                cout<<"----q_lambda0 = " << q_lambda0<<endl;
                // 迭代参数为 lambda0 方法采用牛顿迭代法
                double a = (k0+1)/2;
                double b = 1/(k0-1);
                double c = (k0-1)/(k0+1);
                double d = q_lambda0;

                for(int i =0; i<10000; i++)
                {
                    double funcxn = pow(a,b)*x_n*pow((1-c*x_n*x_n),b)-d;
                    double difffuncxn = pow(a,b)*pow((1-c*x_n*x_n),b)-
                                        2*pow(a,b)*c*x_n*x_n*b*pow((1-c*x_n*x_n),b-1);


                    x_n1 = x_n-funcxn/difffuncxn;
                    cout << "----迭代第" << i+1<<"次 求解λ0 = "<<setprecision(8) << x_n1 <<" 当期迭代偏差="<< setprecision(8)<<x_n1-x_n << endl;
                    //迭代终止条件
                    if(abs(x_n1-x_n)<=0.001)
                    {
                        cout<<"迭代求解结果 lambda_0 = "<<setprecision(8) <<x_n1<<endl;
                        break;
                    }
                    x_n=x_n1;
                }
                lambda_0 = x_n1;
                p_0 = p_0x * pow( 1 - (k0-1)/(k0+1)*lambda_0*lambda_0 , k0/(k0-1));

                T_0 = T_0x * (1 - (k0-1)/(k0+1) * lambda_0 *lambda_0 );

                cout <<"静温T_0 = "<<setprecision(8) << T_0 <<"pc"<<T_0-T_0pre << endl;
                cout<<endl;
            }
            while(abs(T_0-T_0pre)>0.002);
            //
            //        x_n1 为求得解

            // 静叶速度损失系数
            //1
            double A=110,B=1.63;
            if(alpha_1/3.1415926*180+90<=110)
                A =alpha_1/3.1415926*180+90;
            if(sin(0.5)/sin(alpha_1)<1.63)
                B = sin(0.5)/sin(alpha_1);
            double xi_1 = 0.02185 * (0.01065*A*A - 2.295*A+160.5)*(0.1055*B*B-0.3427*B+0.295);
            //2
            double xi_2 = 0.2*1.6/24.42/sin(alpha_1);

            double rho = p_0 / (290.056 * T_0);
            double v = lambda_0 * sqrt( 2*k0 /(k0+1) * 290.056 *T_0x );
            double T_qu = T_0 * 1e-3;
            double mu = (0.229*pow(T_qu,3) - 1.3333*T_qu*T_qu+4.89*T_qu +
                         5.05 - 0.275/60.24096) * 1e-5;
            double R_e = rho*v*0.036964/mu;
            double xi_3 = 2100 / R_e - 0.0021;

            double xi_4 = 0.09 *(1 - ( xi_1 + xi_2 +  xi_3 )) * (8.6667/98.6667)*(8.6667/98.6667);
            double xi_5 = 0.509965 * (xi_1 + xi_2 +  xi_3+xi_4);
            psi = sqrt(1-(xi_1 + xi_2  +  xi_3 + xi_4 + xi_5 ));
            // sigma_s
            sigma_spre = sigma_s;
            sigma_s = getpi(lambda_1/psi,k) / getpi(lambda_1,k);//使用lambda1已确认
            p_1x = p_0x * sigma_s;
            //计算气流角alpha1
            alpha_1pre = alpha_1;
            double Ma = lambda_1 * sqrt(2/((k+1)-(k-1)*lambda_1*lambda_1) );
            alpha_1 = getalpha(Ma,k);

            if(abs(sigma_spre-sigma_s)<=0.01 && abs(alpha_1pre-alpha_1)<=0.01)
            {
                cout<<"偏差="<<sigma_spre-sigma_s<<"  "<<alpha_1pre-alpha_1<<endl;
                break;
            }
            else
            {
                cout<<sigma_spre-sigma_s<<"  "<<alpha_1pre-alpha_1<<endl;
            }

        }
        cout<<"********************************************"<<endl;
        cout<<"T1="<<T_1<<"   T0="<<T_0<<endl;
        cout<<"psi="<<psi<<endl;
        cout<<"sigma_s = " <<sigma_s <<"   alpha1="<<alpha_1*180/3.1415926<<endl;
        cout<<"lambda_0= "<<lambda_0<<" q_lambda0="<<q_lambda0<<endl;
        cout<<"********************************************"<<endl;




        /// 动叶部分
        //c_0 c_1 w_1 w_2后续进行求解
        double u = 354.394;
        double c_0 = lambda_0 *sqrt( (2*k0/(k0+1))*290.056*T_1x );
        double c_1 = lambda_1 *sqrt( (2*k/(k+1))*290.056*T_0x );
        double w_1 = sqrt(c_1*c_1 + u*u - 2*u*c_1*cos(alpha_1));
        //流动损失
        //T_1 = T_0x * (1 - (k-1)/(k+1)*lambda_1*lambda_1); 已求
        double T_1wx = T_1x + (w_1*w_1-c_1*c_1) /(580.112*k/(k-1));
        p_1 = p_1x * pow(1-(k-1)/(k+1)*lambda_1*lambda_1,k/(k-1));
        double p_1wx = p_1 * pow(T_1wx/T_1,k/(k-1));

        //p_2 需要循环在下面进行
        double lambda_1w = w_1 / sqrt( 290.056 * 2* k *T_1wx /(k+1) );
        double q_lambda_1w = pow((k+1)/2,1/(k-1)) * lambda_1w *
                             pow(1-(k-1)/(k+1)*lambda_1w*lambda_1w,1/(k-1));

        double beta_1 = asin(c_1*sin(alpha_1)/w_1);
        //if(beta_1 < 3.1415926/2)
            //beta_1 =3.1415926/2+beta_1;
        double p_2wx = sigma_R*p_1wx;
        double T_2wx = T_1wx;
        while(1)
        {
            do
            {
                T_2pre = T_2;
                double C2 = 1.011554*pow(1.8*T_2pre,7)*1e-25 -
                            1.452677*pow(1.8*T_2pre,6)*1e-21 +
                            7.6215767*1e-18*pow(1.8*T_2pre,5) -
                            1.5128259*1e-14*pow(1.8*T_2pre,4) -
                            6.67178376 * 1e-12 *pow(1.8*T_2pre,3) +
                            6.5519486 *1e-8*pow(1.8*T_2pre,2) -
                            5.1536879 * 1e-5 * 1.8*T_2pre + 0.25020051;

                k2 = C2/(C2-29.56738/437);
                cout<<"----k2 = "<<k2<<endl;
                K2 = sqrt( k2/290.056 * pow(2/(k2+1),(k2+1)/(k2-1)) );
                //
                q_lambda_2w = 0.944448*K/K2* q_lambda_1w * sin(beta_1)/sin(beta_2) *p_1wx/p_2wx ;
                //
                double x_n1=0,x_n=0;
                double a = (k2+1)/2;
                double b = 1/(k2-1);
                double c = (k2-1)/(k2+1);
                double d = q_lambda_2w;

                for(int i =0; i<10000; i++)
                {
                    double funcxn = pow(a,b)*x_n*pow((1-c*x_n*x_n),b)-d;
                    double difffuncxn = pow(a,b)*pow((1-c*x_n*x_n),b)-
                                        2*pow(a,b)*c*x_n*x_n*b*pow((1-c*x_n*x_n),b-1);


                    x_n1 = x_n-funcxn/difffuncxn;
                    cout << "    迭代第" << i+1<<"次 求解λ2w = "<<setprecision(8) << x_n1 <<" 当期迭代偏差="<< setprecision(8)<<x_n1-x_n << endl;
                    //迭代终止条件
                    if(abs(x_n1-x_n)<=0.001)
                    {
                        cout<<"迭代求解结果λ2w = "<<setprecision(8) <<x_n1<<endl;
                        break;
                    }
                    x_n=x_n1;
                }
                lambda_2w = x_n1;
                //
                //p_2 = p_1wx * sigma_R * pow(1-(k2-1)/(k2+1)*lambda_2w*lambda_2w ,k2/(k2-1) ) ;
                p_2 =p_2wx * pow( T_2wx/T_2,k/(1-k) );
                w_2 = psi * sqrt(2 * k2/(k2-1) * 290.056 * T_1wx*
                                 (1 - pow(p_2 / p_1wx,(k2-1)/k2)));

                T_2 = T_2wx - w_2*w_2 / 580.112 / k2 *(k2-1);
                cout<<T_2<<"----"<<T_2pre<<"-----------"<<q_lambda_2w<<endl;
            }
            while(abs(T_2pre-T_2)>0.002);



            ///计算动叶流动损失
            double Ad=110,Bd=1.63;
            if((beta_1+beta_2)/3.1415926*180<=110)
                Ad =(beta_1+beta_2)/3.1415926*180;
            if(sin(beta_1)/sin(beta_2)<1.63)
                Bd = sin(beta_1)/sin(beta_2);
            double xi_1d = 0.02185 * (0.01065*Ad*Ad - 2.295*Ad+160.5)*(0.1055*Bd*Bd-0.3427*Bd+0.295);
            //2
            double xi_2d = 0.2 * 1.2 / 13.937 /sin(beta_2);

            double rho_d = p_1 / 290.056 / T_1;
            double v_d = lambda_2w * sqrt( 2*k/(k+1) * 290.056 * T_1x);
            double T_d = 1e-3 * T_1;
            double mu_d = (0.229*pow(T_d,3)-1.3333*T_d*T_d+4.89*T_d+5.05-0.275/60.24096)*1e-5;
            double R_e_d = rho_d * v_d *0.020536 /mu_d;
            double xi_3d = 2100/R_e_d - 0.0021;

            double C_d=0.15;
            if(beta_1/3.1415926*180<45)
                C_d = 0.09;
            //double param = asin(c_1*sin(alpha_1)/w_1);
            //double xi_4d = C_d*(1-xi_1d - xi_2d-xi_3d)*(1-param/40)*(1-param/40);
            double xi_4d = C_d*(1-xi_1d - xi_2d-xi_3d)*(beta_1/3.1415926*180-40)/40*(beta_1/3.1415926*180-40)/40;


            double xi_5d = 0.398618 * (xi_1d + xi_2d+xi_3d+xi_4d);

            double xi_6d = 0.6 * 12.359 * 1.3/36.275/15.184/sin(beta_2);
            //输出psi_d 动叶损失系数
            double psi_d = sqrt(1-xi_1d-xi_2d-xi_3d-xi_4d-xi_5d-xi_6d);

            ///计算气流角beta2
            sigma_Rpre = sigma_R;
            sigma_R = getpi(lambda_2w/psi_d,k2) / getpi(lambda_2w,k2);

            ///
            p_2wx = sigma_R * p_1wx;
            double Ma2 = lambda_2w*sqrt(2/((k2+1)-(k2-1)*lambda_2w*lambda_2w));
            beta_2pre = beta_2;
            beta_2 = getbeta(Ma2,k2);

            if(abs(beta_2-beta_2pre)<=0.02 && abs(sigma_R-sigma_Rpre)<=0.02)
                break;
        }
        cout<<"*************************************"<<endl;
        cout<<"λ1="<<lambda_1<<"----λ2w="<<lambda_2w<<endl;
        cout<<"\t"<<"\t m="<<m<<endl;
        cout<<"\t"<<"\t c1="<<c_1<<endl;
        cout<<"\t"<<"\t c0="<<c_0<<endl;
        cout<<"\t"<<"\t w1="<<w_1<<endl;
        cout<<"\t"<<"\t w2="<<w_1<<endl;
        //cout<<"\t"<<"\t α0="<<alpha_0<<endl;
        cout<<"\t"<<"\t α1="<<alpha_1<<endl;
        cout<<"\t"<<"\t β1="<<beta_1<<endl;
        cout<<"\t"<<"\t β2="<<beta_2<<endl;
        cout<<"\t"<<"\t λ0="<<lambda_0<<endl;
        cout<<"\t"<<"\t λ1="<<lambda_1<<endl;
        cout<<"\t"<<"\t λ1w="<<lambda_1w<<endl;
        cout<<"\t"<<"\t λ2w="<<lambda_2w<<endl;
        cout<<"\t"<<"\t T0*="<<T_0x<<endl;
        cout<<"\t"<<"\t T0="<<T_0<<endl;
        cout<<"\t"<<"\t T1="<<T_1<<endl;
        cout<<"\t"<<"\t T2="<<T_2<<endl;
        cout<<"\t"<<"\t p0*="<<p_0x<<endl;
        cout<<"\t"<<"\t p1*="<<p_1x<<endl;
        //cout<<"\t"<<"\t p2*="<<p_2*<<endl;
        cout<<"\t"<<"\t p0="<<p_0<<endl;
        cout<<"\t"<<"\t p1="<<p_1<<endl;
        cout<<"\t"<<"\t p2="<<p_2<<endl;
        cout<<"\t"<<"\t p1w*="<<p_1wx<<endl;
        cout<<"\t"<<"\t p2w*="<<p_2wx<<endl;
        //cout<<"m="<<m<<c_1<<c_0<<w_1<<w_2<<alpha_1<<beta_1<<beta_2<<
        cout<<"*************************************"<<endl;
        system("pause");
        if(lambda_1>=1)
            break;
        else
            lambda_1+=0.1;
    }




    ///lambda==0.1
    if(lambda_1 == 1)
    {
        do
        {
            T_1pre = T_1;
            double C = 1.011554*pow(1.8*T_1pre,7)*1e-25 -
                       1.452677*pow(1.8*T_1pre,6)*1e-21 +
                       7.6215767*1e-18*pow(1.8*T_1pre,5) -
                       1.5128259*1e-14*pow(1.8*T_1pre,4) -
                       6.67178376 * 1e-12 *pow(1.8*T_1pre,3) +
                       6.5519486 *1e-8*pow(1.8*T_1pre,2) -
                       5.1536879 * 1e-5 * 1.8*T_1pre + 0.25020051;

            k = C/(C-29.56738/437);
            cout<<"----k = "<<k<<endl;
            //气动函数
            //静叶出口静压
            p_1 = p_1x * pow( ( 1 - (k-1)/(k+1)*pow( lambda_1 , 2 ) ),1/(k-1));
            T_1 = T_1x *( 1 - (k-1) / (k+1)*lambda_1*lambda_1  );
            //
        }
        while(abs(T_1-T_1pre)>0.002);

        cout<<"T_1="<<T_1<<endl;
        K = sqrt( k/290.056 * pow(2/(k+1),(k+1)/(k-1)) );
        q_lambda1 = 1;
        T_0=330;
        do
        {
            T_0pre = T_0;
            double C0 = 1.011554*pow(1.8*T_0pre,7)*1e-25 -
                        1.452677*pow(1.8*T_0pre,6)*1e-21 +
                        7.6215767*1e-18*pow(1.8*T_0pre,5) -
                        1.5128259*1e-14*pow(1.8*T_0pre,4) -
                        6.67178376 * 1e-12 *pow(1.8*T_0pre,3) +
                        6.5519486 *1e-8*pow(1.8*T_0pre,2) -
                        5.1536879 * 1e-5 * 1.8*T_0pre + 0.25020051;

            k0 = C0/(C0-29.56738/437);
            cout<<"----k0 = "<<k0<<endl;
            //气动函数
            q_lambda1 = pow((k0+1)/2,1/(k0-1)) * lambda_1 *
                        pow( ( 1 - (k0-1)/(k0+1)*pow(lambda_1,2) ),1/(k0-1));
            cout<<"----q_lambda1 = " << q_lambda1<<endl;
            q_lambda0 = 1.035813 * q_lambda1 * sin(alpha_1) * K / K0 *sigma_s;
            cout<<"----q_lambda0 = " << q_lambda0<<endl;
            // 迭代参数为 lambda0 方法采用牛顿迭代法
            double a = (k0+1)/2;
            double b = 1/(k0-1);
            double c = (k0-1)/(k0+1);
            double d = q_lambda0;

            double x_n=0,x_n1=0;
            for(int i =0; i<10000; i++)
            {
                double funcxn = pow(a,b)*x_n*pow((1-c*x_n*x_n),b)-d;
                double difffuncxn = pow(a,b)*pow((1-c*x_n*x_n),b)-
                                    2*pow(a,b)*c*x_n*x_n*b*pow((1-c*x_n*x_n),b-1);


                x_n1 = x_n-funcxn/difffuncxn;
                cout << "----迭代第" << i+1<<"次 求解λ0 = "<<setprecision(8) << x_n1 <<" 当期迭代偏差="<< setprecision(8)<<x_n1-x_n << endl;
                //迭代终止条件
                if(abs(x_n1-x_n)<=0.001)
                {
                    cout<<"迭代求解结果 lambda_0 = "<<setprecision(8) <<x_n1<<endl;
                    break;
                }
                x_n=x_n1;
            }
            lambda_0 = x_n1;
            //静叶出口静压
            p_0 = 521276 * pow( ( 1 - (k0-1)/(k0+1)*pow( lambda_1 , 2 ) ),1/(k0-1));
            T_0 = T_1x *( 1 - (k0-1) / (k0+1)*lambda_1*lambda_1  );
            //
        }
        while(abs(T_0-T_0pre)>0.002);
        K0 = sqrt( k0/290.056 * pow(2/(k0+1),(k0+1)/(k0-1)) );
        m = 680.108797 * K* sigma_s *sin(alpha_1);

    }

    ///lambda_1>1
    if(lambda_1>1)
    {
        double lambda_1s[4] = {1.05,1.1,1.15,1.2};
        for(int i=0; i<4; i++)
        {
            lambda_1 = lambda_1s[i];
            sigma_s = getpi(lambda_1/psi,k)/getpi(lambda_1,k);
            double Ma = lambda_1 * sqrt(2/((k+1)-(k-1)*lambda_1*lambda_1) );
            alpha_1 = getalpha(Ma,k);
            p_1 = p_1x *pow((1 - (k-1)/(k+1)*lambda_1*lambda_1),k/(k-1));
            T_1 = T_1x *(1 - (k-1)/(k+1)*lambda_1*lambda_1);

            //
            double u = 354.394;
            double c_0 = lambda_0 *sqrt( (2*k0/(k0+1))*290.056*T_1x );
            double c_1 = lambda_1 *sqrt( (2*k/(k+1))*290.056*T_0x );
            double w_1 = sqrt(c_1*c_1 + u*u - 2*u*c_1*cos(alpha_1));
            //流动损失
            //T_1 = T_0x * (1 - (k-1)/(k+1)*lambda_1*lambda_1); 已求
            double T_1wx = T_1x + (w_1*w_1-c_1*c_1) /(580.112*k/(k-1));
            p_1 = p_1x * pow(1-(k-1)/(k+1)*lambda_1*lambda_1,k/(k-1));
            double p_1wx = p_1 * pow(T_1wx/T_1,(k-1)/k);

            //p_2 需要循环在下面进行
            double lambda_1w = w_1 / sqrt( 290.056 * 2* k *T_1wx /(k+1) );
            double q_lambda_1w = pow((k+1)/2,1/(k-1)) * lambda_1w *
                                 pow(1-(k-1)/(k+1)*lambda_1w*lambda_1w,1/(k-1));

            double beta_1 = asin(c_1*sin(alpha_1)/w_1);
            //if(beta_1 < 3.1415926/2)
            //beta_1 =3.1415926/2+beta_1;
            double p_2wx = sigma_R*p_1wx;
            double T_2wx = T_1wx;
            while(1)
            {
                do
                {
                    T_2pre = T_2;
                    double C2 = 1.011554*pow(1.8*T_2pre,7)*1e-25 -
                                1.452677*pow(1.8*T_2pre,6)*1e-21 +
                                7.6215767*1e-18*pow(1.8*T_2pre,5) -
                                1.5128259*1e-14*pow(1.8*T_2pre,4) -
                                6.67178376 * 1e-12 *pow(1.8*T_2pre,3) +
                                6.5519486 *1e-8*pow(1.8*T_2pre,2) -
                                5.1536879 * 1e-5 * 1.8*T_2pre + 0.25020051;

                    k2 = C2/(C2-29.56738/437);
                    cout<<"----k2 = "<<k2<<endl;
                    K2 = sqrt( k2/290.056 * pow(2/(k2+1),(k2+1)/(k2-1)) );
                    //
                    q_lambda_2w = 0.944448*K/K2* q_lambda_1w * sin(beta_1)/sin(beta_2) *p_1wx/p_2wx ;
                    //
                    double x_n1=0,x_n=0;
                    double a = (k2+1)/2;
                    double b = 1/(k2-1);
                    double c = (k2-1)/(k2+1);
                    double d = q_lambda_2w;

                    for(int i =0; i<10000; i++)
                    {
                        double funcxn = pow(a,b)*x_n*pow((1-c*x_n*x_n),b)-d;
                        double difffuncxn = pow(a,b)*pow((1-c*x_n*x_n),b)-
                                            2*pow(a,b)*c*x_n*x_n*b*pow((1-c*x_n*x_n),b-1);


                        x_n1 = x_n-funcxn/difffuncxn;
                        cout << "    迭代第" << i+1<<"次 求解λ2w = "<<setprecision(8) << x_n1 <<" 当期迭代偏差="<< setprecision(8)<<x_n1-x_n << endl;
                        //迭代终止条件
                        if(abs(x_n1-x_n)<=0.001)
                        {
                            cout<<"迭代求解结果λ2w = "<<setprecision(8) <<x_n1<<endl;
                            break;
                        }
                        x_n=x_n1;
                    }
                    lambda_2w = x_n1;
                    if(lambda_2w>=1)
                        break;
                    //
                    p_2 = p_1wx * sigma_R * pow(1-(k2-1)/(k2+1)*lambda_2w*lambda_2w ,k2/(k2-1) ) ;
                    w_2 = psi * sqrt(2 * k2/(k2-1) * 290.056 * T_1wx*
                                     (1 - pow(p_2 / p_1wx,(k2-1)/k2)));

                    T_2 = T_2wx - w_2*w_2 / 580.112 / k2 *(k2-1);
                    cout<<T_2<<"----"<<T_2pre<<"-----------"<<q_lambda_2w<<endl;
                }
                while(abs(T_2pre-T_2)>0.002);
                ///计算动叶流动损失
                double Ad=110,Bd=1.63;
                if((beta_1+beta_2)/3.1415926*180<=110)
                    Ad =(beta_1+beta_2)/3.1415926*180;
                if(sin(beta_1)/sin(beta_2)<1.63)
                    Bd = sin(beta_1)/sin(beta_2);
                double xi_1d = 0.02185 * (0.01065*Ad*Ad - 2.295*Ad+160.5)*(0.1055*Bd*Bd-0.3427*Bd+0.295);
                //2
                double xi_2d = 0.2 * 1.2 / 13.937 /sin(beta_2);

                double rho_d = p_1 / 290.056 / T_1;
                double v_d = lambda_2w * sqrt( 2*k/(k+1) * 290.056 * T_1x);
                double T_d = 1e-3 * T_1;
                double mu_d = (0.229*pow(T_d,3)-1.3333*T_d*T_d+4.89*T_d+5.05-0.275/60.24096)*1e-5;
                double R_e_d = rho_d * v_d *0.020536 /mu_d;
                double xi_3d = 2100/R_e_d - 0.0021;

                double C_d=0.15;
                if(beta_1/3.1415926*180<45)
                    C_d = 0.09;
                //double param = asin(c_1*sin(alpha_1)/w_1);
                //double xi_4d = C_d*(1-xi_1d - xi_2d-xi_3d)*(1-param/40)*(1-param/40);
                double xi_4d = C_d*(1-xi_1d - xi_2d-xi_3d)*(beta_1/3.1415926*180-40)/40*(beta_1/3.1415926*180-40)/40;


                double xi_5d = 0.398618 * (xi_1d + xi_2d+xi_3d+xi_4d);

                double xi_6d = 0.6 * 12.359 * 1.3/36.275/15.184/sin(beta_2);
                //输出psi_d 动叶损失系数
                double psi_d = sqrt(1-xi_1d-xi_2d-xi_3d-xi_4d-xi_5d-xi_6d);

                ///计算气流角beta2
                sigma_Rpre = sigma_R;
                sigma_R = getpi(lambda_2w/psi_d,k2) / getpi(lambda_2w,k2);

                ///
                p_2wx = sigma_R * p_1wx;
                double Ma2 = lambda_2w*sqrt(2/((k2+1)-(k2-1)*lambda_2w*lambda_2w));
                beta_2pre = beta_2;
                beta_2 = getbeta(Ma2,k2);

                if(abs(beta_2-beta_2pre)<=0.02 && abs(sigma_R-sigma_Rpre)<=0.02)
                    break;
            }

        }
         double p_1wx;
         double T_1wx;
         double T_2wx;
        if(lambda_2w==1)
        {
            if(q_lambda_2w>1)//如果大于1
            {
                //q_lambda_2w==1

            }
            sigma_R = 1;
            beta_2 = 20/180*3.1415926;
            double u = 354.394;
            double c_0 = lambda_0 *sqrt( (2*k0/(k0+1))*290.056*T_1x );
            double c_1 = lambda_1 *sqrt( (2*k/(k+1))*290.056*T_0x );
            double w_1 = sqrt(c_1*c_1 + u*u - 2*u*c_1*cos(alpha_1));
            //流动损失
            //T_1 = T_0x * (1 - (k-1)/(k+1)*lambda_1*lambda_1); 已求
            T_1wx = T_1x + (w_1*w_1-c_1*c_1) /(580.112*k/(k-1));
            p_1 = p_1x * pow(1-(k-1)/(k+1)*lambda_1*lambda_1,k/(k-1));
            p_1wx = p_1 * pow(T_1wx/T_1,(k-1)/k);

            //p_2 需要循环在下面进行
            double lambda_1w = w_1 / sqrt( 290.056 * 2* k *T_1wx /(k+1) );
            double q_lambda_1w = pow((k+1)/2,1/(k-1)) * lambda_1w *
                                 pow(1-(k-1)/(k+1)*lambda_1w*lambda_1w,1/(k-1));

            double beta_1 = asin(c_1*sin(alpha_1)/w_1);
            double p_2wx = sigma_R*p_1wx;
            T_2wx = T_1wx;
            double m_2= K*sigma_R*p_1wx/sqrt(T_2wx)*0.0467481*sin(beta_2);
        }
        if(lambda_2w>1)
        {
            double lambda_2ws[4]={1.05,1.1,1.15,1.2};
            p_2 = p_1wx * sigma_R * pow(1-(k2-1)/(k2+1)*lambda_2w*lambda_2w ,k2/(k2-1) ) ;
            w_2 = psi * sqrt(2 * k2/(k2-1) * 290.056 * T_1wx*
            (1 - pow(p_2 / p_1wx,(k2-1)/k2)));

            T_2 = T_2wx - w_2*w_2 / 580.112 / k2 *(k2-1);
            ///
            cout<<"---"<<endl;

        }


    }
    return 0;
}
double getpi(double lambda,double k)
{
    double pi = pow(1-(k-1)/(k+1)*lambda*lambda,k/(k-1));
    return pi;
}
double getbeta(double Ma2,double k)
{
    if(Ma2<=0.5)
    {
        double r= 34.11367;
        return 3.1415926*r/180;
    }
    if(Ma2 == 1)
    {
        double r= 31.49331;
        return 3.1415926*r/180;
    }
    if( Ma2>0.5 && Ma2 <1)
    {
        double r = 31.49331 + 5.24072*(1-Ma2);
        return 3.1415926*r/180;
    }
    if( Ma2>1 )
    {
        double qMa2 = pow((k+1)/2,1/(k-1))*Ma2*pow(1-(k-1)/(k+1)*Ma2*Ma2,1/(k-1));
        double r= 0.24431+asin(0.518758/qMa2)/3.1415926*180;
        return 3.1415926*r/180;
    }
}
double getalpha(double Ma,double k)
{
    if(Ma<=0.5)
        return 3.1415926*22.926798/180;
    if(Ma == 1)
        return 3.1415926*28.26336/180;
    if( Ma>0.5 && Ma <1)
    {
        double r= 28.26336 + 10.673124 *(Ma - 1);
        return 3.1415926*r/180;
    }

    if( Ma >1 )
    {
        double qMa = pow((k+1)/2,1/(k-1))*Ma*pow(1-(k-1)/(k+1)*Ma*Ma,1/(k-1));
        double r= 7.08456+asin(0.36128/qMa)/3.1415926*180;
        return 3.1415926*r/180;
    }
}
