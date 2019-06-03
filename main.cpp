#include <iostream>
#include<stdlib.h>
#include <math.h>

#include <iomanip> //控制输出精度

using namespace std;
//pi函数
double getpi(double lambda,double k);
double getbeta(double Ma2,double k);//计算beta
double getalpha(double Ma,double k);//计算alpha
double lambda_0,lambda_1,lambda_2w,lambda_2;
double k,K ;
double c_2;
//待输出参数
float alpha_1 = 3.1415926/9,beta_2 = 3.1415926/9,beta_2pre=beta_2,alpha_1pre = alpha_1,beta_1,alpha_2; // pi /(180/20)
float sigma_R = 1,sigma_Rpre = sigma_R,sigma_s = 1,sigma_spre = sigma_s;//循环优化参数
//分别是λ2w大于1 , qλ2w第一次大于1,qλ2w第二次大于1的,qλ2w等于1的标志
bool islambda2w=false,isqlambda2wOnce=false,isqlambda2wTwice=false,isqlambda2weq=false;
double lambda_1foot = 0.1;
double p_1,p_2,w_2,w_1,c_1,c_0;
double T_0,p_0,p_1wx,p_2wx,T_1wx,T_2wx,T_2x,lambda_1w;
double psi,pi_Tx,p_2x,eta_Tx,L_T,X_m,X_L,X_n;
double m ;
//100转速
double u = 354.394,n=16955;

int main()
{
    //
    double A_1 = 0.044145,A_0 = 0.04268845;
    islambda2w = false;
    double T_0x = 1144.9,p_0x=521276;
    //
    double p_1x = sigma_s*p_0x;
    double T_1x = T_0x;
    //
    lambda_1 = 0.5;
    //
    double q_lambda_2w,q_lambda1,q_lambda0,p_1;
    double T_2pre,T_2 = 330;
    double T_1 = 1200  ;//初始静温 待使用参数
    double T_0 = 1200,T_k;

    T_k = 1105;
    double C = 1.011554*pow(1.8*T_k,7)*1e-25 -
    1.452677*pow(1.8*T_k,6)*1e-21 +
    7.6215767*1e-18*pow(1.8*T_k,5) -
    1.5128259*1e-14*pow(1.8*T_k,4) -
    6.67178376 * 1e-12 *pow(1.8*T_k,3) +
    6.5519486 *1e-8*pow(1.8*T_k,2) -
    5.1536879 * 1e-5 * 1.8*T_k + 0.25020051;
    k = C/(C-29.56738/437);
    K = sqrt( k/290.056 * pow(2/(k+1),(k+1)/(k-1)) );

    //lambda 从0.5开始递增 大于1跳出循环
    while(lambda_1<1)
    {
        //qλ1
        q_lambda1 = pow((k+1)/2,1/(k-1)) * lambda_1 *
        pow( ( 1 - (k-1)/(k+1)*pow(lambda_1,2) ),1/(k-1));

        m = 0.0013046658 * K * p_1x* q_lambda1 *sin(alpha_1);

        p_1 = p_1x * pow( ( 1 - (k-1)/(k+1)*pow( lambda_1 , 2 ) ),1/(k-1));
        T_1 = T_1x *( 1 - (k-1) / (k+1)*lambda_1*lambda_1  );

      //静叶 迭代求解sigma_s 和 alpha 的循环过程
        while(1)
        {
            double x_n1=0,x_n=0;

            q_lambda0 = A_1/A_0* q_lambda1 * sin(alpha_1) *sigma_s;

            // 迭代参数为 lambda0 方法采用牛顿迭代法
            double a = (k+1)/2;
            double b = 1/(k-1);
            double c = (k-1)/(k+1);
            double d = q_lambda0;

            for(int i =0; i<10000; i++)
            {
                double funcxn = pow(a,b)*x_n*pow((1-c*x_n*x_n),b)-d;
                double difffuncxn = pow(a,b)*pow((1-c*x_n*x_n),b)-
                2*pow(a,b)*c*x_n*x_n*b*pow((1-c*x_n*x_n),b-1);


                x_n1 = x_n-funcxn/difffuncxn;
                //cout << "----迭代第" << i+1<<"次 求解λ0 = "<<setprecision(8) << x_n1 <<" 当期迭代偏差="<< setprecision(8)<<x_n1-x_n << endl;
                //迭代终止条件
                if(fabs(x_n1-x_n)<=0.001)
                {
                    //cout<<"迭代求解结果 lambda_0 = "<<setprecision(8) <<x_n1<<endl;
                    break;
                }
                x_n=x_n1;
            }
            lambda_0 = x_n1;
            //x_n1 为求得解 牛顿迭代完成


            // 静叶速度损失系数
            double A=110,B=1.63;
            if(alpha_1/3.1415926*180+90<=110)
                A =alpha_1/3.1415926*180+90;
            if(sin(0.5)/sin(alpha_1)<1.63)
                B = sin(0.5)/sin(alpha_1);
            double xi_1 = 0.02185 * (0.01065*A*A - 2.295*A+160.5)*(0.1055*B*B-0.3427*B+0.295);
            //2
            double xi_2 = 0.2*1.6/24.42/sin(alpha_1);

            double rho = p_1 / (290.056 * T_1);
            double v = lambda_0 * sqrt( 2*k /(k+1) * 290.056 *T_0x );
            double T_qu = T_1 * 1e-3;
            double mu = (0.229*pow(T_qu,3) - 1.3333*T_qu*T_qu+4.89*T_qu +
                         5.05 - 0.275/60.24096) * 1e-5;
            double R_e = rho*v*0.036964/mu;
            double xi_3 = 2100 / R_e - 0.0021;

            double xi_4 = 0.09 *(1 - ( xi_1 + xi_2 +  xi_3 )) * (8.6667/98.6667)*(8.6667/98.6667);
            double xi_5 = 0.509965 * (xi_1 + xi_2 +  xi_3+xi_4);
            psi = sqrt(1-(xi_1 + xi_2  +  xi_3 + xi_4 + xi_5 ));
            // sigma_spre，alpha_1pre用与存放前一次的sigma_s和alpha_1的值
            sigma_spre = sigma_s;
            sigma_s = getpi(lambda_1/psi,k) / getpi(lambda_1,k);
            p_1x = p_0x * sigma_s;
            //计算气流角alpha1
            alpha_1pre = alpha_1;
            double Ma = lambda_1 * sqrt(2/((k+1)-(k-1)*lambda_1*lambda_1) );
            alpha_1 = getalpha(Ma,k); //调用函数求解α1
            //当前alpha和前一次alpha相差小于0.01
            //并且当前sigma_spre和前一次sigma_s相差小于0.01
            //跳出迭代循环
            if(fabs(sigma_spre-sigma_s)<=0.01 && fabs(alpha_1pre-alpha_1)<=0.01)
            {
                //cout<<"偏差="<<sigma_spre-sigma_s<<"  "<<alpha_1pre-alpha_1<<endl;
                break;
            }

        }
        p_0 = p_0x * pow( 1 - (k-1)/(k+1)*lambda_0*lambda_0 , k/(k-1));
        T_0 = T_0x * (1 - (k-1)/(k+1) * lambda_0 *lambda_0 );
        cout<<"********************************************"<<endl;
        cout<<"静叶计算完成"<<endl;
        cout<<"T1="<<T_1<<"   T0="<<T_0<<endl;
        cout<<"psi="<<psi<<endl;
        cout<<"m="<<m<<endl;
        //cout<<"sigma_s = " <<sigma_s <<"   alpha1="<<alpha_1*180/3.1415926<<endl;
        //cout<<"lambda_0= "<<lambda_0<<" q_lambda0="<<q_lambda0<<endl;
        cout<<"********************************************"<<endl;

        /// 动叶部分

        c_0 = lambda_0 *sqrt( (2*k/(k+1))*290.056*T_1x );
        c_1 = lambda_1 *sqrt( (2*k/(k+1))*290.056*T_0x );
        double c_1u = c_1 * cos(alpha_1);
        double w_1u = c_1u - u;
        double c_1a = c_1*sin(alpha_1);
        w_1 = sqrt(c_1a*c_1a+w_1u*w_1u);

        T_1wx = T_1 + w_1 *w_1/(2*290.056*k/(k-1));
        p_1wx = p_1 * pow(T_1wx/T_1,k/(k-1));

        lambda_1w = w_1 / sqrt( 290.056 * 2* k *T_1wx /(k+1) );
        double q_lambda_1w = pow((k+1)/2,1/(k-1)) * lambda_1w *
                             pow(1-(k-1)/(k+1)*lambda_1w*lambda_1w,1/(k-1));

        //gbeta_1为求解beta_1的中间参数

        //版本原因atan求解负值会出错
        double gbeta_1 = c_1a/w_1u;
        if(gbeta_1<0)
            beta_1 = - atan(gbeta_1);
        else
        beta_1 = atan(gbeta_1);

        p_2wx = sigma_R*p_1wx;
        T_2wx = T_1wx;

        //
        while(1)
        {
            q_lambda_2w =  m*sqrt(T_2wx)/0.0467418/K/p_2wx/sin(beta_2);
            if(fabs(q_lambda_2w-1.0000)<=0.01) //如果等于1
            {
                isqlambda2weq=true;//=1的标志
                break;//跳出循环
            }
            if(q_lambda_2w>1) //如果大于1
            {
                isqlambda2wTwice = true;//大于1的标志
                break;//跳出循环至
            }


            //牛顿迭代法
            double x_n1=0.5,x_n=0.5;
            double a = (k+1)/2;
            double b = 1/(k-1);
            double c = (k-1)/(k+1);
            double d = q_lambda_2w;
            for(int i =0; i<10000; i++)
            {
                double funcxn = pow(a,b)*x_n*pow((1-c*x_n*x_n),b)-d;
                double difffuncxn = pow(a,b)*pow((1-c*x_n*x_n),b)-
                2*pow(a,b)*c*x_n*x_n*b*pow((1-c*x_n*x_n),b-1);

                x_n1 = x_n-funcxn/difffuncxn;
                //cout << "    迭代第" << i+1<<"次 求解λ2w = "<<setprecision(8) << x_n1 <<" 当期迭代偏差="<< setprecision(8)<<x_n1-x_n << endl;
                //迭代终止条件
                if(fabs(x_n1-x_n)<=0.001)
                {
                    //cout<<"迭代求解结果λ2w = "<<setprecision(8) <<x_n1<<endl;
                    break;
                }
                x_n=x_n1;
                if(x_n>=1)
                {
                    islambda2w=true; //λ2w大于1的标志
                    break;//λ2w大于1跳出牛顿迭代的循环 至222行
                }
            }
            lambda_2w = x_n1;
            //牛顿迭代法至此结束

            if(islambda2w)
            {
                break;//λ2w大于1继续 跳出迭代更新sigma_R和beta2的循环 至280行
            }

            T_2 = (1-(k-1)/(k+1)*lambda_2w*lambda_2w)*T_2wx;
            p_2 = p_2wx * pow( T_2wx/T_2,k/(1-k) );

            w_2 = psi * sqrt(2 * k/(k-1) * 290.056 * T_1wx*
            (1 - pow(p_2 / p_1wx,(k-1)/k)));


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

            ///计算气流角beta2和sigma_R
            sigma_Rpre = sigma_R;
            sigma_R = getpi(lambda_2w/psi_d,k) / getpi(lambda_2w,k);

            p_2wx = sigma_R * p_1wx;//更新p_2wx的值

            double Ma2 = lambda_2w*sqrt(2/((k+1)-(k-1)*lambda_2w*lambda_2w));
            beta_2pre = beta_2;
            beta_2 = getbeta(Ma2,k);

            if(fabs(beta_2-beta_2pre)<=0.01 && fabs(sigma_R-sigma_Rpre)<=0.01)//满足要求跳出循环
                break;
        }

        if(lambda_1>=1||islambda2w||isqlambda2weq) //λ2w大于1或者 qλ2w=1 或者 λ1大于1 跳出 λ1+0.1的循环 至424行
            break;
        else if(isqlambda2wTwice) //如果qλ2w大于1 减小λ1的值 步长减半
        {
            lambda_1foot = lambda_1foot/2;
            lambda_1-=lambda_1foot;
        }else  //满足要求继续运算
        {
            double c_2u = w_2 * cos(beta_2) - u;
            double w_2a = w_2*sin(beta_2);

            c_2 = sqrt(w_2a*w_2a+c_2u*c_2u);
            alpha_2 = acos((w_2*cos(beta_2)-u)/c_2);

            alpha_2 = 3.1415926 - alpha_2 ;
            T_2x = T_2 + c_2*c_2/(2*290.056*k)*(k-1);
            //T_2x = T_2wx + (u*u - 2*w_2*u*cos(beta_2))/580.112/k*(k-1);
            lambda_2 = c_2/sqrt(2*290.056*k/(k+1)*T_2x);
            p_2x = p_2 / pow((1-(k-1)/(k+1)*lambda_2*lambda_2),k/(k-1));
            pi_Tx = p_0x/p_2x;

            eta_Tx = 0.975*(1-T_2x/T_1x)/(1-1/pow(pi_Tx,(k-1)/k));
            L_T = 332085.1144*k/(k-1)*(1-1/pow(pi_Tx,(k-1)/k))*eta_Tx;
            X_m = m*sqrt(T_1x)/p_1x;
            X_L = L_T/T_1x;
            X_n = n/T_1x;

            cout<<"*************************************"<<endl;
            cout<<"λ1="<<lambda_1<<"----λ2w="<<lambda_2w<<endl;
            cout<<"\t"<<"\t c2="<<c_2<<endl;
            cout<<"\t"<<"\t λ2="<<lambda_2<<endl;
            cout<<"\t"<<"\t m="<<m<<endl;
            cout<<"\t"<<"\t c1="<<c_1<<endl;
            cout<<"\t"<<"\t c0="<<c_0<<endl;
            cout<<"\t"<<"\t w1="<<w_1<<endl;
            cout<<"\t"<<"\t w2="<<w_2<<endl;
            cout<<"\t"<<"\t α2="<<alpha_2 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t α1="<<alpha_1 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t β1="<<beta_1 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t β2="<<beta_2 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t λ0="<<lambda_0<<endl;
            cout<<"\t"<<"\t λ1="<<lambda_1<<endl;
            cout<<"\t"<<"\t λ1w="<<lambda_1w<<endl;
            cout<<"\t"<<"\t λ2w="<<lambda_2w<<endl;
            cout<<"\t"<<"\t T0*="<<T_0x<<endl;
            cout<<"\t"<<"\t T0="<<T_0<<endl;
            cout<<"\t"<<"\t T1="<<T_1<<endl;
            cout<<"\t"<<"\t T2="<<T_2<<endl;
            cout<<"\t"<<"\t T2*="<<T_2x<<endl;
            cout<<"\t"<<"\t p0*="<<p_0x<<endl;
            cout<<"\t"<<"\t p1*="<<p_1x<<endl;
            cout<<"\t"<<"\t p2*="<<p_2x<<endl;
            cout<<"\t"<<"\t p0="<<p_0<<endl;
            cout<<"\t"<<"\t p1="<<p_1<<endl;
            cout<<"\t"<<"\t p2="<<p_2<<endl;
            cout<<"\t"<<"\t p1w*="<<p_1wx<<endl;
            cout<<"\t"<<"\t p2w*="<<p_2wx<<endl;
            cout<<"\t"<<"\t πT*="<<pi_Tx<<endl;
            cout<<"\t"<<"\t ηT*="<<eta_Tx<<endl;
            cout<<"\t"<<"\t LT="<<L_T<<endl;
            cout<<"\t"<<"\t Xm*="<<X_m<<endl;
            cout<<"\t"<<"\t XL="<<X_L<<endl;
            cout<<"\t"<<"\t Xn="<<X_n<<endl;
            //cout<<"m="<<m<<c_1<<c_0<<w_1<<w_2<<alpha_1<<beta_1<<beta_2<<
            cout<<"*************************************"<<endl;
            //system("pause");
            lambda_1+=lambda_1foot; //λ1+0.1
        }
    }
    //qlmabda2w==1时处理
    if(isqlambda2weq)
    {
        double lambda_2ws[4]={1.05,1.1,1.15,1.2};
        for (int i=0;i<4;i++)
        {
            lambda_2w = lambda_2ws[i]; //依次给lambda_2w赋值1.05,1.1,1.15,1.2
            p_2 = p_1wx * sigma_R * pow(1-(k-1)/(k+1)*lambda_2w*lambda_2w ,k/(k-1) ) ;
            w_2 = psi * sqrt(2 * k/(k-1) * 290.056 * T_1wx*
            (1 - pow(p_2 / p_1wx,(k-1)/k)));
            T_2 = T_2wx - w_2*w_2 / 580.112 / k *(k-1);

            double c_2u = w_2 * cos(beta_2) - u;
            double w_2a = w_2*sin(beta_2);

            c_2 = sqrt(w_2a*w_2a+c_2u*c_2u);
            alpha_2 = acos((w_2*cos(beta_2)-u)/c_2);

            alpha_2 = 3.1415926 - alpha_2;
            T_2x = T_2 + c_2*c_2/(2*290.056*k)*(k-1);
            //T_2x = T_2wx + (u*u - 2*w_2*u*cos(beta_2))/580.112/k*(k-1);
            lambda_2 = c_2/sqrt(2*290.056*k/(k+1)*T_2x);
            p_2x = p_2 / pow((1-(k-1)/(k+1)*lambda_2*lambda_2),k/(k-1));
            pi_Tx = p_0x/p_2x;

            eta_Tx = 0.975*(1-T_2x/T_1x)/(1-1/pow(pi_Tx,(k-1)/k));
            L_T = 332085.1144*k/(k-1)*(1-1/pow(pi_Tx,(k-1)/k))*eta_Tx;
            X_m = m*sqrt(T_1x)/p_1x;
            X_L = L_T/T_1x;
            X_n = n/T_1x;

            alpha_2 = acos((w_2*cos(beta_2)-u)/c_2);
            alpha_2 = 3.1415926- alpha_2;
            cout<<"*************************************"<<endl;
            cout<<"λ1="<<lambda_1<<"----λ2w="<<lambda_2w<<endl;
            cout<<"\t"<<"\t λ2="<<lambda_2<<endl;
            cout<<"\t"<<"\t c2="<<c_2<<endl;
            cout<<"\t"<<"\t m="<<m<<endl;
            cout<<"\t"<<"\t c1="<<c_1<<endl;
            cout<<"\t"<<"\t c0="<<c_0<<endl;
            cout<<"\t"<<"\t w1="<<w_1<<endl;
            cout<<"\t"<<"\t w2="<<w_2<<endl;
            cout<<"\t"<<"\t α2="<<alpha_2 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t α1="<<alpha_1 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t β1="<<beta_1 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t β2="<<beta_2 / 3.1415926 * 180<<endl;
            cout<<"\t"<<"\t λ0="<<lambda_0<<endl;
            cout<<"\t"<<"\t λ1="<<lambda_1<<endl;
            cout<<"\t"<<"\t λ1w="<<lambda_1w<<endl;
            cout<<"\t"<<"\t λ2w="<<lambda_2w<<endl;
            cout<<"\t"<<"\t T0*="<<T_0x<<endl;
            cout<<"\t"<<"\t T0="<<T_0<<endl;
            cout<<"\t"<<"\t T1="<<T_1<<endl;
            cout<<"\t"<<"\t T2="<<T_2<<endl;
            cout<<"\t"<<"\t T2*="<<T_2x<<endl;
            cout<<"\t"<<"\t p0*="<<p_0x<<endl;
            cout<<"\t"<<"\t p1*="<<p_1x<<endl;
            cout<<"\t"<<"\t p2*="<<p_2x<<endl;
            cout<<"\t"<<"\t p0="<<p_0<<endl;
            cout<<"\t"<<"\t p1="<<p_1<<endl;
            cout<<"\t"<<"\t p2="<<p_2<<endl;
            cout<<"\t"<<"\t p1w*="<<p_1wx<<endl;
            cout<<"\t"<<"\t p2w*="<<p_2wx<<endl;
            cout<<"\t"<<"\t πT*="<<pi_Tx<<endl;
            cout<<"\t"<<"\t ηT*="<<eta_Tx<<endl;
            cout<<"\t"<<"\t LT="<<L_T<<endl;
            cout<<"\t"<<"\t Xm*="<<X_m<<endl;
            cout<<"\t"<<"\t XL="<<X_L<<endl;
            cout<<"\t"<<"\t Xn="<<X_n<<endl;
            //cout<<"m="<<m<<c_1<<c_0<<w_1<<w_2<<alpha_1<<beta_1<<beta_2<<
            cout<<"*************************************"<<endl;
        }

    }

    ///lambda==0.1
    if(lambda_1 == 1)
    {
        //气动函数
        //静叶出口静压
        /*p_1 = p_1x * pow( ( 1 - (k-1)/(k+1)*pow( lambda_1 , 2 ) ),1/(k-1));
        T_1 = T_1x *( 1 - (k-1) / (k+1)*lambda_1*lambda_1  );
        //
        cout<<"T_1="<<T_1<<endl;
        K = sqrt( k/290.056 * pow(2/(k+1),(k+1)/(k-1)) );
        q_lambda1 = 1;
        T_0=330;

        //气动函数
        q_lambda1 = pow((k+1)/2,1/(k-1)) * lambda_1 *
        pow( ( 1 - (k-1)/(k+1)*pow(lambda_1,2) ),1/(k-1));
        cout<<"----q_lambda1 = " << q_lambda1<<endl;
        q_lambda0 = 1.035813 * q_lambda1 * sin(alpha_1) *sigma_s;
        cout<<"----q_lambda0 = " << q_lambda0<<endl;
        // 迭代参数为 lambda0 方法采用牛顿迭代法
        double a = (k+1)/2;
        double b = 1/(k-1);
        double c = (k-1)/(k+1);
        double d = q_lambda0;

        double x_n=0,x_n1=0;
        for(int i =0; i<10000; i++)
        {
            double funcxn = pow(a,b)*x_n*pow((1-c*x_n*x_n),b)-d;
            double difffuncxn = pow(a,b)*pow((1-c*x_n*x_n),b)-
            2*pow(a,b)*c*x_n*x_n*b*pow((1-c*x_n*x_n),b-1);


            x_n1 = x_n-funcxn/difffuncxn;
            //cout << "----迭代第" << i+1<<"次 求解λ0 = "<<setprecision(8) << x_n1 <<" 当期迭代偏差="<< setprecision(8)<<x_n1-x_n << endl;
            //迭代终止条件
            if(fabs(x_n1-x_n)<=0.001)
            {
                break;
            }
            x_n=x_n1;
        }
        lambda_0 = x_n1;
        //静叶出口静压
        p_0 = 521276 * pow( ( 1 - (k-1)/(k+1)*pow( lambda_1 , 2 ) ),1/(k-1));
        T_0 = T_1x *( 1 - (k-1) / (k+1)*lambda_1*lambda_1  );
        //
        m = 680.108797 * K* sigma_s *sin(alpha_1);*/

    }

    ///lambda_1>1
    if(lambda_1>1)
    {
        double lambda_1s[4] = {1.05,1.1,1.15,1.2};
        for(int i=0; i<4; i++)
        {
            lambda_1 = lambda_1s[i];
          /*  sigma_s = getpi(lambda_1/psi,k)/getpi(lambda_1,k);
            double Ma = lambda_1 * sqrt(2/((k+1)-(k-1)*lambda_1*lambda_1) );
            alpha_1 = getalpha(Ma,k);
            p_1 = p_1x *pow((1 - (k-1)/(k+1)*lambda_1*lambda_1),k/(k-1));
            T_1 = T_1x *(1 - (k-1)/(k+1)*lambda_1*lambda_1);

            //
            double u = 354.394;
            double c_0 = lambda_0 *sqrt( (2*k/(k+1))*290.056*T_1x );
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

                //
                q_lambda_2w = 0.944448* q_lambda_1w * sin(beta_1)/sin(beta_2) *p_1wx/p_2wx ;
                //
                double x_n1=0,x_n=0;
                double a = (k+1)/2;
                double b = 1/(k-1);
                double c = (k-1)/(k+1);
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
                p_2 = p_1wx * sigma_R * pow(1-(k-1)/(k+1)*lambda_2w*lambda_2w ,k/(k-1) ) ;
                w_2 = psi * sqrt(2 * k/(k-1) * 290.056 * T_1wx*
                                 (1 - pow(p_2 / p_1wx,(k-1)/k)));

                T_2 = T_2wx - w_2*w_2 / 580.112 / k *(k-1);
                cout<<T_2<<"----"<<T_2pre<<"-----------"<<q_lambda_2w<<endl;

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
                sigma_R = getpi(lambda_2w/psi_d,k) / getpi(lambda_2w,k);

                ///
                p_2wx = sigma_R * p_1wx;
                double Ma2 = lambda_2w*sqrt(2/((k+1)-(k-1)*lambda_2w*lambda_2w));
                beta_2pre = beta_2;
                beta_2 = getbeta(Ma2,k);

                if(abs(beta_2-beta_2pre)<=0.02 && abs(sigma_R-sigma_Rpre)<=0.02)
                    break;
            }*/
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
