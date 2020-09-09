#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<NTL/ZZ.h>
#include<NTL/RR.h>
#include<NTL/vector.h>
#include<NTL/vec_vec.h>
#include<NTL/xdouble.h>
#include<NTL/quad_float.h>
#include<NTL/matrix.h>
#include<NTL/mat_ZZ.h>

using namespace std;
using namespace NTL;

MAT<ZZ> RandomMatrix(long n,long m,ZZ q) //模q的随机矩阵生成
{
  RR y;
  ZZ x;
  Mat<ZZ> A;
  A.SetDims(n,m);
  srand((unsigned)time(NULL));
  for(long i=1;i<=n;i++){
    for(long j=1;j<=m;j++)
    {
      x=rand();
      x=x%q;
      A(i,j)=x;
    }
  }
  return A;
}


ZZ SampleZ(RR,c,RR s,RR t)//拒绝抽样技术
{
  RR left,right,condition,pi,y;
  pi=3.1415926;
  ZZ a,b,x;

  left=c-s*(log(t)/log(2));
  right=c+s*(log(t)/log(2));

  do
  {
    x=rand();
    FloorToZZ(a,right-left);
    x=x%a;
    FloorToZZ(b,left);
    x=x+b;
    ConvPrec(y,x,8);
    condition=exp(-1*pi*(y-c)*(y-c)/(s*s));
  }while ((condition<=0.0)&&(condition>1.0));
  
  return x;
  
}

int main(){
  long n,w,m1,m,k,singal;
  ZZ q,x,bina;
  RR y;
  float k1,t;
  Mat<ZZ> AR,A1,A2,A,R,G,I,RI,VER,e,u,g,h,Bn,V,X,Ax;
  cin >> n;
  cin >> q;
  k1=log(q)/log(2);

  k=ceil(k1);
  w=n*k;
  srand((unsigned)time(NULL));
  m1=rand();
  m1=m1%(w-n-1);
  m1=m1+n+1;
  m=w+m1;
  A1.SetDims(n,m1);
  A.SetDims(n,m1);
  R.SetDims(m1,w);
  G.SetDims(n,w);
  AR.SetDims(n,w);
  A2.SetDims(n,w);
  Bn.SetDims(w,1);
  I.SetDims(w,w);
  RI.SetDims(m,w);
  VER.SetDims(n,w);
  e.SetDims(m,1);
  V.SetDims(m,1);
  u.SetDims(n,1);
  Ax.SetDims(n,1);
  g.SetDims(n,1);
  h.SetDims(n,1);

  //生成随机矩阵A1
  A1=RandomMatrix(n,m1,q);
  //cout <<A1<<'\n';

  //陷门R用拒绝抽样抽取，参照典型的子高斯分布
  for(long i=1;i<=m1;i++)
    {
      for(long j=1;j<=w;j++)
      {
        x=SampleZ(RR(1),RR(15),RR(20));
        x=x%q;
        R(i,j)=x;
      }
    }
    //cout <<R<<'\n';

    for(long i=1;i<=w;i++)
    {
      for(long j=1;j<=w;j++)
      {
        if(i==j)
        {
          I(i,j)=1;
        }
        else
        {
          I(i,j)=0;
        }
      }   
    }
    //cout <<I<<'\n';

    for(long i=1;i<n;i++)    
      {
        for(long j=1;j<=w;j++)
        {
          if((j<i*k)&&(j>(i-1)*k))
          {
            G(i,j)=pow(2,(j-(i-1)));
          }
          else
          {
            G(i,j)=0;
          }
        }
      }
    //cout <<G<<'\n';

    mul(AR,A1,R);
    for(long i=1;i<n;i++)    
      {
        for(long j=1;j<=w;j++)
        {
          AR(i,j)=AR(i,j)%q;
        }
      }
    //cout <<AR<<'\n';

    sub(A2,G,AR);
    for(long i=1;i<n;i++)    
      {
        for(long j=1;j<=w;j++)
        {
          A2(i,j)=A2(i,j)%q;
        }
      }
    //cout <<A2<<'\n';

    for(long i=1;i<=n;i++)
    {
      for(long j=1;j<=m;j++)
        {
          if(j<=m1)
          {
            A(i,j)=A1(i,j);
          }
          else
          {
            A(i,j)=A2(i,j-m1);
          }
        }   
    }
    //cout <<A<<'\n';

    for(long i=1;i<=m;i++)
      {
      for(long j=1;j<=w;j++)
        {
          if(i<=m1)
          {
            RI(i,j)=R(i,j);
          }
          else
          {
            RI(i,j)=I(i-m1,j);
          }
        }   
      }
    //cout <<RI<<'\n';

    mul(VER,A,RI);
    for(long i=1;i<=n;i++)
      {
      for(long j=1;j<=w;j++)
        {
          VER(i,j)=VER(i,j)%q;
        }   
      }
    //cout <<VER<<'\n';

    //盲化 对任意的消息m，经过安全Hash函数g(m)得到固定长度的
    //伪随机数，暂且用随机向量代替

    for(long i=1;i<=m;i++)
    {
      x=SampleZ(RR(1),RR(15),RR(20));
      x=x%q;
      e(i,1)=x;
    }
    //cout <<e<<'\n';

    mul(u,A,e);
    for(long i=1;i<=n;i++)
      {
        u(i,1)=u(i,1)%q;
      }
    //cout <<u<<'\n';

    for(long i=1;i<=n;i++)
      {
        x=SampleZ(RR(1),RR(15),RR(20));
        x=x%q;
        g(i,1)=x;
      }
    //cout <<g<<'\n';

    add(h,g,u);
    for(long i=1;i<=n;i++)
    {
      h(i,1)=h(i,1)%q;
    }
    //cout <<h<<'\n';

    //签名
    for(long i=1;i<=n;i++)
    {
      bina=h(i,1);
      for(long j=1;j<=k;j++){

        DivRem(bina,Bn(j+(i-1)*k,1),bina,ZZ(2));
      }
    }
    //cout <<Bn<<'\n';

    mul(V,RI,Bn);
    for(long i=1;i<=m;i++){
      V(i,1)=V(i,1)%q;
    }
    //cout <<V<<'\n';

    //去盲
    sub(X,V,e);
    for(long i=1;i<m;i++)
    {
      X(i,1)=X(i,1)%q;
    }
    //cout <<X<<'\n';
    
    //验证
    mul(Ax,A,X);
    for(long i=1;i<=n;i++)
    {
      Ax(i,1)=Ax(i,1)%q;
    }
    //cout <<Ax<<'\n';
    //cout <<g<<'\n';

    for(long i=1;i<=n;i++)
    {
      if(Ax(i,1)!=g(i,1))
      {
        signal=0;
        break;
      }
      else
      {
        signal=1;
      }
    } 

    cout<<"keyGen,生成校验矩阵A和陷门矩阵R"<<'\n';
    cout<<"输出校验矩阵A"<<'\n';
    cout<<"A"<<'\n';
    cout<<A<<'\n';
    cout<<"输出陷门矩阵R"<<'\n';
    cout< <"R"<<'\n';
    cout<<R<<'\n';
    cout<<"运行盲化算法"<<'\n';
    cout<<"g(m): g()为 Hash函数，m为任意长度的消息"<<'\n';
    cout<<"g(m)"<<'\n';
    cout<<g<<'\n'; .
    cout<<"输出盲化后的消息h"<<'\n';
    cout<<"h"<<'\n';
    cout<<h<<'\n';
    cout<<"对盲化后的消息h,运行签名算法，输出签名"<<'\n';
    cout<<"v"<<'\n';
    cout<<V<<'\n';
    cout<<"对签名x，运行去盲算法，输出去盲后的签名"<<'\n';
    cout<<"x"<<'\n';
    cout<<X<<'\n';
    cour<<"运行验证算法"<<'\n';
    cout<<"A*x"<<'\n';
    cout<<Ax<<'\n';
    cout<<"g(m)"<<'\n';
    cout<<g<<'\n';
    if(singal==0)
      {
        cout<<"验证失败"<<'\n';
      }
    else
      {
        cout<<"验证成功"<<'\n';
      }
  }
}