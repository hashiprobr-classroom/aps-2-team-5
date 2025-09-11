#include <math.h>
#include "fourier.h"

void normalize(complex s[], int n){
    for(int k=0;k<n;k++){
        s[k].a/=n;s[k].b/=n;}
}
void nft(complex s[],complex t[],int n,int sign){
    for(int k=0;k<n;k++){
        t[k].a=0; t[k].b=0;
        for(int j=0;j<n;j++){
            double x=sign*2*PI*k*j/n;
            double cosx=cos(x),sinx=sin(x);
            t[k].a+=s[j].a*cosx-s[j].b*sinx;
            t[k].b+=s[j].a*sinx+s[j].b*cosx;
        }
    }
}

void nft_forward(complex s[],complex t[],int n){
    nft(s,t,n,-1);}
void nft_inverse(complex t[],complex s[],int n){
    nft(t,s,n,1);normalize(s,n);}

void fft(complex s[],complex t[],int n,int sign){
    if(n==1){t[0]=s[0];return;}
    int m=n/2;

    complex p[m],i[m],tp[m],ti[m];
    for(int x=0;x<m;x++){
        p[x]=s[2*x];i[x]=s[2*x+1];}
    fft(p,tp,m,sign);
    fft(i,ti,m,sign);
    for(int k=0;k<m;k++){
        double ang=sign*-2.0*PI*k/n;
        double ca=cos(ang),sa=sin(ang);
        complex w;w.a=ca;w.b=sa;
        complex mult;
        mult.a=w.a*ti[k].a-w.b*ti[k].b;
        mult.b=w.a*ti[k].b+w.b*ti[k].a;

        t[k].a=tp[k].a+mult.a;
        t[k].b=tp[k].b+mult.b;
        t[k+m].a=tp[k].a-mult.a;
        t[k+m].b=tp[k].b-mult.b;
    }
}

void fft_forward(complex s[],complex t[],int n){fft(s,t,n,-1);}

void fft_inverse(complex t[],complex s[],int n){fft(t,s,n,1);normalize(s,n);}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE],int width,int height){
    complex s[MAX_SIZE],t[MAX_SIZE];

    for(int y=0;y<height;y++){
        for(int x=0;x<width;x++){
            s[x]=matrix[y][x];}
        fft_forward(s,t,width);
        for(int x=0;x<width;x++){
            matrix[y][x]=t[x];}
    }
    for(int x=0;x<width;x++){
        for(int y=0;y<height;y++){
            s[y]=matrix[y][x];}
        fft_forward(s,t,height);
        for(int y=0;y<height;y++){
            matrix[y][x]=t[y];}
    }
}

void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE],int width,int height){
    complex s[MAX_SIZE],t[MAX_SIZE];
    for(int x=0;x<width;x++){
        for(int y=0;y<height;y++){
            s[y]=matrix[y][x];}
        fft_inverse(s,t,height);
        for(int y=0;y<height;y++){
            matrix[y][x]=t[y];}
    }
    for(int y=0;y<height;y++){
        for(int x=0;x<width;x++){
            s[x]=matrix[y][x];}
        fft_inverse(s,t,width);
        for(int x=0;x<width;x++){
            matrix[y][x]=t[x];}
    }
}
