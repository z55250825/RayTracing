#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <GLUT/glut.h>

#define eps 1e-8
#define inf -10000
#define PI 3.1415926535898

using namespace std;

//faraway：定义最远处，如果光线延伸了这么远则停止追踪
const double faraway=1e3;
//PHONG_N：用于镜面反射
const int PHONG_N=10;

//光线追踪最深深度
const int maxdepth=8;

//像平面的高和宽（这里是2048*2048)
const int Len=2048;
const double LEN=Len;

//光线衰减函数的常数
const double c1=0.5,c2=0.045,c3=0.0001;

/*
 *  光线衰减函数，跟距离相关
 */
double reduction(double d)
{
    if (d>faraway)return 0.001;
    return (min(1/(c1+c2*d+c3*d*d),0.95));
}

inline bool bige(double a,double b){if (a-b>-eps)return true;return false;}

//如果a的绝对值在eps内则认为是0，精度误差控制
inline bool Z(double a){if (a>-eps&&a<eps)return true;return false;}

inline double sqr(double n){return n*n;}

int max(int a,int b){if (a<b)return b;else return a;}
double min(double a,double b){if (a<b)return a;else return b;}

/*
 *  三维向量
 */
struct vec
{
    double x,y,z,len,slen;
    vec(){x=y=z=len=slen=0.0;}
    vec(double _x,double _y,double _z){x=_x;y=_y;z=_z;slen=x*x+y*y+z*z;len=sqrt(slen);}
    vec(double _x,double _y,double _z,double l,double sl):x(_x),y(_y),z(_z),len(l),slen(sl){}
    vec(const vec& a):x(a.x),y(a.y),z(a.z),len(a.len),slen(a.slen){}
    
    void set(int _x,int _y,int _z)
    {
        x=_x;y=_y;z=_z;slen=x*x+y*y+z*z;len=sqrt(slen);
    }
    
    friend inline vec operator +(const vec&a,const vec&b)
    {
        return vec(a.x+b.x,a.y+b.y,a.z+b.z);
    }
    
    friend inline vec operator -(const vec&a,const vec&b)
    {
        return vec(a.x-b.x,a.y-b.y,a.z-b.z);
    }
    friend inline vec operator *(const double c,const vec&a)
    {
        return vec(c*a.x,c*a.y,c*a.z);
    }
    
    friend inline vec operator *(const vec&a,const double c)
    {
        return vec(c*a.x,c*a.y,c*a.z);
    }
    
    friend inline vec operator /(const vec&a,const double c)
    {
        return vec(a.x/c,a.y/c,a.z/c);
    }
    
    friend inline double dot(const vec &a,const vec &b){return a.x*b.x+a.y*b.y+a.z*b.z;}
    
    friend inline vec mul(const vec &a,const vec &b){return vec(a.x*b.x,a.y*b.y,a.z*b.z);}
    
    friend inline vec det(const vec &a,const vec &b)
    {
        return vec(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
    }
    
    friend inline vec unit(const vec &a)
    {if (Z(a.len))return vec(0,0,0,0,0);else return vec(a.x/a.len,a.y/a.len,a.z/a.len);}
    
    friend inline bool operator ==(const vec &a,const vec &b)
    {
        return (Z(a.x-b.x)&&Z(a.y-b.y)&&Z(a.z-b.z));
    }
    
    friend inline double dist(const vec &a,const vec &b)
    {
        return (sqrt(sqr(a.x-b.x)+sqr(a.y-b.y)+sqr(a.z-b.z)));
    }
    
    friend inline vec neg(const vec &a){return vec(-a.x,-a.y,-a.z);}
    
    inline vec operator =(const vec &a){x=a.x;y=a.y;z=a.z;len=a.len;slen=a.slen;return *this;}
};

//ZERO常量
const vec ZERO(0,0,0,0,0);
const vec INF(inf,inf,inf);

typedef vec color;

//环境光常量
const color Ienvironment(10.0,10.0,10.0);
//const color Ienvironment(0,0,0);

//调试用，打印出某个vec的坐标
void printvec(const vec &a,const string &s){cout<<s<<a.x<<" "<<a.y<<" "<<a.z<<endl;}

/*
 *  Ray:光线类
 *  st: 光线起点
 *  dir:光线方向（单位化）
 */
struct Ray
{
    vec st,dir;
    Ray(){dir.set(0,0,0);}
    Ray(const vec&St,const vec&Dir):st(St),dir(unit(Dir)){}
    Ray(const Ray&r):st(r.st),dir(r.dir){}
    void set(const vec&St,const vec&Dir){st=St;dir=Dir;}
    
    //求出光线上的某点a的参数坐标t(所有点表示为Q=st+t*dir)
    double getT(const vec &a) const {return (a-st).len;}
};

/*
 *  refraction:简单的折射
 *  ray:入射光线
 *  dots:入射点
 *  N:入射处表面的法线（朝入射所在介质）
 *  n1:入射介质的折射率（空气默认为1.0)
 *  n2:折射介质的折射率。
 */
Ray refraction(const vec &ray,const vec &dots,const vec &N,const double n1,const double n2)
{
    vec L=unit(neg(ray));
    double cosL=dot(N,L);
    vec NN=ZERO;
    if (cosL<-eps)NN=unit(neg(N));else NN=unit(N);
    double sinL=1.0-cosL*cosL;
    double n=n1/n2;
    double ns=n*n;
    if (sinL*ns+eps>=1.0)return Ray(ZERO,ZERO);
    vec T=(n*cosL-sqrt(1-ns*sinL))*NN-n*L;
    return Ray(dots,unit(T));
}

/*
 * shape:所有环境中的物品的父类，所有物体都继承该类
 * （方便实现多态）
 */
struct shape
{
    vec kA;//ambient Reflection Coefficient [0,1]
    vec kD;//diffuse Reflection Coefficient [0,1]
    vec kS;//specular Reflection Coefficient [0,1]
    vec kT;//refraction Reflection Coeeficient [0,1]
    double n;
    shape(){}
    shape(const vec&A,const vec&D,const vec&S,const vec &T,double N)
    {
        kA=A;kD=D;kS=S;kT=T;n=N;
    }
    shape(const shape &a):kA(a.kA),kD(a.kD),kS(a.kS),kT(a.kT),n(a.n){}
    
    /*
     *  反射，所有物体都要实现
     *  ray: 入射光线
     *  dot: 入射点
     *  返回值，反射光线
     */
    virtual Ray reflection(const Ray &ray,const vec &dot){return Ray();}
    
    /*
     *  得到某处的法向量，所有物体都要实现
     *  dots:物体上的指定位置
     *  返回值，该点的法向量
     */
    virtual vec getN(const vec &dots){return vec();}
    
    /*
     * 内部光线追踪，所有物体都要实现
     * 内部单独的光线追踪，区别在于没有环境光，加上不需要特地扫描所有物体求交
     * ray：光线
     * dep：追踪深度
     * 返回值：求出来的光线颜色
     */
    virtual color insideRayTracing(const Ray&ray,int dep){return color();}
    
    /*
     * 判定光线是否相交，所有物体都要实现
     * 判断某条光线是否与该物体相交，如果相交，根据num
     * num==1,求出最近交点放在dots里
     * num==2,求出第二近的交点放在dots里，如果只有一个交点，则把那个
     * 交点放在dots里（为了判阴影所需）
     * 返回值：布尔，如果相交则为true，否则为false
     */
    virtual bool intersect(const Ray &ray,vec &dots,int num){return false;}
};

//核心函数
color RayTracing(const Ray ray,int dep);

/*
 * 球类
 */
struct circle:public shape
{
    //center:球心，dim:半径
    vec center;
    double dim;
    circle(){}
    circle(vec Center,double Dim,const vec&A,const vec&D,const vec&S,const vec&T,double N)
    {
        kA=A;kD=D;kS=S;kT=T;n=N;
        center=Center;
        dim=Dim;
    }
    circle(const circle &c)
    {
        kA=c.kA;kD=c.kD;kS=c.kS;kT=c.kT;n=c.n;
        center=c.center;dim=c.dim;
    }
    
    vec getN(const vec &dots){return unit(dots-center);}
    
    bool intersect(const Ray &ray,vec &dots,int num)
    {
        double tmp=dot(ray.st-center,ray.dir);
        double b=2.0*tmp;
        double c=ray.st.slen+center.slen-2*dot(ray.st,center)-dim*dim;
        double delta=b*b-4.0*c;
        if (delta<=eps)return false;
        double d=sqrt(delta)/2.0;
        double t1=-tmp-d,t2=-tmp+d;
        if (num==1)
        {
            if (t1>eps)
            {
                dots=ray.st+t1*ray.dir;
                return true;
            }
            else
                if (t2>eps)
                {
                    dots=ray.st+t2*ray.dir;
                    return true;
                }
                else
                    return false;
        }
        if (num==2)
        {
            if (t2>eps)
            {
                dots=ray.st+t2*ray.dir;
                return true;
            }
            return false;
        }
        return false;
    }
    
    Ray reflection(const Ray &ray,const vec &dots)
    {
        vec n=getN(dots);
        return Ray(dots,ray.dir+(2.0*(-dot(n,ray.dir)))*n);
    }
    
    //内部反射，其实好像推了一下发现和外部是一样的
    Ray insideReflection(const Ray &ray,const vec &dots)
    {
        vec n=neg(getN(dots));
        return Ray(dots,ray.dir+(2.0*(-dot(n,ray.dir)))*n);
    }
    
    color insideRayTracing(const Ray &ray,int dep)
    {
        if (dep>maxdepth)return ZERO;
        if (ray.dir==ZERO)return ZERO;
        vec Dim=center-ray.st;
        double s=2.0*dot(Dim,ray.dir);
        vec dots=ray.st+s*ray.dir;
        vec N=unit(dots-center);
        color iR=ZERO;
        if (dep<maxdepth)insideRayTracing(insideReflection(ray,dots),dep+1);
        color iT=ZERO;
        if (dep<maxdepth)RayTracing(refraction(ray.dir, dots, neg(N),n,1.0), dep+1);
        //printvec(iR,"iR ");
        //printvec(iT,"iT ");
        return mul(iR,kS)+mul(iT,kT);
    }
};

//平面类
struct plane:public shape
{
    //n:法向量，dis:平面上任意一点到原点的向量与
    //法向量n的点积，实际上就是平面到原点的距离或者
    //距离的相反数。
    
    //平面上所有点到原点的向量与n点积都是固定值dis，用这点
    //可以判断出某个点是否在平面上
    vec n;
    double dis;
    plane(){}
    plane(const vec& N,double Dis,const vec &A,const vec &D,const vec &S,const vec &T)
    {
        n=unit(N);dis=Dis;kA=A;kD=D;kS=S;kT=ZERO;
    }
    plane(const plane &a)
    {
        n=a.n;dis=a.dis;kA=a.kA;kD=a.kD;kS=a.kS;kT=a.kT;
    }
    
    vec getN(const vec &dots){return n;}
    
    bool intersect(const Ray &ray,vec &dots,int num)
    {
        double par=dot(ray.dir,n);
        if (Z(par))return false;
        double t=(dis-dot(ray.st,n))/par;
        if (t<eps||t>faraway)return false;
        dots=ray.st+t*ray.dir;
        return true;
    }
    Ray reflection(const Ray &ray,const vec &dots)
    {
        return Ray(dots,ray.dir+(2.0*(-dot(n,ray.dir)))*n);
    }
};

/*
 * 长方形（其实是四边形也可以）
 */
struct rectangle
{
    //a,b,c,d为逆时针或者顺时针的顶点
    //n为法向量
    vec a,b,c,d,n;
    double dis;
    rectangle(){}
    rectangle(const vec& _a,const vec&_b,const vec&_c,const vec&_d)
    {
        a=_a;b=_b;c=_c;d=_d;
        n=unit(det(a-b,c-b));
        dis=dot(a,n);
    }
    rectangle(const rectangle &a):a(a.a),b(a.b),c(a.c),d(a.d),n(a.n),dis(a.dis){}
    void set(const vec &_a,const vec &_b,const vec &_c,const vec &_d)
    {
        a=_a;b=_b;c=_c;d=_d;
        n=unit(det(a-b,c-b));
        dis=dot(a,n);
    }
    
    //d是否在三角形abc中，其中在bc上也算
    inline bool Init(const vec &a,const vec &b,const vec &c,const vec &d)
    {
        vec P=d-a;
        vec P1=b-a;
        vec P2=c-a;
        double A=P1.slen;
        double B=dot(P1,P2);
        double C=P2.slen;
        double D1=dot(P,P1);
        double D2=dot(P,P2);
        double s=A*C-B*B;
        double u=(D1*C-D2*B)/s;
        double v=(D2*A-D1*B)/s;
        if (u>0.0&&u+eps<1.0&&v>0.0&&v+eps<1.0&&u+v+eps<1.0)return true;
        return false;
    }
    //判断某个点dots是否在该四边形内部
    //如果是返回true,否则返回false
    //点在边ab,bc上也认为是在内部
    //但点在ad,cd上则不认为在内部
    //注意：该函数已经确定了该点在四边形平面上
    bool inIt(const vec &dots)
    {
        vec at=unit(a-dots);
        vec bt=unit(b-dots);
        vec ct=unit(c-dots);
        vec dt=unit(d-dots);
        if (neg(at)==bt)return true;
        if (neg(bt)==ct)return true;
        if (neg(ct)==dt)return false;
        if (neg(dt)==at)return false;
        if (Init(a,b,d,dots))return true;
        if (Init(c,d,b,dots))return true;
        return false;
    }
    
    //判断某条光线ray是否和该四边形相交
    //如果是，将交点存在dots中，返回true
    //否则返回false
    bool intersect(const Ray &ray,vec &dots)
    {
        double par=dot(ray.dir,n);
        //printvec(ray.st,"intersect st ");
        //printvec(ray.dir,"intersect dir ");
        //cout<<"par "<<par<<endl;
        if (Z(par))return false;
        double t=(dis-dot(ray.st,n))/par;
        //cout<<"t "<<t<<endl;
        if (t<eps||t>faraway)return false;
        vec dott=ray.st+t*ray.dir;
        if (inIt(dott))
        {
            dots=dott;
            return true;
        }
        else
            return false;
    }
};

struct cuboid:public shape
{
    rectangle r[6];
    cuboid(){}
    cuboid(const rectangle &a,const vec &N,double h,double Ns,
           const vec&A,const vec&D,const vec&S,const vec&T)
    {
        kA=A;kD=D;kS=S;kT=T;
        vec Nh=h*unit(N);
        n=Ns;
        vec s=det(a.c-a.b,a.a-a.b);
        if (dot(s,N)>0)r[0].set(a.a,a.b,a.c,a.d);
                else   r[0].set(a.a,a.d,a.c,a.b);
        vec an=r[0].a;vec bn=r[0].b;
        vec cn=r[0].c;vec dn=r[0].d;
        vec au=an+Nh;vec bu=bn+Nh;
        vec cu=cn+Nh;vec du=dn+Nh;
        r[1]=rectangle(au,du,cu,bu);
        r[2]=rectangle(an,au,bu,bn);
        r[3]=rectangle(bn,bu,cu,cn);
        r[4]=rectangle(dn,cn,cu,du);
        r[5]=rectangle(an,dn,du,au);
        //print();
    }
    
    void print()
    {
        for (int i=0;i<6;++i)
        {
            printvec(r[i].a,"a ");
            printvec(r[i].b,"b ");
            printvec(r[i].c,"c ");
            printvec(r[i].d,"d ");
            printvec(r[i].n,"n ");
            putchar('\n');
        }
    }
    
    vec getN(const vec &dots)
    {
        for (int i=0;i<6;++i)
        {
            double dis=dot(r[i].n,dots);
            if (Z(dis-r[i].dis))
                if (r[i].inIt(dots))
                    return r[i].n;
        }
        return ZERO;
    }
    
    vec getN(int i)
    {
        return r[i].n;
    }
    
    bool intersect(const Ray &ray,vec &dots,int num)
    {
        double t[7];
        int tots=0;
        vec dott=ZERO;
        for (int i=0;i<6;++i)
            if (r[i].intersect(ray, dott))
            {
                t[tots]=ray.getT(dott);
                tots++;
            }
        if (tots==0)return false;
        sort(t,t+tots);
        int ind=0;
        while (ind<tots&&t[ind]<eps)ind++;
        if (ind==tots)return false;
        if (num==1)
        {
            dots=ray.st+t[ind]*ray.dir;
            return true;
        }
        if (num==2)
        {
            if (ind==tots-1)
                dots=ray.st+t[ind]*ray.dir;
            else
                dots=ray.st+t[ind+1]*ray.dir;
            return true;
        }
        return false;
    }
    
    Ray insideReflection(const Ray &ray,const vec &dots,const vec &N)
    {
        return Ray(dots,ray.dir+(-2.0*dot(N,ray.dir)*N));
    }
    
    Ray reflection(const Ray &ray,const vec &dots)
    {
        for (int i=0;i<6;++i)
        {
            double dis=dot(dots,r[i].n);
            if (Z(r[i].dis-dis))
                if (r[i].inIt(dots))
                    return Ray(dots,ray.dir+(-2.0*dot(r[i].n,ray.dir))*r[i].n);
        }
        return Ray(ZERO,ZERO);
    }
    
    color insideRayTracing(const Ray &ray,int dep)
    {
        if (dep>maxdepth)return ZERO;
        if (ray.dir==ZERO)return ZERO;
        //cout<<dep<<endl;
        //printvec(ray.st,"1 st ");
        //printvec(ray.dir,"1 dir ");
        vec dots=ZERO,dott,N;
        for (int i=0;i<6;++i)
            if (r[i].intersect(ray, dott))
            {
                double t=ray.getT(dott);
                //cout<<"t "<<t<<endl;
                if (t<eps)continue;
                    else
                    {
                        dots=dott;
                        N=neg(r[i].n);
                    }
            }
        if (dots==ZERO)
        {
            //getchar();
            return RayTracing(ray, dep);
        }
        color iR=ZERO;
        if (dep<maxdepth)iR=insideRayTracing(insideReflection(ray, dots, N), dep+1);
        color iT=ZERO;
        if (dep<maxdepth)
            iT=RayTracing(refraction(ray.dir, dots, N, n, 1.0),dep+1);
        //printvec(ray.st,"st ");
        //printvec(ray.dir,"dir ");
        //printvec(iR,"iR ");
        //printvec(iT,"iT ");
        return mul(iR,kS)+mul(iT,kT);
    }
};

//摄影机类
struct camera
{
    /*
     *  eye:摄影机位置
     *  forward:摄影机正对方向的向量（单位化，下同）
     *  right：摄影机的正右端方向的向量
     *  up:摄影机的正上方方向的向量
     *  stdUp:初始构造摄影机类时候初步判断摄影机的上下的向量
     *  （注：right,up都是构造时利用forward和stdUp算出来的）
     *  fov：视角角度
     *  scale：面前1单位的屏幕的视角宽或者高
     */
    vec eye,forward,right,up,stdUp;
    double fov,scale;
    camera(){}
    camera(const vec&e,const vec&f,const vec& StdUp,double fo):eye(e),forward(unit(f)),fov(fo),stdUp(StdUp)
    {
        right=unit(det(forward,stdUp));
        up=unit(det(right,forward));
        scale=2*tan(fov*PI/360.00);
    }
    camera(const camera &c):eye(c.eye),forward(c.forward),stdUp(c.stdUp),fov(c.fov),
    right(c.right),up(c.up),scale(c.scale){}
    
    void set(const vec&e,const vec&f,const vec& StdUp,double fo)
    {
        eye=e;forward=unit(f);fov=fo;stdUp=StdUp;
        right=unit(det(forward,stdUp));
        up=unit(det(right,forward));
        scale=2*tan(fov*PI/360.00);
    }
    
    /*
     * 对于像平面（默认(0,0)~(1,1))(x,y)处的像素所产生的光线
     */
    Ray generateRay(double x,double y)const
    {
        vec r=((x-0.50)*scale)*right;
        vec u=((y-0.50)*scale)*up;
        return Ray(eye,unit(forward+r+u));
    }
}mySpot;

//光源类，pos表示位置，color表示颜色
struct light
{
    vec pos,color;
    light(){}
    light(const vec &p,const vec &c){pos=p;color=c;}
    light(const light &l):pos(l.pos),color(l.color){}
};

light *lightLst[10];
shape *shapeLst[10];
int shapetot=0,lighttot=0;

/*
 * 判断某条光线和某个点dots之间是否有其他物体遮挡
 * ray : 光线
 * dots: 点
 * 返回值：判断ray和dots间是否有物体遮挡住，如果是，返回true，
 * 否则返回false
 */
bool RayDirect(const Ray ray,const vec &dots)
{
    vec dott=ZERO;
    double st=ray.getT(dots);
    //枚举所有物体，求交点，判断是否在中间
    for (int i=0;i<shapetot;++i)
            if (shapeLst[i]->intersect(ray,dott,2))
            {
                double tt=ray.getT(dott);
                if (tt>eps&&tt+eps<=st)return true;
            }
    return false;
}

//光线追踪
color RayTracing(const Ray ray,int dep)
{
    if (ray.dir==ZERO)return ZERO;
    if (dep>maxdepth)return ZERO;
    bool find=false;
    vec dots;
    vec dott;
    double maxt=0.0;
    int maxi=0;
    //求最近交点
    for (int i=0;i<shapetot;++i)
        if (shapeLst[i]->intersect(ray,dott,1))
        {
            if (find==false)
            {
                find=true;
                dots=dott;
                maxt=ray.getT(dots);
                maxi=i;
            }
            else
            {
                double t=ray.getT(dott);
                if (t+eps<faraway&&t+eps<maxt)
                {
                    maxt=t;
                    maxi=i;
                    dots=dott;
                }
            }
        }
    if (find)
    {
        //求环境光影响iL
        vec iL=mul(Ienvironment,shapeLst[maxi]->kA);
        vec N=shapeLst[maxi]->getN(dots);
        vec V=unit(ray.st-dots);
        bool finds=false;
        for (int i=0;i<lighttot;++i)
            if (!RayDirect(Ray(dots,unit(lightLst[i]->pos-dots)), lightLst[i]->pos))
            {
                finds=true;
                vec L=unit(lightLst[i]->pos-dots);
                double d=dist(lightLst[i]->pos,dots);
                double cosSeta=dot(N,L);
                vec iD=ZERO;
                //求漫反射影响iD
                if (cosSeta>eps)
                    iD=cosSeta*mul(lightLst[i]->color,shapeLst[maxi]->kD);
                
                vec H=unit(L+V);
                vec iS=ZERO;
                double cosSeta2=dot(N,H);
                //求镜面反射影响iS
                if (cosSeta2>eps)
                    iS=pow(cosSeta2,PHONG_N)*mul(lightLst[i]->color,shapeLst[maxi]->kS);
                iL=iL+reduction(d)*(iS+iD);
            }
        //求完局部光照影响，继续光线追踪
        if (dep<maxdepth)
        {
            color iR=RayTracing(shapeLst[maxi]->reflection(ray, dots), dep+1);
            color iT=ZERO;
            if (shapeLst[maxi]->kT==ZERO)iT.set(0,0,0);
                else iT=shapeLst[maxi]->insideRayTracing(
                                                         refraction(ray.dir,dots,shapeLst[maxi]->getN(dots),1.0,shapeLst[maxi]->n),
                                                         dep+1);
            return iL+mul(shapeLst[maxi]->kS,iR)+mul(shapeLst[maxi]->kT,iT);
        }
        else
            return iL;
    }
    else return ZERO;
    return ZERO;
}

vec ans[Len][Len];

void readAns()
{
    freopen("output.txt","r",stdin);
    for (int i=0;i<Len;++i)
    {
        for (int j=0;j<Len;++j)
            scanf("%lf %lf %lf ",&ans[i][j].x,&ans[i][j].y,&ans[i][j].z);
    }
    fclose(stdin);
}

double maxx,maxy,maxz;
//画图
void Display()
{
    //设置背景颜色为白色
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    //初始化画布
    glClear(GL_COLOR_BUFFER_BIT);
    
    glPointSize(0.5f);
    
    glBegin(GL_POINTS);
    for (int i=0;i<Len;++i)
        for (int j=0;j<Len;++j)
        {
            glColor3d((double)ans[i][j].x/(double)maxx,
                      (double)ans[i][j].y/(double)maxy,
                      (double)ans[i][j].z/(double)maxz);
            glVertex2d((double)i*2.0/LEN-1, (double)j*2.0/LEN-1);
        }
    glEnd();
    
    glFlush();
}

int main(int argc, char * argv[])
{
    //(-1,1,0) may be problem?
    //摄影机视角 camera(位置vec，正前方vec，粗略的确定正上方的vec，视角大小（角度多少度）)
    //mySpot.set(vec(80,25,5),unit(vec(-1.0,0.0,0.0)),vec(0,0,1),90.0);
    
    mySpot.set(vec(6,-25,0),unit(vec(0,1,0)),vec(0,0,1),90.0);
    //mySpot.set(vec(5,-35,10),unit(vec(0,1.0,0.0)),vec(0,0,1),90.0);
    
    //球a，circle(球心vec,半径,kA,kD,kS,kT,折射率）
    circle a(vec(-10,10,0), 10.0, vec(0.05,0.45,0.1), vec(0.01,0.95,0.01),
                                vec(0.3,0.63,0.3), vec(0.3,0.63,0.3), 1.33);
    circle b(vec(13,10,13),8.0,vec(0.45,0.01,0.01),vec(0.95,0.05,0.05),
                                 vec(0.85,0.5,0.5),vec(0.7,0.3,0.3),1.33);
    
    //长方体k cuboid(底面rectangle(四个顶点)，高度，折射率，kA,kD,kS,kT)
    cuboid k(rectangle(vec(5,10,-10),vec(15,0,-10),vec(25,10,-10),vec(15,20,-10)),
             vec(0,0,1),
             15.0,1.33,
             vec(0.65,0.65,0.25),vec(0.65,0.65,0.05),vec(0.85,0.85,0.4),vec(0.85,0.85,0.3));
    circle e(vec(0,38,10),20.0,vec(0.05,0.05,0.45),vec(0.05,0.05,1.0),vec(0.3,0.3,0.63),vec(0.3,0.3,0.73),1.33);
    circle j(vec(0,-3,-5.0),5.0,vec(0.0,0.0,0.0),vec(0.00,0.00,0.00),vec(0.5,0.5,0.5),
             vec(0.5,0.5,0.5),1.33);
    
    cuboid l(rectangle(vec(-10,38,-10),vec(10,18,-10),vec(30,38,-10),vec(10,58,-10)),
             vec(0,0,1),
             40.0,1.2,
             vec(0.25,0.75,0.75),vec(0.2,0.75,0.75),vec(0.2,0.5,0.5),vec(0.07,0.7,0.7));
    //放在xy平面以下10.5单位的平面（底部墙壁）
    plane d(vec(0,0,1),-10.5,vec(0.6,0.6,0.6),vec(0.5,0.5,0.5),vec(0.5,0.5,0.5),vec(0,0,0));
    //放在xz平面向后60单位的平面（红色）（前部墙壁）
    plane f(vec(0,-1,0),-60.0,vec(0.3,0.3,0.3),vec(0.6,0.0,0.0),vec(0.0,0.0,0.0),vec(0,0,0));
    //放在yz平面向左20.5cm的平面（左墙壁）
    plane g(vec(1,0,0),-20.5,vec(0.1,0.1,0.1),vec(0.2,0.2,0.2),vec(0.5,0.5,0.5),vec(0,0,0));
    
    //三个光源 light(位置vec,光的颜色vec)
    light c(vec(10,5,25),vec(255,255,255));
    light h(vec(-10,15,25),vec(255,255,255));
    light i(vec(0,30,40),vec(255,255,255));
    light m(vec(28,10,45),vec(255,255,255));
    //将需要绘制的物体加入到shapeLst指针数组即可
    shapeLst[shapetot++]=&a;
    shapeLst[shapetot++]=&b;
    //shapeLst[shapetot++]=&e;
    shapeLst[shapetot++]=&j;
    shapeLst[shapetot++]=&k;
    shapeLst[shapetot++]=&l;
    shapeLst[shapetot++]=&d;
    shapeLst[shapetot++]=&f;
    shapeLst[shapetot++]=&g;
    
    //将需要的光源加入到lightLst指针数组即可
    lightLst[lighttot++]=&c;
    lightLst[lighttot++]=&h;
    lightLst[lighttot++]=&i;
    lightLst[lighttot++]=&m;
    
    //枚举像素，由camera产生追踪光线
    //maxx=maxy=maxz=0.0;
    maxx=maxy=maxz=255.0;
    //color cs=shapeLst[2]->insideRayTracing(Ray(vec(16.5394,1.53944,-9.63859),vec(-0.853694,0.502995,-0.134916)),7);
    //printvec(cs," ");
    //return 0;
    for (int ii=0;ii<Len;++ii)
    {
        double si=(double)(ii)/LEN;
        for (int jj=0;jj<Len;++jj)
        {
            double sj=(double)(jj)/LEN;
            ans[ii][jj]=RayTracing(
                                 mySpot.generateRay(si,sj),
                                 1
                                );
            ans[ii][jj].x=min(ans[ii][jj].x,255.0);
            ans[ii][jj].y=min(ans[ii][jj].y,255.0);
            ans[ii][jj].z=min(ans[ii][jj].z,255.0);
            /*maxx=max(ans[ii][jj].x,maxx);
            maxy=max(ans[ii][jj].y,maxy);
            maxz=max(ans[ii][jj].z,maxz);*/
        }
    }
    
    //将结果打印到output.txt里面，方便以后直接输出
    //printAns();*/
    //readAns();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB);
    glutInitWindowPosition(200, 200);
    glutInitWindowSize(700, 700);
    glutCreateWindow("RayTracing");
    glutDisplayFunc(&Display);
    glutMainLoop();/*
     */
    return 0;
}
