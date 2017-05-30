#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <GLUT/glut.h>

#define eps 1e-8
#define PI 3.1415926535898

using namespace std;

const int PHONG_N=10;
const int maxdepth=8;

inline bool bige(double a,double b){if (a-b>-eps)return true;return false;}
inline bool Z(double a){if (a>-eps&&a<eps)return true;return false;}

inline double sqr(double n){return n*n;}

int max(int a,int b){if (a<b)return b;else return a;}
double min(double a,double b){if (a<b)return a;else return b;}

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
    
    friend inline vec unit(const vec &a){return vec(a.x/a.len,a.y/a.len,a.z/a.len);}
    
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

const vec ZERO(0,0,0,0,0);

typedef vec color;

const color Ienvironment(10.0,10.0,10.0);
//const color Ienvironment(0,0,0);

void printvec(const vec &a,const string &s){cout<<s<<a.x<<" "<<a.y<<" "<<a.z<<endl;}

struct Ray
{
    vec st,dir;
    Ray(){dir.set(0,0,0);}
    Ray(const vec&St,const vec&Dir):st(St),dir(unit(Dir)){}
    Ray(const Ray&r):st(r.st),dir(r.dir){}
    void set(const vec&St,const vec&Dir){st=St;dir=Dir;}
    double getT(const vec &a) const {return (a-st).len;}
};

/*
 *  refraction:简单的折射
 *  ray:
 */
Ray refraction(const vec &ray,const vec &dots,const vec &N,const double n1,const double n2)
{
    vec L=neg(ray);
    double cosL=dot(N,L);
    double sinL=1.0-cosL*cosL;
    double n=n1/n2;
    double ns=n*n;
    if (sinL*ns>1.0)return Ray(ZERO,ZERO);
    vec T=(n*cosL-sqrt(1-ns*sinL))*N-n*L;
    return Ray(dots,unit(T));
}

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
    virtual Ray reflection(const Ray &ray,const vec &dot){return Ray();}
    //Ray refraction(const Ray &ray,const vec &dot){return Ray();}
    virtual vec getN(const vec &dots){return vec();}
    virtual color insideRayTracing(const Ray&ray,int dep){return color();}
    virtual bool intersect(const Ray &ray,vec &dots,int num){return false;}
};

color RayTracing(const Ray ray,int dep);

struct circle:public shape
{
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
        color iR=insideRayTracing(insideReflection(ray,dots),dep+1);
        color iT=RayTracing(refraction(ray.dir, dots, neg(N),n,1.0), dep+1);
        return mul(iR,kS)+mul(iT,kT);
    }
};

struct plane:public shape
{
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
        if (par==0)return false;
        double t=(dis-dot(ray.st,n))/par;
        if (t<eps)return false;
        dots=ray.st+t*ray.dir;
        return true;
    }
    Ray reflection(const Ray &ray,const vec &dots)
    {
        return Ray(dots,ray.dir+(2.0*(-dot(n,ray.dir)))*n);
    }
};

struct rectangle
{
    vec a,b,c,d,n;
    double dis;
    rectangle(){}
    rectangle(const vec& _a,const vec&_b,const vec&_c,const vec&_d):a(_a),b(_b),c(_c),d(_d)
    {
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
    
    bool inIt(const vec &dots)
    {
        double ra=0.0;
        vec at=unit(a-dots);
        vec bt=unit(b-dots);
        vec ct=unit(c-dots);
        vec dt=unit(d-dots);
        if (neg(at)==bt)return true;
        if (neg(bt)==ct)return true;
        if (neg(ct)==dt)return false;
        if (neg(dt)==at)return false;
        double ab=dot(at,bt);
        double bc=dot(bt,ct);
        double cd=dot(ct,dt);
        double da=dot(dt,at);
        ra+=acos(ab);
        ra+=acos(bc);
        ra+=acos(cd);
        ra+=acos(da);
        if (Z(ra-2.0*PI))
            return true;
        else
            return false;
    }
    
    bool intersect(const Ray &ray,vec &dots)
    {
        double par=dot(ray.dir,n);
        if (par==0)return false;
        double t=(dis-dot(ray.st,n))/par;
        if (t<eps)return false;
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
    }
    
    void print()
    {
        for (int i=0;i<6;++i)
        {
            printvec(r[i].a,"a ");
            printvec(r[i].b,"b ");
            printvec(r[i].c,"c ");
            printvec(r[i].d,"d ");
            putchar('\n');
        }
    }
    
    vec getN(const vec &dots)
    {
        for (int i=0;i<6;++i)
            if (r[i].inIt(dots))
                return r[i].n;
        return ZERO;
    }
    
    vec getN(int i)
    {
        return r[i].n;
    }
    
    bool intersect(const Ray &ray,vec &dots,int num)
    {
        double t[2];
        int tot=0;
        vec dott=ZERO;
        for (int i=0;i<6;++i)
            if (r[i].intersect(ray, dott))
            {
                t[tot]=ray.getT(dott);
                tot++;
            }
        if (tot==0)return false;
        if (num==1)
        {
            if (tot==1)
            {
                dots=ray.st+t[0]*ray.dir;
                return true;
            }
            else
            {
                if (t[0]+eps<t[1])
                {
                    dots=ray.st+t[0]*ray.dir;
                    return true;
                }
                else
                {
                    dots=ray.st+t[1]*ray.dir;
                    return true;
                }
            }
        }
        if (num==2)
        {
            if (tot==1)
            {
                dots=ray.st+t[0]*ray.dir;
                return true;
            }
            else
            {
                if (t[0]+eps<t[1])
                {
                    dots=ray.st+t[1]*ray.dir;
                    return true;
                }
                else
                {
                    dots=ray.st+t[0]*ray.dir;
                    return true;
                }
            }
        }
        return false;
    }
    
    Ray reflection(const Ray &ray,const vec &dots)
    {
        for (int i=0;i<6;++i)
            if (r[i].inIt(dots))
                return Ray(dots,ray.dir+(-2.0*dot(r[i].n,ray.dir))*r[i].n);
        return Ray(ZERO,ZERO);
    }
    
    color insideRayTracing(const Ray &ray,int dep)
    {
        if (dep>maxdepth)return ZERO;
        if (ray.dir==ZERO)return ZERO;
        vec dots,dott,N;
        for (int i=0;i<6;++i)
            if (r[i].intersect(ray, dott))
            {
                double t=ray.getT(dott);
                if (t<eps)continue;
                    else
                    {
                        dots=dott;
                        N=r[i].n;
                    }
            }
        color iR=insideRayTracing(reflection(ray, dots), dep+1);
        color iT=RayTracing(refraction(ray.dir, dots, neg(N), n, 1.0),dep+1);
        return mul(iR,kS)+mul(iT,kT);
    }
};

struct camera
{
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
    /*Ray generateRay(double x,double y)
    {
        vec r=((x-0.50)*scale)*right;
        vec u=((y-0.50)*scale)*up;
        return Ray(eye,unit(forward+r+u));
    }*/
    Ray generateRay(double x,double y)const
    {
        vec r=((x-0.50)*scale)*right;
        vec u=((y-0.50)*scale)*up;
        return Ray(eye,unit(forward+r+u));
    }
}mySpot;

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

const double c1=0.5,c2=0.04,c3=0.0001;

double reduction(double d)
{
    return (min(1/(c1+c2*d+c3*d*d),0.95));
}

bool RayDirect(const Ray ray,const vec &dots)
{
    vec dott=ZERO;
    double st=ray.getT(dots);
    for (int i=0;i<shapetot;++i)
            if (shapeLst[i]->intersect(ray,dott,2))
            {
                double tt=ray.getT(dott);
                if (tt>eps&&tt+eps<=st)return true;
            }
    return false;
}

color RayTracing(const Ray ray,int dep)
{
    if (ray.dir==ZERO)return ZERO;
    bool find=false;
    vec dots;
    vec dott;
    double maxt=0.0;
    int maxi=0;
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
                if (t+eps<maxt)
                {
                    maxt=t;
                    maxi=i;
                    dots=dott;
                }
            }
        }
    if (find)
    {
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
                if (cosSeta>eps)
                    iD=cosSeta*mul(lightLst[i]->color,shapeLst[maxi]->kD);
                
                vec H=unit(L+V);
                vec iS=ZERO;
                double cosSeta2=dot(N,H);
                if (cosSeta2>eps)
                    iS=pow(cosSeta2,PHONG_N)*mul(lightLst[i]->color,shapeLst[maxi]->kS);
                iL=iL+reduction(d)*(iS+iD);
            }
        if (dep<maxdepth)
        {
            color iR=RayTracing(shapeLst[maxi]->reflection(ray, dots), dep+1);
            color iT;
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

const int Len=2048;
const double LEN=Len;
const double units=1/LEN;
vec ans[Len][Len];
int maxx,maxy,maxz;

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
            glColor3d((double)ans[i][j].x/(double)255.0,
                      (double)ans[i][j].y/(double)255.0,
                      (double)ans[i][j].z/(double)255.0);
            glVertex2d((double)i*2.0/LEN-1, (double)j*2.0/LEN-1);
        }
    glEnd();
    
    glFlush();
}

void renderDepth(const camera &cam)
{
    for (int i=0;i<Len;++i)
    {
        double si=1-(double)(i)/LEN;
        for (int j=0;j<Len;++j)
        {
            double sj=(double)(j)/LEN;
            Ray ray=cam.generateRay(si,sj);
            vec dots;
            dots.set(0,0,0);
            shapeLst[0]->intersect(ray, dots,1);
            if (dots==ZERO)
            {
                ans[i][j]=ZERO;
            }
            else
            {
                double depth=(dots-cam.eye).len;
                ans[i][j]=vec(depth,depth,depth);
                maxx=max(maxx,depth);
                maxy=maxz=maxx;
            }
        }
    }
}

int main(int argc, char * argv[])
{
    mySpot.set(vec(0,0,25),unit(vec(0,1.0,1.0)),vec(0,1,0),90.0);
    circle a(vec(-10,10,0), 10.0, vec(0.05,0.45,0.1), vec(0.01,0.95,0.01),
                                vec(0.3,0.63,0.3), vec(0.3,0.63,0.3), 1.33);
    circle b(vec(13,10,13),8.0,vec(0.45,0.01,0.01),vec(0.95,0.05,0.05),
                                 vec(0.85,0.5,0.5),vec(0.7,0.3,0.3),1.33);
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
    plane d(vec(0,0,1),-10.1,vec(0.6,0.6,0.6),vec(0.5,0.5,0.5),vec(0.5,0.5,0.5),vec(0,0,0));
    plane f(vec(0,-1,0),-60.0,vec(0.3,0.3,0.3),vec(0.6,0.0,0.0),vec(0.0,0.0,0.0),vec(0,0,0));
    plane g(vec(1,0,0),-20.5,vec(0.1,0.1,0.1),vec(0.2,0.2,0.2),vec(0.5,0.5,0.5),vec(0,0,0));
    light c(vec(10,5,25),vec(255,255,255));
    light h(vec(-10,15,25),vec(255,255,255));
    light i(vec(0,30,40),vec(255,255,255));
    shapeLst[shapetot++]=&a;
    shapeLst[shapetot++]=&b;
    //shapeLst[shapetot++]=&e;
    //shapeLst[shapetot++]=&j;
    shapeLst[shapetot++]=&k;
    shapeLst[shapetot++]=&l;
    shapeLst[shapetot++]=&d;
    shapeLst[shapetot++]=&f;
    shapeLst[shapetot++]=&g;
    lightLst[lighttot++]=&c;
    lightLst[lighttot++]=&h;
    //lightLst[lighttot++]=&i;
    maxx=maxy=maxz=0;
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
        }
    }
    /*
    for (int i=0;i<Len;++i)
    {
        for (int j=0;j<Len;++j)
            cout<<ans[i][j].x<<"."<<ans[i][j].y<<"."<<ans[i][j].z<<" ";
        cout<<endl;
    }*/
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB);
    glutInitWindowPosition(200, 200);
    glutInitWindowSize(700, 700);
    glutCreateWindow("???");
    glutDisplayFunc(&Display);
    glutMainLoop();/*
     */
    return 0;
}
