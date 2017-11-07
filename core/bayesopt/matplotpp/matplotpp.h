/****************************************************************************
Copyright (c) 2011 Yuichi Katori All Rights Reserved
License: Gnu Public license (GPL) v3
Author: Yuichi Katori (yuichi.katori@gmail.com)
Project:MATPLOT++ (MATLAB-like plotting tool in C++).
Version:0.3.13

Modified: Ruben Martinez-Cantin (2013)
    - Fixed std namespace issue
    - Fixed bugs
****************************************************************************/

//FreeGLUT does not need this. But Nate's GLUT for Windows does.
#if defined(_WIN32) || defined(_WIN64)
#    include <windows.h>
#endif

//Using default GLUT in Mac OS
#ifdef __APPLE__
#    include <GLUT/glut.h>
#else
#    include <GL/glut.h>
#endif

#include <vector>
#include <deque>
#include <string>
#include <valarray>
#include <iostream>
#include <cmath>
#include <ctime>
#include "gl2ps.h"

#define PI 3.14159265358979323846264

typedef std::vector<double> dvec;
typedef std::vector< std::vector<double> > dmat;
typedef std::vector< std::vector<float> > tcvec;
typedef std::vector< std::vector< std::vector<float> > > tcmat;

inline dvec linspace(double min,double max,int n){
    dvec a;
    if(n<1){n=1;}
    a.resize(n);
    for(int i=0;i<n;++i){a[i]=min+(max-min)*i/(n-1);}
    return a;
};

inline std::valarray<double> valinspace(double min,double max,int n){
    std::valarray<double> a; 
    a.resize(n);
    for(int i=0;i<n;++i){a[i]=min+(max-min)*i/(n-1);}
    return a;
};

inline int glutCreateWindow(int left,int top,int width,int height,char c[]){
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
    glutInitWindowPosition( left, top );
    glutInitWindowSize( width, height );   
    return glutCreateWindow(c);
};
inline int glutCreateWindow(int left,int top,int width,int height){
    char c[]="adf";
    return glutCreateWindow(left,top,width,height,
			    c);
};

class Figure{///
 public:
    int id;
    //int Status;// 0:minimized, 1:default position, 2:maximized 
    int Position[4];//left top width height
    int Visible;
    std::vector<int> Children;

    void add_child(int i);
    Figure(int id_){
	id=id_;
	//Status=1;
	//Position[0]=100;
	//Position[1]=100;
	//Position[2]=800;
	//Position[3]=800;
	Visible=1;
    };
};


class Layer{///
 public:
    int id;
    int Visible;
    std::string layername;
    std::vector<int> Children;
    Layer(int id_);
    void add_child(int i);
};

class Axes{///
 private:
    
 protected:
    
 public:
    int id;

    float cta,phi;  // controled by mouse
    float cta0,phi0;// default value or specified by command line
    // Mouse 
    double XMouse,YMouse;
    int Mouse;

    double xmin,xmax,ymin,ymax,zmin,zmax;
    int num_child;

    void reset();
    void config();
    int ID();
    int selected();
    void selected(int i);
    void add_child(int i);    
    dvec make_tick(double min,double max);
    
    int View;// 0:2D, 1:3D

    tcvec ColorMap;// for colorbar

    // Matlab variables //
    // styles
    int Box;//0:Off, 1:On
    std::string GridLineStyle;
    float LineWidth;
    std::string TickDir;// {in} | out
    //string TickDirMode;
    //TickLength
    int Visible;//0:Off, 1:On
    int XGrid,YGrid,ZGrid;// {0:Off}, 1:On

    // General Information 
    int Parent;
    std::vector<int> Children;
    int Selected;
    float Position[4];//left bottom width height
    float Viewport3d[4];//left bottom width height

    //Scale
    std::string XAxisLocation;
    std::string YAxisLocation;

    //string XDir,YDir,ZDir;

    double XLim[2],YLim[2],ZLim[2];//plot range
    int XLimMode,YLimMode,ZLimMode;//0:Auto 1:Manual

    int XScale,YScale,ZScale;// linear | log

    dvec XTick,YTick,ZTick;
    std::string XTickMode,YTickMode,ZTickMode;
    int TickLabel;// 0:Off, {1:On}
    //View
    float CameraPosition[3];
    float CameraTarget[3];
    float CameraUpVector[3];

    // Label
    std::string Title;
    std::string XLabel,YLabel,ZLabel;

    double CLim[2];

    Axes(int id_){
	id=id_;
	Selected=0;
	Position[0]=0.13;
	Position[1]=0.11;
	Position[2]=0.775;
	Position[3]=0.815;
    
	Viewport3d[0]=0.0;
	Viewport3d[1]=0.0;
	Viewport3d[2]=1.0;
	Viewport3d[3]=1.0;
	
	Mouse=1;
	View=0;
	Visible=1;
	Box=1;
	Children.clear();
	
	cta0=30; 
	phi0=30;
	cta=cta0; 
	phi=cta0;
	
	CameraPosition[0]=1; CameraPosition[1]=1; CameraPosition[2]=1;
	CameraTarget[0]=0.;  CameraTarget[1]=0; CameraTarget[2]=0;
	CameraUpVector[0]=0; CameraUpVector[1]=0; CameraUpVector[2]=1;
	
	LineWidth=1;
	
	GridLineStyle=":";
	XGrid=0;
	YGrid=0;
	ZGrid=0;
	
	XLim[0]=0;    XLim[1]=10;
	YLim[0]=0;    YLim[1]=10;
	ZLim[0]=0;    ZLim[1]=10;
	
	XLimMode=0; YLimMode=0; ZLimMode=0;
	XAxisLocation="bottom";//top | bottom
	YAxisLocation="left";// left | right
	XScale=0;// {0:linear} | 1:log
	YScale=0;
	ZScale=0;
    
	TickLabel=1;
	TickDir="in";
	
	xmin= 1e99;    xmax=-1e99;
	ymin= 1e99;    ymax=-1e99;
	zmin= 1e99;    zmax=-1e99;
	
	CLim[0]=0;CLim[1]=0;

	num_child=0;
	//MakeTick();
	//Parent=i_figure;
    };
};


class Line{///
 public:
    int id;
    int Errorbar;

    void reset();
    void color(float r,float g,float b);

    // Matlab oriented variables //

    dvec XData,YData,ZData;
    dvec YPData,YMData;
    //dmat XData,YData,ZData;
    //dmat EData,UData,LData;
    //dmat VData,WData;

    std::string Color;
    std::string LineStyle;// {-} | - - | : | -. | none
    float  LineWidth;
    std::string Marker;// {none}
    float  MarkerSize;
    std::string MarkerEdgeColor;
    std::string MarkerFaceColor;

    int Clipping;
    //string EraseMode;
    int SelectionHighlight;
    int Visible;

    // General Information 
    int Parent;
    int Selected;
    
    Line(int id_){	
	id=id_;

	Color="b";
	LineStyle="-";
	LineWidth=0.5;

	Marker="none";
	MarkerSize=6;
	MarkerEdgeColor="k";
	MarkerFaceColor="w";

	Errorbar=0;
	//PlotStyle=0;
    }   
};
class Surface{///
 public:
    int type;
    int id;
    std::string ColorMap;

    //dvec XData,YData;
    dmat XData,YData,ZData,CDataIndex;
    dvec V,UserData;
    tcmat CData;
    
    std::string FaceColor;//ColorSpec    | none | {flat} 
    std::string EdgeColor;//ColorSpec{k} | none | flat
    
    std::string LineStyle;// {-} | - - | : | -. | none
    float  LineWidth;
    std::string Marker;// {none}
    float  MarkerSize;
    std::string MarkerEdgeColor;
    std::string MarkerFaceColor;

    int Parent;

    int NContour;

    Surface(int id_){
	id=id_;

	ColorMap="Gray";
	//Shading="faceted";
	FaceColor="flat"; 
	EdgeColor="b";
	LineStyle="-";
	LineWidth=0.5;
	NContour=10;
	V.clear();
	
    }
    void get(){
	std::cout <<"FaceColor: "<< FaceColor <<std::endl;
	std::cout <<"EdgeColor: "<< EdgeColor <<std::endl;
	std::cout <<"LineStyle: "<< LineStyle <<std::endl;
	std::cout <<"LineWidth: "<< LineWidth <<std::endl;
    }
};
//Note: 
// [m,n] = size(Z), length(x) = n length(y) = m, (x(j),y(i),Z(i,j))

class Patch{///
 public:
    int id;
    int type;
    std::vector< std::vector<int> > Faces;
    dmat Vertices;
    dmat FaceVertexCData;
    dmat XData,YData,ZData;

    //dvec ICVec;
    //dmat ICMat;    
    //tcmat CData;
    tcvec CData;

    std::string EdgeColor,FaceColor;//{ColorSpec}|none|flat|interp 

    std::string LineStyle; // {-} | - - | : | -. | none
    float  LineWidth;

    Patch(int id_){	
	id=id_;
	type=0;
	LineWidth=1;
	FaceColor="r"; 
	EdgeColor="k";
	LineStyle="-";
    }
};
//Note: XData[iv][if]

class Text{///  
 public:
    int id;
    std::string String;
    float Position[3];
    int Parent;
    int type;//0:axis 1:figure
    Text(int id_);
};



const int tFigure=1;
const int tAxes=2;
const int tLine=3;
const int tSurface=4;
const int tText=5;
const int tLayer=6;
const int tPatch=7;

/// contour
struct ContourPoint{
    double x,y;
    int xj,yi;
    int xy;
    int done;
};
dmat contourc(dvec x, dvec y, dmat Z, dvec v);


class MatPlot{///
 private:
    int is_debug1;
    int is_debug2;
    
    tcvec cmap;//TODO move to class Figure
    
    int mode;//0:initialization 1:configuration
    int init_level;// initialization level of objects
    int hObj;// handle number of current object

    int time_layer_clicked,time_layer_clicked_last;

    // Events //
    int window_w,window_h;
    float xButtonDown,yButtonDown;// last clicked mouse position
    float ctaButtonDown,phiButtonDown;
    int xPassive,yPassive;

    // pointer to current objects //
    Figure *cf;
    Layer *cfr;
    Axes *ca;
    Line *cl;
    Surface *cs;
    Patch *cp;
    Text *ct;    

    // objects containers //
    std::vector< Figure > vFigure;
    std::vector< Layer > vLayer;
    std::vector< Axes > vAxes; 
    std::vector< Line > vLine;
    std::vector< Surface > vSurface;
    std::vector< Patch > vPatch;
    std::vector< Text > vText;

    // objects counter //
    int iFigure;
    int iLayer;
    int iAxes;
    int iLine;
    int iSurface;
    int iPatch;
    int iText;

    // Selected object //
    int iAxesSelected;
    
    // coordinate transform  //
    // figure coordination
    float ctx2(double x);
    float cty2(double y);
    // axes coordination
    float ctx(double x);
    float cty(double y);
    float ct3x(double x);
    float ct3y(double y);
    float ct3z(double z);
    
    int figure();

    // display //
    void display_figure();    
    void display_layer();
    void display_layer2();
    
    void display_axes();
    void display_axes_2d();
    void display_axes_3d();
    void display_axes_colorbar();
    
    void display_line();
    
    void display_surface();
    void display_surface_2d();
    void display_surface_3d();
    void display_pcolor();
    void display_contour();

    void display_patch();
    void display_patch_2d();
    void display_patch_3d();


    void display_bar();

    void display_text();

    // mouse //
    void Layer_mouse(int button, int state, int x, int y );
    void Axes_mouse(int button, int state, int x, int y );
    void Axes_motion(int x, int y );


    void surface_config();
    void line_config();
    void patch_config();
    tcvec Index2TrueColor(dvec IC);

 public:
    
    MatPlot();
    ~MatPlot();

    void virtual DISPLAY(){};

    void inline debug1(){is_debug1=1;}
    void inline debug2(){is_debug2=1;}

    // GLUT Callback Functions ///
    void display();
    void reshape(int w, int h);
    void mouse(int button, int state, int x, int y );
    void motion(int x, int y );
    void passivemotion(int x,int y);
    void keyboard(unsigned char key, int x, int y);
    
    // Layer ///
    int layer();
    //int layer(string s);    
    int layer(std::string s,int Visible);
    int frame(std::string s,int Visible);// do not use

    // Axes ///
    
    int axes();
    int gca();
    int subplot(int m,int n,int p);
    
    int colorbar();

    void axis(double xMin,double xMax,double yMin,double yMax);
    void axis(double xMin,double xMax,double yMin,double yMax,double zMin,double zMax);

    void axis(std::string s);
    void axis(int s);
    
    void grid(std::string s);
    void grid(int s);

    void ticklabel(std::string s);
    void ticklabel(int s);

    void title(std::string s);
    void xlabel(std::string s);
    void ylabel(std::string s);
    void zlabel(std::string s);

    //void xlim(double min,double max);
    //void xlim(string s);

    //void legend(string s,int N);

    //int plotyy

    void mouse_capture(double *xmouse,double *ymouse);   
    
    // set, General Object Handling ///
    void set(std::string v);
    void set(float v);  
    void set(std::string p,float v);
    void set(std::string p,std::string v);
    void set(int h,std::string p,std::string v);
    void set(int h,std::string p,float v);      
    int gco();
    
    // Line ///
    
    int begin();//do not use
    void end();//do not use
    void vertex(double x,double y);
    void vertex(double x,double y,double z);

    int line();
    int line(dvec x,dvec y);
    int line(dvec x,dvec y,dvec z);
    //line(X,Y)
    //line(X,Y,Z)

    int plot(dvec y);
    int plot(dvec x,dvec y);    
    //int plot(dmat Y);
    //int plot(dvec x,dmat Y);
    //int plot(dmat X,dmat Y);
    int plot(std::valarray<double> x,std::valarray<double> y);
    
    int plot3(dvec x,dvec y,dvec z);
    //int plot3(dvec X,dvec Y,dvec Z);

    int semilogx(dvec x,dvec y);
    int semilogy(dvec x,dvec y);
    //int loglog(dvec y);
    int loglog(dvec x,dvec y);    
    
    //int polar(dvec theta,dvec rho);

    void vertex(double x,double y,double ep,double em);
    int errorbar(dvec x,dvec y,dvec e);
    int errorbar(dvec x,dvec y,dvec ep,dvec em);

    //int quiver(U,V);
    //int quiver(X,Y,U,V);
    
    //int scatter(X,Y,S,C)
    //int scatter(X,Y,S)
    //int scatter(X,Y)

    // Surface, Contour ///
    dmat peaks(int n);
    //dmat peaks(int m,int n);
    //dmat peaks(int m,int n,string type);

    int surface();
    int surface(dmat Z);
    int surface(dmat Z,dmat C);
    int surface(dmat Z,tcmat C); //!!   
    int surface(dvec x,dvec y,dmat Z);
    int surface(dvec x,dvec y,dmat Z,dmat C);
    int surface(dvec x,dvec y,dmat Z,tcmat C);//!!
    int surface(dmat X,dmat Y,dmat Z);
    int surface(dmat X,dmat Y,dmat Z,dmat C);
    int surface(dmat X,dmat Y,dmat Z,tcmat C);//!!
    
    int pcolor(dmat C);
    int pcolor(tcmat C);
    int pcolor(dvec x,dvec y,dmat C);
    int pcolor(dvec x,dvec y,tcmat C);
    int pcolor(dmat X,dmat Y,dmat C);
    int pcolor(dmat X,dmat Y,tcmat C);

    int contour(dmat Z);
    int contour(dmat Z,int n);
    int contour(dmat Z,dvec v);
    int contour(dvec x, dvec y, dmat Z);
    int contour(dvec x, dvec y, dmat Z,int n);
    int contour(dvec x, dvec y, dmat Z,dvec v);
    //int contour(dmat X, dmat Y, dmat Z);
    //int contour(dmat X, dmat Y, dmat Z,int n);
    //int contour(dmat X, dmat Y, dmat Z,dvec v);

    //int mesh(dmat Z);
    //int mesh(dmat Z,dmat C);
    //int mesh(dmat Z,tcmat C);    
    int mesh(dvec x, dvec y, dmat Z);
    //int mesh(dvec x, dvec y, dmat Z,dmat C);
    //int mesh(dvec x, dvec y, dmat Z,tcmat C);    
    //int mesh(dmat X,dmat Y,dmat Z);
    //int mesh(dmat X,dmat Y,dmat Z,dmat C);
    //int mesh(dmat X,dmat Y,dmat Z,tcmat C);
    // meshc()
    // meshz()

    int surf(dvec x, dvec y, dmat Z);

    // Patch ///

    int patch();
    int patch(dmat X,dmat Y);
    int patch(dmat X,dmat Y,dvec C);
    int patch(dmat X,dmat Y,tcvec C);    
    int patch(dmat X,dmat Y,dmat Z);//!!
    int patch(dmat X,dmat Y,dmat Z,dvec C);//!!    
    int patch(dmat X,dmat Y,dmat Z,tcvec C);//!!
    //int patch(dmat X,dmat Y,tcmat C);
    //int patch(dmat X,dmat Y,dmat Z,tcmat C);

    int bar(dvec y);
    int bar(dvec y,float width);
    int bar(dvec x,dvec y);
    int bar(dvec x,dvec y,float width);

    //int bar(Y)
    //int bar(Y,float width);
    //int bar(Y,string style);
    //int bar(Y,float width,string style);

    //int bar(x,Y)
    //int bar(x,Y,float width);
    //int bar(x,Y,string style);
    //int bar(x,Y,float width,string style);
      
    //int hist(y);
    //int hist(y,x);
    //int hist(y,nbins);

    //int pie(dvec x);
    //int pie(dvec x, vector<int> Explode);

    // Text ///
    //TODO: more fonts    
    int text();
    int text(double x,double y,std::string s);
    void set_font(char font_[],int size);
    void ptext(float x,float y,std::string s);
    void ptext3(float x,float y,float z,std::string s);
    void ptext3c(float x,float y,float z,std::string s);

    // Colors ///
    void color(float r,float g,float b);
    std::vector<float> colormap(std::string c,float t);
    void colormap(std::string c);
    void colormap(tcvec c);

    void gray();
    void jet();
    void hsv();
    void hot();
    void cool();
    void spring();
    void summer();
    void autumn();
    void winter();

    std::vector<float> map2color(double x,double xmin,double xmax);
    
    void Shading(std::string c);
    void shading(std::string c);
    std::vector<float> ColorSpec2RGB(std::string c);
    std::string rgb2colorspec(std::vector<float> rgb);

    // print ///
    void print();
    void print(std::string filename, std::string title, GLint format);

};

