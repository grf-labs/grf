/****************************************************************************
Copyright (c) 2011 Yuichi Katori All Rights Reserved
License: Gnu Public license (GPL) v3
Author: Yuichi Katori (yuichi.katori@gmail.com)
Project:MATPLOT++ (MATLAB-like plotting tool in C++).
Version:0.3.13

TODO: This file only works in Debug mode. Revise!

****************************************************************************/
using namespace std;
#include "matplotpp.h"

#if defined(_WIN32) || defined(_WIN64)
#  define fmax max
#  define fmin min
#  pragma warning (disable:4996)
#  define snprintf sprintf_s
#endif

/// Figure
void Figure::add_child(int i){Children.push_back(i);}
/// Axes
void Axes::reset(){
    num_child=0;
    xmin=1e99; xmax=-1e99;
    ymin=1e99; ymax=-1e99;
    zmin=1e99; zmax=-1e99;	
}
void Axes::config(){
    float extent=0,extent_linear=0.03;
    if((XLimMode==0)&&(xmax>xmin)){
	if(XScale==0){extent=extent_linear;}
	if(XScale==1){extent=0;}
	XLim[0]=xmin-extent*(xmax-xmin); 
	XLim[1]=xmax+extent*(xmax-xmin);
    }
    if((YLimMode==0)&&(ymax>ymin)){
	if(YScale==0){extent=extent_linear;}
	if(YScale==1){extent=0;}
	YLim[0]=ymin-extent*(ymax-ymin); 
	YLim[1]=ymax+extent*(ymax-ymin);
    }
    if((ZLimMode==0)&&(zmax>zmin)){
	ZLim[0]=zmin-extent*(zmax-zmin); 
	ZLim[1]=zmax+extent*(zmax-zmin);
    }	
    //printf("Z: %d,%f,%f\n",ZLimMode,ZLim[0],ZLim[1]);
    //if(num_child){Visible=1;}else{Visible=0;}
    
    XTick=make_tick(XLim[0],XLim[1]);
    YTick=make_tick(YLim[0],YLim[1]);
    ZTick=make_tick(ZLim[0],ZLim[1]);
}
    
int Axes::ID(){return id;}
int Axes::selected(){return Selected;}
void Axes::selected(int i){Selected=i;}
void Axes::add_child(int i){Children.push_back(i);}
    
dvec Axes::make_tick(double min,double max){
	int i,j;
	double dg;
	double x,y;
	int z;
	x=fabs(max-min);
	z=(int)log10(x);
	y=x/pow((double)10,(double)z);
	dg=1*pow((double)10,(double)z);
	if(y<5){dg=0.5*pow((double)10,(double)z);}
	if(y<2){dg=0.2*pow((double)10,(double)z);}

	double min0;
	min0=min-fmod(min,dg);j=0;

	dvec tick;
	tick.clear();
	if(max>min){ i=-2; while(max>=min0+dg*i){ if(min<=min0+dg*i){ tick.push_back(min0+dg*i); j++; } i+=1;} }
	if(max<min){ i=-2; while(max<=min0-dg*i){ if(min>=min0-dg*i){ tick.push_back(min0-dg*i); j++; } i+=1;} }
	return tick;
}


/// Line

void Line::reset(){
    XData.clear();
    YData.clear();
    ZData.clear();
    YPData.clear();
    YMData.clear();
}
void Line::color(float r,float g,float b){
    //Color[0]=r;
    //Color[1]=g;
    //Color[2]=b;
}
/// Surface

/// Text
Text::Text(int id_){
    id=id_;
    type=0;
}


/// Layer
Layer::Layer(int id_){ 
    id=id_;
    Visible=1;
    Children.clear();	
}
/// Patch

void Layer::add_child(int i){Children.push_back(i);}

//// MatPlot //
MatPlot::MatPlot(){
    is_debug1=0;
    is_debug2=0;

    if(is_debug1){cout<<"init()..."<<endl;}

    mode=0;
    init_level=0;

    iAxesSelected=-1;

    colormap("Jet");

    xPassive=100;
    yPassive=0;
    
    if(is_debug1){cout<<"init()...done"<<endl;}
}
MatPlot::~MatPlot(){
    vFigure.clear();
    vAxes.clear();
    vLine.clear();	
}
//void MatPlot::virtual DISPLAY(){};


// display ///

void MatPlot::display(){

    int is_done=0;
    int mode_next=-1;
    while(!is_done){
	if(mode==0){// Allocation of objects

	    if(is_debug1){cout<<"============================= allocation ..."<<endl;}

	    iFigure =0;
	    iLayer  =0; 	
	    iAxes   =0;
	    iLine   =0;
	    iSurface=0;
	    iPatch  =0;
	    iText   =0;	    
    
	    if(init_level==0){vFigure.clear();}
	    if(init_level<=1){vLayer.clear();}
	    if(init_level<=2){vAxes.clear();}
	    if(init_level<=3){vLine.clear();}
	    if(init_level<=3){vSurface.clear();}
	    if(init_level<=3){vText.clear();}
	    if(init_level<=3){vPatch.clear();}
	    
	    //default objects
	    figure();
	    layer();
	    axes();

	    DISPLAY();

	    if(is_debug1){cout<<"done"<<endl;}
    
	    mode_next=1;
	}
	if(mode==1){// Configuration of objects
	    if(is_debug1){cout<<"============================= configulation ..."<<endl;}
	    iFigure =0;
	    iLayer  =0; 	
	    iAxes   =0;
	    iLine   =0;
	    iSurface=0;
	    iPatch  =0;
	    iText   =0;

	    //default objects
	    figure();
	    layer();
	    axes();

	    DISPLAY();
	
	    mode_next=2;

	    // checking number of objects. if it is inconsistent, then objects are reallocated.
	    if(vFigure.size() !=iFigure ){ mode_next=0; if(is_debug1){cout<<"not match in Figure"<<endl;} }
	    if(vAxes.size()   !=iAxes   ){ mode_next=0; if(is_debug1){cout<<"not match in Axes"<<endl;}}
	    if(vLine.size()   !=iLine   ){ mode_next=0; if(is_debug1){cout<<"not match in Line"<<endl;}}
	    if(vSurface.size()!=iSurface){ mode_next=0; if(is_debug1){cout<<"not match in Surface"<<endl;}}
	    if(vText.size()   !=iText   ){ mode_next=0; if(is_debug1){cout<<"not match in Text"<<endl;}}
	    if(vLayer.size()  !=iLayer  ){ mode_next=0; if(is_debug1){cout<<"not match in Layer"<<endl;}}

	    for(int i=0;i<vAxes.size();++i){ vAxes[i].config(); }

	    if(is_debug1){cout<<"done"<<endl;}
	}
	if(mode==2){// display

	    if(is_debug1){cout<<"============================= display ..."<<endl;}	

	    display_figure();

	    if(is_debug1){cout<<" done "<<endl;}

	    mode_next=1;
	    is_done=1;
	}
	mode=mode_next;
    }	
}

void MatPlot::color(float r,float g,float b){
    if(mode){
	//cl->Color[0]=r;
	//cl->Color[1]=g;
	//cl->Color[2]=b;
    }
}
// coordinate transform  ///

// figure coordination
float MatPlot::ctx2(double x){
    double t;
    if(ca->XScale==0){//linear
	return ca->Position[0] +ca->Position[2]*( (x-ca->XLim[0])/(ca->XLim[1]-ca->XLim[0]) );
    }
    if(ca->XScale==1){//log
	t=( log10(x) - log10(ca->XLim[0]) )/( log10(ca->XLim[1])-log10(ca->XLim[0]) );
	if(x<=0){t=-1;}
	return ca->Position[0] +ca->Position[2]*t;
    }
    
}
float MatPlot::cty2(double y){ 
    if(ca->YScale==0){//linear
	return ca->Position[1] +ca->Position[3]*( (y-ca->YLim[0])/(ca->YLim[1]-ca->YLim[0]) );
    }
    if(ca->YScale==1){//log
	return ca->Position[1] +ca->Position[3]*( log10(y) - log10(ca->YLim[0]) )/( log10(ca->YLim[1])-log10(ca->YLim[0]) );
    }
}
// axes coordination
float MatPlot::ctx(double x){
    if(ca->XScale==0){//linear
	return (x-ca->XLim[0])/(ca->XLim[1]-ca->XLim[0]);
    }
    if(ca->XScale==1){//log
	return ( log10(x) - log10(ca->XLim[0]) )/( log10(ca->XLim[1])-log10(ca->XLim[0]) );
    }
}
float MatPlot::cty(double y){
    if(ca->YScale==0){//linear
	return (y-ca->YLim[0])/(ca->YLim[1]-ca->YLim[0]);
    }
    if(ca->YScale==1){//log
	return ( log10(y) - log10(ca->YLim[0]) )/( log10(ca->YLim[1])-log10(ca->YLim[0]) );
    }
}
float MatPlot::ct3x(double x){	
    return -1+2*(x-ca->XLim[0])/(ca->XLim[1]-ca->XLim[0]);
}
float MatPlot::ct3y(double y){	
    return -1+2*(y-ca->YLim[0])/(ca->YLim[1]-ca->YLim[0]);
}
float MatPlot::ct3z(double z){	
    return -1+2*(z-ca->ZLim[0])/(ca->ZLim[1]-ca->ZLim[0]);
}
// set ///
/// set(v)
void MatPlot::set(string v){
    int h=gco();
    int tObj=h%100;
    int iObj=h/100;
    //if( (tObj==tLine) && (iObj<vLine.size()) ){
    if( v=="k" ){ set(h,"COLOR","k"); }
    if( v=="r" ){ set(h,"COLOR","r"); }
    if( v=="b" ){ set(h,"COLOR","b"); }
    if( v=="g" ){ set(h,"COLOR","g"); }
    if( v=="c" ){ set(h,"COLOR","c"); }
    if( v=="m" ){ set(h,"COLOR","m"); }
    if( v=="y" ){ set(h,"COLOR","y"); }
    if( v=="w" ){ set(h,"COLOR","w"); }

    if( v=="dr" ){ set(h,"COLOR","dr"); }
    if( v=="db" ){ set(h,"COLOR","db"); }
    if( v=="dg" ){ set(h,"COLOR","dg"); }
    if( v=="dc" ){ set(h,"COLOR","dc"); }
    if( v=="dm" ){ set(h,"COLOR","dm"); }
    if( v=="dy" ){ set(h,"COLOR","dy"); }

    if( v=="lr" ){ set(h,"COLOR","lr"); }
    if( v=="lb" ){ set(h,"COLOR","lb"); }
    if( v=="lg" ){ set(h,"COLOR","lg"); }
    if( v=="lc" ){ set(h,"COLOR","lc"); }
    if( v=="lm" ){ set(h,"COLOR","lm"); }
    if( v=="ly" ){ set(h,"COLOR","ly"); }

    if( v=="ur" ){ set(h,"COLOR","ur"); }
    if( v=="ub" ){ set(h,"COLOR","ub"); }
    if( v=="ug" ){ set(h,"COLOR","ug"); }
    if( v=="uy" ){ set(h,"COLOR","uy"); }
    if( v=="uc" ){ set(h,"COLOR","uc"); }
    if( v=="up" ){ set(h,"COLOR","up"); }
    if( v=="uo" ){ set(h,"COLOR","uo"); }
    if( v=="um" ){ set(h,"COLOR","um"); }
    if( v=="ubr" ){set(h,"COLOR","ubr"); }

    if( v=="-"  ){ set(h,"LineStyle","-");   set(h,"Marker","none"); }
    if( v=="- -"){ set(h,"LineStyle","- -"); set(h,"Marker","none"); }
    if( v==":"  ){ set(h,"LineStyle",":");   set(h,"Marker","none"); }
    if( v=="-." ){ set(h,"LineStyle","-.");  set(h,"Marker","none"); }
    
    if( v=="." ){ set(h,"Marker","."); set(h,"LineStyle","none"); }
    if( v=="+" ){ set(h,"Marker","+"); set(h,"LineStyle","none"); }
    if( v=="x" ){ set(h,"Marker","x"); set(h,"LineStyle","none"); }
    if( v=="d" ){ set(h,"Marker","d"); set(h,"LineStyle","none"); }
    if( v=="^" ){ set(h,"Marker","^"); set(h,"LineStyle","none"); }
    if( v=="v" ){ set(h,"Marker","v"); set(h,"LineStyle","none"); }
    if( v=="o" ){ set(h,"Marker","o"); set(h,"LineStyle","none"); }
    if( v=="*" ){ set(h,"Marker","*"); set(h,"LineStyle","none"); }
    if( v=="s" ){ set(h,"Marker","s"); set(h,"LineStyle","none"); }
    if( v==">" ){ set(h,"Marker",">"); set(h,"LineStyle","none"); }
    if( v=="<" ){ set(h,"Marker","<"); set(h,"LineStyle","none"); }
    if( v=="p" ){ set(h,"Marker","p"); set(h,"LineStyle","none"); }
    if( v=="h" ){ set(h,"Marker","h"); set(h,"LineStyle","none"); }

	//}
    
}
void MatPlot::set(float v){
    int h=gco();
    int tObj=h%100;
    int iObj=h/100;
    if( (tObj==tLine) && (iObj<vLine.size()) ){
	vLine[iObj].LineWidth=v;
	vLine[iObj].MarkerSize=v;	
    }
}
/// set(p,v),set(h,p,v)
void MatPlot::set(string p,string v){
    int h=gco();
    set(h,p,v);
}
void MatPlot::set(string p,float v){
    int h=gco();
    set(h,p,v);
}
void MatPlot::set(int h,string p,string v){
    int tObj=h%100;
    int iObj=h/100;
    // Axes //
    if( (tObj==tAxes) && (iObj<vAxes.size()) ){
	if(p=="TickDir"){ vAxes[iObj].TickDir=v; }
	//if(p==""){ vLine[iObj].=v; }	
    }
    // Line //
    if( (tObj==tLine) && (iObj<vLine.size()) ){

	if(p=="COLOR"){ vLine[iObj].Color=v; vLine[iObj].MarkerEdgeColor=v; }
	if(p=="Color"){ vLine[iObj].Color=v; }
	if(p=="Marker"){ vLine[iObj].Marker=v; }
	if(p=="LineStyle"){ vLine[iObj].LineStyle=v; }
	if(p=="MarkerEdgeColor"){ vLine[iObj].MarkerEdgeColor=v; }
	if(p=="MarkerFaceColor"){ vLine[iObj].MarkerFaceColor=v; }
	//if(p==""){ vLine[iObj].=v; }	
    }
    // Surface //
    if( (tObj==tSurface) && (iObj<vSurface.size()) ){
	if(p=="COLOR"){ vSurface[iObj].EdgeColor=v; }
	if(p=="LineStyle"){ vSurface[iObj].LineStyle=v; }
	if(p=="EdgeColor"){ vSurface[iObj].EdgeColor=v; }
	if(p=="FaceColor"){ vSurface[iObj].FaceColor=v; }
	//if(p==""){ vLine[iObj].=v; }	
    }    
    // Patch //
    if( (tObj==tPatch) && (iObj<vPatch.size()) ){
	//if(p==""){ vPatch[iObj].=v; }	
	if(p=="COLOR"){ vPatch[iObj].EdgeColor=v; }
	if(p=="LineStyle"){ vPatch[iObj].LineStyle=v; }
	if(p=="EdgeColor"){ vPatch[iObj].EdgeColor=v; }
	if(p=="FaceColor"){ vPatch[iObj].FaceColor=v; }

    }

}
void MatPlot::set(int h,string p,float v){
    int tObj=h%100;
    int iObj=h/100;
    // Line //
    if( (tObj==tLine) && (iObj<vLine.size()) ){
	if(p=="LineWidth"){ vLine[iObj].LineWidth=v; }	
	if(p=="MarkerSize"){ vLine[iObj].MarkerSize=v; }	
	
    }
    // Surface //
    if( (tObj==tSurface) && (iObj<vSurface.size()) ){
	if(p=="LineWidth"){ vSurface[iObj].LineWidth=v; }
	//if(p==""){ vLine[iObj].=v; }	
    }
    // Patch //
    if( (tObj==tPatch) && (iObj<vPatch.size()) ){
	if(p=="LineWidth"){ vPatch[iObj].LineWidth=v; }
	//if(p==""){ vLine[iObj].=v; }	
    }
    
}
int MatPlot::gco(){
    return hObj;
}
//void MatPlot::set_line_width(float linewidth){
//    cl->LineWidth=linewidth;
//}

// Events ///
void MatPlot::reshape(int w, int h){
    window_w=w;
    window_h=h;
    if(mode==2){
	glViewport(0,0,w,h);
    }
    //cout <<"window size: "<< w <<" "<<h<<endl;
}
void MatPlot::mouse(int button, int state, int x, int y ){
    Layer_mouse(button,state,x,y);
    Axes_mouse(button,state,x,y);	
    display();
}
void MatPlot::motion(int x, int y ){
    Axes_motion(x,y);
    //glutPostRedisplay(); 
}
void MatPlot::passivemotion(int x,int y){
    xPassive=x;
    yPassive=y;
    //cout <<"Passive: "<<x<<" "<<y<<endl;
}
void MatPlot::keyboard(unsigned char key, int x, int y){
    switch(key){
    case 'q':
	exit(0);
	break;
    case 'p':
	print();
	break;
    }
}
// print ///
void MatPlot::print(){
    FILE *fp;
    int state = GL2PS_OVERFLOW, buffsize = 0;
    
    fp = fopen("out.eps", "wb");
    printf("Writing 'out.eps'... ");
    while(state == GL2PS_OVERFLOW){
	buffsize += 2024*2024;
	gl2psBeginPage("test", "gl2ps", NULL, GL2PS_EPS, GL2PS_SIMPLE_SORT, 
		       GL2PS_USE_CURRENT_VIEWPORT, 
		       GL_RGBA, 0, NULL, 0, 0, 0, buffsize, fp, "out.eps");
	display();
	state = gl2psEndPage();
    }
    fclose(fp);
    printf("Done!\n");
}

void MatPlot::print(std::string filename, std::string title, GLint format){
    FILE *fp;
    int state = GL2PS_OVERFLOW, buffsize = 0;
    
    fp = fopen(filename.c_str(), "wb");
    cout << "Writing '" << filename << "'... ";
    while(state == GL2PS_OVERFLOW){
	buffsize += 2024*2024;
	gl2psBeginPage(title.c_str(), "", NULL, format, GL2PS_SIMPLE_SORT, 
		       GL2PS_USE_CURRENT_VIEWPORT, 
		       GL_RGBA, 0, NULL, 0, 0, 0, buffsize, fp, filename.c_str());
	display();
	state = gl2psEndPage();
    }
    fclose(fp);
    cout << "Done!" << endl;
}


// Figure ///

int MatPlot::figure(){
    int h=iFigure*100 + tFigure; hObj=h;
    if(is_debug1){printf("mode: %d handle: %4d Figure\n",mode,h);}
    if(mode==0){	    
	vFigure.push_back(Figure(h));
    }
    if(mode==1){
	
    }
    if(iFigure<vFigure.size()){cf=&vFigure[iFigure];}
    iFigure++;
    return h;
}

void MatPlot::display_figure(){

    if(is_debug1){printf("mode: %d handle: %4d Figure %d layers\n",
			 mode,cf->id, cf->Children.size());}

    int tObj;// type of child object
    int iObj;// index of child object
    if(cf->Visible){
	glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_NEVER);

	glClearColor(1, 1, 1, 0.);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glViewport(0,0, (int)(window_w),(int)(window_h));
	glLoadIdentity();
	//gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
	glOrtho(0, 1, 0, 1, -1, 3);
	
	for(int i=0;i<cf->Children.size();++i){
	    tObj=cf->Children[i]%100;
	    iObj=cf->Children[i]/100;
	    //printf("Obj t:%3d i:%3d \n",tObj,iObj);
	    
	    if(tObj==tLayer){
		cfr=&vLayer[iObj];
		if(cfr->Visible){ display_layer();}
	    }
	}
	/*
	  for(int i1=0;i1<cf->Children.size();++i1){
	  tObj1=cf->Children[i1]%100;
		iObj1=cf->Children[i1]/100;		
		if(tObj1==tLayer){
		cfr=&vLayer[iObj1];
		if(cfr->Visible){ display_layer();}
		
		for(int i2=0;i2<cfr->Children.size();++i2){
		tObj2=cfr->Children[i2]%100;
		iObj2=cfr->Children[i2]/100;
		if(tObj2==tAxes){
		ca=&vAxes[iObj2];
		if(ca->Children.size()){display_axes();}
		}
		}
		}
		}
	*/
	    	
	//Visible(1);

	display_layer2();
	
	glFlush();
	glViewport(0,0,window_w,window_h);
	
    }
    glutSwapBuffers();
}
    

// Layer ///
int MatPlot::layer(){
    int h=iLayer*100 + tLayer;hObj=h;
    if(is_debug1){printf("mode: %d handle: %4d Layer\n",mode,h);}
    
    if(mode==0){	    
	vLayer.push_back(Layer(h));
	cf->add_child(h);
    }
    if(mode==1){
	
    }
    if(iLayer<vLayer.size()){cfr=&vLayer[iLayer];}
    iLayer++;

    axes();//default axes

    return h;
}

void MatPlot::display_layer(){
    if(is_debug1){printf("mode: %d handle: %4d Layer ( %s ) \n",
			 mode,cfr->id,cfr->layername.c_str());}

    int tObj;// type  of child object
    int iObj;// index of child object	
    for(int i=0;i<cfr->Children.size();++i){
	tObj=cfr->Children[i]%100;
	iObj=cfr->Children[i]/100;
	if(tObj==tAxes){
	    ca=&vAxes[iObj];
	    display_axes();
	}
    }//i
}
void MatPlot::display_layer2(){
    int l,t,w,h,r;
    string s;
    l=1;
    t=1;
    w=20;//button_width;
    h=20;//button_height;
    r=3;

    if(xPassive<25){
	

	glViewport(0,0, (int)(window_w),(int)(window_h));
	glLoadIdentity();
	gluOrtho2D( 0.0, (int)(window_w), (int)(window_h),0 );
	
	glDisable(GL_LINE_STIPPLE); 
	gl2psDisable(GL2PS_LINE_STIPPLE);

	glLineWidth(2);
	glColor3d(0,0,1);
	
	for(int j=0;j<vLayer.size();++j){
	    //glBegin(GL_LINE_STRIP);
	    //glVertex2d(l   ,h*j );
	    //glVertex2d(l   ,h*j+h );
	    //glVertex2d(l+w ,h*j+h );
	    //glVertex2d(l+w ,h*j );
	    //glEnd();
	    if(vLayer[j].Visible==1){// [v]
		    
		glBegin(GL_LINE_STRIP);
		glVertex2d(l+r   ,h*j+r );
		glVertex2d(l+r   ,h*j+h-r );
		glVertex2d(l+w-r ,h*j+h-r );
		glVertex2d(l+w-r ,h*j+r );
		glVertex2d(l+r   ,h*j+r );
		glEnd();
		
		glBegin(GL_LINE_STRIP);
		glVertex2d(l+9 ,h*j+5 );
		glVertex2d(l+8 ,h*j+15 );
		glVertex2d(l+15  ,h*j+7 );
		glEnd();
		//ptext(5,h*j+h-3,"[v]");
	    }
	    if(vLayer[j].Visible==0){// []
		glBegin(GL_LINE_STRIP);
		glVertex2d(l+r   ,h*j+r );
		glVertex2d(l+r   ,h*j+h-r );
		glVertex2d(l+w-r ,h*j+h-r );
		glVertex2d(l+w-r ,h*j+r );
		glVertex2d(l+r   ,h*j+r );
		glEnd();
		
		//ptext(5,h*j+h-3,"[  ]");
	    }
	    //ptext(22,h*j+h-3,vLayer[j].layername);
	    glColor3f(0,0,1);
	    glRasterPos2f( 22, h*j+h-6 );
	    s=vLayer[j].layername;
	    gl2psText(s.c_str(), "Arial", 12);
	    for(int i=0;i<s.size();++i){
		glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );
	    }
	}
    }
}

/// events (mouse)
void MatPlot::Layer_mouse(int button, int state, int x, int y ){
    int l,t,w,h;
    int is;
    double dt;
    l=1;
    t=1;
    w=20;//button_width;
    h=20;//button_height;
    //cout <<"button:"<< button<<" state:"<<state<<" x:"<<x<<" y:"<<y<<endl;
    if ( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) {// Left Click
	for(int j=0;j<vLayer.size();++j){
	    if((l<x)&&(x<w)&&(h*j<y)&&(y<h*j+h)){
		
		//visibility
		if(vLayer[j].Visible==1){is=0;}
		if(vLayer[j].Visible==0){is=1;}
		vLayer[j].Visible=is;
		
		//double click to focus layer
		time_layer_clicked=clock();
		dt=(float)(time_layer_clicked - time_layer_clicked_last)/CLOCKS_PER_SEC;
		if(dt < 0.05 ){
		    for(int k=0;k<vLayer.size();++k){ vLayer[k].Visible=0; }
		    vLayer[j].Visible=1;
		    //cout <<"!!"<<endl;
		}
		//cout <<"click: "<<time_layer_clicked<<" "<<CLOCKS_PER_SEC<<" "<<dt<<endl;
		time_layer_clicked_last=time_layer_clicked;

		//cout<<"pushed:"<<j<<" Visible:"<<vLayer[j].Visible<<" "<<is<<endl;
		for(int i=0;i<vAxes.size();++i){
		    if(vAxes[i].Parent=vLayer[j].id){vAxes[i].Visible=is;}
		}
	    }
	}
    }	
}
int MatPlot::layer(string s,int Visible){
    int h=layer();
    if(mode==0){
	cfr->Visible=Visible;
	cfr->layername=s;
    }
    return h;
}
int MatPlot::frame(string s,int Visible){
    int h=layer();
    if(mode==0){
	cfr->Visible=Visible;
	cfr->layername=s;
    }
    return h;
}

// Axes ///
int MatPlot::gca(){
    return ca->id;	
}
int MatPlot::axes(){
    int h=iAxes*100 + tAxes; hObj=h;
    //int h=handle(iAxes,tAxes);
    if(is_debug1){printf("mode: %d handle: %4d Axes \n",mode,h);}
    
    if(mode==0){
	cfr->add_child(h);
	vAxes.push_back(Axes(h));
	}
    if(mode==1){
	if(iAxes<vAxes.size()){ vAxes[iAxes].reset(); }
    }
    if(iAxes<vAxes.size()){ca=&vAxes[iAxes];}
    ca->Visible=cfr->Visible;
    iAxes++;
    return h;
}
void MatPlot::display_axes(){
    if(is_debug1){printf("mode: %d handle: %4d Axes  #Children %d\n",
			 mode,ca->id,ca->Children.size());}
	
	if(ca->Children.size()){
	    if(ca->View==0){//2D
		display_axes_2d();
	    }
	    if(ca->View==1){//2D
		display_axes_3d();
	    }	
	}
	if(ca->View==2){//colorbar
	    display_axes_colorbar();
	}

	// childlen //
	int tObj;// type  of child object
	int iObj;// index of child object
	for(int i=0;i<ca->Children.size();++i){
	    tObj=ca->Children[i]%100;
	    iObj=ca->Children[i]/100;
	    if(tObj==tLine){
		cl=&vLine[iObj];
		display_line();
	    }
	    if(tObj==tSurface){
		cs=&vSurface[iObj];
		display_surface();
	    }
	    if(tObj==tText){
		ct=&vText[iObj];
		display_text();
	    }
	    if(tObj==tPatch){
		cp=&vPatch[iObj];
		display_patch();
	    }
	}//i
    }
/// display
void MatPlot::display_axes_2d(){
	
	
    char ctmp[100];
    float l,b,w,h;//left,bottom,width,height
    float r=0.01;

    l=ca->Position[0];
    b=ca->Position[1];
    w=ca->Position[2];
    h=ca->Position[3];
    
    // Viewport Figure (VpF) for drawing axes
    glViewport(0,0, (int)(window_w),(int)(window_h));
    glLoadIdentity();
    gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
    
    glDisable(GL_LINE_STIPPLE); 
    gl2psDisable(GL2PS_LINE_STIPPLE);
    
    float x_axis_location=1,y_axis_location=1;
    if(ca->XAxisLocation=="top"  ){x_axis_location=1;}else{x_axis_location=-1;}
    if(ca->YAxisLocation=="right"){y_axis_location=1;}else{y_axis_location=-1;}
 
    int char_w=6,char_h=12;
    float offset=0.01;
    int num_char=4;
   
    int gridlinestyle;

    int tickdir=1;//1:in, -1:out
    
    if(ca->Box){
	// box //  
	glLineWidth(ca->LineWidth);
	gl2psLineWidth(ca->LineWidth);
	if(ca->Selected){
	    glLineWidth( 2* ca->LineWidth);
	    gl2psLineWidth(ca->LineWidth);
	}
	glColor3f(0,0,0);
	glBegin(GL_LINE_LOOP);	    
	glVertex2d(l,  b);
	glVertex2d(l+w,b);
	glVertex2d(l+w,b+h);
	glVertex2d(l,  b+h);
	glEnd();

	// mouse capture //
	if(ca->Selected){
	    sprintf(ctmp,"Mouse:(%f,%f)",ca->XMouse,ca->YMouse);
	    ptext( l,b+h+r,ctmp );
	}

	// Grid //
	gridlinestyle=3;
	if(ca->GridLineStyle=="-"  ){gridlinestyle=1;}
	if(ca->GridLineStyle=="- -"){gridlinestyle=2;}
	if(ca->GridLineStyle==":"  ){gridlinestyle=3;}
	if(ca->GridLineStyle=="-." ){gridlinestyle=4;}

	if(ca->XGrid){
	    glLineWidth(ca->LineWidth);
	    gl2psLineWidth(ca->LineWidth);
	    for(int i=0;i<ca->XTick.size();++i){

		//cout <<"grid "<<gridlinestyle<<" "<<ca->XTick[i]<<endl;

		if(gridlinestyle==1){// -
		    glDisable(GL_LINE_STIPPLE); 
		    gl2psDisable(GL2PS_LINE_STIPPLE);
		}
		if(gridlinestyle==2){//- -
		    glEnable(GL_LINE_STIPPLE);
		    glLineStipple(1, 0xF0F0);
		    gl2psEnable(GL2PS_LINE_STIPPLE);
		}
		if(gridlinestyle==3){//:		    
		    glEnable(GL_LINE_STIPPLE);
		    glLineStipple(1, 0xCCCC);
		    gl2psEnable(GL2PS_LINE_STIPPLE);
		}
		if(gridlinestyle==4){//-.
		    glEnable(GL_LINE_STIPPLE);	   
		    glLineStipple(1, 0x087F);
		    gl2psEnable(GL2PS_LINE_STIPPLE);
		}
		glBegin(GL_LINE_STRIP);
		glVertex2d( ctx2(ca->XTick[i]),b );
		glVertex2d( ctx2(ca->XTick[i]),b+h );//!! TODO
		glEnd();		
	    }
	}
	if(ca->YGrid){	    
	    for(int i=0;i<ca->XTick.size();++i){
		if(gridlinestyle==1){// -
		    glDisable(GL_LINE_STIPPLE); 
		    gl2psDisable(GL2PS_LINE_STIPPLE);
		}
		if(gridlinestyle==2){//- -
		    glEnable(GL_LINE_STIPPLE);
		    glLineStipple(1, 0xF0F0);
		    gl2psEnable(GL2PS_LINE_STIPPLE);
		}
		if(gridlinestyle==3){//:
		    glEnable(GL_LINE_STIPPLE);
		    glLineStipple(1, 0xCCCC);
		    gl2psEnable(GL2PS_LINE_STIPPLE);
		}
		if(gridlinestyle==4){//-.
		    glEnable(GL_LINE_STIPPLE);	   
		    glLineStipple(1, 0x087F);
		    gl2psEnable(GL2PS_LINE_STIPPLE);
		}
		glBegin(GL_LINE_STRIP);
		glVertex2d( l,  cty2(ca->YTick[i]) );
		glVertex2d( l+w,cty2(ca->YTick[i]) );
		glEnd();
	    }
	}

	// Ticks // 
	if(ca->TickDir=="in"){tickdir=1;}
	if(ca->TickDir=="out"){tickdir=-1;}
	
	glDisable(GL_LINE_STIPPLE); 
	gl2psDisable(GL2PS_LINE_STIPPLE);
	//TODO precise adjustment of tick location
	// x tick 
	for(int i=0;i<ca->XTick.size();++i){
	    glBegin(GL_LINE_STRIP);
	    glVertex2d( ctx2(ca->XTick[i]),b );
	    glVertex2d( ctx2(ca->XTick[i]),b+tickdir*0.01 );//b-0.02*h
	    glEnd();
	}
	// x tick label
	if(ca->TickLabel){
	for(int i=0;i<ca->XTick.size();++i){
	    sprintf(ctmp,"%4.1f",ca->XTick[i]);
	    //ptext( ctx2(ca->XTick[i])-0.02, b-0.025,ctmp );//b-0.05*h
	    ptext( ctx2(ca->XTick[i])-(float)num_char*char_w/window_w/2.0, 
		   b-offset-1.0*char_h/window_h,ctmp );//b-0.05*h
	}
	}
	// y tick 
	for(int i=0;i<ca->YTick.size();++i){
	    glBegin(GL_LINE_STRIP);
	    glVertex2d( l,             cty2(ca->YTick[i]) );
	    glVertex2d( l+tickdir*0.01,cty2(ca->YTick[i]) );		
	    glEnd();
	}
	// y tick label
	if(ca->TickLabel){
	for(int i=0;i<ca->YTick.size();++i){
	    sprintf(ctmp,"%4.1f",ca->YTick[i]);
	    //ptext( l-0.05,cty2(ca->YTick[i])-0.0,ctmp );
	    ptext( l-(float)num_char*char_w/window_w-offset, 
		   cty2(ca->YTick[i])-0.5*char_h/window_h,ctmp );
	}
	}
    }//Box
        
    //Title
    num_char=ca->Title.length();    
    ptext( l+w/2.0-(float)num_char*char_w/window_w/2.0, 
	   b+h+offset,
	   ca->Title );

    //XLabel
    num_char=ca->XLabel.length();    
    ptext( l+w/2.0-(float)num_char*char_w/window_w/2.0, 
	   b-offset-2.0*char_h/window_h,
	   ca->XLabel );

    //YLabel
    num_char=ca->YLabel.length();   
    ptext( l,   b+h+offset,
	   ca->YLabel );

    // Viewport Axes (VpA) for drawing lines and surfaces
    glViewport((int)(ca->Position[0]*window_w ),
	       (int)(ca->Position[1]*window_h ),
	       (int)(ca->Position[2]*window_w ),
	       (int)(ca->Position[3]*window_h ));
    glLoadIdentity();
    gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
}
void MatPlot::display_axes_3d(){
	
    char ctmp[100];
    float l,b,w,h;//left,bottom,width,height
    float r=0.01;
    
    l=ca->Position[0];
    b=ca->Position[1];
    w=ca->Position[2];
    h=ca->Position[3];
    
    // Viewport Axes 
    /*
    glViewport((int)(ca->Position[0]*window_w ),
	       (int)(ca->Position[1]*window_h ),
	       (int)(ca->Position[2]*window_w ),
	       (int)(ca->Position[3]*window_h ));
    */

    glViewport((int)(ca->Viewport3d[0]*window_w ),
	       (int)(ca->Viewport3d[1]*window_h ),
	       (int)(ca->Viewport3d[2]*window_w ),
	       (int)(ca->Viewport3d[3]*window_h ));
    
    glLoadIdentity();
    //glOrtho(-1.7, 1.7, -1.7, 1.7, -1.5, 3);
    glOrtho(-1.8, 1.8, -1.8, 1.8, -1.5, 3);
    
    gluLookAt(cos(ca->cta*PI/180)*cos(ca->phi*PI/180), 
	      sin(ca->cta*PI/180)*cos(ca->phi*PI/180), 
	      sin(ca->phi*PI/180),
	      //gluLookAt(CameraPosition[0],CameraPosition[1],CameraPosition[2],
	      ca->CameraTarget[0],  ca->CameraTarget[1],  ca->CameraTarget[2],
	      ca->CameraUpVector[0],ca->CameraUpVector[1],ca->CameraUpVector[2]);
    
    /*
    // x,y,z axes for test
    glColor3f(1,0,0);
    glBegin(GL_LINE_STRIP);
    glVertex3d(0,0,0); glVertex3d(1,0,0);
    glEnd();
    
    glColor3f(0,1,0);
    glBegin(GL_LINE_STRIP);
    glVertex3d(0,0,0); glVertex3d(0,1,0);
    glEnd();
    
    glColor3f(0,0,1);
    glBegin(GL_LINE_STRIP);
    glVertex3d(0,0,0); glVertex3d(0,0,1);
    glEnd();
    
    glColor3f(0,0,0);
    glBegin(GL_LINE_LOOP);
    glVertex3d(-1,-1,0); glVertex3d(1,-1,0); glVertex3d(1,1,0); glVertex3d(-1,1,0);
    glEnd();
    */
   
    int char_w=6,char_h=12;
    float offset=0.01;
    int num_char=4;
   
    if(ca->Box){
	// tick 
	float cta0;
	float r1=1.05;//tick width
	float r2=1.2;
	float r3=1.4;
	int signx,signy;
	cta0=ca->cta;cta0=fmod(ca->cta,360);
	if((  0<=cta0)&&(cta0< 90)){signx= 1;signy= 1;}
	if(( 90<=cta0)&&(cta0<190)){signx=-1;signy= 1;}
	if((180<=cta0)&&(cta0<270)){signx=-1;signy=-1;}
	if((270<=cta0)&&(cta0<360)){signx= 1;signy=-1;}

	glColor3f(0,0,0);

	// axes //
	// x
	glBegin(GL_LINE_STRIP);
	glVertex3d(-1,signy,-1); 
	glVertex3d( 1,signy,-1);
	glEnd();
	// y
	glBegin(GL_LINE_STRIP);
	glVertex3d(signx,-1,-1); 
	glVertex3d(signx, 1,-1);
	glEnd();
	// z
	glBegin(GL_LINE_STRIP);
	glVertex3d(signy,-signx,-1); 
	glVertex3d(signy,-signx, 1);
	glEnd();

	// Tick //
	//x
	for(int i=0;i<ca->XTick.size();++i){
	    glBegin(GL_LINE_STRIP);
	    glVertex3d( ct3x(ca->XTick[i]),signy   ,-1 );
	    glVertex3d( ct3x(ca->XTick[i]),signy*r1,-1 );
	    glEnd();
	}
	// y
	for(int i=0;i<ca->YTick.size();++i){
	    glBegin(GL_LINE_STRIP);
	    glVertex3d( signx   ,ct3y(ca->YTick[i]),-1 );
	    glVertex3d( signx*r1,ct3y(ca->YTick[i]),-1 );
	    glEnd();
	}
	// z
	for(int i=0;i<ca->YTick.size();++i){
	    glBegin(GL_LINE_STRIP);
	    glVertex3d( signy   ,-signx   ,ct3z(ca->ZTick[i]) );
	    glVertex3d( signy*r1,-signx,ct3z(ca->ZTick[i]) );
	    glEnd();
	}
	// Tick Label //
	if(ca->TickLabel){
	    //x
	    for(int i=0;i<ca->XTick.size();++i){
		sprintf(ctmp,"%4.1f",ca->XTick[i]);
		ptext3c( ct3x(ca->XTick[i]),signy*r2 ,-1,ctmp );
	    }
	    // y	
	    for(int i=0;i<ca->YTick.size();++i){
		sprintf(ctmp,"%4.1f",ca->YTick[i]);
		ptext3c( signx*r2,ct3y(ca->YTick[i]),-1,ctmp );
	    }
	    // z
	    for(int i=0;i<ca->ZTick.size();++i){
		sprintf(ctmp,"%4.1f",ca->ZTick[i]);
		ptext3c( signy*r2,-signx,ct3z(ca->ZTick[i]),ctmp );
	    }
	}
	// xyz Label //
	ptext3c(0,signy*r3,-1,"x");
	ptext3c(signx*r3,0,-1,"y");
	ptext3c(signy*r3,-signx,0,"z");

    }//box
}

/// colorbar
void MatPlot::display_axes_colorbar(){
    char ctmp[100];
	float l,b,w,h;//left,bottom,width,height
	float r=0.01;

	l=ca->Position[0];
	b=ca->Position[1];
	w=ca->Position[2];
	h=ca->Position[3];
	    
	// Viewport Figure (VpF) for drawing axes
	glViewport(0,0, (int)(window_w),(int)(window_h));
	glLoadIdentity();
	gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
	
	glDisable(GL_LINE_STIPPLE); 
	gl2psDisable(GL2PS_LINE_STIPPLE);
	
	if(ca->Box){
	    // box 	    
	    glLineWidth(ca->LineWidth);
	    gl2psLineWidth(ca->LineWidth);
	    glColor3f(0,0,0);
	    glBegin(GL_LINE_LOOP);	    
	    glVertex2d(l,  b);
	    glVertex2d(l+w,b);
	    glVertex2d(l+w,b+h);
	    glVertex2d(l,  b+h);
	    glEnd();	    
	    
	    // z tick 
	    for(int i=0;i<ca->ZTick.size();++i){
		glBegin(GL_LINE_STRIP);
		glVertex2d( l+w,       cty2(ca->ZTick[i]) );
		glVertex2d( l+w+0.01,  cty2(ca->ZTick[i]) );		
		glEnd();
	    }
	    // z tick number
	    for(int i=0;i<ca->ZTick.size();++i){
		sprintf(ctmp,"%4.1f",ca->ZTick[i]);
		ptext( l+w+0.01,cty2(ca->ZTick[i]),ctmp );
	    }
	}//Box

	vector<float> rgb;
	int n=cmap.size();
	for(int i=0;i<n;++i){
	    rgb=ca->ColorMap[i];
	    glColor3f(rgb[0],rgb[1],rgb[2]);
	    
	    glBegin(GL_QUADS); 
	    glVertex2d(l  ,b+h*i/n);
	    glVertex2d(l+w,b+h*i/n);
	    glVertex2d(l+w,b+h*(i+1)/n);
	    glVertex2d(l  ,b+h*(i+1)/n);
	    glEnd();	    
	}

    }
/// events (mouse, motion)
void MatPlot::Axes_mouse(int button, int state, int x, int y ){

	float X,Y;
	double rx,ry,mx,my;//mouse
	X=(float)          x /window_w;
	Y=(float)(window_h-y)/window_h;
	float l,b,w,h;

	if ( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) {// Left Click
	    xButtonDown=x;
	    yButtonDown=y;
	
	    //cout <<"window w h"<<window_w << " "<<window_h<<endl;
	    if(is_debug2){cout <<"Left click: "<< X <<" "<< Y <<endl; }
	    
	    // mouse capture axes //
	    if(iAxesSelected>=0){

		ca=&vAxes[iAxesSelected];
		if(ca->Visible){
		    l=ca->Position[0];
		    b=ca->Position[1];
		    w=ca->Position[2];
		    h=ca->Position[3];
		    
		    if( (l<=X)&&(X<=l+w)&&(b<=Y)&&(Y<=b+h)&&(ca->Mouse==1) ){
			rx=(X-l)/w;
			ry=(Y-b)/h;
			mx=rx*(ca->XLim[1]-ca->XLim[0])+ca->XLim[0];
			my=ry*(ca->YLim[1]-ca->YLim[0])+ca->YLim[0];
			ca->XMouse = mx;
			ca->YMouse = my;
			if(is_debug2){cout <<"mouse capture: "<< mx <<" "<< my <<endl; }
		    }
		}		
	    }
	    
	    // axes select //
	    iAxesSelected=-1;
	    
	    for(int i=0;i<vAxes.size();++i){
		ca=&vAxes[i];
		ca->Selected=0;
		if(ca->Visible){
		    l=ca->Position[0];
		    b=ca->Position[1];
		    w=ca->Position[2];
		    h=ca->Position[3];
		    
		    if( (l<=X)&&(X<=l+w)&&(b<=Y)&&(Y<=b+h) ){
			iAxesSelected=i;
			ca->Selected=1;
			if(ca->View==1){//3D			
			    ctaButtonDown = ca->cta;
			    phiButtonDown = ca->phi;
			}
			//cout << "axes "<< i << " is selected "<< vAxes[i].selected()  << endl;
			//cout <<"(cta,phi) = "<< ctaButtonDown <<" "<< phiButtonDown << endl;		
		    }
		}
	    }

	    if(is_debug2){ cout <<"selected axes"<< iAxesSelected<<endl; }

	}// left click

    }
void MatPlot::Axes_motion(int x, int y ){
	float phi,cta;
	for(int i=1;i<vAxes.size();++i){
	    ca=&vAxes[i];
	    if((ca->Selected)&&(ca->View==1)){
		cta = ctaButtonDown - (float) (x - xButtonDown)*1;
		phi = phiButtonDown + (float) (y - yButtonDown)*1;
		if(phi>= 90){phi= 90;}
		if(phi<=-90){phi=-90;}
		if(cta>360){cta+=-360;}
		if(cta<  0){cta+= 360;}
		
		ca->phi = phi;
		ca->cta = cta;
		//cout <<"( phi,cta ) = ( "<< vAxes[i].phi <<","<< vAxes[i].cta <<" )"<<endl;
	    }
	}

	//float phi,cta;
	//phi = phiButtonDown +(float) (y - yButtonDown)*0.01;
	//cta = ctaButtonDown -(float) (x - xButtonDown)*0.01;
	//cout <<"( dx,dy ) = ( "<< x-xButtonDown <<","<<y-yButtonDown <<" )"<<endl;
	//cout <<"( phi,cta ) = ( "<< phi <<","<< cta <<" )"<<endl;
	
	//glutPostRedisplay(); 
    }


/// subplot
int MatPlot::subplot(int m,int n,int p){	
    int h=axes();
    int ix,iy;
    ix=(p-1)%n;
    iy=(m-1)-(p-1)/n;
    ca->Position[0]=(ix+0.13)/n;
    ca->Position[1]=(iy+0.11)/m;
    ca->Position[2]=0.775/n;
    ca->Position[3]=0.815/m;

    ca->Viewport3d[0]=1.0*ix/n;
    ca->Viewport3d[1]=1.0*iy/m;
    ca->Viewport3d[2]=1.0/n;
    ca->Viewport3d[3]=1.0/m;

    return h;
}
/// colorbar
int MatPlot::colorbar(){
    float l,b,w,h;
    l=ca->Position[0];
    b=ca->Position[1];
    w=ca->Position[2];
    h=ca->Position[3];
    float zmin,zmax;
    zmin=ca->ZLim[0];
    zmax=ca->ZLim[1];
    
    // TODO use in 3D
    
    int hh=axes();
    ca->ColorMap=cmap;
    ca->View=2;
    ca->Position[0]=l+w+w*0.01;
    ca->Position[1]=b;
    ca->Position[2]=w*0.05;
    ca->Position[3]=h;
    ca->ZLim[0]=zmin;
    ca->ZLim[1]=zmax;
    ca->YLim[0]=zmin;
    ca->YLim[1]=zmax;
    return hh;
}
/// axis
void MatPlot::axis(double xMin,double xMax,double yMin,double yMax){	
    if(xMin!=xMax){
	ca->XLim[0]=xMin;
	ca->XLim[1]=xMax;
	ca->XLimMode=1;
    }
    if(yMin!=yMax){
	ca->YLim[0]=yMin;
	ca->YLim[1]=yMax;		
	ca->YLimMode=1;
    }
    ca->View=0;//2D
}

void MatPlot::axis(double xMin,double xMax,double yMin,double yMax,double zMin,double zMax){
    ca->XLim[0]=xMin; ca->XLim[1]=xMax;
    ca->YLim[0]=yMin; ca->YLim[1]=yMax;
    ca->ZLim[0]=zMin; ca->ZLim[1]=zMax;
    ca->XLimMode=1;
    ca->YLimMode=1;
    ca->ZLimMode=1;
    ca->View=1;//3D
}
void MatPlot::axis(string s){
    if( s=="on"  ){ca->Box=1;}
    if( s=="off" ){ca->Box=0;}    
}
void MatPlot::axis(int s){
    if(s){ca->Box=1;}
    else{ ca->Box=0;}    
}
void MatPlot::grid(string s){
    if( s=="on" ){ca->XGrid=1; ca->YGrid=1; ca->ZGrid=1;}
    if( s=="off"){ca->XGrid=0; ca->YGrid=0; ca->ZGrid=0;}    
}
void MatPlot::grid(int s){
    if(s){ca->XGrid=1; ca->YGrid=1; ca->ZGrid=1;}
    else{ca->XGrid=0; ca->YGrid=0; ca->ZGrid=0;}    
}
void MatPlot::ticklabel(int s){
    if(s){ca->TickLabel=1;}
    else{ ca->TickLabel=0;}    
}
void MatPlot::title(string s){
    ca->Title=s;
}
void MatPlot::xlabel(string s){
    ca->XLabel=s;
}
void MatPlot::ylabel(string s){
    ca->YLabel=s;
}
void MatPlot::mouse_capture(double *xmouse,double *ymouse){
    ca->Mouse=1;
    //ca->XMouse = xmouse;
    //ca->YMouse = ymouse;
}
/// Fmax Fmin
double Fmax(dvec x){
    double max=-1e99;
    for(int i=0;i<x.size();++i){
	if(x[i]>max){max=x[i];}
    }
    return max;
}
double Fmin(dvec x){
    double min=1e99;
    for(int i=0;i<x.size();++i){
	if(x[i]<min){min=x[i];}
    }
    return min;
}
double Fmax(dmat x){
    double max=-1e99;
    for(int i=0;i<x.size();++i){
	for(int j=0;j<x[i].size();++j){
	    if(x[i][j] > max){max=x[i][j];}
	}
    }
    return max;
}
double Fmin(dmat x){
    double min=1e99;
    for(int i=0;i<x.size();++i){
	for(int j=0;j<x[i].size();++j){
	    if(x[i][j] < min){min=x[i][j];}
	}
    }
    return min;
}
// Line ///
/// line
int MatPlot::line(){
    int h=iLine*100 + tLine; hObj=h;
    if(is_debug1){printf("mode: %d handle: %4d Line\n",mode,h);}
    
    if(mode==0){	    
	ca->add_child(h);	    	    
	vLine.push_back(Line(h));
	//al.Parent=gca();
    }
    if(mode==1){
	if(iLine<vLine.size()){ vLine[iLine].reset(); }
    }
    if(iLine<vLine.size()){cl=&vLine[iLine];}
    iLine++;    
    return h;
}
void MatPlot::line_config(){
    int n;
    double t; 
    n=cl->XData.size();
    if(ca->XScale==0){//linear
	for(int i=0;i<n;++i){
	    t=cl->XData[i];	    
	    if(ca->xmin>t){ca->xmin=t;}
	    if(ca->xmax<t){ca->xmax=t;}
	}
    }
    if(ca->XScale==1){//log
	for(int i=0;i<n;++i){
	    t=cl->XData[i];	    
	    if((ca->xmin>t)&&(t>0)){ca->xmin=t;}
	    if(ca->xmax<t){ca->xmax=t;}
	}
    }
    double y;
    n=cl->YData.size();
    for(int i=0;i<n;++i){
	y=cl->YData[i];	    
	if(ca->ymin>y){ca->ymin=y;}
	if(ca->ymax<y){ca->ymax=y;}
    }
    double z;
    n=cl->ZData.size();
    for(int i=0;i<n;++i){
	z=cl->ZData[i];  
	if(ca->zmin>z){ca->zmin=z;}
	if(ca->zmax<z){ca->zmax=z;}
    }
}
int MatPlot::line(dvec x,dvec y){
    int h=line();
    if(cfr->Visible){
	cl->XData=x;
	cl->YData=y;
	line_config();
    }
    return h;
}
int MatPlot::line(dvec x,dvec y,dvec z){
    int h=line();
    if(cfr->Visible){
	cl->XData=x;
	cl->YData=y; 
	cl->ZData=z;
	line_config();
    }
    return h;
}
/// vertex
int MatPlot::begin(){ return line(); }
void MatPlot::end(){}
void  MatPlot::vertex(double x,double y){
    if(cfr->Visible){
	if(ca->xmin>x){ca->xmin=x;}
	if(ca->xmax<x){ca->xmax=x;}
	if(ca->ymin>y){ca->ymin=y;}
	if(ca->ymax<y){ca->ymax=y;}	    
	cl->XData.push_back(x);
	cl->YData.push_back(y);
    }
}
/// plot, semilogx, semilogy, loglog
int MatPlot::plot(dvec y){
    int n=y.size();
    dvec x;
    x.resize(n);    
    for(int i=0;i<n;++i){x[i]=1.0*i/(n-1);}
    return line(x,y);
}
int MatPlot::plot(dvec x,dvec y){	
    return line(x,y);
}
int MatPlot::plot(valarray<double> x,valarray<double> y){	
    dvec xx,yy;
    for(int i=0;i<x.size();++i){xx.push_back(x[i]);}
    for(int i=0;i<y.size();++i){yy.push_back(y[i]);}    
    return line(xx,yy);
}
int MatPlot::semilogx(dvec x,dvec y){
    ca->XScale=1;
    int h=line();
    if(cfr->Visible){
	cl->XData=x;
	cl->YData=y; 
	line_config();
    }
    return h;
}
int MatPlot::semilogy(dvec x,dvec y){
    ca->YScale=1;    
    int h=line();
    if(cfr->Visible){
	cl->XData=x;
	cl->YData=y; 
	line_config();
    }
    return h;
}
int MatPlot::loglog(dvec x,dvec y){
    ca->XScale=1;
    ca->YScale=1;
    int h=line();
    if(cfr->Visible){
	cl->XData=x;
	cl->YData=y; 
	line_config();
    }
    return h;
}
/// errorbar
void MatPlot::vertex(double x,double y,double ep,double em){//for errorbar
	if(ca->xmin>x){ca->xmin=x;}
	if(ca->xmax<x){ca->xmax=x;}
	if(ca->ymin>y+ep){ca->ymin=y+ep;}
	if(ca->ymax<y-em){ca->ymax=y-em;}
	cl->XData.push_back(x);
	cl->YData.push_back(y);
	cl->YPData.push_back(ep);
	cl->YMData.push_back(em);	
    }
int MatPlot::errorbar(dvec x,dvec y,dvec e){	
	begin();	
	for(int i=0;i<x.size();++i){ vertex(x[i],y[i],e[i],e[i]); }
	end();
	cl->Errorbar=1;
	return 0;
    }
int MatPlot::errorbar(dvec x,dvec y,dvec ep,dvec em){
	begin();
	for(int i=0;i<x.size();++i){ vertex(x[i],y[i],ep[i],em[i]); }
	end();
	cl->Errorbar=1;
	return 0;
    }
/// 3D line
void MatPlot::vertex(double x,double y,double z){
    if(cfr->Visible){
	if(ca->xmin>x){ca->xmin=x;}
	if(ca->xmax<x){ca->xmax=x;}
	if(ca->ymin>y){ca->ymin=y;}
	if(ca->ymax<y){ca->ymax=y;}
	if(ca->zmin>z){ca->zmin=z;}
	if(ca->zmax<z){ca->zmax=z;}
	cl->XData.push_back(x);
	cl->YData.push_back(y);
	cl->ZData.push_back(z);
    }
    //if(cl->LineStyle){//Line
    //glVertex3d(ct3x(x),ct3y(y),ct3z(z));
    //}
}    
int MatPlot::plot3(dvec x,dvec y,dvec z){
    ca->View=1;
    begin();
    for(int i=0;i<x.size();++i){ vertex(x[i],y[i],z[i]); }
    end();
    return 0;
}

/// display_line
void MatPlot::display_line(){

    if(is_debug1){printf("mode: %d handle: %4d Line \n",
			 mode,cl->id);}
    
    float xx,yy;// transformed coordination 
    float r;//marker size
    float rx,ry;
    vector<float> rgb=ColorSpec2RGB(cl->Color);
    glColor3f(rgb[0],rgb[1],rgb[2]);
    
    glLineWidth(cl->LineWidth);
    glPointSize(cl->LineWidth);
    gl2psLineWidth(cl->LineWidth);
    gl2psPointSize(cl->LineWidth);
    // 2D //
    if(ca->View==0){
	if(cl->LineStyle !="none" ){// Line //
	    
	    if(cl->LineStyle=="-"){
		glDisable(GL_LINE_STIPPLE); 
		gl2psDisable(GL2PS_LINE_STIPPLE);
	    }
	    if(cl->LineStyle=="- -"){
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1, 0xF0F0);
		gl2psEnable(GL2PS_LINE_STIPPLE);
	    }
	    if(cl->LineStyle==":"){
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1, 0xCCCC);
		gl2psEnable(GL2PS_LINE_STIPPLE);
	    }
	    if(cl->LineStyle=="-."){
		glEnable(GL_LINE_STIPPLE);	   
		glLineStipple(1, 0x087F);
		gl2psEnable(GL2PS_LINE_STIPPLE);
	    }
	    if(cl->XData.size()){
		glBegin(GL_LINE_STRIP);
		for(int i=0;i<cl->XData.size();++i){
		    //printf("i:%d %f %f\n",i,xx,yy);
		    xx=ctx(cl->XData[i]);
		    yy=cty(cl->YData[i]);
		    glVertex2d(xx,yy);
		}
		glEnd();
	    }
	}
	 
	if(cl->Marker != "none"){// Marker //

	    r=cl->MarkerSize/500.0;
	    rx=cl->MarkerSize/window_w;
	    ry=cl->MarkerSize/window_h;


	    glDisable(GL_LINE_STIPPLE); 
	    gl2psDisable(GL2PS_LINE_STIPPLE);

	    for(int i=0;i<cl->XData.size();++i){
		xx=ctx(cl->XData[i]);
		yy=cty(cl->YData[i]);		
		
		if(cl->Marker=="."){//.
		    glPointSize(cl->LineWidth);
		    glBegin(GL_POINTS);
		    glVertex2d(xx,yy);
		    glEnd();
		}		
		if(cl->Marker=="+"){//+
		    glBegin(GL_LINE_STRIP); 
		    glVertex2d(xx-rx,yy);
		    glVertex2d(xx+rx,yy);
		    glEnd();
		    glBegin(GL_LINE_STRIP); 
		    glVertex2d(xx,yy-ry);
		    glVertex2d(xx,yy+ry);
		    glEnd();
		}
		if(cl->Marker=="x"){//x
		    glBegin(GL_LINE_STRIP);
		    glVertex2d(xx-rx,yy-ry);
		    glVertex2d(xx+rx,yy+ry);
		    glEnd();
		    glBegin(GL_LINE_STRIP);
		    glVertex2d(xx+rx,yy-ry);
		    glVertex2d(xx-rx,yy+ry);
		    glEnd();
		}
		if(cl->Marker=="d"){//d diamond
		    glBegin(GL_LINE_LOOP);
		    glVertex2d(xx,   yy+ry);
		    glVertex2d(xx+rx,yy);
		    glVertex2d(xx,   yy-ry);
		    glVertex2d(xx-rx,yy);
		    glEnd();
		}
		if(cl->Marker=="^"){//^
		    glBegin(GL_LINE_LOOP);
		    glVertex2d(xx,   yy+ry);
		    glVertex2d(xx+rx,yy-ry);
		    glVertex2d(xx-rx,yy-ry);
		    glEnd();
		}
		if(cl->Marker=="v"){//v
		    glBegin(GL_LINE_LOOP);
		    glVertex2d(xx,   yy-ry);
		    glVertex2d(xx+rx,yy+ry);
		    glVertex2d(xx-rx,yy+ry);
		    glEnd();
		}
		if(cl->Marker=="o"){//o
		    glBegin(GL_LINE_LOOP);
		    for(int i=0;i<20;i++){
			glVertex2d(xx+rx*cos(2*PI*(double)i/(double)(20)),
				   yy+ry*sin(2*PI*(double)i/(double)(20)));
		    }
		    glEnd();
		}
		if(cl->Marker=="*"){//*
		    glBegin(GL_LINE_STRIP); 
		    glVertex2d(xx-rx,yy);
		    glVertex2d(xx+rx,yy);
		    glEnd();
		    glBegin(GL_LINE_STRIP); 
		    glVertex2d(xx,yy-ry);
		    glVertex2d(xx,yy+ry);
		    glEnd();
		    glBegin(GL_LINE_STRIP); 
		    glVertex2d(xx-rx,yy-ry);
		    glVertex2d(xx+rx,yy+ry);
		    glEnd();
		    glBegin(GL_LINE_STRIP); 
		    glVertex2d(xx+rx,yy-ry);
		    glVertex2d(xx-rx,yy+ry);
		    glEnd();

		}
		if(cl->Marker=="s"){//s :squire
		    glBegin(GL_LINE_LOOP);
		    glVertex2d(xx-rx,yy-ry);
		    glVertex2d(xx-rx,yy+ry);
		    glVertex2d(xx+rx,yy+ry);
		    glVertex2d(xx+rx,yy-ry);
		    glEnd();
		}
		//TODO use switch
	    }//i
	}// Marker

	
	if(cl->Errorbar){// Errorbar //
	    float xx,yy,yyp,yym;// transformed coordination 
	
	    glDisable(GL_LINE_STIPPLE); 
	    gl2psDisable(GL2PS_LINE_STIPPLE);
	    //r=cl->MarkerSize/500;

	    for(int i=0;i<cl->XData.size();++i){
		xx=ctx(cl->XData[i]);
		yy=cty(cl->YData[i]);
		yyp=cty(cl->YData[i] + cl->YPData[i]);
		yym=cty(cl->YData[i] - cl->YMData[i]);

		glBegin(GL_LINE_STRIP);
		glVertex2d(xx,yyp);
		glVertex2d(xx,yym);
		glEnd();
		glBegin(GL_LINE_STRIP);
		glVertex2d(xx-rx,yy);
		glVertex2d(xx+rx,yy);
		glEnd();
		glBegin(GL_LINE_STRIP);
		glVertex2d(xx-rx,yyp);
		glVertex2d(xx+rx,yyp);
		glEnd();
		glBegin(GL_LINE_STRIP);
		glVertex2d(xx-rx,yym);
		glVertex2d(xx+rx,yym);
		glEnd();
	    }//i
	}//Errorbar
	//TODO:selection of error bar type 
	}//2D

	// 3D //
	if(ca->View==1){
	    glBegin(GL_LINE_STRIP);
	    for(int i=0;i<cl->XData.size();++i){
		glVertex3d(ct3x(cl->XData[i]),
			   ct3y(cl->YData[i]),
			   ct3z(cl->ZData[i]));
	    }
	    glEnd();
	}
    }

// Surface ///
int MatPlot::surface(){
    int h=iSurface*100 + tSurface; hObj=h;
    if(is_debug1){printf("mode: %d handle: %4d Surface\n",mode,h);}

    if(mode==0){
	ca->add_child(h);	   
	vSurface.push_back(Surface(h));
	//as.Parent=gca();
    }
    if(mode==1){
	
    }
    if(iSurface<vSurface.size()){cs=&vSurface[iSurface];}
    iSurface++;
    return h;    
}

void MatPlot::surface_config(){

    // check data size    
    int nzi,nzj;
    nzi=cs->ZData.size();
    if(nzi){nzj=cs->ZData[0].size(); }
    
    int nci,ncj;
    nci=cs->CDataIndex.size();
    if(nci){ncj=cs->CDataIndex[0].size();}

    // generate x and y data
    int nx=0,ny=0;
    if(nzi){ny=nzi; nx=nzj;}
    if(nci){ny=nci; nx=ncj;}
    
    //printf("%s %s:%d: %d %d %d %d \n", __func__, __FILE__, __LINE__,nzi,nci,nx,ny);
    if(cs->XData.size()==0){
	cs->XData.resize(1); 
	cs->XData[0]=linspace(1.0,(double)nx,nx);
    }
    if(cs->YData.size()==0){
	cs->YData.resize(1); 
	cs->YData[0]=linspace(1.0,(double)ny,ny);
    }
    
    // config data range
    ca->xmax=fmax(Fmax(cs->XData),ca->xmax);
    ca->xmin=fmin(Fmin(cs->XData),ca->xmin);
    ca->ymax=fmax(Fmax(cs->YData),ca->ymax);
    ca->ymin=fmin(Fmin(cs->YData),ca->ymin);
    ca->zmax=fmax(Fmax(cs->ZData),ca->zmax );
    ca->zmin=fmin(Fmin(cs->ZData),ca->zmin );

    // set CLim 
    //note: first called surface effects CLim
    if(ca->CLim[0]==ca->CLim[1]){
	ca->CLim[0]=Fmin(cs->CDataIndex);
	ca->CLim[1]=Fmax(cs->CDataIndex);
    }

    // CData !!!
    if( (cs->CData.size()==0) && (cs->CDataIndex.size()) ){
	vector<float> rgb;
	//tcmat cdata(ny,nx);
	tcmat cdata(ny);
	
	for(int i=0;i<ny;++i){
	    cdata[i].resize(nx);
	    for(int j=0;j<nx;++j){
		rgb=map2color(cs->CDataIndex[i][j],ca->CLim[0],ca->CLim[1]);
		cdata[i][j]=rgb;
	    }
	}
	cs->CData=cdata;
    }

    // contour plot
    if(cs->V.size()==0){
	if(cs->NContour<1){
	    cs->NContour=10;		    
	}
	cs->V=linspace(Fmin(cs->ZData),Fmax(cs->ZData),cs->NContour);	
    }
}

/// create surface
int MatPlot::surface(dmat Z){    
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->ZData=Z;
    cs->CDataIndex=Z;
    cs->CData.clear();
    surface_config();
    return h;
}
int MatPlot::surface(dmat Z,dmat C){    
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->ZData=Z;
    cs->CDataIndex=C;
    cs->CData.clear();
    surface_config();
    return h;
}
int MatPlot::surface(dmat Z,tcmat C){    
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->ZData=Z;
    cs->CDataIndex.clear();
    cs->CData=C;
    surface_config();
    return h;
}
int MatPlot::surface(dvec x,dvec y,dmat Z){    
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;
    cs->ZData=Z;
    cs->CDataIndex=Z;
    cs->CData.clear();
    surface_config();
    return h;
}
int MatPlot::surface(dvec x,dvec y,dmat Z,dmat C){    
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;
    cs->ZData=Z;
    cs->CDataIndex=C;
    cs->CData.clear();
    surface_config();
    return h;
}
int MatPlot::surface(dvec x,dvec y,dmat Z,tcmat C){    
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;
    cs->ZData=Z;
    cs->CDataIndex.clear();
    cs->CData=C;
    surface_config();
    return h;
}

int MatPlot::surface(dmat X,dmat Y,dmat Z){
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData=X;
    cs->YData=Y;
    cs->ZData=Z;
    cs->CDataIndex=Z;
    cs->CData.clear();
    surface_config();
    return h;
}
int MatPlot::surface(dmat X,dmat Y,dmat Z,dmat C){
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData=X;
    cs->YData=Y;
    cs->ZData=Z;
    cs->CDataIndex=C;
    cs->CData.clear();

    surface_config();
    return h;
}
int MatPlot::surface(dmat X,dmat Y,dmat Z,tcmat C){
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData=X;
    cs->YData=Y;
    cs->ZData=Z;
    cs->CDataIndex.clear();
    cs->CData=C;
    
    surface_config();
    return h;
}
/// surf
int MatPlot::surf(dvec x, dvec y, dmat Z){
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;
    cs->ZData=Z;
    cs->CDataIndex=Z;
    cs->CData.clear();    
    cs->EdgeColor="k";
    cs->FaceColor="flat";
    surface_config();
    return h;
}
/// create pcolor
int MatPlot::pcolor(dmat C){
    
    int h; h=surface();

    cs->type=0;
    cs->XData.clear();
    cs->YData.clear();
    cs->ZData.clear();
    cs->CDataIndex=C;
    cs->CData.clear();

    surface_config();

    return h;
}
int MatPlot::pcolor(tcmat C){
    int h=surface();
    cs->type=0;
    cs->XData.clear();
    cs->YData.clear();
    cs->ZData.clear();
    cs->CDataIndex.clear();
    cs->CData=C;    
    surface_config();
    return h;
}
int MatPlot::pcolor(dvec x, dvec y, dmat C){
    int h=surface();
    cs->type=0;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;
    cs->ZData.clear();
    cs->CDataIndex=C;
    cs->CData.clear();
    surface_config();
    return h;
}
int MatPlot::pcolor(dvec x, dvec y, tcmat C){
    int h=surface();
    cs->type=0;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;
    cs->ZData.clear();
    cs->CDataIndex.clear();
    cs->CData=C;
    surface_config();
    return h;
}
int MatPlot::pcolor(dmat X,dmat Y,dmat C){
    int h=surface();
    cs->type=0; 
    cs->XData=X;
    cs->YData=Y;	
    cs->ZData.clear();
    cs->CDataIndex=C;
    cs->CData.clear();    
    surface_config();
    return h;
}
int MatPlot::pcolor(dmat X,dmat Y,tcmat C){
    int h=surface();
    cs->type=0; 
    cs->XData=X;
    cs->YData=Y;	
    cs->ZData.clear();
    cs->CDataIndex.clear();
    cs->CData=C;    
    surface_config();
    return h;
}

/// mesh
int MatPlot::mesh(dvec x, dvec y, dmat Z){
    int h=surface();
    ca->View=1;
    cs->type=1;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;
    cs->ZData=Z;
    cs->CDataIndex=Z;
    cs->CData.clear();    
    cs->EdgeColor="k";
    cs->FaceColor="w";
    surface_config();
    return h;
}

/// contour
int MatPlot::contour(dmat Z){
    int h=surface();
    cs->type=3;
    cs->XData.clear();
    cs->YData.clear();
    cs->ZData=Z;
    cs->NContour=0;
    cs->V.clear();      

    surface_config();
    return h;
}
int MatPlot::contour(dmat Z,int n){
    int h=surface();
    cs->type=3;
    cs->XData.clear();
    cs->YData.clear();
    cs->ZData=Z;
    cs->NContour=n; 
    cs->V.clear();
   
    surface_config();
    return h;
}
int MatPlot::contour(dmat Z, dvec v){
    int h=surface();
    cs->type=3;
    cs->XData.clear();
    cs->YData.clear();
    cs->ZData=Z;
    cs->NContour=v.size();
    cs->V=v;
    
    surface_config();
    return h;
}
int MatPlot::contour(dvec x,dvec y,dmat Z){
    int h=surface();
    cs->type=3;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;	
    cs->ZData=Z;
    cs->NContour=0;
    cs->V.clear();      

    surface_config();
    return h;
}
int MatPlot::contour(dvec x,dvec y,dmat Z,int n){
    int h=surface();
    cs->type=3;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;	
    cs->ZData=Z;
    cs->NContour=n; 
    cs->V.clear();
   
    surface_config();
    return h;
}
int MatPlot::contour(dvec x, dvec y, dmat Z, dvec v){
    int h=surface();
    cs->type=3;
    cs->XData.resize(1); cs->XData[0]=x;
    cs->YData.resize(1); cs->YData[0]=y;	
    cs->ZData=Z;
    cs->NContour=v.size();
    cs->V=v;
    
    surface_config();
    return h;
}
/// display 
void MatPlot::display_surface(){
    if(is_debug1){printf("mode: %d handle: %4d Surface type %d, Z.size %d\n",
			 mode,cs->id,cs->type,cs->ZData.size());}
    
    //if(cs->type==0){ display_pcolor(); }
    if(cs->type==0){ display_surface_2d(); }
    if(cs->type==1){ display_surface_3d(); }
    if(cs->type==2){ display_surface_3d(); }
    if(cs->type==3){ display_contour(); }
    
}
void MatPlot::display_surface_2d(){
    int nxi,nxj,nyi,nyj,nzi,nzj;
    vector<float> rgb;
    nzi=cs->ZData.size(); if(nzi){nzj=cs->ZData[0].size();}    
    nxi=cs->XData.size(); if(nxi){nxj=cs->XData[0].size();}
    nyi=cs->YData.size(); if(nyi){nyj=cs->YData[0].size();}

    //printf("%s %s:%d  (%d %d) (%d %d) (%d %d) \n", __func__, __FILE__, __LINE__,nxi,nxj,nyi,nyj,nzi,nzj);

    // (Z) // do not use
    if(nxi==0){
      //printf("%s %s:%d\n", __func__, __FILE__, __LINE__);
	// Edge
	if(cs->EdgeColor != "none"){	    
	    glLineWidth(cs->LineWidth);
	    rgb=ColorSpec2RGB(cs->EdgeColor);
	    glColor3d(rgb[0],rgb[1],rgb[2]);
	    
	    for(int i=0;i<nzi;++i){		
		glBegin(GL_LINE_STRIP); //TODO add more style		    
		for(int j=0;j<nzj;++j){ 
		    glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j)/(nzj-1)),
			       cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i)/(nzi-1)));
		}
		glEnd();
	    }
	    /*
	    for(int j=0;j<nzj;++j){ 	    		
		glBegin(GL_LINE_STRIP); 		    
		for(int i=0;i<nzi;++i){
		    glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j)/(nzj-1)),
			       cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i)/(nzi-1)));
		}
		glEnd();
	    }
	    */
	}
	// Face 
	if(cs->FaceColor != "none"){
	    for(int i=0;i<nzi-1;++i){
		for(int j=0;j<nzj-1;++j){
		    
		    rgb=ColorSpec2RGB(cs->FaceColor);		
		    if(cs->FaceColor=="flat"){ rgb=cs->CData[i][j]; }
		    glColor3f(rgb[0],rgb[1],rgb[2]);

		    glBegin(GL_QUADS);
		    glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j  )/(float)(nzj-1)),
			       cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i  )/(float)(nzi-1)));
		    glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j+1)/(float)(nzj-1)),
			       cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i  )/(float)(nzi-1)));
		    glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j+1)/(float)(nzj-1)),
			       cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i+1)/(float)(nzi-1)));
		    glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j  )/(float)(nzj-1)),
			       cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i+1)/(float)(nzi-1)));
		    glEnd();
		}		
	    }
	}//Face
    }

    // (x,y,Z) //
    if(nxi==1){
	//printf("%s %s:%d  %d %d \n", __func__, __FILE__, __LINE__,nxj,nyj);
	// Edge
	if(cs->EdgeColor != "none"){	    
	    
	    glLineWidth(cs->LineWidth);
	    rgb=ColorSpec2RGB(cs->EdgeColor);
	    glColor3d(rgb[0],rgb[1],rgb[2]);
	    
	    for(int i=0;i<nyj;++i){	
		glBegin(GL_LINE_STRIP); //TODO add more style		    
		//for(int j=0;j<nxj;++j){
		    //glVertex2d(ctx(cs->XData[0][j]),cty(cs->YData[0][i]));
		//}
		glVertex2d( ctx(cs->XData[0][    0]),cty(cs->YData[0][i]) );
		glVertex2d( ctx(cs->XData[0][nxj-1]),cty(cs->YData[0][i]) );
		
		glEnd();
	    }
	    for(int j=0;j<nxj;++j){ 	    		
		glBegin(GL_LINE_STRIP); 		    
		//for(int i=0;i<nxj;++i){
		//glVertex2d(ctx(cs->XData[0][j]), cty(cs->YData[0][i]));
		//}
		glVertex2d(ctx(cs->XData[0][j]),cty(cs->YData[0][0]));
		glVertex2d(ctx(cs->XData[0][j]),cty(cs->YData[0][nyj-1]));
		glEnd();
	    }
	}  
	// Face
	if(cs->FaceColor != "none"){
	    for(int i=0;i<nyj-1;++i){	    	     	    
		for(int j=0;j<nxj-1;++j){ 
		    // color 
		    rgb=ColorSpec2RGB(cs->FaceColor);
		    if(cs->FaceColor=="flat"){ rgb=cs->CData[i][j]; }
		    glColor3f(rgb[0],rgb[1],rgb[2]);

		    glBegin(GL_QUADS);
		    glVertex2d(ctx(cs->XData[0][j]),
			       cty(cs->YData[0][i]) ); 
		    glVertex2d(ctx(cs->XData[0][j]),
			       cty(cs->YData[0][i+1]) ); 
		    glVertex2d(ctx(cs->XData[0][j+1]),
			       cty(cs->YData[0][i+1]) );
		    glVertex2d(ctx(cs->XData[0][j+1]),
			       cty(cs->YData[0][i]) );		
		    glEnd();
		}
	    }
	}
    }//nxi==1
    
    // (X,Y,C) //
    if(nxi>1){
	// Edge
	//printf("%s %s:%d\n", __func__, __FILE__, __LINE__);
	if(cs->EdgeColor != "none"){
	    glLineWidth(cs->LineWidth);
	    rgb=ColorSpec2RGB(cs->EdgeColor);
	    glColor3d(rgb[0],rgb[1],rgb[2]);
	    
	    for(int i=0;i<nxi;++i){		
		glBegin(GL_LINE_STRIP); //TODO add more style		    
		for(int j=0;j<nxj;++j){ 
		    glVertex2d(ctx(cs->XData[i][j]),
			       cty(cs->YData[i][j]));
		}
		glEnd();
	    }
	    for(int j=0;j<nxi;++j){ 	    		
		glBegin(GL_LINE_STRIP); 		    
		for(int i=0;i<nxj;++i){
		    glVertex2d(ctx(cs->XData[i][j]),
			       cty(cs->YData[i][j]));
		}
		glEnd();
	    }
	}  
	// Face
	if(cs->FaceColor != "none"){
	    for(int i=0;i<nxi-1;++i){	    	     	    
		for(int j=0;j<nxj-1;++j){ 
		    // color 
		    rgb=ColorSpec2RGB(cs->FaceColor);		
		    if(cs->FaceColor=="flat"){rgb=cs->CData[i][j]; }
		    glColor3f(rgb[0],rgb[1],rgb[2]);
		    
		    glBegin(GL_QUADS);
		    glVertex2d(ctx(cs->XData[i  ][j]),
			       cty(cs->YData[i  ][j]) ); 
		    glVertex2d(ctx(cs->XData[i  ][j+1]),
			       cty(cs->YData[i  ][j+1]) ); 
		    glVertex2d(ctx(cs->XData[i+1][j+1]),
			       cty(cs->YData[i+1][j+1]) );
		    glVertex2d(ctx(cs->XData[i+1][j]),
			       cty(cs->YData[i+1][j]) );		
		    glEnd();
		}
	    }
	}
    }
}

/// display 3d
void MatPlot::display_surface_3d(){
    vector<float> rgb;
    int ny=cs->ZData.size();
    int nx=cs->ZData[0].size();
    
    if(cs->XData.size()==1){
	//(x,y,Z);
	//Face
	if(cs->FaceColor != "none"){
	    for(int i=0;i<ny-1;++i){	    	     	    
		for(int j=0;j<nx-1;++j){ 
		    rgb=ColorSpec2RGB(cs->FaceColor);		
		    glColor3d(rgb[0], rgb[1], rgb[2]);
		    if(cs->FaceColor=="flat"){
			rgb=cs->CData[i][j];
			glColor3f(rgb[0],rgb[1],rgb[2]);
		    }
		    glBegin(GL_TRIANGLES);
		    glVertex3d(ct3x(cs->XData[0][j]),
			       ct3y(cs->YData[0][i]),
			       ct3z(cs->ZData[i][j]) ); 
		    glVertex3d(ct3x(cs->XData[0][j]),
			       ct3y(cs->YData[0][i+1]),
			       ct3z(cs->ZData[i+1][j]) ); 
		    glVertex3d(ct3x(cs->XData[0][j+1]),
			       ct3y(cs->YData[0][i+1]),
			       ct3z(cs->ZData[i+1][j+1]) ); 
		    glEnd();
		    glBegin(GL_TRIANGLES);
		    glVertex3d(ct3x(cs->XData[0][j]),
			       ct3y(cs->YData[0][i]),
			       ct3z(cs->ZData[i][j]) ); 
		    glVertex3d(ct3x(cs->XData[0][j+1]),
			       ct3y(cs->YData[0][i]),
			       ct3z(cs->ZData[i][j+1]) ); 
		    glVertex3d(ct3x(cs->XData[0][j+1]),
			       ct3y(cs->YData[0][i+1]),
			       ct3z(cs->ZData[i+1][j+1]) ); 
		    glEnd();
		}
	    }
	}
	// Edge
	if(cs->EdgeColor != "none"){
	    glLineWidth(cs->LineWidth);
	    rgb=ColorSpec2RGB(cs->EdgeColor);
	    glColor3d(rgb[0], rgb[1], rgb[2]);
	    
	    for(int i=0;i<ny;++i){		
		glBegin(GL_LINE_STRIP); 		    
		for(int j=0;j<nx;++j){ 
		    glVertex3d(ct3x(cs->XData[0][j]),
			       ct3y(cs->YData[0][i]),
			       ct3z(cs->ZData[i][j]) ); 
		}
		glEnd();
	    }
	    for(int j=0;j<nx;++j){ 	    		
		glBegin(GL_LINE_STRIP); 		    
		for(int i=0;i<ny;++i){
		    glVertex3d(ct3x(cs->XData[0][j]),
			       ct3y(cs->YData[0][i]),
			       ct3z(cs->ZData[i][j]) ); 
		}
		glEnd();
	    }
	}    
    }// (x,y,Z)
    // (X,Y,Z) //
    if(cs->XData.size()>1){

    //Face
    if(cs->FaceColor != "none"){
	for(int i=0;i<ny-1;++i){	    	     	    
	    for(int j=0;j<nx-1;++j){ 
		rgb=ColorSpec2RGB(cs->FaceColor);		
		glColor3d(rgb[0], rgb[1], rgb[2]);
		if(cs->FaceColor=="flat"){
		    rgb=cs->CData[i][j];
		    glColor3f(rgb[0],rgb[1],rgb[2]);
		}

		glBegin(GL_TRIANGLES);
		glVertex3d(ct3x(cs->XData[i][j]),
			   ct3y(cs->YData[i][j]),
			   ct3z(cs->ZData[i][j]) ); 
		glVertex3d(ct3x(cs->XData[i][j+1]),
			   ct3y(cs->YData[i][j+1]),
			   ct3z(cs->ZData[i][j+1]) ); 
		glVertex3d(ct3x(cs->XData[i+1][j+1]),
			   ct3y(cs->YData[i+1][j+1]),
			   ct3z(cs->ZData[i+1][j+1]) ); 
		glEnd();

		glBegin(GL_TRIANGLES);
		glVertex3d(ct3x(cs->XData[i][j]),
			   ct3y(cs->YData[i][j]),
			   ct3z(cs->ZData[i][j]) ); 
		glVertex3d(ct3x(cs->XData[i+1][j]),
			   ct3y(cs->YData[i+1][j]),
			   ct3z(cs->ZData[i+1][j]) ); 
		glVertex3d(ct3x(cs->XData[i+1][j+1]),
			   ct3y(cs->YData[i+1][j+1]),
			   ct3z(cs->ZData[i+1][j+1]) ); 
		glEnd();
	    }
	}
    }
    // Edge
    if(cs->EdgeColor != "none"){
	glLineWidth(cs->LineWidth);
	rgb=ColorSpec2RGB(cs->EdgeColor);
	glColor3d(rgb[0], rgb[1], rgb[2]);
	
	for(int i=0;i<ny;++i){		
	    glBegin(GL_LINE_STRIP); 		    
	    for(int j=0;j<nx;++j){ 
		glVertex3d(ct3x(cs->XData[i][j]),
			   ct3y(cs->YData[i][j]),
			   ct3z(cs->ZData[i][j]) ); 
	    }
	    glEnd();
	}
	for(int j=0;j<nx;++j){ 	    		
	    glBegin(GL_LINE_STRIP); 		    
	    for(int i=0;i<ny;++i){
		glVertex3d(ct3x(cs->XData[i][j]),
			   ct3y(cs->YData[i][j]),
			   ct3z(cs->ZData[i][j]) ); 
	    }
	    glEnd();
	}
    }    
    }//(X,Y,Z)
}
/// dispaly contour
dmat contourc(dvec x, dvec y, dmat Z, dvec v){
    //Z(i,j), x(j),y(i)
    double x0,y0,z0;    
    int ny=Z.size();
    int nx=Z[0].size();
    dmat C;
    ContourPoint c;
    vector<ContourPoint> vc;
    deque<ContourPoint> ac;
    C.resize(2);
    //z0=0.1;
    //v.clear();v.push_back(0.2);v.push_back(0.1);
    int is_print=0;

    

    for(int iv=0;iv<v.size();++iv){
	z0=v[iv];

	// find contour points
	vc.clear();
	for(int i=0;i<ny;++i){
	    for(int j=0;j<nx;++j){
		if( (j<nx-1)&&( (Z[i][j+1]-z0)*(Z[i][j]-z0)<0 ) ){
		    x0 = x[j]+(x[j+1]-x[j])*(z0-Z[i][j])/( Z[i  ][j+1]-Z[i][j] );
		    c.xj=j;c.yi=i;c.xy=0;c.done=0;
		    c.x=x0;c.y=y[i];
		    vc.push_back(c);		
		}
		if( (i<ny-1)&&( (Z[i+1][j]-z0)*(Z[i][j]-z0)<0 ) ){
		    y0 = y[i]+(y[i+1]-y[i])*(z0-Z[i][j])/( Z[i+1][j  ]-Z[i][j] );
		    c.xj=j;c.yi=i;c.xy=1;c.done=0;
		    c.x=x[j];c.y=y0;
		    vc.push_back(c);
		}
	    }
	} 
	if(is_print){
	    printf("vc.size %d\n",vc.size());
	    for(int k=0;k<vc.size();++k){
		printf("vc: %2d : %2d %2d :%f %f\n",k,vc[k].xj,vc[k].yi,vc[k].x,vc[k].y);
	    }
	}
	// sort contour points    
	int is,is_connect=0;	
	int m,k,l,kk;
	int mode,mode_next;

	k=0;
	mode=0;
	while( mode<5 ){

	    if(mode==0){// set initial contour point
		ac.clear();
		is=0; m=0;		
		while( !is && (m<vc.size()) ){
		    if(!vc[m].done){ is=1; kk=m; }
		    m++;
		}
		if(is){		    
		    vc[kk].done=2; 
		    c=vc[kk];
		    ac.push_back(vc[kk]); 
		    mode_next=1;
		}else{
		    mode_next=5;
		}
		
	    }
	    if( (mode==1)||(mode==3) ){//search next contour point
		is=0;
		m =0;
		while( !is && (m<vc.size()) ){
		    is=0;
		    if( (!vc[m].done) || ((vc[m].done==2)&&(ac.size()>2)) ){		
			if( (c.xy==0) && (vc[m].xy==0) && (vc[m].xj==c.xj  ) && (vc[m].yi==c.yi-1) ){ is=1; }
			if( (c.xy==0) && (vc[m].xy==0) && (vc[m].xj==c.xj  ) && (vc[m].yi==c.yi+1) ){ is=2; }
			if( (c.xy==0) && (vc[m].xy==1) && (vc[m].xj==c.xj  ) && (vc[m].yi==c.yi  ) ){ is=3; }
			if( (c.xy==0) && (vc[m].xy==1) && (vc[m].xj==c.xj+1) && (vc[m].yi==c.yi  ) ){ is=4; }
			if( (c.xy==0) && (vc[m].xy==1) && (vc[m].xj==c.xj  ) && (vc[m].yi==c.yi-1) ){ is=5; }
			if( (c.xy==0) && (vc[m].xy==1) && (vc[m].xj==c.xj+1) && (vc[m].yi==c.yi-1) ){ is=6; }
			if( (c.xy==1) && (vc[m].xy==1) && (vc[m].xj==c.xj+1) && (vc[m].yi==c.yi  ) ){ is=7; }
			if( (c.xy==1) && (vc[m].xy==1) && (vc[m].xj==c.xj-1) && (vc[m].yi==c.yi  ) ){ is=8; }
			if( (c.xy==1) && (vc[m].xy==0) && (vc[m].xj==c.xj  ) && (vc[m].yi==c.yi  ) ){ is=9; }
			if( (c.xy==1) && (vc[m].xy==0) && (vc[m].xj==c.xj  ) && (vc[m].yi==c.yi+1) ){ is=10; }
			if( (c.xy==1) && (vc[m].xy==0) && (vc[m].xj==c.xj-1) && (vc[m].yi==c.yi  ) ){ is=11; }
			if( (c.xy==1) && (vc[m].xy==0) && (vc[m].xj==c.xj-1) && (vc[m].yi==c.yi+1) ){ is=12; }
		    }
		    if(is){kk=m;}
		    m++;
		}
		if(is){		    
		    vc[kk].done=1; 
		    c=vc[kk];
		}
		if(mode==1){
		    if(is){		    
			ac.push_back(vc[kk]); 
			mode_next=1;
		    }else{
			mode_next=2;
		    }
		}
		if(mode==3){
		    if(is){		    
			ac.push_front(vc[kk]); 
			mode_next=3;
		    }else{
			mode_next=4;
		    }
		}
	    }
	    if(mode==2){// set first terminal
		c=ac[0];
		mode_next=3;
	    }
	    if(mode==4){// move the accumlated points 
		if(ac.size()){
		    C[0].push_back(z0);
		    C[1].push_back(ac.size());
		    for(int i=0;i<ac.size();++i){
			C[0].push_back( ac[i].x );
			C[1].push_back( ac[i].y );
		    }
		}
		mode_next=0;
	    }

	    mode=mode_next;
	}
	
    }//iv
    // print
    if(is_print){
	for(int i=0;i<C[0].size();++i){
	    printf("C: %3d %f %f\n",i,C[0][i],C[1][i]);
	}
	cout <<"done"<<endl;
    }
    return C;
}


void MatPlot::display_contour(){
    //vector<float> rgb;
    //int ny=cs->ZData.size();
    //int nx=cs->ZData[0].size();
    vector< vector< double > > C;
    double xx,yy;

    C=contourc(cs->XData[0],cs->YData[0],cs->ZData,cs->V);

    glDisable(GL_LINE_STIPPLE);
    gl2psDisable(GL2PS_LINE_STIPPLE);

    //glLineWidth(cl->LineWidth);
    glLineWidth(2);
    gl2psLineWidth(2);

    //rgb=cs->CData[i][j];
    //glColor3f(rgb[0],rgb[1],rgb[2]);
    glColor3f(0,0,0);

    //TODO: adjust line color and properties
    int k=0,nk;
    for(int i=0;i<C[0].size();++i){
	if(k==0){
	    nk=(int)C[1][i];
	    glBegin(GL_LINE_STRIP);
	}else{
	    if(k<=nk){
		xx=ctx( C[0][i] );
		yy=cty( C[1][i] );
		glVertex2d(xx,yy);
	    }
	}
	k++;
	if(k>nk){
	    k=0;
	    glEnd();
	}
    }   
}

dmat MatPlot::peaks(int n){

    float x1=1,y1=0;
    float x2=-1,y2=1;
    float x3=-0.5,y3=-1;
    float sr1,sr2,sr3;
    float sigma=0.4;
    float a1=1,a2=0.5,a3=0.3;
    double x,y;
    //vector< vector< double > > Z(n,n);
    dmat Z(n,dvec(n));
    for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
	x=-2.0+4.0*j/(n-1);
	y=-2.0+4.0*i/(n-1);
	sr1=(x-x1)*(x-x1)+(y-y1)*(y-y1);
	sr2=(x-x2)*(x-x2)+(y-y2)*(y-y2);
	sr3=(x-x3)*(x-x3)+(y-y3)*(y-y3);
	Z[i][j]=a1/(1+sr1/sigma) +a2/(1+sr2/sigma) +a3/(1+sr3/sigma) ;
    }
    }
    return Z;
}


// Patch ///
int MatPlot::patch(){
    int h=iPatch*100 + tPatch; hObj=h;
    if(is_debug1){printf("mode: %d handle: %4d Patch\n",mode,h);}
    
    if(mode==0){	    
	ca->add_child(h);	   
	vPatch.push_back(Patch(h));
	//as.Parent=gca();
    }
    if(mode==1){
	
    }
    if(iPatch<vPatch.size()){cp=&vPatch[iPatch];}
    iPatch++;    
    return h;
}
void MatPlot::patch_config(){
    ca->xmax=fmax( Fmax(cp->XData), ca->xmax);
    ca->xmin=fmin( Fmin(cp->XData), ca->xmin);
    ca->ymax=fmax( Fmax(cp->YData), ca->ymax);
    ca->ymin=fmin( Fmin(cp->YData), ca->ymin);
    ca->zmax=fmax( Fmax(cp->ZData), ca->zmax );
    ca->zmin=fmin( Fmin(cp->ZData), ca->zmin );
}
tcvec MatPlot::Index2TrueColor(dvec IC){
    if(ca->CLim[0]==ca->CLim[1]){
	ca->CLim[0]=fmin( Fmin(IC), Fmin(IC) );
	ca->CLim[1]=fmax( Fmax(IC), Fmax(IC) );	
    }
    vector<float> rgb;
    tcvec tc;	
    for(int j=0;j<IC.size();++j){
	rgb=map2color(IC[j],ca->CLim[0],ca->CLim[1]);
	tc.push_back(rgb);
    }
    return tc;    
}
/// patch
int MatPlot::patch(dmat X,dmat Y){
    // Single color
    int h=patch();
    cp->type=0;

    cp->XData=X;
    cp->YData=Y;
    cp->ZData.clear();
    cp->CData.clear();

    patch_config();
    return h;
}
int MatPlot::patch(dmat X,dmat Y,dvec C){
    // One color per face with index color 
    int h=patch();
    cp->type=0;
    cp->XData=X;
    cp->YData=Y;
    cp->ZData.clear();
    cp->CData=Index2TrueColor(C);

    patch_config();
    return h;
}
int MatPlot::patch(dmat X,dmat Y,tcvec C){
    // One color per face with true color
    int h=patch();
    cp->type=0;
    cp->XData=X;
    cp->YData=Y;
    cp->ZData.clear();
    cp->CData=C;

    patch_config();
    return h;
}
int MatPlot::patch(dmat X,dmat Y,dmat Z){
    // Single color
    int h=patch(); 
    ca->View=1;
    cp->type=1;
    cp->XData=X;
    cp->YData=Y;
    cp->ZData=Z;
    cp->CData.clear();

    patch_config();
    return h;
}
int MatPlot::patch(dmat X,dmat Y,dmat Z,dvec C){
    // One color per face    
    int h=patch();
    ca->View=1;
    cp->type=1;

    cp->XData=X;
    cp->YData=Y;
    cp->ZData=Z;
    cp->CData=Index2TrueColor(C);

    patch_config();
    return h;
}
int MatPlot::patch(dmat X,dmat Y,dmat Z,tcvec C){
    // One color per face    
    int h=patch();    
    ca->View=1;
    cp->type=1;

    cp->XData=X;
    cp->YData=Y;
    cp->ZData=Z;
    cp->CData=C;

    patch_config();
    return h;
}
/// bar

int MatPlot::bar(dvec y){ 
    dvec x;
    x.resize(y.size());
    for(int i=0;i<y.size();++i){
	x[i]=1.0*(1+i);
    }
    return bar(x,y,0.8);
}
int MatPlot::bar(dvec y,float width){ 
    dvec x;
    x.resize(y.size());
    for(int i=0;i<y.size();++i){
	x[i]=1.0*(1+i);
    }
    return bar(x,y,width);
}
int MatPlot::bar(dvec x,dvec y){ 
    return bar(x,y,0.8);
}
int MatPlot::bar(dvec x,dvec y,float width){    
    int h=patch();
    cp->type=0;
    cp->XData.clear();
    cp->YData.clear();
    cp->ZData.clear();
    //cp->CData.clear();
    
    double wx=width*(Fmax(x)-Fmin(x))/x.size();

    dvec X(4),Y(4);
    for(int i=0;i<x.size();++i){
	X[0]=x[i]-wx/2.0; Y[0]=0;
	X[1]=x[i]+wx/2.0; Y[1]=0;
	X[2]=x[i]+wx/2.0; Y[2]=y[i];
	X[3]=x[i]-wx/2.0; Y[3]=y[i];
	cp->XData.push_back(X);
	cp->YData.push_back(Y);
    }
    patch_config();
    return h;	    
}

/// display

void MatPlot::display_patch(){
    if(is_debug1){printf("display:%4d Patch :type %d, Faces.size %d\n",
			 cp->id,cp->type,cp->Faces.size());}

    if(cp->type==0){ display_patch_2d(); }
    if(cp->type==1){ display_patch_3d(); }
}

void MatPlot::display_patch_2d(){
    // FaceVertex & CData //
    
    vector<float> v(3);
    vector<int> f(3);
    float x,y;
    for(int i=0;i<cp->Faces.size();++i){
	f=cp->Faces[i];
	glBegin(GL_TRIANGLES);	
	x=ctx(cp->Vertices[f[0]][0]);
	y=cty(cp->Vertices[f[0]][1]);
	glVertex2d(x,y);
	x=ctx(cp->Vertices[f[1]][0]);
	y=cty(cp->Vertices[f[1]][1]);
	glVertex2d(x,y);
	x=ctx(cp->Vertices[f[2]][0]);
	y=cty(cp->Vertices[f[2]][1]);
	glVertex2d(x,y);	
	glEnd();
    }
    

    // XYZ Data //
    int nf,nv;//number of faces and vertex
    nf=cp->XData.size();
    vector<float> rgb;

    for(int i=0;i<nf;++i){
	nv=cp->XData[i].size();

	// Edge
	if(cp->EdgeColor!="none"){
	    
	    glLineWidth(cp->LineWidth);
	    gl2psLineWidth(cp->LineWidth);
	    
	    rgb=ColorSpec2RGB(cp->EdgeColor);
	    glColor3f(rgb[0],rgb[1],rgb[2]); 	    
		
	    glBegin(GL_LINE_LOOP);
	    for(int iv=0;iv<nv;++iv){
		glVertex2d( ctx(cp->XData[i][iv]),
			    cty(cp->YData[i][iv]) );
	    }
	    glEnd();
		
	}
	// Face
	if(cp->FaceColor!="none"){
	    
	    rgb=ColorSpec2RGB(cp->FaceColor);
	    glColor3f(rgb[0],rgb[1],rgb[2]); 	    
	    
	    if(cp->CData.size()){
		rgb=cp->CData[i];
		glColor3d(rgb[0],rgb[1],rgb[2]);
	    }
	    
	    glBegin(GL_POLYGON);
	    for(int iv=0;iv<nv;++iv){
		glVertex2d( ctx(cp->XData[i][iv]),
			    cty(cp->YData[i][iv]) );
	    }
	    glEnd();
		
	}
    }
    

}
void MatPlot::display_patch_3d(){
    // FaceVertex & CData //
    // XYZ Data //
    int nf,nv;//number of faces and vertex
    nf=cp->XData.size();
    vector<float> rgb;

    for(int i=0;i<nf;++i){
	nv=cp->XData[i].size();

	// Edge
	if(cp->EdgeColor!="none"){
	    
	    glLineWidth(cp->LineWidth);
	    gl2psLineWidth(cp->LineWidth);
	    
	    rgb=ColorSpec2RGB(cp->EdgeColor);
	    glColor3f(rgb[0],rgb[1],rgb[2]); 	    
		
	    glBegin(GL_LINE_LOOP);
	    for(int iv=0;iv<nv;++iv){
		glVertex3d( ct3x(cp->XData[i][iv]),
			    ct3y(cp->YData[i][iv]), 
			    ct3z(cp->ZData[i][iv]) );
	    }
	    glEnd();
		
	}
	// Face
	if(cp->FaceColor!="none"){
	    
	    rgb=ColorSpec2RGB(cp->FaceColor);
	    glColor3f(rgb[0],rgb[1],rgb[2]); 	    
	    
	    if(cp->CData.size()){
		rgb=cp->CData[i];
		glColor3d(rgb[0],rgb[1],rgb[2]);
	    }
	    
	    glBegin(GL_POLYGON);
	    for(int iv=0;iv<nv;++iv){
		glVertex3d( ct3x(cp->XData[i][iv]),
			    ct3y(cp->YData[i][iv]), 
			    ct3z(cp->ZData[i][iv]) 			    
			    );
	    }
	    glEnd();
		
	}
    }
}
// Text ///
//TODO more fonts
/// text 
int MatPlot::text(){
    int h=iText*100 + tText; hObj=h;
    if(is_debug1){printf("mode: %d handle: %4d Text\n",mode,h);}
    
    if(mode==0){
	vText.push_back(Text(h));
	ca->add_child(h);
	//at.Parent=gca();
    }
    if(mode==1){
    }
    if(iText<vText.size()){ct=&vText[iText];}
    iText++;    
    return h;
}
void MatPlot::display_text(){
    if(is_debug1){printf("display:%4d Text\n",ct->id);}    
    glColor3f(0,0,0);
    glRasterPos2f( ctx(ct->Position[0]), cty(ct->Position[1]) );
    gl2psText(ct->String.c_str(), "Arial", 12);		
    for(int i=0; i<(int)ct->String.size(); ++i ){
	glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, ct->String[i] );
    }
}
int MatPlot::text(double x,double y,string s){
    // text on current axes	
    //if(is_debug1){cout<<"mode:"<<mode<<" Text (text):"<<i_text<<endl;}
    int h=text();    
    ct->Position[0]=x;
    ct->Position[1]=y;
    ct->String=s;
    return h;
}
void MatPlot::set_font(char font_[],int size){
    //font=font_;
    //font_size=size;
}
/// ptext
void MatPlot::ptext(float x,float y,string s){
    // Viewport Figure 
    glViewport(0,0, (int)(window_w),(int)(window_h));
    glLoadIdentity();
    gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
    
    glColor3f(0,0,0);
    glRasterPos2f( x, y );
    gl2psText(s.c_str(), "Arial", 12);
    //gl2psText(test, "Times-Roman", 24);
    
    for(int i=0; i<(int)s.size(); i++ ){
	glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );				
    }
}    
void MatPlot::ptext3(float x,float y,float z,string s){
    glColor3f(0,0,0);
    glRasterPos3d(x,y,z);
    gl2psText(s.c_str(), "Arial", 12);    
    for(int i=0; i<(int)s.size(); i++ ){
	glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );				
    }
}  
void MatPlot::ptext3c(float x,float y,float z,string s){
    int char_w=6,char_h=12;
    int num_char=s.length();
    glColor3f(0,0,0);
    glRasterPos3d(x,y,z);
    glBitmap(0,0,0,0,-char_w*num_char/2,-char_h/2,NULL);
    gl2psText(s.c_str(), "Arial", 12);
    
    for(int i=0; i<(int)s.size(); i++ ){
	glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );				
    }
} 


/// Color ///

void MatPlot::Shading(string c){
    shading(c);
}
void MatPlot::shading(string c){

    int tObj=hObj%100;
    int iObj=hObj/100;
    if( (tObj==tSurface) && (iObj<vSurface.size()) ){
	if(c=="faceted"){ 
	    vSurface[iObj].EdgeColor="k"; 
	}
	if(c=="flat"){ 
	    vSurface[iObj].EdgeColor="none";
	}
    }
}

/// colormap
    
vector<float> MatPlot::colormap(string c,float t){
    
    vector<float> rgb(3);
    if(t>1){t=1;}
    if(t<0){t=0;}
    
    if( c=="Gray" ){
	rgb[0]=t;
	rgb[1]=t;
	rgb[2]=t;
	return rgb;
    }
    if( c=="HSV" ){
	t*=6;
	if(  0<=t && t<=1.0){rgb[0] = 1;        rgb[1] = t;        rgb[2] = 0;}
	if(1.0<=t && t<=2.0){rgb[0] = 1-(t-1);  rgb[1] = 1;        rgb[2] = 0;}
	if(2.0<=t && t<=3.0){rgb[0] = 0;        rgb[1] = 1;        rgb[2] = t-2;}
	if(3.0<=t && t<=4.0){rgb[0] = 0;        rgb[1] = 1-(t-3);  rgb[2] = 1;}
	if(4.0<=t && t<=5.0){rgb[0] = t-4;      rgb[1] = 0;        rgb[2] = 1;}
	if(5.0<=t && t<=6.0){rgb[0] = 1;        rgb[1] = 0;        rgb[2] = 1-(t-5);}
	
	return rgb;	
    }
    if( c=="Jet" ){
	t*=8;
	if(  0<=t && t<=1.0){rgb[0] = 0;        rgb[1] = 0;          rgb[2] = 0.5+0.5*t;}
	if(1.0<=t && t<=3.0){rgb[0] = 0;        rgb[1] = 0.5*(t-1);  rgb[2] = 1;}
	if(3.0<=t && t<=5.0){rgb[0] = 0.5*(t-3);rgb[1] = 1;          rgb[2] = 1-0.5*(t-3);}
	if(5.0<=t && t<=7.0){rgb[0] = 1;        rgb[1] = 1-0.5*(t-5);rgb[2] = 0;}
	if(7.0<=t && t<=8.0){rgb[0] = 1-0.5*(t-7);rgb[1] = 0;        rgb[2] = 0;}
	return rgb;	
    }
    if( c=="Hot" ){
	t*=3;
	if(  0<=t && t<=1.0){rgb[0] = t; rgb[1] = 0;   rgb[2] = 0;}
	if(1.0<=t && t<=2.0){rgb[0] = 1; rgb[1] = t-1; rgb[2] = 0;}
	if(2.0<=t && t<=3.0){rgb[0] = 1; rgb[1] = 1;   rgb[2] = t-2;}
	return rgb;
    }
    
    if( c=="Cool" ){
	rgb[0]=t;
	rgb[1]=1-t;
	rgb[2]=1;
	return rgb;
    }
    if( c=="Spring" ){// Magenta - Yellow
	rgb[0]=1;
	rgb[1]=t;
	rgb[2]=1-t;
	return rgb;
    }
    if( c=="Summer" ){// Green Yellow
	rgb[0]=t;
	rgb[1]=1;
	rgb[2]=0;
	return rgb;
    }
    if( c=="Autumn" ){
	rgb[0]=1;
	rgb[1]=t;
	rgb[2]=0;
	return rgb;
    }
    if( c=="Winter" ){
	rgb[0]=0;
	rgb[1]=t;
	rgb[2]=1-t;
	return rgb;
    }
    if( c=="Bone" ){
	rgb[0]=t;
	//rgb[1]=t;if(t<0.8){rgb[1]=t;}
	rgb[2]=t;
	return rgb;
    }
}
void MatPlot::colormap(string c){
    //if(is_debug1){printf("colormap %s \n",c.c_str());}
    int n;
    float t;
    vector<float> rgb(3);

    cmap.clear();
    n=64;
    for(int i=0;i<n;++i){
	rgb=colormap(c,(float)i/(n-1));		
	cmap.push_back(rgb);
    }    
    if(mode==1){
	//cs->iColorMap=iColorMap;
	cs->ColorMap=c;
    }    
}
void MatPlot::colormap(vector<vector<float> > c){
    cmap=c;
}
void MatPlot::gray(){colormap("Gray");};
void MatPlot::jet(){ colormap("Jet"); }
void MatPlot::hsv(){ colormap("HSV"); }
void MatPlot::hot(){ colormap("Hot"); }
void MatPlot::cool(){colormap("Cool");}  
void MatPlot::spring(){colormap("Spring");}   
void MatPlot::summer(){colormap("Summer");}     
void MatPlot::autumn(){colormap("Autumn");}
void MatPlot::winter(){colormap("Winter");}     

vector<float> MatPlot::map2color(double x,double xmin,double xmax){
    int n=cmap.size();
    float normx;
    vector<float> rgb(3);
    
    normx=(x-xmin)/(xmax-xmin);
    if(x>xmax){normx=1;}
    if(x<xmin){normx=0;}
    rgb=cmap[(int)(normx*(n-1))];
    //cout << "c: "<<(int)(normx*n) <<endl;
    //cout << "rgb: "<<rgb[0]<<" "<<endl;
    return rgb;
}

vector<float> MatPlot::ColorSpec2RGB(string c){
    float r,g,b;
    //line
    if( c=="k" ){r=0;g=0;b=0;}// black
    if( c=="r" ){r=1;g=0;b=0;}// red
    if( c=="b" ){r=0;g=0;b=1;}// blue
    if( c=="g" ){r=0;g=1;b=0;}// green	    
    if( c=="c" ){r=0;g=1;b=1;}// cyan
    if( c=="m" ){r=1;g=0;b=1;}// magenta
    if( c=="y" ){r=1;g=1;b=0;}// yellow
    if( c=="w" ){r=1;g=1;b=1;}// white

    //dark color
    float h=0.6;
    if( c=="dr" ){r=h;g=0;b=0;}// red
    if( c=="db" ){r=0;g=0;b=h;}// blue
    if( c=="dg" ){r=0;g=h;b=0;}// green	    
    if( c=="dc" ){r=0;g=h;b=h;}// cyan
    if( c=="dm" ){r=h;g=0;b=h;}// magenta
    if( c=="dy" ){r=h;g=h;b=0;}// yellow
    
    //light color
    h=0.5;
    if( c=="lr" ){r=1;g=h;b=h;}// red
    if( c=="lb" ){r=h;g=h;b=1;}// blue
    if( c=="lg" ){r=h;g=1;b=h;}// green	    
    if( c=="lc" ){r=h;g=1;b=1;}// cyan
    if( c=="lm" ){r=1;g=h;b=1;}// magenta
    if( c=="ly" ){r=1;g=1;b=h;}// yellow
    
    //universal color
    if( c=="ur" ){r=1;   g=0.2; b=0;  }//red
    if( c=="ub" ){r=0;   g=0.25;b=1;  }//blue
    if( c=="ug" ){r=0.2; g=0.6; b=0.4;}//green
    if( c=="uy" ){r=1;   g=1;   b=1;  }//yellow
    if( c=="uc" ){r=0.4; g=0.8; b=1;  }//sky blue
    if( c=="up" ){r=1;   g=0.6; b=0.6;}//pink
    if( c=="uo" ){r=1;   g=0.6; b=0;  }//orange
    if( c=="um" ){r=0.6; g=0;   b=0.4;}//perple
    if( c=="ubr"){r=0.4; g=0.2; b=0;  }//brown

    char line[1000],*tp;
    char d[]=" ,:\t\n[]";//delimiter
    /*
    if(c.size()){
	if( (c[0]=='[') && (c[c.size()-1]==']') ){
	    sprintf(line,c.c_str());
	    tp=strtok(line,d); if(tp){r=atof(tp);}
	    tp=strtok(NULL,d); if(tp){g=atof(tp);}
	    tp=strtok(NULL,d); if(tp){b=atof(tp);}
	}
    }
    */
    //cout <<"c,r,g,b: "<< c <<" "<<r<<" "<<g <<" "<<b<<endl;
    vector<float> out(3);
    out[0]=r;
    out[1]=g;
    out[2]=b;
    return out;
       
}
string MatPlot::rgb2colorspec(vector<float> rgb){
    char c[100];
    sprintf(c,"[%f %f %f]",rgb[0],rgb[1],rgb[2]);
    string s=c;
    return s;
}


