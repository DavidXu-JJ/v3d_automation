/*
 * Copyright (c)2006-2010  Hanchuan Peng (Janelia Farm, Howard Hughes Medical Institute).
 * All rights reserved.
 */


/************
                                            ********* LICENSE NOTICE ************

This folder contains all source codes for the V3D project, which is subject to the following conditions if you want to use it.

You will ***have to agree*** the following terms, *before* downloading/using/running/editing/changing any portion of codes in this package.

1. This package is free for non-profit research, but needs a special license for any commercial purpose. Please contact Hanchuan Peng for details.

2. You agree to appropriately cite this work in your related studies and publications.

Peng, H., Ruan, Z., Long, F., Simpson, J.H., and Myers, E.W. (2010) “V3D enables real-time 3D visualization and quantitative analysis of large-scale biological image data sets,” Nature Biotechnology, Vol. 28, No. 4, pp. 348-353, DOI: 10.1038/nbt.1612. ( http://penglab.janelia.org/papersall/docpdf/2010_NBT_V3D.pdf )

Peng, H, Ruan, Z., Atasoy, D., and Sternson, S. (2010) “Automatic reconstruction of 3D neuron structures using a graph-augmented deformable model,” Bioinformatics, Vol. 26, pp. i38-i46, 2010. ( http://penglab.janelia.org/papersall/docpdf/2010_Bioinfo_GD_ISMB2010.pdf )

3. This software is provided by the copyright holders (Hanchuan Peng), Howard Hughes Medical Institute, Janelia Farm Research Campus, and contributors "as is" and any express or implied warranties, including, but not limited to, any implied warranties of merchantability, non-infringement, or fitness for a particular purpose are disclaimed. In no event shall the copyright owner, Howard Hughes Medical Institute, Janelia Farm Research Campus, or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; reasonable royalties; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.

4. Neither the name of the Howard Hughes Medical Institute, Janelia Farm Research Campus, nor Hanchuan Peng, may be used to endorse or promote products derived from this software without specific prior written permission.

*************/




/****************************************************************************
**
V3D main program

By Hanchuan Peng

July 21, 2006
060924: v3d v1.2

Last update: 2008-04-25: try to add command line based utilities
Last update: 2010-04-12: add a global try-catch to catch all exceptions
Last update: 2010-11-19: change some of the help info
Last update: 2011-04-19: fix some potential problem of null mainWin pointer
Last update: 2011-08-25: remove some uncalled old code, and adjust the inconsistent return values of the main function

****************************************************************************/

#define COMPILE_TO_COMMANDLINE 1

#include "../3drenderer/v3dr_common.h"

#include "v3d_compile_constraints.h"

#include <QApplication>
#include <QFile>

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
#include <queue>

#include "mainwindow.h"
#include "v3d_application.h"
#include "customdebug.h"

#include <string>
using namespace std;

#include "v3d_core.h"

#include "v3d_commandlineparser.h"
#include "v3d_version_info.h"
#include "../plugin_loader/v3d_plugin_loader.h"


void printHelp_v3d();
void printHelp_align();
void printHelp_straight();
void printHelp_trace();

V3dApplication* V3dApplication::theApp = 0;

void printHelp_v3d()
{
    cout<<endl<<"Vaa3D: a 3D image visualization and analysis platform developed by Hanchuan Peng and colleagues."<<endl;
    cout<<endl<<"Usage: v3d -h -M moduleCode [all other options specific to different modules]"<<endl;

    cout<<"    -h/H         help information."<<endl;

    cout<<"    -i <file>                    open single or multiple image (.tif/.tiff, .lsm, .mrc, .raw/.v3draw) / object (.ano, .apo, .swc, .marker) files"<<endl;
    cout<<"    -o <file>                    indicates single or multiple outputs"<<endl;
    cout<<"    -x <plugin_dll_full_path or unique partial path>    a string indicates the full path or a unique partial path of a dll (for a plugin) to be launched."<<endl;
    cout<<"    -m <menu_name>               a string indicates which menu of a plugin will be called."<<endl;
    cout<<"    -f <function_name>           a string indicates which function of a plugin will be called."<<endl;
    cout<<"    -p <parameters>              a string indicates parameters that plugin function use"<<endl;
    cout<<"    -pf <configuration>          a string read from configuration file indicates parameters that plugin function use"<<endl;

    cout<<"    -v                           force to open a 3d viewer when loading an image, otherwise use the default v3d global setting (from \"Adjust Preference\")"<<endl;
    cout<<"    -na                          open NeuronAnnotator work-mode directly"<<endl;
    cout<<"    -cmd  [headless command-line arguments, intended for compute grid use. Try \'-cmd -h\' for more information on this option]"<<endl;

    //added by Hanchuan Peng, 20120217
    V3dApplication* app = V3dApplication::getInstance();
    if (!app) return;
    MainWindow* mainWin=app->getMainWindow();
    if (!mainWin) return;
    QStringList existingPluginsList = mainWin->pluginLoader->getPluginNameList();
    if (existingPluginsList.size()>0)
        cout << endl << "Found [" << existingPluginsList.size() << "] plugins"<<endl;
    for (int i=0;i<existingPluginsList.size();i++)
        cout << "#" << i+1 << "          " << existingPluginsList.at(i).toLocal8Bit().constData() << endl;

    return;
}

QString generate_apo_name(const std::string & path,const CellAPO &centerAPO){
    return QString::fromStdString(path+"/"+std::to_string(int(centerAPO.x))+".000_"+std::to_string(int(centerAPO.y))+".000_"+std::to_string(int(centerAPO.z))+".000.apo");
}

QString generate_marker_name(const std::string & path,const ImageMarker &startPoint){
    return QString::fromStdString(path+"/"+std::to_string(int(startPoint.x))+".000_"+std::to_string(int(startPoint.y))+".000_"+std::to_string(int(startPoint.z))+".000.marker");
}

QString generate_eswc_name(const std::string & path,const CellAPO &centerAPO){
    return QString::fromStdString(path+"/"+std::to_string(int(centerAPO.x))+".000_"+std::to_string(int(centerAPO.y))+".000_"+std::to_string(int(centerAPO.z))+".000.eswc");
}

QString generate_swc_name(const std::string & path,const CellAPO &centerAPO){
    return QString::fromStdString(path+"/"+std::to_string(int(centerAPO.x))+".000_"+std::to_string(int(centerAPO.y))+".000_"+std::to_string(int(centerAPO.z))+".000.swc");
}

//struct node{
//    int x,y,z,marker_x,marker_y,marker_z;
//    node(int x_,int y_,int z_,int marker_x_,int marker_y_,int marker_z_):x(x_),y(y_),z(z_),marker_x(marker_x_),marker_y(marker_y_),marker_z(marker_z_){}
//    bool operator == (const node & cmp) const{
//        return (x==cmp.x)&&(y==cmp.y)&&(z==cmp.z)&&(marker_x==cmp.marker_x)&&(marker_y==cmp.marker_y)&&(marker_z==cmp.marker_z);
//    }
//    bool operator < (const node & cmp) const{
//        return x<cmp.x;
//    }
//};

struct bbox_extend{
    int center_id;
    CellAPO centerAPO;
    ImageMarker startPoint;
    int centerSWC;
    bbox_extend(int _id_,CellAPO _centerAPO_,ImageMarker _startPoint_,int _centerSWC_):center_id(_id_),centerAPO(_centerAPO_),startPoint(_startPoint_),centerSWC(_centerSWC_){}
};

QString Res_Path="G:/18454/RES(26298x35000x11041)";
QString Vaa3d_App_Path="C:/3.603c";
QString Work_Dir="G:/work_dir_debug";
QString Answer_File="G:/18454_answer/whole_image.eswc";
QMap<int,QVector<int> > Answer_Graph;
QMap<int,NeuronSWC> Answer_Map;
NeuronTree Ans_Tree;
QMap<int,bool> Ans_used;

int cnt=0;

QProcess p;

int dx[6]={1,-1,0,0,0,0};
int dy[6]={0,0,1,-1,0,0};
int dz[6]={0,0,0,0,1,-1};
//QMap<node,bool> vis;
QVector<NeuronSWC> used_swc;

int X=14530,Y=10693,Z=3124;
double close_distance=4;
double identical_threshold=300;   //(undetermined)
double seg_identical_threshold_mean=50;
double seg_identical_threshold_mx=100;
double Border_Threshold=50;

double distance_square(const NeuronSWC & point_a,const NeuronSWC & point_b){
    return (point_a.x-point_b.x)*(point_a.x-point_b.x)+(point_a.y-point_b.y)*(point_a.y-point_b.y)+(point_a.z-point_b.z)*(point_a.z-point_b.z);
}

bool between(const int & mid,const int & left,const int & right){   //double?
    return mid>=left&&mid<=right;
}

bool if_finish(){
    int cnt=0;
    for(const NeuronSWC & i:Ans_Tree.listNeuron){
        if(Ans_used.count(i.n)){
            ++cnt;
        }
    }
    return cnt>=Ans_Tree.listNeuron.size();
}
int unused_id(){
    int id=-1;
    double mxmn=-1e8;
    for(const NeuronSWC & i:Ans_Tree.listNeuron){
        if(!Ans_used.count(i.n)){
            double now_mn=1e8;
            for(auto it=Ans_used.begin();it!=Ans_used.end();++it){
                now_mn=std::min(now_mn,distance_square(Answer_Map[it.key()],i));
            }
            if(mxmn<now_mn){
                mxmn=now_mn;
                id=i.n;
            }
        }
    }
    return id;
}

//int unused_id(){
//    int id=-1;
//    for(const NeuronSWC & i:Ans_Tree.listNeuron){
//        if(!Ans_used.count(i.n)){
//            id=i.n;
//        }
//    }
//    return id;
//}

bool if_need_extend(const CellAPO & centerAPO,const int & blocksize){
    const int & X=centerAPO.x;
    const int & Y=centerAPO.y;
    const int & Z=centerAPO.z;
    for(const NeuronSWC & swc:Ans_Tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            if(!Ans_used.count(swc.n)){
                return true;
            }
        }
    }
    return false;
}

void update_Ans_used(const NeuronTree & updateTree){
    for(const NeuronSWC & swc1:Ans_Tree.listNeuron){
        for(const NeuronSWC & swc2:updateTree.listNeuron)
        if(swc1.n==swc2.n){
            Ans_used[swc1.n]=true;
        }
    }
}

void update_Ans_used(const NeuronSWC & update){
    Ans_used[update.n]=true;
}

NeuronTree Get_Ans_In_BBox(const CellAPO & centerAPO,const int & blocksize){
    NeuronTree ret;
    QMap<int,bool> Inbox;
    const int & X=centerAPO.x;
    const int & Y=centerAPO.y;
    const int & Z=centerAPO.z;
    for(const NeuronSWC & swc:Ans_Tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            Inbox[swc.n]=true;
            NeuronSWC add=swc;
            add.x=swc.x-X+blocksize/2;
            add.y=swc.y-Y+blocksize/2;
            add.z=swc.z-Z+blocksize/2;
            ret.listNeuron.push_back(add);
        }
    }
    for(NeuronSWC & i:ret.listNeuron){
        if(!Inbox.count(i.parent)){
            i.parent=-1;
        }
    }
    return ret;
}


double distance_XYZ(const XYZ & p1,const XYZ & p2){
    XYZ vec=p2-p1;
    return sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
}
double distance_unit(const V_NeuronSWC_unit & p1,const V_NeuronSWC_unit & p2){
    XYZ vec=XYZ(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
    return sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
}
double dot_mul(const XYZ & mid,const XYZ & p1,const XYZ & p2){  //dot multiply mid->p1 and mid-p2
    XYZ diff1=p1-mid,diff2=p2-mid;
    return diff1.x*diff2.x+diff1.y*diff2.y+diff1.z*diff2.z;
}
double dot_mul(const XYZ & p1,const XYZ & p2){                  //dot multiply vector p1 and vector p2
    return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
}
XYZ projection(const XYZ & cal_point,const XYZ & p1,const XYZ & p2){
    double dm=dot_mul(p1,cal_point,p2);
    return p1+(p2-p1)*(dm/distance_XYZ(p1,p2))/distance_XYZ(p1,p2);
}
double dot_to_line(const XYZ & cal_point,const XYZ & p1,const XYZ & p2){
    XYZ proj=projection(cal_point,p1,p2);
    double dm=dot_mul(proj-p1,proj-p2);
    if(dm<0) return distance_XYZ(cal_point,proj);
    return std::min(distance_XYZ(cal_point,p1),distance_XYZ(cal_point,p2) );
}
double Distance_Unit_To_Seg(const V_NeuronSWC_unit & p, const V_NeuronSWC & s){  //传参XYZ可以优化
    double mn=1e8;
    XYZ check_point;
    check_point.x=p.x;check_point.y=p.y;check_point.z=p.z;
    QMap<int,V_NeuronSWC_unit> mp;
    for(int i=0;i<s.row.size()-1;++i){
        mp[s.row[i].n]=s.row[i];
//        qDebug()<<s.row[i].n<<" "<<s.row[i].parent<<" "<<s.row[i].nchild<<" "<<s.row[i].nodeinseg_id<<" "<<s.row[i].seg_id;
    }
    for(int i=0;i<s.row.size()-1;++i){
        XYZ p1,p2;
        int parent_id=s.row[i].parent;
        if(parent_id>0){
            p1.x=s.row[i].x;p1.y=s.row[i].y;p1.z=s.row[i].z;
            p2.x=mp[parent_id].x;p2.y=mp[parent_id].y;p2.z=mp[parent_id].z;
            mn=std::min(mn,dot_to_line(check_point,p1,p2));
        }
    }

    return mn;
}
double Distance_Unit_To_Tree(const V_NeuronSWC_unit & p,const V_NeuronSWC_list & Check_Tree){
    if(Check_Tree.seg.empty()) return 1e-8;
    double mn=1e8;
    for(const V_NeuronSWC & Seg:Check_Tree.seg){
        mn=std::min(mn,Distance_Unit_To_Seg(p,Seg));
    }
    return mn;
}
double Distance_Unit_To_Tree(const NeuronSWC & p,const V_NeuronSWC_list & Check_Tree){
    V_NeuronSWC_unit u;
    u.x=p.x;u.y=p.y;u.z=p.z;
    return Distance_Unit_To_Tree(u,Check_Tree);
}
double Tree_Length(const V_NeuronSWC_list & tree){
    double ret=0;
    for(const V_NeuronSWC & seg:tree.seg){
        for(int i=0;i<seg.row.size()-1;++i){
            ret+=distance_unit(seg.row[i],seg.row[i+1]);
        }
    }
    return ret;
}
bool Check_Seg_Identical(const V_NeuronSWC & Check_Seg,const V_NeuronSWC_list & Answer_Tree){
    double mean=0;
    double mxx=-1e8;
    int num=0;
    for(const V_NeuronSWC_unit & Check_Point:Check_Seg.row){
         double mx=-1e8;
         for(const V_NeuronSWC & Answer_Seg:Answer_Tree.seg){
            mx=std::max(mx,Distance_Unit_To_Seg(Check_Point,Answer_Seg));
            mxx=std::max(mxx,mx);
         }
         ++num;
         mean+=mx;
    }
    mean/=num;
    return mean<seg_identical_threshold_mean && mxx<seg_identical_threshold_mx;
}
bool Check_Tree_Identical(const V_NeuronSWC_list & Check_Tree,const V_NeuronSWC_list & Answer_Tree){
    if(Answer_Tree.seg.empty()) return false;
    bool flag=true;
    for(const V_NeuronSWC & Check_Seg:Check_Tree.seg){
        flag&=Check_Seg_Identical(Check_Seg,Answer_Tree);
        if(!flag) break;
    }
    return flag;
}
bool Check_Tree_Identical(const NeuronTree & Check_Tree,const V_NeuronSWC_list & Answer_Tree){
    NeuronTree cpy=Check_Tree;
    V_NeuronSWC_list ls=NeuronTree__2__V_NeuronSWC_list(cpy);
    return Check_Tree_Identical(ls,Answer_Tree);
}
bool Check_Tree_Length(const V_NeuronSWC_list & Check_Tree,const V_NeuronSWC_list & Answer_Tree){
    double len1=Tree_Length(Check_Tree),len2=Tree_Length(Answer_Tree);
    return std::abs(len1-len2)/std::max(len1,len2)<0.15;
}
double Distance_Point_To_Border(const NeuronSWC & point,const int & blocksize){
    return std::min({point.x,abs(blocksize-point.x),point.y,abs(blocksize-point.y),point.z,abs(blocksize-point.z)});
}
double Vector_Angle(const XYZ & a,const XYZ & b){
    static const double PI=std::acos(-1);
    return std::acos(dot_mul(a,b)/distance_XYZ(XYZ(0.0),a)/distance_XYZ(XYZ(0.0),b))*(180/PI);
}
bool has_same_vector(const V_NeuronSWC & Seg,const QVector<XYZ> & v){
    XYZ v1;
    v1.x=Seg.row.front().x;
    v1.y=Seg.row.front().y;
    v1.z=Seg.row.front().z;
    XYZ v2;
    v2.x=Seg.row.back().x;
    v2.y=Seg.row.back().y;
    v2.z=Seg.row.back().z;
    XYZ seg=v1-v2;
    for(const XYZ & vec:v){
        if(Vector_Angle(seg,vec)<30||Vector_Angle(seg,vec)>150){
            return true;
        }
    }
    return false;
}
void Drop_NeuronSWC(const QString &path,const QList<NeuronSWC> & output){
    QFile file(path);
    qDebug()<<file.open(QIODevice::WriteOnly|QIODevice::Text);

    QTextStream out(&file);
    for(const NeuronSWC & swc:output){
        out<<swc.n<<" "<<swc.type<<" "<<swc.x<<" "<<swc.y<<" "<<swc.z<<" "<<swc.r<<" "<<swc.pn<<"\n";
    }
    file.close();
}
void crop_ans_swc(const QString& input_file,const int &X,const int &Y,const int &Z,const int &blocksize,const QString& output_file){
    NeuronTree ans_tree=readSWC_file(input_file);

    QList<NeuronSWC> output;
    for(NeuronSWC & swc:ans_tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            output.push_back(swc);
        }
    }

    QFile file(output_file);
    qDebug()<<file.open(QIODevice::WriteOnly|QIODevice::Text);

    QTextStream out(&file);
    for(NeuronSWC & swc:output){
        out<<swc.n<<" "<<swc.type<<" "<<swc.x-X+blocksize/2<<" "<<swc.y-Y+blocksize/2<<" "<<swc.z-Z+blocksize/2<<" "<<swc.r<<" "<<swc.pn<<"\n";
    }
    file.close();
}
QVector<std::pair<NeuronSWC,QQueue<XYZ> > > Find_Border(const NeuronTree & App2_Tree,const int & blocksize,const QMap<int,QVector<int> > & son, const QMap<int,NeuronSWC> & mp){

    const int compacity=5;      //decide how many vectors to memorize
    QVector<std::pair<int,QQueue<XYZ> > > border;
    QQueue<std::pair<int,QQueue<XYZ> > > q;
    q.push_back(std::make_pair(son[-1][0],QQueue<XYZ>() ) );
    //BFS to find the border
    while(!q.empty()){
        std::pair<int,QQueue<XYZ> > now=q.front();
        q.pop_front();
        int now_id=now.first;
        const NeuronSWC & now_swc=mp[now_id];
        for(const int & to:son[now_id]){
            QQueue nex_queue=now.second;
            const NeuronSWC & nex_swc=mp[to];
            if(nex_queue.size()>=compacity){
                nex_queue.pop_front();
            }
            nex_queue.push_back(XYZ(nex_swc.x-now_swc.x,nex_swc.y-now_swc.y,nex_swc.z-now_swc.z));
            if(son[to].empty()){    //if want to cut more, edit here
                border.push_back(std::make_pair(to,nex_queue));
            }else{
                q.push_back(std::make_pair(to,nex_queue));
            }
        }
    }

    QVector<std::pair<NeuronSWC,QQueue<XYZ> > > ret;
    for(const std::pair<int,QQueue<XYZ> >  & i:border){
        const int & id=i.first;
        const XYZ & point=mp[id];

        double mn=std::min({point.x,abs(blocksize-point.x),point.y,abs(blocksize-point.y),point.z,abs(blocksize-point.z)});

        if(mn>Border_Threshold) continue;

        //save Border_Point and its previous vector
        ret.push_back(std::make_pair(mp[id],i.second));

    }
    return ret;
}
std::vector<int> Judge_Direction(const XYZ & vec){
    int weight[3][3]={{0,1,1},{1,0,1},{1,1,0}};
    double mn=1e10;
    int pick= -1;
    //decide to expand in which axis
    //pick=0 equals weight[i][0]=0, vec is more close to axis x
    //pick=2 equals weight[i][1]=0, vec is more close to axis y
    //pick=4 equals weight[i][2]=0, vec is more close to axis z

    std::vector<int> direction={0,1,2};
    std::sort(direction.begin(),direction.end(),[&](int a,int b){
        return weight[a][0]*vec.x*vec.x+weight[a][1]*vec.y*vec.y+weight[a][2]*vec.z*vec.z < weight[b][0]*vec.x*vec.x+weight[b][1]*vec.y*vec.y+weight[b][2]*vec.z*vec.z;
    });
    for(int i=0;i<3;++i){
        direction[i]*=2;
    }
    for(int i=0;i<3;++i){
        if(direction[i]==0 && vec.x<0)
            ++direction[i];
        if(direction[i]==2 && vec.y<0)
            ++direction[i];
        if(direction[i]==4 && vec.z<0)
            ++direction[i];
    }
    for(int i=2;i>0;--i){
        direction.push_back(direction[i]+((direction[i]&1)?-1:1 ) );
    }

    return direction;

    //decide to expand on positive axis or negative axis
    //pick=0 positive axis x
    //pick=1 negative axis x
    //pick=2 positive axis y
    //pick=3 negative axis y
    //pick=4 positive axis z
    //pick=5 negative axis z
}
int Find_Nearest_Id(const ImageMarker & startPoint,const QVector<NeuronSWC> & Points){
    int mn=1e8;
    int pick_id=-1;
    for(const NeuronSWC & i:Points){
        double dis=distance_XYZ(XYZ(startPoint),XYZ(i));
        if(mn>dis){
            mn=dis;
            pick_id=i.n;
        }
    }
    return pick_id;
}
NeuronTree Get_Answer_Tree(const int & depth,const CellAPO & centerAPO,const ImageMarker & startPoint,const int & blocksize,const int & centerSWC){
    const double eps=1e-8;

    const int & X=centerAPO.x;
    const int & Y=centerAPO.y;
    const int & Z=centerAPO.z;

    QVector<NeuronSWC> In_Block_Points;
    QMap<int,int> weight;
    for(const NeuronSWC & swc:Ans_Tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            NeuronSWC Add_Point=swc;
            Add_Point.x=Add_Point.x-X+blocksize/2;
            Add_Point.y=Add_Point.y-Y+blocksize/2;
            Add_Point.z=Add_Point.z-Z+blocksize/2;
            In_Block_Points.push_back(Add_Point);
            ++weight[Add_Point.n];
            ++weight[Add_Point.pn];
        }
    }

    QMap<int,NeuronSWC> mp;
    QMap<int,QVector<int> > G;
    for(const NeuronSWC & i:In_Block_Points){
        if(i.n>0)
            mp[i.n]=i;
    }
    for(const NeuronSWC & i:In_Block_Points){
        if(mp.count(i.n) && mp.count(i.pn)){
            G[i.n].push_back(i.pn);
            G[i.pn].push_back(i.n);
        }
    }

    if(depth==0){
        NeuronTree Return_Tree;
        NeuronSWC StartPoint_In_Answer=mp[centerSWC];
        StartPoint_In_Answer.n=1;

        int amount=0;
        QQueue<std::pair<int,int> > q;
        QMap<int, bool> vis;
        vis[-1]=true;
        q.push_back(std::make_pair(-1,centerSWC));    //-1 is new, pick_id is old
        while(!q.empty()){
            int prev=q.front().first;   //new
            int now=q.front().second;   //old
            q.pop_front();
            vis[now]=true;

            NeuronSWC Add_Point=mp[now];
            if(distance_XYZ(XYZ(Add_Point), XYZ(0.0))<eps) continue;
            ++amount;
            Add_Point.n = amount;
            Add_Point.pn= prev;
            Return_Tree.listNeuron.push_back(Add_Point);

            for(const int & i:G[now]){
                if(!vis[i])
                    q.push_back(std::make_pair(amount,i));
            }
        }
        return Return_Tree;
    }

    if(In_Block_Points.empty()) return NeuronTree();

    QVector<int> Border_Points;
    for(auto i=weight.begin();i!=weight.end();++i){
        if(i.value()==1){
            Border_Points.push_back(i.key());
        }
    }

    XYZ Marker_XYZ=XYZ(startPoint);

    int pick_id;
    double mn=1e10;
    for(const int & i:Border_Points){
        double dis=distance_XYZ(Marker_XYZ,XYZ(mp[i]));
        if(mn>dis){
            mn=dis;
            pick_id=i;
        }
    }

    //这个参数我认为有点底层，不太应该给用户调整，
    //表示如果答案在这个bbox的边界点中，距离本次DFS距离最小的Marker距离如果大于threshold，就舍弃，认为这个本次的bbox中没有以Marker为起点的树
    //之前生成的图为了保证尽可能搜得广，所以在blocksize为256的时候，直接取500，实际合理的threshold没有测试过，姑且取了blocksize/4
//    const int threshold=1.171875*blocksize;
//    if(mn>threshold) return NeuronTree();

    NeuronTree Return_Tree;
    NeuronSWC StartPoint_In_Answer=mp[pick_id];
    StartPoint_In_Answer.n=1;
    StartPoint_In_Answer.pn=-1;

    int amount=0;
    QQueue<std::pair<int,int> > q;
    QMap<int, bool> vis;
    vis[-1]=true;
    q.push_back(std::make_pair(-1,pick_id));    //-1 is new, pick_id is old
    while(!q.empty()){
        int prev=q.front().first;   //new
        int now=q.front().second;   //old
        q.pop_front();
        vis[now]=true;

        NeuronSWC Add_Point=mp[now];
        if(distance_XYZ(XYZ(Add_Point), XYZ(0.0))<eps) continue;
        ++amount;
        Add_Point.n = amount;
        Add_Point.pn= prev;
        Return_Tree.listNeuron.push_back(Add_Point);

        for(const int & i:G[now]){
            if(!vis[i])
                q.push_back(std::make_pair(amount,i));
        }
    }

    return Return_Tree;
}

QVector<NeuronSWC>  Find_Extend_Marker(const CellAPO & centerAPO,const ImageMarker &Absolute_Marker,const int & blocksize){
    QVector<NeuronSWC> Extend_Marker;
    const int & X=centerAPO.x;
    const int & Y=centerAPO.y;
    const int & Z=centerAPO.z;
    for(const NeuronSWC & swc:Ans_Tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            for(const int & to:Answer_Graph[swc.n]){
                if(to<=0) continue;
                const NeuronSWC & check = Answer_Map[to];
                if(!between(check.x,X-blocksize/2,X+blocksize/2)||!between(check.y,Y-blocksize/2,Y+blocksize/2)||!between(check.z,Z-blocksize/2,Z+blocksize/2)){
                    NeuronSWC Add_Point=swc;
                    Extend_Marker.push_back(Add_Point);
                    break;
                }
            }
        }
    }
    int delete_id=-1;
    double mn=1e8;
    for(int i=0;i<Extend_Marker.size();++i){
        double dis=distance_XYZ(XYZ(Extend_Marker[i]),XYZ(Absolute_Marker));
        if(mn>dis){
            mn=dis;
            delete_id=i;
        }
    }
    if(!Extend_Marker.empty() && delete_id<Extend_Marker.size() && delete_id!=-1)
      Extend_Marker.erase(Extend_Marker.begin()+delete_id);
    return Extend_Marker;
}

bool Make_Dir(const QString & dirname){
    QDir dir(dirname);
    if(!dir.exists()){
        bool if_success = dir.mkdir(dirname);
        if(!if_success){
            qDebug()<<"mkdir "<<dirname<< " fails";
            return false;
        }
    }
    qDebug()<<dirname<<" is already";
    return true;
}

void init(const QString & Answer_File){

    while(!Make_Dir(Work_Dir+QString("/APOFile")));
    while(!Make_Dir(Work_Dir+QString("/MarkerFile")));
    while(!Make_Dir(Work_Dir+QString("/SwcFile")));
    while(!Make_Dir(Work_Dir+QString("/testV3draw")));
    while(!Make_Dir(Work_Dir+QString("/EswcFile")));

    Ans_Tree=readSWC_file(Answer_File);
    for(const NeuronSWC & i:Ans_Tree.listNeuron){
        if(i.type==1){
            X=i.x;
            Y=i.y;
            Z=i.z;
        }
        Answer_Map[i.n]=i;
        if(i.pn>0 && i.n>0){
            Answer_Graph[i.n].push_back(i.pn);
            Answer_Graph[i.pn].push_back(i.n);
        }
    }
}

NeuronTree Vanilla_App2(const CellAPO & centerAPO,const ImageMarker & startPoint,const int & blocksize,const QString & rawFileName){
    //drop marker file(ok)
    QList<ImageMarker> List_Marker_Write;
    List_Marker_Write.push_back(startPoint);
    //Marker is changed to absolute location
    //location transform has been checked(ok)
    ImageMarker Absolute_Marker=startPoint;
    Absolute_Marker.x+=centerAPO.x-blocksize/2;
    Absolute_Marker.y+=centerAPO.y-blocksize/2;
    Absolute_Marker.z+=centerAPO.z-blocksize/2;
    QString Marker_File_Name=generate_marker_name(Work_Dir.toStdString()+"/MarkerFile",Absolute_Marker);	//make file in ./MarkerFile/xxx.000_xxx.000_xxx.000.marker
    writeMarker_file(Marker_File_Name,List_Marker_Write);

  //app2
   QString App2_Eswc_File_Name=generate_eswc_name(Work_Dir.toStdString()+"/SwcFile",centerAPO);
   qDebug()<<"app2:"
          <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"), QStringList()
              <<"/x"<<Vaa3d_App_Path+QString("/plugins/neuron_tracing/Vaa3D_Neuron2/vn2.dll")
              <<"/f"<<"app2"<<"/p"<<Marker_File_Name<<QString::number(0)<<QString::number(-1)
              <<"/i"<< QString(Work_Dir+QString("/testV3draw/thres_")+rawFileName)<<"/o"<<App2_Eswc_File_Name);

   return readSWC_file(App2_Eswc_File_Name);
}

NeuronTree Find_Valid_App2_Tree(const CellAPO & centerAPO,ImageMarker & startPoint,const int & blocksize,const QString & rawFileName){
    const int range=1;
    QVector<int> d={0};
    for(int i=1;i<=range;++i){
        d.push_back(2*i);
        d.push_back(-2*i);
    }
    NeuronTree ret;
    std::map<std::tuple<int,int,int>,bool> mp;
    if(startPoint.x==blocksize/2 && startPoint.y==blocksize/2 &&startPoint.z==blocksize/2){
        for(int i=0;i<d.size();++i){
            for(int j=0;j<d.size();++j){
                for(int k=0;k<d.size();++k){
                    if(mp.count(std::tuple<int,int,int> (d[i],d[j],d[k]))) continue;
                    mp[std::tuple<int,int,int> (d[i],d[j],d[k])]=true;
                    ImageMarker new_marker=startPoint;
                    new_marker.x+=d[i];
                    new_marker.y+=d[j];
                    new_marker.z+=d[k];
                    ret=Vanilla_App2(centerAPO,new_marker,blocksize,rawFileName);
                    if(!ret.listNeuron.empty()){
                        startPoint=new_marker;
                        break;
                    }
                }
                if(!ret.listNeuron.empty()){
                    break;
                }
            }
            if(!ret.listNeuron.empty()){
                break;
            }
        }
        return ret;
    }
    if(startPoint.x!=blocksize/2){
        for(int i=0;i<d.size();++i){
            for(int j=0;j<=i;++j){
                const int & d1=d[i];
                const int & d2=d[j];
                ImageMarker new_marker=startPoint;
                new_marker.y+=d1;
                new_marker.z+=d2;
                ret=Vanilla_App2(centerAPO,new_marker,blocksize,rawFileName);
                if(!ret.listNeuron.empty()){
                    startPoint=new_marker;
                    break;
                }
                if(i==j) continue;
                new_marker=startPoint;
                new_marker.z+=d1;
                new_marker.y+=d2;
                ret=Vanilla_App2(centerAPO,new_marker,blocksize,rawFileName);
                if(!ret.listNeuron.empty()){
                    startPoint=new_marker;
                    break;
                }
            }
            if(!ret.listNeuron.empty()){
                break;
            }
        }
    }else{
        if(startPoint.y!=blocksize/2){
            for(int i=0;i<d.size();++i){
                for(int j=0;j<=i;++j){
                    const int & d1=d[i];
                    const int & d2=d[j];
                    ImageMarker new_marker=startPoint;
                    new_marker.x+=d1;
                    new_marker.z+=d2;
                    ret=Vanilla_App2(centerAPO,new_marker,blocksize,rawFileName);
                    if(!ret.listNeuron.empty()){
                        startPoint=new_marker;
                        break;
                    }
                    if(i==j) continue;
                    new_marker=startPoint;
                    new_marker.z+=d1;
                    new_marker.x+=d2;
                    ret=Vanilla_App2(centerAPO,new_marker,blocksize,rawFileName);
                    if(!ret.listNeuron.empty()){
                        startPoint=new_marker;
                        break;
                    }
                }
                if(!ret.listNeuron.empty()){
                    break;
                }
            }
        }else{
            if(startPoint.z!=blocksize/2){
                for(int i=0;i<d.size();++i){
                    for(int j=0;j<=i;++j){
                        const int & d1=d[i];
                        const int & d2=d[j];
                        ImageMarker new_marker=startPoint;
                        new_marker.x+=d1;
                        new_marker.y+=d2;
                        ret=Vanilla_App2(centerAPO,new_marker,blocksize,rawFileName);
                        if(!ret.listNeuron.empty()){
                            startPoint=new_marker;
                            break;
                        }
                        if(i==j) continue;
                        new_marker=startPoint;
                        new_marker.y+=d1;
                        new_marker.x+=d2;
                        ret=Vanilla_App2(centerAPO,new_marker,blocksize,rawFileName);
                        if(!ret.listNeuron.empty()){
                            startPoint=new_marker;
                            break;
                        }
                    }
                    if(!ret.listNeuron.empty()){
                        break;
                    }
                }
            }
        }
    }

   return ret;
}

unsigned Get_Marker_Idx(const int & x,const int & y,const int & z,const int & blocksize){
    return z*blocksize*blocksize+y*blocksize+x;
}

unsigned Get_Marker_Idx(const ImageMarker & startPoint,const int & blocksize){
    return Get_Marker_Idx(startPoint.x,startPoint.y,startPoint.z,blocksize);
}

struct marker_with_intensity{
    ImageMarker marker;
    unsigned char intensity;
    marker_with_intensity(const ImageMarker & marker_,const unsigned char & intensity_):marker(marker_),intensity(intensity_){}
    bool operator < (const marker_with_intensity & cmp) const{
        return intensity<cmp.intensity;
    }
};

std::pair<QVector<ImageMarker>,QVector<ImageMarker> > Get_Valid_Marker(const QString & v3draw_File,const ImageMarker & startPoint,const int & blocksize){
    Image4DSimple p4dImage;
    p4dImage.loadImage(v3draw_File.toStdString().c_str(),true);
    unsigned char * indata1d = p4dImage.getRawDataAtChannel(0);

    std::pair<QVector<ImageMarker>,QVector<ImageMarker> > ret;

    const int range=3;
    QVector<int> d={0};
    for(int i=1;i<=range;++i){
        d.push_back(i);
        d.push_back(-i);
    }
    std::priority_queue<marker_with_intensity> q;
    std::map<std::tuple<int,int,int>,bool> mp;

    for(int i=0;i<d.size();++i){
        for(int j=0;j<d.size();++j){
            for(int k=0;k<d.size();++k){
                if(mp.count(std::tuple<int,int,int> (d[i],d[j],d[k]))) continue;
                mp[std::tuple<int,int,int> (d[i],d[j],d[k])]=true;
                ImageMarker new_marker=startPoint;
                new_marker.x+=d[i];
                new_marker.y+=d[j];
                new_marker.z+=d[k];
                if(new_marker.x<=0||new_marker.x>=blocksize||new_marker.y<=0||new_marker.y>=blocksize||new_marker.z<=0||new_marker.z>=blocksize) continue;
                unsigned idx=Get_Marker_Idx(new_marker,blocksize);
                q.push(marker_with_intensity(new_marker,indata1d[idx]));
            }
        }
    }
    QVector<ImageMarker> ret1;
    for(int i=0;i<5;++i){
        ret1.push_back(q.top().marker);
        q.pop();
    }

    if(startPoint.x==blocksize/2 && startPoint.y==blocksize/2 &&startPoint.z==blocksize/2){
        return std::make_pair(ret1,QVector<ImageMarker> () );
    }

//    q=std::priority_queue<marker_with_intensity>();
//    if(startPoint.x!=blocksize/2){
//        for(int i=0;i<d.size();++i){
//            for(int j=0;j<=i;++j){
//                const int & d1=d[i];
//                const int & d2=d[j];
//                ImageMarker new_marker=startPoint;
//                new_marker.y+=d1;
//                new_marker.z+=d2;
//                unsigned idx=Get_Marker_Idx(new_marker,blocksize);
//                q.push(marker_with_intensity(new_marker,indata1d[idx]));
//                if(i==j) continue;
//                new_marker=startPoint;
//                new_marker.z+=d1;
//                new_marker.y+=d2;
//                idx=Get_Marker_Idx(new_marker,blocksize);
//                q.push(marker_with_intensity(new_marker,indata1d[idx]));

//            }
//        }
//    }else{
//        if(startPoint.y!=blocksize/2){
//            for(int i=0;i<d.size();++i){
//                for(int j=0;j<=i;++j){
//                    const int & d1=d[i];
//                    const int & d2=d[j];
//                    ImageMarker new_marker=startPoint;
//                    new_marker.x+=d1;
//                    new_marker.z+=d2;
//                    unsigned idx=Get_Marker_Idx(new_marker,blocksize);
//                    q.push(marker_with_intensity(new_marker,indata1d[idx]));
//                    if(i==j) continue;
//                    new_marker=startPoint;
//                    new_marker.z+=d1;
//                    new_marker.x+=d2;
//                    idx=Get_Marker_Idx(new_marker,blocksize);
//                    q.push(marker_with_intensity(new_marker,indata1d[idx]));

//                }
//            }
//        }else{
//            if(startPoint.z!=blocksize/2){
//                for(int i=0;i<d.size();++i){
//                    for(int j=0;j<=i;++j){
//                        const int & d1=d[i];
//                        const int & d2=d[j];
//                        ImageMarker new_marker=startPoint;
//                        new_marker.x+=d1;
//                        new_marker.y+=d2;
//                        unsigned idx=Get_Marker_Idx(new_marker,blocksize);
//                        q.push(marker_with_intensity(new_marker,indata1d[idx]));
//                        if(i==j) continue;
//                        new_marker=startPoint;
//                        new_marker.y+=d1;
//                        new_marker.x+=d2;
//                        idx=Get_Marker_Idx(new_marker,blocksize);
//                        q.push(marker_with_intensity(new_marker,indata1d[idx]));

//                    }
//                }
//            }
//        }
//    }
//    for(int i=0;i<5;++i){
//        ret1.push_back(q.top().marker);
//        q.pop();
//    }

    const int low=4;
    const int high=10;
    d.clear();
    for(int i=low;i<=high;++i){
        d.push_back(i);
        d.push_back(-i);
    }
    q=std::priority_queue<marker_with_intensity>();
    if(startPoint.x!=blocksize/2){
        for(int i=0;i<d.size();++i){
            for(int j=0;j<=i;++j){
                const int & d1=d[i];
                const int & d2=d[j];
                ImageMarker new_marker=startPoint;
                new_marker.y+=d1;
                new_marker.z+=d2;
                unsigned idx=Get_Marker_Idx(new_marker,blocksize);
                q.push(marker_with_intensity(new_marker,indata1d[idx]));
                if(i==j) continue;
                new_marker=startPoint;
                new_marker.z+=d1;
                new_marker.y+=d2;
                idx=Get_Marker_Idx(new_marker,blocksize);
                q.push(marker_with_intensity(new_marker,indata1d[idx]));

            }
        }
    }else{
        if(startPoint.y!=blocksize/2){
            for(int i=0;i<d.size();++i){
                for(int j=0;j<=i;++j){
                    const int & d1=d[i];
                    const int & d2=d[j];
                    ImageMarker new_marker=startPoint;
                    new_marker.x+=d1;
                    new_marker.z+=d2;
                    unsigned idx=Get_Marker_Idx(new_marker,blocksize);
                    q.push(marker_with_intensity(new_marker,indata1d[idx]));
                    if(i==j) continue;
                    new_marker=startPoint;
                    new_marker.z+=d1;
                    new_marker.x+=d2;
                    idx=Get_Marker_Idx(new_marker,blocksize);
                    q.push(marker_with_intensity(new_marker,indata1d[idx]));

                }
            }
        }else{
            if(startPoint.z!=blocksize/2){
                for(int i=0;i<d.size();++i){
                    for(int j=0;j<=i;++j){
                        const int & d1=d[i];
                        const int & d2=d[j];
                        ImageMarker new_marker=startPoint;
                        new_marker.x+=d1;
                        new_marker.y+=d2;
                        unsigned idx=Get_Marker_Idx(new_marker,blocksize);
                        q.push(marker_with_intensity(new_marker,indata1d[idx]));
                        if(i==j) continue;
                        new_marker=startPoint;
                        new_marker.y+=d1;
                        new_marker.x+=d2;
                        idx=Get_Marker_Idx(new_marker,blocksize);
                        q.push(marker_with_intensity(new_marker,indata1d[idx]));

                    }
                }
            }
        }
    }
    QVector<ImageMarker> ret2;
    for(int i=0;i<5;++i){
        ret2.push_back(q.top().marker);
        q.pop();
    }
    ret.first=ret1;
    ret.second=ret2;
    return ret;
}

bool if_complicated(const CellAPO & centerAPO,const int & blocksize){
    int st_id=-1;
    QMap<int,bool> Inbox;
    const int & X=centerAPO.x;
    const int & Y=centerAPO.y;
    const int & Z=centerAPO.z;
    for(const NeuronSWC & swc:Ans_Tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            Inbox[swc.n]=true;
            st_id=swc.n;
        }
    }
    if(st_id==-1) return true;
    QQueue<int> q;
    q.push_back(st_id);
    QMap<int,bool> vis;
    while(!q.empty()){
        int now=q.front();
        q.pop_front();
        vis[now]=true;
        for(const int & to:Answer_Graph[now]){
            if(!vis[to]){
                q.push_back(to);
            }
        }
    }
    for(auto i=Inbox.begin();i!=Inbox.end();++i){
        if(!vis.count(i.key())){
            return true;
        }
    }
    return false;
}

NeuronTree Get_Everything_Inbox(const CellAPO & centerAPO,const int & blocksize){
    NeuronTree ret;
    for(const NeuronSWC & swc:Ans_Tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            ret.listNeuron.push_back(swc);
        }
    }
    return ret;
}

void App2_non_recursive_DFS(const int & Start_x,const int & Start_y,const int & Start_z,const int & blocksize,const QString & File_Name){

    init(Answer_File);

    CellAPO init_centerAPO;
    init_centerAPO.x=Start_x;
    init_centerAPO.y=Start_y;
    init_centerAPO.z=Start_z;

    ImageMarker init_startPoint;
    init_startPoint.x=blocksize/2;
    init_startPoint.y=blocksize/2;
    init_startPoint.z=blocksize/2;

    V_NeuronSWC_list App2_Generate;

    int amount=0;
    int not_change=0;
    int max_not_change=20;
    unsigned long long v3draw_size=-1;
    QQueue<bbox_extend> bbox_queue;
    QMap<int,bool> has_extend;
    int rt_id=1;
    for(const NeuronSWC & i:Ans_Tree.listNeuron){
        if(i.type==1){
            rt_id=i.n;
            break;
        }
    }
    bbox_queue.push_back(bbox_extend(0,init_centerAPO,init_startPoint,rt_id));

    while(!bbox_queue.empty()||!if_finish()){
        if(not_change==max_not_change){
            not_change=0;
            bbox_queue.clear();
            max_not_change=5;
        }
        if(bbox_queue.empty()){
            used_swc.clear();
//            vis.clear();
            has_extend.clear();
            amount=0;
            int id=unused_id();
            CellAPO centerAPO;
            const NeuronSWC & swc=Answer_Map[id];
            update_Ans_used(swc);
            centerAPO.x=int(swc.x+0.5);
            centerAPO.y=int(swc.y+0.5);
            centerAPO.z=int(swc.z+0.5);
            ImageMarker startPoint;
            startPoint.x=blocksize/2;
            startPoint.y=blocksize/2;
            startPoint.z=blocksize/2;
            bbox_queue.push_back(bbox_extend(0,centerAPO,startPoint,swc.n));
            continue;
        }
        bbox_extend now=bbox_queue.back();
        bbox_queue.pop_back();
        CellAPO centerAPO=now.centerAPO;
        ImageMarker startPoint=now.startPoint;
        int center_id=now.center_id;
        int centerSWC=now.centerSWC;
        if(has_extend.count(center_id)) continue;
        if(if_finish()) continue;
        if(!if_need_extend(centerAPO,blocksize)) continue;

        qDebug()<<"start DFS,center_id="<<center_id;

//        node now_node=node(centerAPO.x,centerAPO.y,centerAPO.z,startPoint.x,startPoint.y,startPoint.z);

        NeuronTree App2_Tree=NeuronTree();
        bool use_answer=false;
        NeuronTree Answer_Tree=NeuronTree();
        V_NeuronSWC_list V_Answer_Tree;


        ImageMarker Absolute_Marker=startPoint;
        Absolute_Marker.x+=centerAPO.x-blocksize/2;
        Absolute_Marker.y+=centerAPO.y-blocksize/2;
        Absolute_Marker.z+=centerAPO.z-blocksize/2;

        //if the square of distance between used_swc and new_border_point is lower than identical_threshold, don't search it
        QVector<NeuronSWC> Border_Points=Find_Extend_Marker(centerAPO,Absolute_Marker,blocksize);

        if(Border_Points.size()>=3){
            use_answer=true;
            Answer_Tree=Get_Ans_In_BBox(centerAPO,blocksize);
            V_Answer_Tree=NeuronTree__2__V_NeuronSWC_list(Answer_Tree);
        }
        else{
            //drop apo file(ok)
            QList<CellAPO> List_APO_Write;
            List_APO_Write.push_back(centerAPO);
            QString APO_File_Name=generate_apo_name(Work_Dir.toStdString()+"/APOFile",centerAPO);	//make file in ./APOFile/xxx.000_xxx.000_xxx.000.apo
            writeAPO_file(APO_File_Name,List_APO_Write);


            QString rawFileName=QString("%1.000_%2.000_%3.000.v3draw").arg(centerAPO.x).arg(centerAPO.y).arg(centerAPO.z);

            //crop3D(ok)
            qDebug()<<"crop3D:"
               <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"),QStringList()
                   <<"/x"<<Vaa3d_App_Path+QString("/plugins/image_geometry/crop3d_image_series/cropped3DImageSeries.dll")
                   <<"/f"<<"cropTerafly"<<"/i"<<Res_Path<<APO_File_Name<<Work_Dir+QString("/testV3draw/")
                   <<"/p"<<QString::number(blocksize)<<QString::number(blocksize)<<QString::number(blocksize));

            QFile crop_file(Work_Dir+QString("/testV3draw/")+rawFileName);
            if(v3draw_size==-1){
                v3draw_size=crop_file.size();
            }
            else {
                if(crop_file.size()!=v3draw_size){
                    NeuronTree Ans_In_BBox=Get_Ans_In_BBox(centerAPO,blocksize);
                    V_NeuronSWC_list Segments=NeuronTree__2__V_NeuronSWC_list(Ans_In_BBox);
                    for(const V_NeuronSWC & Seg:Segments.seg){
                        V_NeuronSWC Add_Seg=Seg;
                        QMap<int,bool> redundant;
                        for(int i=0; i<Add_Seg.row.size();++i){
                            //(ok)
                            V_NeuronSWC_unit & swc=Add_Seg.row[i];
                            swc.x+=centerAPO.x-blocksize/2;
                            swc.y+=centerAPO.y-blocksize/2;
                            swc.z+=centerAPO.z-blocksize/2;
                            if(App2_Generate.seg.empty()) continue;
                            if(Distance_Unit_To_Tree(swc,App2_Generate)<close_distance){
                                redundant[i]=true;
                            }
                        }

                        for(int i=0;i<Add_Seg.row.size();++i){
                            if(!redundant.count(i)){
                                V_NeuronSWC Not_Redundant;
                                int now=i;
                                int amount=0;
                                while(now<Add_Seg.row.size() && !redundant.count(now)){
                                    V_NeuronSWC_unit unit=Add_Seg.row[now];
                                    ++amount;
                                    unit.type=3;
                                    unit.n=amount;
                                    unit.parent=amount+1;
                                    Not_Redundant.row.push_back(unit);
                                    ++now;
                                }
                                Not_Redundant.row.back().parent=-1;
                                i=now;
                                if(Not_Redundant.row.size()<=2) continue;
                                App2_Generate.append(Not_Redundant);

                            }
                        }
                    }
                    update_Ans_used(Ans_In_BBox);
                    QFile::remove(Work_Dir+QString("/testV3draw/")+rawFileName);
                    continue;
                }
            }

            //ada(ok)
            qDebug()<<"ada:"
               <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"),QStringList()
                   <<"/x"<<Vaa3d_App_Path+QString("/plugins/image_thresholding/Simple_Adaptive_Thresholding/ada_threshold.dll")
                   <<"/f"<<"adath"<<"/i"<<Work_Dir+QString("/testV3draw/")+rawFileName<<"/o"<<QString(Work_Dir+QString("/testV3draw/thres_")+rawFileName));

            ImageMarker startPointBackup=startPoint;
//            Answer_Tree = Get_Answer_Tree(center_id,centerAPO,startPoint,blocksize,centerSWC);
            Answer_Tree = Get_Ans_In_BBox(centerAPO,blocksize);
            if(Answer_Tree.listNeuron.empty()) continue;
            V_Answer_Tree = NeuronTree__2__V_NeuronSWC_list(Answer_Tree);

            if(center_id!=0){
                QString v3draw=Work_Dir+QString("/testV3draw/thres_")+rawFileName;
                std::pair<QVector<ImageMarker>,QVector<ImageMarker> > possible=Get_Valid_Marker(v3draw,startPoint,blocksize);
                for(const ImageMarker & i:possible.first){
                    startPoint=i;
                    App2_Tree=Vanilla_App2(centerAPO,startPoint,blocksize,rawFileName);
                    if(!App2_Tree.listNeuron.empty()){
                        V_NeuronSWC_list V_App2_Tree = NeuronTree__2__V_NeuronSWC_list(App2_Tree);
                        if(!Check_Tree_Identical(V_App2_Tree,V_Answer_Tree)){
                            App2_Tree=NeuronTree();
                            continue;
                        }
                        break;
                    }
                }
                if(App2_Tree.listNeuron.empty()){
                    for(const ImageMarker & i:possible.second){
                        startPoint=i;
                        App2_Tree=Vanilla_App2(centerAPO,startPoint,blocksize,rawFileName);
                        if(!App2_Tree.listNeuron.empty()){
                            V_NeuronSWC_list V_App2_Tree = NeuronTree__2__V_NeuronSWC_list(App2_Tree);
                            if(!Check_Tree_Identical(V_App2_Tree,V_Answer_Tree)){
                                App2_Tree=NeuronTree();
                                continue;
                            }
                            break;
                        }
                    }
                }
    //            if(App2_Tree.listNeuron.empty()){
    //                 App2_Tree=Find_Valid_App2_Tree(centerAPO,startPoint,blocksize,v3draw);
    //            }

               QFile::remove(Work_Dir+QString("/testV3draw/")+rawFileName);
               QFile::remove(Work_Dir+QString("/testV3draw/thres_")+rawFileName);
            }
            else {
                QString v3draw=Work_Dir+QString("/testV3draw/thres_")+rawFileName;
                std::pair<QVector<ImageMarker>,QVector<ImageMarker> > possible=Get_Valid_Marker(v3draw,startPoint,blocksize);
                for(const ImageMarker & i:possible.first){
                    startPoint=i;
                    App2_Tree=Vanilla_App2(centerAPO,startPoint,blocksize,rawFileName);
                    if(!App2_Tree.listNeuron.empty()){
                        V_NeuronSWC_list V_App2_Tree = NeuronTree__2__V_NeuronSWC_list(App2_Tree);
                        if(!Check_Tree_Identical(V_App2_Tree,V_Answer_Tree)){
                            App2_Tree=NeuronTree();
                            continue;
                        }
                        break;
                    }
                }
                if(App2_Tree.listNeuron.empty()){
                    startPoint=startPointBackup;
                    App2_Tree=Vanilla_App2(centerAPO,startPoint,blocksize,rawFileName);
                }

                QFile::remove(Work_Dir+QString("/testV3draw/")+rawFileName);
                QFile::remove(Work_Dir+QString("/testV3draw/thres_")+rawFileName);
            }

            if(App2_Tree.listNeuron.empty())    startPoint=startPointBackup;

        }

       if(App2_Tree.listNeuron.empty()){
           use_answer=true;
       }

       QMap<int,QVector<int> > son;    //record son
       QMap<int,NeuronSWC> mp;         //map id to single_swc
       QVector<std::pair<NeuronSWC,QQueue<XYZ> > > Border_Point_Vector;
       if(!use_answer){
           //find_border_point
           for(const NeuronSWC & swc:App2_Tree.listNeuron){
               son[swc.pn].push_back(swc.n);
               mp[swc.n]=swc;
           }
           Border_Point_Vector=Find_Border(App2_Tree,blocksize,son,mp);

           //we think the furthest border point is not accurate
           //cut the border with one short line, we don't want to add these points to the answer
           //App2_Tree will be cut, Border_Point_Vector will not
           //so delete them(ok)
           QVector<int> need_cut;
           for(int i=App2_Tree.listNeuron.size()-1;i>=0;i--){
               for(const std::pair<NeuronSWC,QQueue<XYZ> > & unsecure_point_vector:Border_Point_Vector){
                   const NeuronSWC & unsecure_point=unsecure_point_vector.first;
                   if(App2_Tree.listNeuron[i].n==unsecure_point.n){
                       need_cut.push_back(i);
                       break;
                   }
               }
           }
           //the delete operation has been checked in venilla c++(ok)
           for(const int & i:need_cut){
               App2_Tree.listNeuron.erase(App2_Tree.listNeuron.begin()+i);
           }

           if(center_id>0){
               V_NeuronSWC_list V_App2_Tree = NeuronTree__2__V_NeuronSWC_list(App2_Tree);
               //app2 is not accurate
               if((V_Answer_Tree.seg.size()==1?3:V_Answer_Tree.seg.size()*2)<V_App2_Tree.seg.size()||!Check_Tree_Identical(V_App2_Tree,V_Answer_Tree)||!Check_Tree_Identical(V_Answer_Tree,V_App2_Tree)){
                   use_answer=true;
                   App2_Tree=Answer_Tree;
               }
           }
       }
       else{
           App2_Tree=Answer_Tree;
       }

       V_NeuronSWC_list Segments=NeuronTree__2__V_NeuronSWC_list(App2_Tree);
       if(!use_answer){
           QVector<XYZ> Ans_Vec;
           for(const V_NeuronSWC & Seg:V_Answer_Tree.seg){
               if(Seg.row.size()>1){
                   XYZ v1;
                   v1.x=Seg.row.front().x;
                   v1.y=Seg.row.front().y;
                   v1.z=Seg.row.front().z;
                   XYZ v2;
                   v2.x=Seg.row.back().x;
                   v2.y=Seg.row.back().y;
                   v2.z=Seg.row.back().z;
                   Ans_Vec.push_back(v1-v2);
               }
           }
         for(const V_NeuronSWC & Seg:Segments.seg){
             if(!has_same_vector(Seg,Ans_Vec)) {
                 use_answer=true;
                 App2_Tree=Answer_Tree;
                 Segments=NeuronTree__2__V_NeuronSWC_list(App2_Tree);
                 break;
             }
         }
       }

       //save the answer(ok)
       bool if_change=false;
       for(const V_NeuronSWC & Seg:Segments.seg){
           V_NeuronSWC Add_Seg=Seg;
           QMap<int,bool> redundant;
           for(int i=0; i<Add_Seg.row.size();++i){
               //(ok)
               V_NeuronSWC_unit & swc=Add_Seg.row[i];
               swc.x+=centerAPO.x-blocksize/2;
               swc.y+=centerAPO.y-blocksize/2;
               swc.z+=centerAPO.z-blocksize/2;
               if(App2_Generate.seg.empty()) continue;
               if(Distance_Unit_To_Tree(swc,App2_Generate)<close_distance){
                   redundant[i]=true;
               }
           }

           for(int i=0;i<Add_Seg.row.size();++i){
               if(!redundant.count(i)){
                   V_NeuronSWC Not_Redundant;
                   int now=i;
                   int amount=0;
                   while(now<Add_Seg.row.size() && !redundant.count(now)){
                       V_NeuronSWC_unit unit=Add_Seg.row[now];
                       ++amount;
                       if(use_answer) unit.type=3;
                       else unit.type=2;
                       unit.n=amount;
                       unit.parent=amount+1;
                       Not_Redundant.row.push_back(unit);
                       ++now;
                   }
                   Not_Redundant.row.back().parent=-1;
                   i=now;
                   if(Not_Redundant.row.size()<=2) continue;
                   if_change=true;
                   App2_Generate.append(Not_Redundant);

               }
           }
       }

       if(!if_change)  ++not_change;
       else not_change=0;


       NeuronTree output=V_NeuronSWC_list__2__NeuronTree(App2_Generate);
       writeSWC_file(Work_Dir+"/EswcFile/"+QString(QString::fromStdString(std::to_string(++cnt)))+".eswc",output);


       //check if border is identical to last marker
       NeuronSWC Start_Marker_Location;
       //marker has been changed to absolute location
       Start_Marker_Location.x=Absolute_Marker.x;
       Start_Marker_Location.y=Absolute_Marker.y;
       Start_Marker_Location.z=Absolute_Marker.z;
       used_swc.push_back(Start_Marker_Location);

//       if(Border_Points.empty()) continue;
//       if(vis.count(now_node)) continue;

       update_Ans_used(Answer_Tree);
//       if(use_answer && center_id==0){
//           update_Ans_used(Get_Ans_In_BBox(centerAPO,blocksize));
//       }

//       vis[now_node]=true;


       for(const NeuronSWC & Border_Point:Border_Points){//absolute
           bool used=false;
           for(const NeuronSWC & swc:used_swc){
               if(distance_square(swc,Border_Point)<identical_threshold){
                   used=true;
               }
           }
           if(used) continue;

           XYZ vc=(XYZ(Border_Point)-XYZ(Absolute_Marker));

           std::vector<int> direction=Judge_Direction(XYZ(Border_Point)-XYZ(Absolute_Marker));

           //calculate the new center of next round DFS
           CellAPO New_Point;
           int offset=blocksize/2-2;
           //the location of Absolute_Border_Location is float, round this float to int
           New_Point.x=int(Border_Point.x+0.5);
           New_Point.y=int(Border_Point.y+0.5);
           New_Point.z=int(Border_Point.z+0.5);

           //int dx[6]={1,-1,0,0,0,0};
           //int dy[6]={0,0,1,-1,0,0};
           //int dz[6]={0,0,0,0,1,-1};
           //direction=0 positive axis x       dx= 1, dy= 0, dz= 0
           //direction=1 negative axis x       dx=-1, dy= 0, dz= 0
           //direction=2 positive axis y       dx= 0, dy= 1, dz= 0
           //direction=3 negative axis y       dx= 0, dy=-1, dz= 0
           //direction=4 positive axis z       dx= 0, dy= 0, dz= 1
           //direction=5 negative axis z       dx= 0, dy= 0, dz=-1
           //move the center of next bounding_box forward in the accordingly direction, make the Border_Point locates in the center of area(面心)

           ++amount;
           for(int i=0;i<3;++i){
               CellAPO New_Point_Offset=New_Point;
               New_Point_Offset.x+=dx[direction[i]]*offset;
               New_Point_Offset.y+=dy[direction[i]]*offset;
               New_Point_Offset.z+=dz[direction[i]]*offset;
               //calculate the marker in the next bounding_box
               ImageMarker New_Marker;
               //Border_Point's relative position to the center of next bounding_box decrease or increase in the opposite direction
               New_Marker.x=blocksize/2-dx[direction[i]]*offset;
               New_Marker.y=blocksize/2-dy[direction[i]]*offset;
               New_Marker.z=blocksize/2-dz[direction[i]]*offset;

               int new_id=amount;
               has_extend[center_id]=true;
               bbox_queue.push_back(bbox_extend(new_id,New_Point_Offset,New_Marker,-1));
           }

       }
    }

    NeuronTree output=V_NeuronSWC_list__2__NeuronTree(App2_Generate);
    writeSWC_file(Work_Dir+QString("/")+File_Name,output);
}



int main(int argc, char **argv)
{
    int blocksize=128;
    QString Output_File_Name="whole_image.eswc";
    for(int i=0;i<argc;++i){
        if(i==1){
            Res_Path=QString::fromStdString((string(argv[i])));
        }
        if(i==2){
            Vaa3d_App_Path=QString::fromStdString((string(argv[i])));
        }
        if(i==3){
            Work_Dir=QString::fromStdString((string(argv[i])));
        }
        if(i==4){
            Answer_File=QString::fromStdString((string(argv[i])));
        }
    }

    App2_non_recursive_DFS(X,Y,Z,blocksize,Output_File_Name);


    qDebug()<<"finish";

}
