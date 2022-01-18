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

struct node{
    int x,y,z,marker_x,marker_y,marker_z;
    node(int x_,int y_,int z_,int marker_x_,int marker_y_,int marker_z_):x(x_),y(y_),z(z_),marker_x(marker_x_),marker_y(marker_y_),marker_z(marker_z_){}
    bool operator == (const node & cmp) const{
        return (x==cmp.x)&&(y==cmp.y)&&(z==cmp.z)&&(marker_x==cmp.marker_x)&&(marker_y==cmp.marker_y)&&(marker_z==cmp.marker_z);
    }
    bool operator < (const node & cmp) const{
        return x<cmp.x;
    }
};

struct bbox_extend{
    int center_id;
    CellAPO centerAPO;
    ImageMarker startPoint;
    bbox_extend(int _id_,CellAPO _centerAPO_,ImageMarker _startPoint_):center_id(_id_),centerAPO(_centerAPO_),startPoint(_startPoint_){}
};

QString Res_Path="G:/18454/RES(26298x35000x11041)";
QString Vaa3d_App_Path="C:/3.603c";
QString Work_Dir="G:/work_dir_new";
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
QMap<node,bool> vis;
QVector<NeuronSWC> used_swc;

int X=14530,Y=10693,Z=3124;
double close_distance=4;
double identical_threshold=300;   //(undetermined)
double seg_identical_threshold=500;
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


double distance_XYZ(const XYZ & p1,const XYZ & p2){
    XYZ vec=p2-p1;
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
    return p1+(p2-p1)*(dm/distance_XYZ(p1,p2));
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
bool Check_Seg_Identical(const V_NeuronSWC & Check_Seg,const V_NeuronSWC_list & Answer_Tree){
    double mx=-1e8;
    for(const V_NeuronSWC_unit & Check_Point:Check_Seg.row){
         for(const V_NeuronSWC & Answer_Seg:Answer_Tree.seg){
            mx=std::max(mx,Distance_Unit_To_Seg(Check_Point,Answer_Seg));
        }
    }
    qDebug()<<"check_mx:"<<mx;
    return mx<seg_identical_threshold;
}
bool Check_Tree_Identical(const V_NeuronSWC_list & Check_Tree,const V_NeuronSWC_list & Answer_Tree){
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
double Distance_Point_To_Border(const NeuronSWC & point,const int & blocksize){
    return std::min({point.x,abs(blocksize-point.x),point.y,abs(blocksize-point.y),point.z,abs(blocksize-point.z)});
}
double Vector_Angle(const XYZ & a,const XYZ & b){
    static const double PI=std::acos(-1);
    return std::acos(dot_mul(a,b)/distance_XYZ(XYZ(0.0),a)/distance_XYZ(XYZ(0.0),b))*(180/PI);
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
NeuronTree Get_Answer_Tree(const int & depth,const CellAPO & centerAPO,const ImageMarker & startPoint,const int & blocksize){
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
        int root_id=-1;
        for(const NeuronSWC & i:In_Block_Points){
            if(i.pn==-1){
                root_id=i.n;
                break;
            }
        }
        if(root_id==-1){
            root_id=Find_Nearest_Id(startPoint,In_Block_Points);
        }
        NeuronSWC StartPoint_In_Answer=mp[root_id];
        StartPoint_In_Answer.n=1;

        int amount=0;
        QQueue<std::pair<int,int> > q;
        QMap<int, bool> vis;
        vis[-1]=true;
        q.push_back(std::make_pair(-1,root_id));    //-1 is new, pick_id is old
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
//        QString output_file="G:/test.swc";
//        QFile file(output_file);
//        qDebug()<<file.open(QIODevice::WriteOnly|QIODevice::Text);

//        QTextStream out(&file);
//        for(NeuronSWC & swc:Return_Tree.listNeuron){
//            out<<swc.n<<" "<<swc.type<<" "<<swc.x<<" "<<swc.y<<" "<<swc.z<<" "<<swc.r<<" "<<swc.pn<<"\n";
//        }
//        file.close();
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
    const int threshold=1.0*300/256*blocksize;
    if(mn>threshold) return NeuronTree();

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

//    QString output_file=generate_swc_name("G:/ans_crop",centerAPO);
//    QFile file(output_file);
//    qDebug()<<file.open(QIODevice::WriteOnly|QIODevice::Text);

//    QTextStream out(&file);
//    for(NeuronSWC & swc:Return_Tree.listNeuron){
//        out<<swc.n<<" "<<swc.type<<" "<<swc.x<<" "<<swc.y<<" "<<swc.z<<" "<<swc.r<<" "<<swc.pn<<"\n";
//    }
//    file.close();
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
//                    Add_Point.x=Add_Point.x-X+blocksize/2;
//                    Add_Point.y=Add_Point.y-Y+blocksize/2;
//                    Add_Point.z=Add_Point.z-Z+blocksize/2;
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


  //app2(ok)
   QString App2_Eswc_File_Name=generate_eswc_name(Work_Dir.toStdString()+"/SwcFile",centerAPO);
   qDebug()<<"app2:"
          <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"), QStringList()
              <<"/x"<<Vaa3d_App_Path+QString("/plugins/neuron_tracing/Vaa3D_Neuron2/vn2.dll")
              <<"/f"<<"app2"<<"/p"<<Marker_File_Name<<QString::number(0)<<QString::number(-1)
              <<"/i"<< QString(Work_Dir+QString("/testV3draw/thres_")+rawFileName)<<"/o"<<App2_Eswc_File_Name);

   return readSWC_file(App2_Eswc_File_Name);
}

NeuronTree Find_Valid_App2_Tree(const CellAPO & centerAPO,ImageMarker & startPoint,const int & blocksize,const QString & rawFileName){
    const int range=3;
    QVector<int> d={0};
    for(int i=1;i<=range;++i){
        d.push_back(i);
        d.push_back(-i);
    }
    NeuronTree ret;
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

   QFile::remove(Work_Dir+QString("/testV3draw/")+rawFileName);
   QFile::remove(Work_Dir+QString("/testV3draw/thres_")+rawFileName);

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
    QQueue<bbox_extend> bbox_queue;
    QMap<int,bool> has_extend;
    bbox_queue.push_back(bbox_extend(0,init_centerAPO,init_startPoint));

    while(!bbox_queue.empty()){
        bbox_extend now=bbox_queue.back();
        bbox_queue.pop_back();
        CellAPO centerAPO=now.centerAPO;
        ImageMarker startPoint=now.startPoint;
        int center_id=now.center_id;
        if(has_extend.count(center_id)) continue;
        if(if_finish()) continue;
        if(!if_need_extend(centerAPO,blocksize)) continue;

        qDebug()<<"start DFS,center_id="<<center_id;

        node now_node=node(centerAPO.x,centerAPO.y,centerAPO.z,startPoint.x,startPoint.y,startPoint.z);
        if(vis.count(now_node)) continue;

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

        //ada(ok)
        qDebug()<<"ada:"
           <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"),QStringList()
               <<"/x"<<Vaa3d_App_Path+QString("/plugins/image_thresholding/Simple_Adaptive_Thresholding/ada_threshold.dll")
               <<"/f"<<"adath"<<"/i"<<Work_Dir+QString("/testV3draw/")+rawFileName<<"/o"<<QString(Work_Dir+QString("/testV3draw/thres_")+rawFileName));

        ImageMarker startPointBackup=startPoint;

        NeuronTree App2_Tree;
        if(center_id!=0)   App2_Tree=Find_Valid_App2_Tree(centerAPO,startPoint,blocksize,rawFileName);
        else {
            App2_Tree=Vanilla_App2(centerAPO,startPoint,blocksize,rawFileName);
            QFile::remove(Work_Dir+QString("/testV3draw/")+rawFileName);
            QFile::remove(Work_Dir+QString("/testV3draw/thres_")+rawFileName);
        }

        if(App2_Tree.listNeuron.empty())    startPoint=startPointBackup;

        ImageMarker Absolute_Marker=startPoint;
        Absolute_Marker.x+=centerAPO.x-blocksize/2;
        Absolute_Marker.y+=centerAPO.y-blocksize/2;
        Absolute_Marker.z+=centerAPO.z-blocksize/2;



       vis[now_node]=true;


       //readfile(ok)
       qDebug()<<"readfile";
       NeuronTree Answer_Tree = Get_Answer_Tree(center_id,centerAPO,startPoint,blocksize);
       update_Ans_used(Answer_Tree);
       bool use_answer_find_border=false;
       if(Answer_Tree.listNeuron.empty()) continue;
       if(App2_Tree.listNeuron.empty()){
           if(Answer_Tree.listNeuron.empty()){
               continue;
           } else{
               use_answer_find_border=true;
           }
       }
       qDebug()<<"readfile ok";

       QMap<int,QVector<int> > son;    //record son
       QMap<int,NeuronSWC> mp;         //map id to single_swc
       QVector<std::pair<NeuronSWC,QQueue<XYZ> > > Border_Point_Vector;
       if(!use_answer_find_border){
           //find_border_point
           for(const NeuronSWC & swc:App2_Tree.listNeuron){
               son[swc.pn].push_back(swc.n);
               mp[swc.n]=swc;
           }
           qDebug()<<"Find_Border";
           Border_Point_Vector=Find_Border(App2_Tree,blocksize,son,mp);
           qDebug()<<"Find_Border ok";

           //we think the furthest border point is not accurate
           //cut the border with one short line, we don't want to add these points to the answer
           //App2_Tree will be cut, Border_Point_Vector will not
           //so delete them(ok)
           qDebug()<<"cut";
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
           qDebug()<<"cut ok";
           //(new)
           if(center_id>0){
               if(Answer_Tree.listNeuron.empty()) continue;
               V_NeuronSWC_list V_Answer_Tree = NeuronTree__2__V_NeuronSWC_list(Answer_Tree);
               V_NeuronSWC_list V_App2_Tree = NeuronTree__2__V_NeuronSWC_list(App2_Tree);
               //app2 is not accurate
               if(!Check_Tree_Identical(V_Answer_Tree,V_App2_Tree)){
                   use_answer_find_border=true;
                   App2_Tree=Answer_Tree;
               }
           }
       }
       else{
           App2_Tree=Answer_Tree;
       }

       //save the answer(ok)
       qDebug()<<"start save";
       V_NeuronSWC_list Segments=NeuronTree__2__V_NeuronSWC_list(App2_Tree);
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
               qDebug()<<"Distance_Unit_To_Tree(swc,App2_Generate)"<<Distance_Unit_To_Tree(swc,App2_Generate);
               if(Distance_Unit_To_Tree(swc,App2_Generate)<close_distance){
                   redundant[i]=true;
               }
           }

           for(int i=0;i<Add_Seg.row.size();++i){
               qDebug()<<Add_Seg.row.size();
               if(!redundant.count(i)){
                   V_NeuronSWC Not_Redundant;
                   int now=i;
                   int amount=0;
                   while(now<Add_Seg.row.size() && !redundant.count(now)){
                       V_NeuronSWC_unit unit=Add_Seg.row[now];
                       ++amount;
                       if(use_answer_find_border) unit.type=3;
                       else unit.type=2;
                       unit.n=amount;
                       unit.parent=amount+1;
                       Not_Redundant.row.push_back(unit);
                       qDebug()<<unit.n<<" "<<unit.parent;
                       ++now;
                   }
                   Not_Redundant.row.back().parent=-1;
                   i=now;
                   if(Not_Redundant.row.size()<=2) continue;
                   App2_Generate.append(Not_Redundant);

               }
           }
       }

       NeuronTree output=V_NeuronSWC_list__2__NeuronTree(App2_Generate);
       writeSWC_file(Work_Dir+"/EswcFile/"+QString(QString::fromStdString(std::to_string(++cnt)))+".eswc",output);

       qDebug()<<"save ok";

       //if the square of distance between used_swc and new_border_point is lower than identical_threshold, don't search it
       qDebug()<<"use_answer_find_border:"<<use_answer_find_border;
       if(use_answer_find_border||Border_Point_Vector.empty()){
           QVector<NeuronSWC> Border_Points=Find_Extend_Marker(centerAPO,Absolute_Marker,blocksize);
           qDebug()<<"Find_Extend_Marker ok";

           //check if border is identical to last marker
           NeuronSWC Start_Marker_Location;
           //marker has been changed to absolute location
           Start_Marker_Location.x=Absolute_Marker.x;
           Start_Marker_Location.y=Absolute_Marker.y;
           Start_Marker_Location.z=Absolute_Marker.z;
           used_swc.push_back(Start_Marker_Location);

           qDebug()<<"Border_Point.size()"<<Border_Points.size();

           if(Border_Points.empty()) continue;

           for(const NeuronSWC & Border_Point:Border_Points){//absolute
               bool used=false;
               for(const NeuronSWC & swc:used_swc){
                   if(distance_square(swc,Border_Point)<identical_threshold){
                       used=true;
                   }
               }
               if(used) continue;

               XYZ vc=(XYZ(Border_Point)-XYZ(Absolute_Marker));

               qDebug()<<vc.x<<" "<<vc.y<<" "<<vc.z;

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

                   qDebug()<<"New_Point_Offset:"<<New_Point_Offset.x<<" "<<New_Point_Offset.y<<" "<<New_Point_Offset.z;
                   qDebug()<<"New_Marker:"<<New_Marker.x<<" "<<New_Marker.y<<" "<<New_Marker.z;


                   int new_id=amount;
                   has_extend[center_id]=true;
                   bbox_queue.push_back(bbox_extend(new_id,New_Point_Offset,New_Marker));
               }

           }
       }
       else{
           //treat with the border point
           for(std::pair<NeuronSWC,QQueue<XYZ> > & unsecure:Border_Point_Vector){

               //we think the father of the furthest inaccurate is real Border_Point
               //        const NeuronSWC & Border_Point=mp[unsecure.first.pn];

               const NeuronSWC & Border_Point=unsecure.first;


               //at first, judge if this point as been searched before
               //decide whether to search or not
               //Border_Point_Location transform checked(ok)
               NeuronSWC Absolute_Border_Location;
               Absolute_Border_Location.x=centerAPO.x-blocksize/2+Border_Point.x;
               Absolute_Border_Location.y=centerAPO.y-blocksize/2+Border_Point.y;
               Absolute_Border_Location.z=centerAPO.z-blocksize/2+Border_Point.z;
               //this operation succeed in avoiding redundancy to some extent
               bool used=false;
               for(const NeuronSWC & swc:used_swc){
                   if(distance_square(swc,Absolute_Border_Location)<identical_threshold){
                       used=true;
                   }
               }
               if(used) continue;

               //judge the direction to expand(undetermined)
               QQueue<XYZ> & Vectors=unsecure.second;
               //        Vectors.pop_back(); //don't take the last inaccuracy line into consideration
               //take the mean(ok)
               XYZ mean=XYZ(0.0,0.0,0.0);
               for(const XYZ & i:Vectors){
                   mean=mean+i;
               }
               int sz=Vectors.size();
               mean=mean/XYZ(sz);
               //judge direction by mean(recommend to read)
               std::vector<int> direction=Judge_Direction(mean);

               //calculate the new center of next round DFS
               CellAPO New_Point;
               int offset=blocksize/2-2;
               //the location of Absolute_Border_Location is float, round this float to int
               New_Point.x=int(Absolute_Border_Location.x+0.5);
               New_Point.y=int(Absolute_Border_Location.y+0.5);
               New_Point.z=int(Absolute_Border_Location.z+0.5);

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
               for(int i=0;i<5;++i){
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
                   qDebug()<<"New_Point_Offset:"<<New_Point_Offset.x<<" "<<New_Point_Offset.y<<" "<<New_Point_Offset.z;
                   qDebug()<<"New_Marker:"<<New_Marker.x<<" "<<New_Marker.y<<" "<<New_Marker.z;


                   int new_id=amount;
                   has_extend[center_id]=true;
                   bbox_queue.push_back(bbox_extend(new_id,New_Point_Offset,New_Marker));
               }
           }
       }
    }

    NeuronTree output=V_NeuronSWC_list__2__NeuronTree(App2_Generate);
    qDebug()<<"output.listNeuron.size():"<<output.listNeuron.size();
    writeSWC_file(Work_Dir+QString("/")+File_Name,output);
}

int main(int argc, char **argv)
{
    int blocksize=256;
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
        if(i==5){
            X=atoi(argv[i]);
        }
        if(i==6){
            Y=atoi(argv[i]);
        }
        if(i==7){
            Z=atoi(argv[i]);
        }
        if(i==8){
            blocksize=atoi(argv[i]);
        }
        if(i==9){
            close_distance=atof(argv[i]);
        }
        if(i==10){
            identical_threshold=atof(argv[i]);
        }
        if(i==11){
            seg_identical_threshold=atof(argv[i]);
        }
        if(i==12){
            Border_Threshold=atof(argv[i]);
        }
        if(i==13){
            Output_File_Name=QString::fromStdString((string(argv[i])));
        }
    }

    App2_non_recursive_DFS(X,Y,Z,blocksize,Output_File_Name);


    qDebug()<<"finish";

}
