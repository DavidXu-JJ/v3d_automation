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

Peng, H., Ruan, Z., Long, F., Simpson, J.H., and Myers, E.W. (2010) ÒV3D enables real-time 3D visualization and quantitative analysis of large-scale biological image data sets,Ó Nature Biotechnology, Vol. 28, No. 4, pp. 348-353, DOI: 10.1038/nbt.1612. ( http://penglab.janelia.org/papersall/docpdf/2010_NBT_V3D.pdf )

Peng, H, Ruan, Z., Atasoy, D., and Sternson, S. (2010) ÒAutomatic reconstruction of 3D neuron structures using a graph-augmented deformable model,Ó Bioinformatics, Vol. 26, pp. i38-i46, 2010. ( http://penglab.janelia.org/papersall/docpdf/2010_Bioinfo_GD_ISMB2010.pdf )

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

V3dApplication* V3dApplication::theApp=0;

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

struct bbox_extend{
    int center_id;
    CellAPO centerAPO;
    ImageMarker startPoint;
    int origin_id;
    bool pre_use_answer;
    NeuronTree pretree;
    bbox_extend(const int & _id_,const CellAPO & _centerAPO_,const ImageMarker & _startPoint_,const int  & _origin_,const bool & _pre_use_answer,const NeuronTree & _pretree_):
        center_id(_id_),centerAPO(_centerAPO_),startPoint(_startPoint_),origin_id(_origin_),pre_use_answer(_pre_use_answer),pretree(_pretree_){}
};

const double eps=1e-7;

bool double_eq(const double & a,const double & b){
    return std::abs(a-b)<eps;
}

bool XYZ_eq(const XYZ & p1,const XYZ & p2){
    return double_eq(p1.x,p2.x) && double_eq(p1.y,p2.y) && double_eq(p1.z,p2.z);
}

struct XYZ_Pair{
    XYZ p1;
    XYZ p2;
    XYZ_Pair(const XYZ & _p1,const XYZ & _p2):p1(_p1),p2(_p2){}
    bool operator == (const XYZ_Pair & cmp) const {
        return (XYZ_eq(p1,cmp.p1) && XYZ_eq(p2,cmp.p2))||(XYZ_eq(p1,cmp.p2) && XYZ_eq(p2,cmp.p1));
    }
    bool operator < (const XYZ_Pair & cmp) const {
        return p1.x<cmp.p1.x;
    }
};

void Delete_Redundant_Seg(V_NeuronSWC_list & segments){
    V_NeuronSWC_list lines;

    for(auto it=segments.seg.begin();it!=segments.seg.end();++it)
    {
        V_NeuronSWC & segment=*it;
        for(int i=1;i<segment.row.size();++i)
        {
            V_NeuronSWC newSeg;
            newSeg.append(segment.row.at(i-1));
            newSeg.append(segment.row.at(i));
            newSeg.row[0].n=1;newSeg.row[0].parent=2;
            newSeg.row[1].n=2;newSeg.row[1].parent=-1;
            lines.append(newSeg);
        }
    }

    static QVector<XYZ_Pair> vc;
    static int low=0;
    std::vector<int> del;
    for(int i=low;i<lines.seg.size();++i){
        const V_NeuronSWC_unit & u1=lines.seg[i].row[0];
        const V_NeuronSWC_unit & u2=lines.seg[i].row[1];
        XYZ p1,p2;
        p1.x=u1.x;
        p1.y=u1.y;
        p1.z=u1.z;
        p2.x=u2.x;
        p2.y=u2.y;
        p2.z=u2.z;
        XYZ_Pair pr=XYZ_Pair(p1,p2);
        bool if_del=false;
        for(const XYZ_Pair & cmp:vc){
            if(cmp==pr){
                if_del=true;
                break;
            }
        }
        if(if_del){
            del.push_back(i);
        }else{
            vc.push_back(pr);
        }
    }
    for(int i=del.size()-1;i>=0;--i){
        lines.seg.erase(lines.seg.begin()+del[i]);
    }
    low=lines.seg.size();

    segments=lines;
}

QString Res_Path="G:/18454/RES(26298x35000x11041)";
QString Vaa3d_App_Path="C:/3.603c";
QString Work_Dir="G:/demo";
QString Answer_File="G:/18454_answer/whole_image.eswc";
QMap<int,QVector<int> > Answer_Graph;
QMap<int,NeuronSWC> Answer_Map;
NeuronTree Ans_Tree;
QMap<int,bool> Ans_used;
QVector<NeuronSWC> Center_Set;

int cnt=0;
int app2_success=0;

QProcess p;

int dx[6]={1,-1,0,0,0,0};
int dy[6]={0,0,1,-1,0,0};
int dz[6]={0,0,0,0,1,-1};
//QMap<node,bool> vis;
QVector<NeuronSWC> used_swc;

int X=14530,Y=10693,Z=3124;
double close_distance=4;
double identical_threshold=300;   //(undetermined)
double seg_identical_threshold_mean=30;
double seg_identical_threshold_mx=50;
double Border_Threshold=50;
double length_indentical_rate=0.25;

double distance_square(const NeuronSWC & point_a,const NeuronSWC & point_b){
    return (point_a.x-point_b.x)*(point_a.x-point_b.x)+(point_a.y-point_b.y)*(point_a.y-point_b.y)+(point_a.z-point_b.z)*(point_a.z-point_b.z);
}

bool between(const double & mid,const double & left,const double & right){
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
            for(const NeuronSWC & swc:Center_Set){
                now_mn=std::min(now_mn,distance_square(swc,i));
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

bool if_need_extend(const NeuronTree & Tree){

    for(const NeuronSWC & swc:Tree.listNeuron){
        if(!Ans_used.count(swc.timestamp)){
            return true;
        }
    }
    return false;
}

void update_Ans_used(const NeuronTree & updateTree){
    for(const NeuronSWC & swc:updateTree.listNeuron)
        Ans_used[swc.timestamp]=true;
}
void update_Ans_used(const QVector<int> & Answer_Tree_Idx){
    for(const int & i:Answer_Tree_Idx){
        Ans_used[i]=true;
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
bool In_Box(const CellAPO & centerAPO,const int & blocksize,const NeuronSWC & point){
    return between(point.x,centerAPO.x-blocksize/2,centerAPO.x+blocksize/2)&&between(point.y,centerAPO.y-blocksize/2,centerAPO.y+blocksize/2)&&between(point.z,centerAPO.z-blocksize/2,centerAPO.z+blocksize/2);
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
double Distance_Unit_To_Seg(const V_NeuronSWC_unit & p, const V_NeuronSWC & s){  //??XYZ????
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
    return std::abs(len1-len2)/std::max(len1,len2)<length_indentical_rate;
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
NeuronTree Get_Answer_Tree(const int & origin_id,const CellAPO & centerAPO,const int & blocksize,QVector<int> & Answer_Tree_Idx){
    NeuronTree Return_Tree;
    int amount=0;
    QQueue<std::pair<int,int> > q;
    QMap<int, bool> vis;
    q.push_back(std::make_pair(-1,origin_id));    //-1 is new, pick_id is old
    while(!q.empty()){
        int prev=q.front().first;   //new
        int now=q.front().second;   //old
        q.pop_front();

        if(!In_Box(centerAPO,blocksize,Answer_Map[now])) continue;
        Answer_Tree_Idx.push_back(now);
        vis[now]=true;

        NeuronSWC Add_Point=Answer_Map[now];
        ++amount;
        Add_Point.n = amount;
        Add_Point.pn= prev;
        Add_Point.timestamp=now;
        Return_Tree.listNeuron.push_back(Add_Point);

        for(const int & i:Answer_Graph[now]){
            if(i>0 && !vis.count(i))
                q.push_back(std::make_pair(amount,i));
        }
    }
    return Return_Tree;
}

QVector<NeuronSWC>  Find_Extend_Marker(const int & center_id,const CellAPO & centerAPO, const int & blocksize, const NeuronTree & Answer_Tree){

    QVector<NeuronSWC> Extend_Marker;
    XYZ st;

    for(const NeuronSWC & swc:Answer_Tree.listNeuron){
        if(swc.pn==-1){
            st=swc;
        }
        if(In_Box(centerAPO,blocksize,swc)){
            for(const int & to:Answer_Graph[swc.timestamp]){
                if(to<=0) continue;
                const NeuronSWC & check = Answer_Map[to];
                if(!In_Box(centerAPO, blocksize, check)){
                    NeuronSWC Add_Point=swc;
                    Extend_Marker.push_back(Add_Point);
                    break;
                }
            }
        }
    }
    if(center_id!=0){
        int delete_id=-1;
        double mn=1e8;
        for(int i=0;i<Extend_Marker.size();++i){
            double dis=distance_XYZ(XYZ(Extend_Marker[i]),st);
            if(mn>dis){
                mn=dis;
                delete_id=i;
            }
        }
        if(!Extend_Marker.empty() && delete_id<Extend_Marker.size() && delete_id!=-1)
          Extend_Marker.erase(Extend_Marker.begin()+delete_id);
    }
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

void init(const QString & Answer_File,int & soma_id){

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
            soma_id=i.n;
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

void Fill_Vacant_Edge(V_NeuronSWC_list & Generate){
    NeuronTree Generate_Tree=V_NeuronSWC_list__2__NeuronTree(Generate);
    QMap<int,QMap<int,bool>> Edge;
    QMap<int, NeuronSWC> mp,ans_mp;
    for(const NeuronSWC & swc:Generate_Tree.listNeuron){
        if(swc.type==3){
            mp[swc.n]=swc;
            ans_mp[swc.timestamp]=swc;
        }
    }
    for(const NeuronSWC & swc:Generate_Tree.listNeuron){
        if(swc.type==3 && mp[swc.pn].type==3){
            Edge[swc.timestamp][mp[swc.pn].timestamp]=true;
        }
    }
    for(auto it=ans_mp.begin();it!=ans_mp.end();++it){
        int from = it.key();
        for(const int & to:Answer_Graph[from]){
            if(ans_mp.count(to) && (!Edge[from].count(to) && !Edge[to].count(from))){
                V_NeuronSWC Add;
                V_NeuronSWC_unit p1,p2;
                const NeuronSWC & swc1=ans_mp[from];
                const NeuronSWC & swc2=ans_mp[to];
                p1.x=swc1.x;
                p1.y=swc1.y;
                p1.z=swc1.z;
                p1.type=4;
                p1.n=1;
                p1.parent=2;
                p2.x=swc2.x;
                p2.y=swc2.y;
                p2.z=swc2.z;
                p2.type=4;
                p2.n=2;
                p2.parent=-1;
                Add.append(p1);
                Add.append(p2);
                Generate.append(Add);
            }
        }
    }
    return;
}

NeuronTree Complicated_App2(const int & center_id,const CellAPO & centerAPO,ImageMarker & startPoint,const int & blocksize,unsigned long long  & v3draw_size,bool & use_answer,const V_NeuronSWC_list & V_Answer_Tree){
    while(!Make_Dir(Work_Dir+QString("/APOFile")));
    while(!Make_Dir(Work_Dir+QString("/MarkerFile")));
    while(!Make_Dir(Work_Dir+QString("/SwcFile")));
    while(!Make_Dir(Work_Dir+QString("/testV3draw")));
    while(!Make_Dir(Work_Dir+QString("/EswcFile")));
    NeuronTree App2_Tree;
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

    //crop file size problem
    QFile crop_file(Work_Dir+QString("/testV3draw/")+rawFileName);
    if(v3draw_size == -1){
        v3draw_size=crop_file.size();
    }
    else{
        if(v3draw_size != crop_file.size()){
            use_answer=true;
            return NeuronTree();
        }
    }


    //ada(ok)
    qDebug()<<"ada:"
        <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"),QStringList()
                <<"/x"<<Vaa3d_App_Path+QString("/plugins/image_thresholding/Simple_Adaptive_Thresholding/ada_threshold.dll")
                <<"/f"<<"adath"<<"/i"<<Work_Dir+QString("/testV3draw/")+rawFileName<<"/o"<<QString(Work_Dir+QString("/testV3draw/thres_")+rawFileName));

    ImageMarker startPointBackup=startPoint;
    if(center_id!=0){
        QString v3draw=Work_Dir+QString("/testV3draw/thres_")+rawFileName;
        std::pair<QVector<ImageMarker>,QVector<ImageMarker> > possible=Get_Valid_Marker(v3draw,startPoint,blocksize);
        for(const ImageMarker & i:possible.first){
            startPoint=i;
            App2_Tree=Vanilla_App2(centerAPO,startPoint,blocksize,rawFileName);
            if(!App2_Tree.listNeuron.empty()){
                for(NeuronSWC & swc:App2_Tree.listNeuron){
                    swc.x+=centerAPO.x-blocksize/2;
                    swc.y+=centerAPO.y-blocksize/2;
                    swc.z+=centerAPO.z-blocksize/2;
                }
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
                    for(NeuronSWC & swc:App2_Tree.listNeuron){
                        swc.x+=centerAPO.x-blocksize/2;
                        swc.y+=centerAPO.y-blocksize/2;
                        swc.z+=centerAPO.z-blocksize/2;
                    }
                    V_NeuronSWC_list V_App2_Tree = NeuronTree__2__V_NeuronSWC_list(App2_Tree);
                    if(!Check_Tree_Identical(V_App2_Tree,V_Answer_Tree)){
                        App2_Tree=NeuronTree();
                        continue;
                    }
                    break;
                }
            }
        }

    }
    else {
        QString v3draw=Work_Dir+QString("/testV3draw/thres_")+rawFileName;
        std::pair<QVector<ImageMarker>,QVector<ImageMarker> > possible=Get_Valid_Marker(v3draw,startPoint,blocksize);
        for(const ImageMarker & i:possible.first){
            startPoint=i;
            App2_Tree=Vanilla_App2(centerAPO,startPoint,blocksize,rawFileName);
            if(!App2_Tree.listNeuron.empty()){
                for(NeuronSWC & swc:App2_Tree.listNeuron){
                    swc.x+=centerAPO.x-blocksize/2;
                    swc.y+=centerAPO.y-blocksize/2;
                    swc.z+=centerAPO.z-blocksize/2;
                }
                V_NeuronSWC_list V_App2_Tree = NeuronTree__2__V_NeuronSWC_list(App2_Tree);
                if(!Check_Tree_Identical(V_App2_Tree,V_Answer_Tree)){
                    App2_Tree=NeuronTree();
                    continue;
                }
                break;
            }
        }
    }

    if(App2_Tree.listNeuron.empty())    startPoint=startPointBackup,use_answer=true;

    QDir dir = QDir(Work_Dir+QString("/APOFile"));
    dir.removeRecursively();
    dir = QDir(Work_Dir+QString("/MarkerFile"));
    dir.removeRecursively();
    dir = QDir(Work_Dir+QString("/SwcFile"));
    dir.removeRecursively();
    dir = QDir(Work_Dir+QString("/testV3draw"));
    dir.removeRecursively();
    dir = QDir(Work_Dir+QString("/EswcFile"));
    dir.removeRecursively();
//    QFile::remove(Work_Dir+QString("/testV3draw/")+rawFileName);
//    QFile::remove(Work_Dir+QString("/testV3draw/thres_")+rawFileName);

    return App2_Tree;
}

V_NeuronSWC Link_Tree(const NeuronTree & t1,const NeuronTree & t2,const bool & pre_use_answer,const bool & use_answer){
    NeuronSWC pick1 ,pick2;
    double mn=1e8;
    for(const NeuronSWC & p1:t1.listNeuron){
        for(const NeuronSWC & p2:t2.listNeuron){
            double dis=distance_XYZ(p1,p2);
            if(mn>dis){
                pick1=p1;
                pick2=p2;
                mn=dis;
            }
        }
    }
    V_NeuronSWC ret;

    V_NeuronSWC_unit u1;
    u1.x=pick1.x;
    u1.y=pick1.y;
    u1.z=pick1.z;
    u1.n=1;
    u1.parent=2;
    V_NeuronSWC_unit u2;
    u2.x=pick2.x;
    u2.y=pick2.y;
    u2.z=pick2.z;
    u2.n=2;
    u2.parent=-1;
    if(pre_use_answer==true && use_answer == true){
        u1.type=2;
        u2.type=2;
    }
    else {
        u1.type=3;
        u1.type=3;
    }
    ret.row.push_back(u1);
    ret.row.push_back(u2);
    return ret;
}

double Max_Length(const NeuronTree & tree, const int & color){
    QVector<double> len;
    QVector<XYZ> head,tail;
    for(int i=0;i<tree.listNeuron.size();++i){
        if(tree.listNeuron[i].type!=color) continue;
        if(tree.listNeuron[i].parent==-1) continue;
        int pn = tree.listNeuron[i].pn;
        XYZ p1=tree.listNeuron[i], p2=tree.listNeuron[pn-1];
        bool insert=false;
        if(head.empty()){
            insert=true;
        }
        if(!insert){
           bool trigger=false;
           for(int j=0;j<head.length();++j){
               const XYZ & cmp1=head[j];
               const XYZ & cmp2=tail[j];
               if(XYZ_eq(cmp1,p1)){
                   head[j]=p2;
                   len[j]+=distance_XYZ(p1,p2);
                   trigger=true;
                   break;
               }
               if(XYZ_eq(cmp1,p2)){
                   head[j]=p1;
                   len[j]+=distance_XYZ(p1,p2);
                   trigger=true;
                   break;
               }
               if(XYZ_eq(cmp2,p1)){
                   tail[j]=p2;
                   len[j]+=distance_XYZ(p1,p2);
                   trigger=true;
                   break;
               }
               if(XYZ_eq(cmp2,p2)){
                   tail[j]=p1;
                   len[j]+=distance_XYZ(p1,p2);
                   trigger=true;
                   break;
               }
           }
           if(!trigger){
              insert=true;
           }
        }

        if(insert){
            head.push_back(p1);
            tail.push_back(p2);
            len.push_back(distance_XYZ(p1, p2));
        }
    }
    double mx = -1e8;
    for(int i=0;i<len.length();++i){
       mx=std::max(mx,len[i]);
    }
    return mx;
}

void App2_non_recursive_DFS(const int & blocksize,const QString & File_Name){

    int soma_id=-1;
    init(Answer_File,soma_id);

    CellAPO init_centerAPO;
    init_centerAPO.x=X;
    init_centerAPO.y=Y;
    init_centerAPO.z=Z;

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
    bbox_queue.push_back(bbox_extend(0,init_centerAPO,init_startPoint,soma_id,false,NeuronTree()));

    while(!bbox_queue.empty()||!if_finish()){
        if(not_change>=max_not_change){
            not_change=0;
            bbox_queue.clear();
            max_not_change=5;
        }
        if(bbox_queue.empty()){
            used_swc.clear();
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
            bbox_queue.push_back(bbox_extend(0,centerAPO,startPoint,swc.n,false,NeuronTree()));
            continue;
        }
        bbox_extend now=bbox_queue.back();
        bbox_queue.pop_back();
        CellAPO centerAPO=now.centerAPO;
        ImageMarker startPoint=now.startPoint;
        int center_id=now.center_id;
        int origin_id=now.origin_id;
        bool pre_use_answer=now.pre_use_answer;
        NeuronTree pretree=now.pretree;
        if(center_id==0){
            Center_Set.push_back(Answer_Map[origin_id]);
        }
        if(has_extend.count(center_id)) continue;
        if(if_finish()) {
            bbox_queue.clear();
            continue;
        }
        QVector<int> Answer_Tree_Idx;
        NeuronTree Answer_Tree=Get_Answer_Tree(origin_id,centerAPO,blocksize,Answer_Tree_Idx);
        if(Answer_Tree.listNeuron.empty()) continue;
        if(!if_need_extend(Answer_Tree)) continue;
        V_NeuronSWC_list V_Answer_Tree = NeuronTree__2__V_NeuronSWC_list(Answer_Tree);

        qDebug()<<"start DFS,center_id="<<center_id;

        QVector<NeuronSWC> Border_Points= Find_Extend_Marker(center_id,centerAPO,blocksize,Answer_Tree);

        NeuronTree App2_Tree=NeuronTree();
        bool use_answer=false;


        ImageMarker Absolute_Marker=startPoint;
        Absolute_Marker.x+=centerAPO.x-blocksize/2;
        Absolute_Marker.y+=centerAPO.y-blocksize/2;
        Absolute_Marker.z+=centerAPO.z-blocksize/2;


        if(Border_Points.size()>=2){
            use_answer=true;
        }
        else{
            App2_Tree=Complicated_App2(center_id,centerAPO,startPoint,blocksize,v3draw_size,use_answer,V_Answer_Tree);
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
               if((V_Answer_Tree.seg.size()==1?3:V_Answer_Tree.seg.size()*2)<V_App2_Tree.seg.size()||!Check_Tree_Identical(V_App2_Tree,V_Answer_Tree)||!Check_Tree_Identical(V_Answer_Tree,V_App2_Tree)||!Check_Tree_Length(V_Answer_Tree,V_App2_Tree)){
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

       if(center_id>0 && !pretree.listNeuron.empty() && !App2_Tree.listNeuron.empty()){
           Segments.seg.push_back(Link_Tree(pretree,App2_Tree,pre_use_answer,use_answer));
           App2_Tree=V_NeuronSWC_list__2__NeuronTree(Segments);
       }

       //save the answer(ok)
       bool if_change=false;
       QVector<V_NeuronSWC> Ready_To_Add;
       for(const V_NeuronSWC & Seg:Segments.seg){
           V_NeuronSWC Add_Seg=Seg;
           QMap<int,bool> redundant;
           for(int i=0; i<Add_Seg.row.size();++i){
               //(ok)
               V_NeuronSWC_unit & swc=Add_Seg.row[i];
//               swc.x+=centerAPO.x-blocksize/2;
//               swc.y+=centerAPO.y-blocksize/2;
//               swc.z+=centerAPO.z-blocksize/2;
               if(App2_Generate.seg.empty()) continue;
//               if(Distance_Unit_To_Tree(swc,App2_Generate)<close_distance){
//                   redundant[i]=true;
//               }
           }

           for(int i=0;i<Add_Seg.row.size();++i){
               if(!redundant.count(i)){
                   V_NeuronSWC Not_Redundant;
                   int now=i;
                   int amount=0;
                   while(now<Add_Seg.row.size() && !redundant.count(now)){
                       V_NeuronSWC_unit unit=Add_Seg.row[now];
                       ++amount;
                       if(use_answer) {
                           unit.type=3;
                           unit.timestamp=unit.n;
                       }
                       else unit.type=2;
                       unit.n=amount;
                       unit.parent=amount+1;
                       Not_Redundant.row.push_back(unit);
                       ++now;
                   }
                   Not_Redundant.row.back().parent=-1;
                   i=now;
//                   if(Not_Redundant.row.size()<=2) continue;
                   if_change=true;
                   Ready_To_Add.append(Not_Redundant);
               }
           }
       }
//       if(center_id!=0){
//           QVector<V_NeuronSWC_unit> head,tail;
//           for(const V_NeuronSWC & seg:Ready_To_Add){
//               head.push_back(seg.row[0]);
//               tail.push_back(seg.row[seg.row.size()-1]);
//           }
//           bool is_head=true;
//           int idx=-1;
//           double mn=1e8;
//           for(int i=0;i<head.size();++i){
//                XYZ hd;
//                hd.x=head[i].x;
//                hd.y=head[i].y;
//                hd.z=head[i].z;
//                double dis=distance_XYZ(origin,hd);
//                if(mn>dis){
//                    is_head=true;
//                    mn=dis;
//                    idx=i;
//                }
//                XYZ tl;
//                tl.x=tail[i].x;
//                tl.y=tail[i].y;
//                tl.z=tail[i].z;
//                dis=distance_XYZ(origin,tl);
//                if(mn>dis){
//                    is_head=false;
//                    mn=dis;
//                    idx=i;
//                }
//           }
//           if(idx!=-1){
//               if(is_head){
//                    V_NeuronSWC new_line;
//                    V_NeuronSWC_unit ori;
//                    ori.x=origin.x;
//                    ori.y=origin.y;
//                    ori.z=origin.z;
//                    ori.n=1;
//                    ori.parent=2;
//                    new_line.append(ori);
//                    int amount=1;
//                    for(const V_NeuronSWC_unit & unit:Ready_To_Add[idx].row){
//                        V_NeuronSWC_unit add=unit;
//                        ++amount;
//                        add.n=amount;
//                        add.parent=amount+1;
//                        new_line.append(add);
//                    }
//                    new_line.row[new_line.row.size()-1].parent=-1;
//                    Ready_To_Add[idx]=new_line;
//               }else{
//                   int last=Ready_To_Add[idx].row.size()-1;
//                   Ready_To_Add[idx].row[last].parent=Ready_To_Add[idx].row[last].n+1;
//                   V_NeuronSWC_unit ori;
//                   ori.x=origin.x;
//                   ori.y=origin.y;
//                   ori.z=origin.z;
//                   ori.n=Ready_To_Add[idx].row[last].parent;
//                   ori.parent=-1;
//                   Ready_To_Add[idx].append(ori);
//               }
//           }
//       }
       for(V_NeuronSWC & seg:Ready_To_Add){
           App2_Generate.append(seg);
       }

       update_Ans_used(Answer_Tree_Idx);

       if(!if_change)  ++not_change;
       else not_change=0;


       Delete_Redundant_Seg(App2_Generate);
       NeuronTree output=V_NeuronSWC_list__2__NeuronTree(App2_Generate);
       if(!use_answer) ++app2_success;
//       writeSWC_file(Work_Dir+"/EswcFile/"+QString(QString::fromStdString(std::to_string(++cnt)))+".eswc",output);


       //check if border is identical to last marker
       NeuronSWC Start_Marker_Location;
       //marker has been changed to absolute location
       Start_Marker_Location.x=Absolute_Marker.x;
       Start_Marker_Location.y=Absolute_Marker.y;
       Start_Marker_Location.z=Absolute_Marker.z;
       used_swc.push_back(Start_Marker_Location);

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
           int offset=blocksize/2-4;
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
           //move the center of next bounding_box forward in the accordingly direction, make the Border_Point locates in the center of area(??)

           ++amount;
           for(int i=0;i<1;++i){
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
               bbox_queue.push_back(bbox_extend(new_id,New_Point_Offset,New_Marker,Border_Point.timestamp,use_answer,App2_Tree));
           }

       }
    }

    Delete_Redundant_Seg(App2_Generate);

    NeuronTree output=V_NeuronSWC_list__2__NeuronTree(App2_Generate);
    writeSWC_file(Work_Dir+QString("/")+File_Name,output);

    QFile file1(Work_Dir+QString("/App2_Success"));
    qDebug()<<file1.open(QIODevice::WriteOnly|QIODevice::Text);

    QTextStream out1(&file1);
    out1<<app2_success;
    file1.close();

    QFile file2(Work_Dir+QString("/Max_App2_Length"));
    qDebug()<<file2.open(QIODevice::WriteOnly|QIODevice::Text);

    QTextStream out2(&file2);
    out2<<Max_Length(output,2);
    file2.close();
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
        if(i==5){
            seg_identical_threshold_mean=std::atof(argv[i]);
        }
        if(i==6){
            seg_identical_threshold_mx=std::atof(argv[i]);
        }
        if(i==7){
            length_indentical_rate=std::atof(argv[i]);
        }
    }

    App2_non_recursive_DFS(blocksize,Output_File_Name);

    qDebug()<<"finish";

}
