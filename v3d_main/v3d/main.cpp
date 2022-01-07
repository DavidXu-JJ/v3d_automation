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

double distance_square(const NeuronSWC & point_a,const NeuronSWC & point_b){
    return (point_a.x-point_b.x)*(point_a.x-point_b.x)+(point_a.y-point_b.y)*(point_a.y-point_b.y)+(point_a.z-point_b.z)*(point_a.z-point_b.z);
}

bool between(const int & mid,const int & left,const int & right){   //double?
    return mid>=left&&mid<=right;
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

    const double Border_Threshold=6;
    QVector<std::pair<NeuronSWC,QQueue<XYZ> > > ret;
    for(const std::pair<int,QQueue<XYZ> >  & i:border){
        const int & id=i.first;
        const XYZ & point=mp[id];

        double mn=std::min({point.x,abs(blocksize-point.x),point.y,abs(blocksize-point.y),point.z,abs(blocksize-point.z)});

        if(mn>Border_Threshold) continue;

        //save Border_Point and its previous vector
        ret.push_back(std::make_pair(mp[id],i.second));

    }
    qDebug()<<"Border_Finish";
    return ret;
}

int Judge_Direction(const XYZ & vec){
    int weight[3][3]={{0,1,1},{1,0,1},{1,1,0}};
    double mn=1e10;
    int pick= -1;
    //decide to expand in which axis
    //pick=0 equals weight[i][0]=0, vec is more close to axis x
    //pick=2 equals weight[i][1]=0, vec is more close to axis y
    //pick=4 equals weight[i][2]=0, vec is more close to axis z
    for(int i=0;i<3;++i){
        double value=weight[i][0]*vec.x*vec.x+weight[i][1]*vec.y*vec.y+weight[i][2]*vec.z*vec.z;
        if(value<mn){
            mn=value;
            pick=2*i;
        }
    }

    //decide to expand on positive axis or negative axis
    if(pick==0) return pick+(vec.x<0);
    if(pick==2) return pick+(vec.y<0);
    return pick+(vec.z<0);
    //pick=0 positive axis x
    //pick=1 negative axis x
    //pick=2 positive axis y
    //pick=3 negative axis y
    //pick=4 positive axis z
    //pick=5 negative axis z
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


QString Res_Path="G:/18454/RES(26298x35000x11041)";
QString Vaa3d_App_Path="C:/3.603c";
QString Work_Dir="G:/work_dir_debug";

QProcess p;

int dx[6]={1,-1,0,0,0,0};
int dy[6]={0,0,1,-1,0,0};
int dz[6]={0,0,0,0,1,-1};
QMap<node,bool> vis;
QVector<NeuronSWC> used_swc;

const int X=14530,Y=10693,Z=3124;

void DFS(const CellAPO & centerAPO,const ImageMarker & startPoint,const int &blocksize,V_NeuronSWC_list & App2_Generate){
//    if(dep==3) return;

    node now=node(centerAPO.x,centerAPO.y,centerAPO.z,startPoint.x,startPoint.y,startPoint.z);
    if(vis.count(now)) return;

    //drop apo file(ok)
    QList<CellAPO> List_APO_Write;
    List_APO_Write.push_back(centerAPO);
    QString APO_File_Name=generate_apo_name(Work_Dir.toStdString()+"/APOFile",centerAPO);	//make file in ./APOFile/xxx.000_xxx.000_xxx.000.apo
    writeAPO_file(APO_File_Name,List_APO_Write);

    //drop marker file(ok)
    QList<ImageMarker> List_Marker_Write;
    List_Marker_Write.push_back(startPoint);
    QString Marker_File_Name=generate_marker_name(Work_Dir.toStdString()+"/MarkerFile",startPoint);	//make file in ./MarkerFile/xxx.000_xxx.000_xxx.000.marker
    writeMarker_file(Marker_File_Name,List_Marker_Write);

    //crop3D(ok)
    qDebug()<<"crop3D:"
              <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"),QStringList()
              <<"/x"<<Vaa3d_App_Path+QString("/plugins/image_geometry/crop3d_image_series/cropped3DImageSeries.dll")
              <<"/f"<<"cropTerafly"<<"/i"<<Res_Path<<APO_File_Name<<Work_Dir+QString("/testV3draw/")
              <<"/p"<<QString::number(blocksize)<<QString::number(blocksize)<<QString::number(blocksize));

    QString rawFileName=QString("%1.000_%2.000_%3.000.v3draw").arg(centerAPO.x).arg(centerAPO.y).arg(centerAPO.z);

    //ada(ok)
    qDebug()<<"ada:"
            <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"),QStringList()
            <<"/x"<<Vaa3d_App_Path+QString("/plugins/image_thresholding/Simple_Adaptive_Thresholding/ada_threshold.dll")
            <<"/f"<<"adath"<<"/i"<<Work_Dir+QString("/testV3draw/")+rawFileName<<"/o"<<QString(Work_Dir+QString("/testV3draw/thres_")+rawFileName));

    //app2(ok)
    QString App2_Eswc_File_Name=generate_eswc_name(Work_Dir.toStdString()+"/SwcFile",centerAPO);
    qDebug()<<App2_Eswc_File_Name;
    qDebug()<<"app2:"
          <<p.execute(Vaa3d_App_Path+QString("/vaa3d_msvc.exe"), QStringList()
          <<"/x"<<Vaa3d_App_Path+QString("/plugins/neuron_tracing/Vaa3D_Neuron2/vn2.dll")
          <<"/f"<<"app2"<<"/p"<<Marker_File_Name<<QString::number(0)<<QString::number(-1)
          <<"/i"<< QString(Work_Dir+QString("/testV3draw/thres_")+rawFileName)<<"/o"<<App2_Eswc_File_Name);

    vis[now]=true;

    //readfile(ok)
    NeuronTree App2_Tree=readSWC_file(App2_Eswc_File_Name);
    if(App2_Tree.listNeuron.empty()) return;



    //find_border_point
    QMap<int,QVector<int> > son;    //record son
    QMap<int,NeuronSWC> mp;         //map id to single_swc
    for(const NeuronSWC & swc:App2_Tree.listNeuron){
        son[swc.pn].push_back(swc.n);
        mp[swc.n]=swc;
    }
    QVector<std::pair<NeuronSWC,QQueue<XYZ> > > Border_Point_Vector=Find_Border(App2_Tree,blocksize,son,mp);

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

    //save the answer(ok)
    V_NeuronSWC_list Segments=NeuronTree__2__V_NeuronSWC_list(App2_Tree);
    for(const V_NeuronSWC & Seg:Segments.seg){
        V_NeuronSWC Add_Seg=Seg;
        for(V_NeuronSWC_unit & swc:Add_Seg.row){
            //(ok)
            swc.x+=centerAPO.x-blocksize/2;
            swc.y+=centerAPO.y-blocksize/2;
            swc.z+=centerAPO.z-blocksize/2;
        }
        App2_Generate.append(Add_Seg);
    }

    //check if border is identical to last marker
    NeuronSWC Start_Marker_Location;
    //location transform has been checked(ok)
    Start_Marker_Location.x=centerAPO.x-blocksize/2+startPoint.x;
    Start_Marker_Location.y=centerAPO.y-blocksize/2+startPoint.y;
    Start_Marker_Location.z=centerAPO.z-blocksize/2+startPoint.z;
    used_swc.push_back(Start_Marker_Location);

    //if the square of distance between used_swc and new_border_point is lower than identical_threshold, don't search it
    const double identical_threshold=600;   //(undetermined)

    //treat with the border point
    for(std::pair<NeuronSWC,QQueue<XYZ> > & unsecure:Border_Point_Vector){

        //we think the father of the furthest inaccurate is real Border_Point
        const NeuronSWC & Border_Point=mp[unsecure.first.pn];

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
        Vectors.pop_back(); //don't take the last inaccuracy line into consideration
        //take the mean(ok)
        XYZ mean=XYZ(0.0,0.0,0.0);
        for(const XYZ & i:Vectors){
            mean=mean+i;
        }
        int sz=Vectors.size();
        mean=mean/XYZ(sz);
        //judge direction by mean(recommend to read)
        int direction=Judge_Direction(mean);

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
        New_Point.x+=dx[direction]*offset;
        New_Point.y+=dy[direction]*offset;
        New_Point.z+=dz[direction]*offset;
        //calculate the marker in the next bounding_box
        ImageMarker New_Marker;
        //Border_Point's relative position to the center of next bounding_box decrease or increase in the opposite direction
        New_Marker.x=blocksize/2-dx[direction]*offset;
        New_Marker.y=blocksize/2-dy[direction]*offset;
        New_Marker.z=blocksize/2-dz[direction]*offset;

        DFS(New_Point,New_Marker,blocksize,App2_Generate);
    }

}

void App2_DFS(const int & Start_x,const int & Start_y,const int & Start_z,const int & blocksize,const QString & File_Name){
    CellAPO centerAPO;
    centerAPO.x=Start_x;
    centerAPO.y=Start_y;
    centerAPO.z=Start_z;

    ImageMarker startPoint;
    startPoint.x=blocksize/2;
    startPoint.y=blocksize/2;
    startPoint.z=blocksize/2;

    V_NeuronSWC_list App2_Generate;
    DFS(centerAPO,startPoint,blocksize,App2_Generate);

    NeuronTree output=V_NeuronSWC_list__2__NeuronTree(App2_Generate);

    writeSWC_file(Work_Dir+QString("/whole_image.eswc"),output);
}


int main(int argc, char **argv)
{


    const int blocksize=128;
    App2_DFS(X,Y,Z,blocksize,"./");


    qDebug()<<"finish";


}
