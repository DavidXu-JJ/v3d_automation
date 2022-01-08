

## Geometry

```c++
double distance_square(const NeuronSWC & point_a,const NeuronSWC & point_b){
    return (point_a.x-point_b.x)*(point_a.x-point_b.x)+(point_a.y-point_b.y)*(point_a.y-point_b.y)+(point_a.z-point_b.z)*(point_a.z-point_b.z);
}

bool between(const int & mid,const int & left,const int & right){   //double?
    return mid>=left&&mid<=right;
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
    const double threshold=10;
    double mx=-1e8;
    for(const V_NeuronSWC_unit & Check_Point:Check_Seg.row){
         for(const V_NeuronSWC & Answer_Seg:Answer_Tree.seg){
            mx=std::max(mx,Distance_Unit_To_Seg(Check_Point,Answer_Seg));
        }
    }
    return mx<threshold;
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

```

## File relative

```c++
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


NeuronTree Get_Answer_Tree(const CellAPO & centerAPO,const ImageMarker & startPoint,const int & blocksize,const QString & Answer_File){
    NeuronTree Ans_Tree=readSWC_file(Answer_File);

    QMap<int,NeuronSWC> mp;
    QMap<int,QVector<int> > G;
    for(const NeuronSWC & i:Ans_Tree.listNeuron){
        mp[i.n]=i;
        G[i.n].push_back(i.pn);
        G[i.pn].push_back(i.n);
    }
    const int & X=centerAPO.x;
    const int & Y=centerAPO.y;
    const int & Z=centerAPO.z;

    QList<NeuronSWC> In_Block_Points;
    QMap<int,int> weight;
    for(NeuronSWC & swc:Ans_Tree.listNeuron){
        if(between(swc.x,X-blocksize/2,X+blocksize/2)&&between(swc.y,Y-blocksize/2,Y+blocksize/2)&&between(swc.z,Z-blocksize/2,Z+blocksize/2)){
            swc.x=swc.x-X+blocksize/2;
            swc.y=swc.y-Y+blocksize/2;
            swc.z=swc.z-Z+blocksize/2;
            In_Block_Points.push_back(swc);
            ++weight[swc.n];
            ++weight[swc.pn];
        }
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

    const int threshold=20;
    if(mn>threshold) return NeuronTree();

    NeuronTree Return_Tree;
    NeuronSWC StartPoint_In_Answer=mp[pick_id];
    StartPoint_In_Answer.n=1;
    StartPoint_In_Answer.pn=-1;
    Return_Tree.listNeuron.push_back(StartPoint_In_Answer);

    int amount=1;
    QQueue<std::pair<int,int> > q;
    q.push_back(std::make_pair(-1,pick_id));    //-1 is new, pick_id is old
    while(!q.empty()){
        int prev=q.front().first;   //new
        int now=q.front().second;   //old
        q.pop_front();

        NeuronSWC Add_Point=mp[now];
        ++amount;
        Add_Point.n = amount;
        Add_Point.pn= prev;
        Return_Tree.listNeuron.push_back(Add_Point);

        for(const int & i:G[now]){
            q.push_back(std::make_pair(amount,i));
        }
    }

    return Return_Tree;
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

```

## Border and Direction

```c++

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

    return direction;

    //decide to expand on positive axis or negative axis
    //pick=0 positive axis x
    //pick=1 negative axis x
    //pick=2 positive axis y
    //pick=3 negative axis y
    //pick=4 positive axis z
    //pick=5 negative axis z
}

```