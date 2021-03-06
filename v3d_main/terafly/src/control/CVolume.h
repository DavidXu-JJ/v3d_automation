//------------------------------------------------------------------------------------------------
// Copyright (c) 2012  Alessandro Bria and Giulio Iannello (University Campus Bio-Medico of Rome).
// All rights reserved.
//------------------------------------------------------------------------------------------------

/*******************************************************************************************************************************************************************************************
*    LICENSE NOTICE
********************************************************************************************************************************************************************************************
*    By downloading/using/running/editing/changing any portion of codes in this package you agree to this license. If you do not agree to this license, do not download/use/run/edit/change
*    this code.
********************************************************************************************************************************************************************************************
*    1. This material is free for non-profit research, but needs a special license for any commercial purpose. Please contact Alessandro Bria at a.bria@unicas.it or Giulio Iannello at
*       g.iannello@unicampus.it for further details.
*    2. You agree to appropriately cite this work in your related studies and publications.
*
*       Bria, A., et al., (2012) "Stitching Terabyte-sized 3D Images Acquired in Confocal Ultramicroscopy", Proceedings of the 9th IEEE International Symposium on Biomedical Imaging.
*       Bria, A., Iannello, G., "A Tool for Fast 3D Automatic Stitching of Teravoxel-sized Datasets", submitted on July 2012 to IEEE Transactions on Information Technology in Biomedicine.
*
*    3. This material is provided by  the copyright holders (Alessandro Bria  and  Giulio Iannello),  University Campus Bio-Medico and contributors "as is" and any express or implied war-
*       ranties, including, but  not limited to,  any implied warranties  of merchantability,  non-infringement, or fitness for a particular purpose are  disclaimed. In no event shall the
*       copyright owners, University Campus Bio-Medico, or contributors be liable for any direct, indirect, incidental, special, exemplary, or  consequential  damages  (including, but not
*       limited to, procurement of substitute goods or services; loss of use, data, or profits;reasonable royalties; or business interruption) however caused  and on any theory of liabil-
*       ity, whether in contract, strict liability, or tort  (including negligence or otherwise) arising in any way out of the use of this software,  even if advised of the possibility of
*       such damage.
*    4. Neither the name of University  Campus Bio-Medico of Rome, nor Alessandro Bria and Giulio Iannello, may be used to endorse or  promote products  derived from this software without
*       specific prior written permission.
********************************************************************************************************************************************************************************************/

#ifndef CLOADSUBVOLUME_H
#define CLOADSUBVOLUME_H

#include "CViewer.h"

//#include <QThread>
#include <string>
#include "CPlugin.h"
#include "CImport.h"
//#include "CViewer.h"

#include <QThread>


class terafly::CVolume : public QThread
{
    Q_OBJECT

    private:

        /*********************************************************************************
        * Singleton design pattern: this class can have one instance only,  which must be
        * instantiated by calling static method "istance(...)"
        **********************************************************************************/
        static CVolume* uniqueInstance;
        CVolume() : QThread()
        {
         /**/tf::debug(tf::LEV1, 0, __itm__current__function__);

            reset();
        }

        //automatically called when current thread is started
        void run();

        //members
        //I suspect the "int" below might be problematic in the long run, but I leave them as is at this moment. I would use V3D_LONG instead.  by PHC 20131029
        int voiResIndex;                                            //volume of interest resolution index
        int voiV0,voiV1,voiH0,voiH1,voiD0,voiD1, voiT0, voiT1;      //volume of interest coordinates
        QWidget* source;                                            //the object that requested the VOI
        int streamingSteps;                                         //
        int cur_t;                                                  // current time frame selected (it is loaded and shown before the other frames)

    public:

        /*********************************************************************************
        * Singleton design pattern: this class can have one instance only,  which must be
        * instantiated by calling static method "istance(...)"
        **********************************************************************************/
        static CVolume* instance()
        {
            if (uniqueInstance == 0)
                uniqueInstance = new CVolume();
            return uniqueInstance;
        }
        static void uninstance();
        ~CVolume();

        //GET and SET methods
//        void initBuffer(tf::uint8* data, int xDim, int yDim, int zDim, int cDim=1, int tDim=1) throw (tf::RuntimeException)
//        {
//            /**/tf::debug(tf::LEV1, strprintf("xDim = %d, yDim=%d, zDim = %d, cDim=%d, tDim = %d",
//                                                xDim,       yDim,    zDim,      cDim,    tDim).c_str(), __itm__current__function__);

//            tf::uint64 size = xDim;
//            size *= yDim;
//            size *= zDim;
//            size *= cDim;
//            size *= tDim;
//            if(buffer)
//                delete[] buffer;

//            try{ buffer = new tf::uint8[size]; }
//            catch(...){ throw tf::RuntimeException("in CVolume::initBuffer(): cannot allocate memory");}

//            for(tf::uint8 *buf_p = buffer, *data_p = data; buf_p - buffer < size; buf_p++, data_p++)
//                *buf_p = *data_p;
//        }
//        tf::uint8* getBuffer(){return buffer;}

        void setStreamingSteps(int nsteps){streamingSteps = nsteps;}
        int getStreamingSteps(){return streamingSteps;}
        void reset()
        {
            /**/tf::debug(tf::LEV1, 0, __itm__current__function__);

            voiResIndex = -1;
            //buffer = 0;
            voiV0 = voiV1 = voiH0 = voiH1 = voiD0 = voiD1 = voiT0 = voiT1 = -1;
            source = 0;
            streamingSteps = 1;
            cur_t = -1;
        }
        int getVoiV0(){return voiV0;}
        int getVoiV1(){return voiV1;}
        int getVoiH0(){return voiH0;}
        int getVoiH1(){return voiH1;}
        int getVoiD0(){return voiD0;}
        int getVoiD1(){return voiD1;}
        int getVoiT0(){return voiT0;}
        int getVoiT1(){return voiT1;}
        int getVoiResIndex(){return voiResIndex;}
        void setVoi(QWidget* _sourceObject, int _voiResIndex, int _V0, int _V1, int _H0, int _H1, int _D0, int _D1, int _T0, int _T1) 
        {
            /**/tf::debug(tf::LEV1, strprintf("_voiResIndex = %d, _V0 = %d, _V1=%d, _H0 = %d, _H1=%d, _D0 = %d, _D1=%d, _T0 = %d, _T1=%d",
                                                _voiResIndex, _V0, _V1, _H0, _H1, _D0, _D1, _T0, _T1).c_str(), __itm__current__function__);

            source = _sourceObject;
            voiResIndex = _voiResIndex;
            iim::VirtualVolume* volume = CImport::instance()->getVolume(voiResIndex);

            //---- Alessandro 2013-08-06: reestabilished automatic VOI adjustement. This way, get methods return the actual VOI instead of the virtual one.
            //add () by PHC, 20131029
            voiV0 = (_V0 >=0)                   ? _V0 : 0;
            voiV1 = (_V1 <= volume->getDIM_V()) ? _V1 : volume->getDIM_V();
            voiH0 = (_H0 >=0)                   ? _H0 : 0;
            voiH1 = (_H1 <= volume->getDIM_H()) ? _H1 : volume->getDIM_H();
            voiD0 = (_D0 >=0)                   ? _D0 : 0;
            voiD1 = (_D1 <= volume->getDIM_D()) ? _D1 : volume->getDIM_D();
            voiT0 = (_T0 >=0)                   ? _T0 : 0;
            voiT1 = (_T1 <  volume->getDIM_T()) ? _T1 : volume->getDIM_T()-1;

            //qDebug()<<"testing bb ... "<<volume->getDIM_V()<<volume->getDIM_H()<<volume->getDIM_D();
            //qDebug()<<"testing voi ... "<<voiV1 - voiV0 << voiH1 - voiH0 << voiD1 - voiD0 << voiT1 - voiT0;

            //---- Alessandro 2013-09-03: added check to detect invalid VOI
            if(voiV1 - voiV0 <= 0 || voiH1 - voiH0 <= 0 || voiD1 - voiD0 <= 0 || voiT1 - voiT0 < 0)
                throw tf::RuntimeException("Invalid VOI selected");
        }
        void setSource(QWidget* _sourceObject){source =_sourceObject;}

        void setVoiT(int _T0, int _T1, int _cur_t = -1) 
        {
            /**/tf::debug(tf::LEV1, strprintf("asked to set [%d, %d] and _cur_t = %d",_T0, _T1, _cur_t).c_str(), __itm__current__function__);
            iim::VirtualVolume* volume = CImport::instance()->getVolume(voiResIndex);

            // correct [voiT0, voiT1]
            voiT0 = (_T0 >=0)                   ? _T0 : 0;
            voiT1 = (_T1 <  volume->getDIM_T()) ? _T1 : volume->getDIM_T()-1;

            // invalidate cur_t if out of range
            if(_cur_t < voiT0 || _cur_t > voiT1)
                cur_t = -1;
            else
                cur_t = _cur_t;

            /**/tf::debug(tf::LEV1, strprintf("but set [%d, %d] and cur_t = %d", voiT0, voiT1, cur_t).c_str(), __itm__current__function__);

            if(voiT1 - voiT0 < 0)
                throw tf::RuntimeException(tf::strprintf("Invalid VOI selected along T: [%d,%d]", voiT0, voiT1).c_str());
        }

        template<typename T>
        static inline T scaleCoord(T coord, int srcRes, int dstRes, iim::axis dir, bool round) 
        {
            #ifdef terafly_enable_debug_max_level
            /**/tf::debug(tf::LEV3, tf::strprintf("coord = %s, srcRes = %d, dstRes = %d, dir = %d, round = %s", tf::num2str(coord).c_str(), srcRes, dstRes, dir, round ? "true": "false").c_str(), __itm__current__function__);
            #endif

            // checks
            if(srcRes < 0 || srcRes >= CImport::instance()->getResolutions())
                throw tf::RuntimeException("scaleCoord(): invalid source resolution");
            if(dstRes < 0 || dstRes >= CImport::instance()->getResolutions())
                throw tf::RuntimeException("scaleCoord(): invalid destination resolution");

            // special case: image spaces match
            if(srcRes == dstRes)
            {
                #ifdef terafly_enable_debug_max_level
                /**/tf::debug(tf::LEV3, tf::strprintf("image spaces match, return %s", tf::num2str(coord).c_str()).c_str(), __itm__current__function__);
                #endif
                return coord;
            }
            // special case: boundary coordinate
            else if(coord >= CImport::instance()->getVolume(srcRes)->getDIM(dir))
            {
                #ifdef terafly_enable_debug_max_level
                /**/tf::debug(tf::LEV3, tf::strprintf("boundary coordinate, return %d", CImport::instance()->getVolume(dstRes)->getDIM(dir)).c_str(), __itm__current__function__);
                #endif
                return CImport::instance()->getVolume(dstRes)->getDIM(dir);
            }
            // special case: 2D image
            else if(CImport::instance()->getVolume(srcRes)->getDIM(dir) == 1)
            {
                #ifdef terafly_enable_debug_max_level
                /**/tf::debug(tf::LEV3, tf::strprintf("srcRes is 2D, return %d", coord ? CImport::instance()->getVolume(dstRes)->getDIM(dir) : 0).c_str(), __itm__current__function__);
                #endif
                return coord ? CImport::instance()->getVolume(dstRes)->getDIM(dir) : 0;
            }
            // normal case: scale mapping
            else
            {
                float rescale = CImport::instance()->getRescaleFactor(dstRes, srcRes, dir);

                #ifdef terafly_enable_debug_max_level
                /**/tf::debug(tf::LEV3, tf::strprintf("normal case, rescale = %f, return %s", rescale, tf::num2str(round? tf::round(coord*rescale) : (coord*rescale)).c_str()).c_str(), __itm__current__function__);
                #endif

                return round? tf::round(coord*rescale) : coord*rescale;
            }
        }

//        static inline int scaleVCoord(int coord, int srcRes, int dstRes) throw (RuntimeException)
//        {
//            #ifdef terafly_enable_debug_max_level
//            /**/tf::debug(tf::LEV3, strprintf("coord = %d, srcRes = %d, dstRes = %d", coord, srcRes, dstRes).c_str(), __itm__current__function__);
//            #endif

//            //checks
//            if(srcRes < 0 || srcRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid source resolution of coordinate mapping operation");
//            if(dstRes < 0 || dstRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid destination resolution of coordinate mapping operation");

//            //computation
//            if(srcRes == dstRes)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes = dstRes, return %d", coord).c_str(), __itm__current__function__);
//                #endif
//                return coord;
//            }
//            // --- Alessandro added on August 8th, 2013: useful for exact boundary conversion ---
//            else if(coord == CImport::instance()->getVolume(srcRes)->getDIM_V())
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("coord = srcRes->DIM_V(), return %d", CImport::instance()->getVolume(dstRes)->getDIM_V()).c_str(), __itm__current__function__);
//                #endif
//                return CImport::instance()->getVolume(dstRes)->getDIM_V();
//            }
//            else if(CImport::instance()->getVolume(srcRes)->getDIM_V() == 1)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes is 2D, return %d", coord ? CImport::instance()->getVolume(dstRes)->getDIM_V() : 0).c_str(), __itm__current__function__);
//                #endif
//                return coord ? CImport::instance()->getVolume(dstRes)->getDIM_V() : 0;
//            }
//            else
//            {
//                float ratio = (CImport::instance()->getVolume(dstRes)->getDIM_V()-1.0f)/(CImport::instance()->getVolume(srcRes)->getDIM_V()-1.0f);
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("final else, return %d", static_cast<int>(coord*ratio + 0.5f)).c_str(), __itm__current__function__);
//                #endif
//                return static_cast<int>(coord*ratio + 0.5f);
//            }
//        }
//        static inline int scaleHCoord(int coord, int srcRes, int dstRes) throw (RuntimeException)
//        {
//            #ifdef terafly_enable_debug_max_level
//            /**/tf::debug(tf::LEV3, strprintf("coord = %d, srcRes = %d, dstRes = %d", coord, srcRes, dstRes).c_str(), __itm__current__function__);
//            #endif

//            //checks
//            if(srcRes < 0 || srcRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid source resolution of coordinate mapping operation");
//            if(dstRes < 0 || dstRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid destination resolution of coordinate mapping operation");

//            //computation
//            if(srcRes == dstRes)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes = dstRes, return %d", coord).c_str(), __itm__current__function__);
//                #endif
//                return coord;
//            }
//            else if(coord == CImport::instance()->getVolume(srcRes)->getDIM_H())
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("coord = srcRes->DIM_H(), return %d", CImport::instance()->getVolume(dstRes)->getDIM_H()).c_str(), __itm__current__function__);
//                #endif
//                return CImport::instance()->getVolume(dstRes)->getDIM_H();
//            }
//            else if(CImport::instance()->getVolume(srcRes)->getDIM_H() == 1)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes is 2D, return %d", coord ? CImport::instance()->getVolume(dstRes)->getDIM_H() : 0).c_str(), __itm__current__function__);
//                #endif
//                return coord ? CImport::instance()->getVolume(dstRes)->getDIM_H() : 0;
//            }
//            else
//            {
//                float ratio = (CImport::instance()->getVolume(dstRes)->getDIM_H()-1.0f)/(CImport::instance()->getVolume(srcRes)->getDIM_H()-1.0f);
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("final else, return %d", static_cast<int>(coord*ratio + 0.5f)).c_str(), __itm__current__function__);
//                #endif
//                return static_cast<int>(coord*ratio + 0.5f);
//            }
//        }
//        static inline int scaleDCoord(int coord, int srcRes, int dstRes) throw (RuntimeException)
//        {
//            #ifdef terafly_enable_debug_max_level
//            /**/tf::debug(tf::LEV3, strprintf("coord = %d, srcRes = %d, dstRes = %d, srcRes->getDIM_D() = %d, dstRes->getDIM_D() = %d", coord, srcRes, dstRes, CImport::instance()->getVolume(srcRes)->getDIM_D(), CImport::instance()->getVolume(dstRes)->getDIM_D()).c_str(), __itm__current__function__);
//            #endif

//            //checks
//            if(srcRes < 0 || srcRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid source resolution of coordinate mapping operation");
//            if(dstRes < 0 || dstRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid destination resolution of coordinate mapping operation");

//            //computation
//            if(srcRes == dstRes)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes = dstRes, return %d", coord).c_str(), __itm__current__function__);
//                #endif
//                return coord;
//            }
//            else if(coord == CImport::instance()->getVolume(srcRes)->getDIM_D())
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("coord = srcRes->DIM_D(), return %d", CImport::instance()->getVolume(dstRes)->getDIM_D()).c_str(), __itm__current__function__);
//                #endif
//                return CImport::instance()->getVolume(dstRes)->getDIM_D();
//            }
//            else if(CImport::instance()->getVolume(srcRes)->getDIM_D() == 1)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes is 2D, return %d", coord ? CImport::instance()->getVolume(dstRes)->getDIM_D() : 0).c_str(), __itm__current__function__);
//                #endif
//                return coord ? CImport::instance()->getVolume(dstRes)->getDIM_D() : 0;
//            }
//            else
//            {
//                float ratio = (CImport::instance()->getVolume(dstRes)->getDIM_D()-1.0f)/(CImport::instance()->getVolume(srcRes)->getDIM_D()-1.0f);
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("final else, return %d", static_cast<int>(coord*ratio + 0.5f)).c_str(), __itm__current__function__);
//                #endif
//                return static_cast<int>(coord*ratio + 0.5f);
//            }
//        }
//        static inline float scaleVCoord(float coord, int srcRes, int dstRes) throw (RuntimeException)
//        {
//            #ifdef terafly_enable_debug_max_level
//            /**/tf::debug(tf::LEV3, strprintf("coord = %.3f, srcRes = %d, dstRes = %d", coord, srcRes, dstRes).c_str(), __itm__current__function__);
//            #endif

//            //checks
//            if(srcRes < 0 || srcRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid source resolution of coordinate mapping operation");
//            if(dstRes < 0 || dstRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid destination resolution of coordinate mapping operation");

//            //computation
//            if(srcRes == dstRes)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes = dstRes, return %f", coord).c_str(), __itm__current__function__);
//                #endif
//                return coord;
//            }
//            else if(coord == static_cast<float>(CImport::instance()->getVolume(srcRes)->getDIM_V()))
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("coord = srcRes->DIM_V(), return %d", CImport::instance()->getVolume(dstRes)->getDIM_V()).c_str(), __itm__current__function__);
//                #endif
//                return static_cast<float>(CImport::instance()->getVolume(dstRes)->getDIM_V());
//            }
//            else if(CImport::instance()->getVolume(srcRes)->getDIM_V() == 1)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes is 2D, return %d", coord ? CImport::instance()->getVolume(dstRes)->getDIM_V() : 0).c_str(), __itm__current__function__);
//                #endif
//                return coord ? static_cast<float>(CImport::instance()->getVolume(dstRes)->getDIM_V()) : 0;
//            }
//            else
//            {
//                float ratio = (CImport::instance()->getVolume(dstRes)->getDIM_V()-1.0f)/(CImport::instance()->getVolume(srcRes)->getDIM_V()-1.0f);
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("final else, return %f", coord*ratio).c_str(), __itm__current__function__);
//                #endif
//                return coord*ratio;
//            }
//        }
//        static inline float scaleHCoord(float coord, int srcRes, int dstRes) throw (RuntimeException)
//        {
//            #ifdef terafly_enable_debug_max_level
//            /**/tf::debug(tf::LEV3, strprintf("coord = %.3f, srcRes = %d, dstRes = %d", coord, srcRes, dstRes).c_str(), __itm__current__function__);
//            #endif

//            //checks
//            if(srcRes < 0 || srcRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid source resolution of coordinate mapping operation");
//            if(dstRes < 0 || dstRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid destination resolution of coordinate mapping operation");

//            //computation
//            if(srcRes == dstRes)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes = dstRes, return %f", coord).c_str(), __itm__current__function__);
//                #endif
//                return coord;
//            }
//            else if(coord == static_cast<float>(CImport::instance()->getVolume(srcRes)->getDIM_H()))
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("coord = srcRes->DIM_H(), return %d", CImport::instance()->getVolume(dstRes)->getDIM_H()).c_str(), __itm__current__function__);
//                #endif
//                return static_cast<float>(CImport::instance()->getVolume(dstRes)->getDIM_H());
//            }
//            else if(CImport::instance()->getVolume(srcRes)->getDIM_H() == 1)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes is 2D, return %d", coord ? CImport::instance()->getVolume(dstRes)->getDIM_H() : 0).c_str(), __itm__current__function__);
//                #endif
//                return coord ? static_cast<float>(CImport::instance()->getVolume(dstRes)->getDIM_H()) : 0;
//            }
//            else
//            {
//                float ratio = (CImport::instance()->getVolume(dstRes)->getDIM_H()-1.0f)/(CImport::instance()->getVolume(srcRes)->getDIM_H()-1.0f);
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("final else, return %f", coord*ratio).c_str(), __itm__current__function__);
//                #endif
//                return coord*ratio;
//            }
//        }
//        static inline float scaleDCoord(float coord, int srcRes, int dstRes) throw (RuntimeException)
//        {
//            #ifdef terafly_enable_debug_max_level
//            /**/tf::debug(tf::LEV3, strprintf("coord = %.3f, srcRes = %d, dstRes = %d", coord, srcRes, dstRes).c_str(), __itm__current__function__);
//            #endif

//            //checks
//            if(srcRes < 0 || srcRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid source resolution of coordinate mapping operation");
//            if(dstRes < 0 || dstRes >= CImport::instance()->getResolutions())
//                throw RuntimeException("Invalid destination resolution of coordinate mapping operation");

//            //computation
//            if(srcRes == dstRes)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes = dstRes, return %f", coord).c_str(), __itm__current__function__);
//                #endif
//                return coord;
//            }
//            else if(coord == static_cast<float>(CImport::instance()->getVolume(srcRes)->getDIM_D()))
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("coord = srcRes->DIM_D(), return %d", CImport::instance()->getVolume(dstRes)->getDIM_D()).c_str(), __itm__current__function__);
//                #endif
//                return static_cast<float>(CImport::instance()->getVolume(dstRes)->getDIM_D());
//            }
//            else if(CImport::instance()->getVolume(srcRes)->getDIM_D() == 1)
//            {
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("srcRes is 2D, return %d", coord ? CImport::instance()->getVolume(dstRes)->getDIM_D() : 0).c_str(), __itm__current__function__);
//                #endif
//                return coord ? static_cast<float>(CImport::instance()->getVolume(dstRes)->getDIM_D()) : 0;
//            }
//            else
//            {
//                float ratio = (CImport::instance()->getVolume(dstRes)->getDIM_D()-1.0f)/(CImport::instance()->getVolume(srcRes)->getDIM_D()-1.0f);
//                #ifdef terafly_enable_debug_max_level
//                /**/tf::debug(tf::LEV3, strprintf("final else, return %f", coord*ratio).c_str(), __itm__current__function__);
//                #endif
//                return coord*ratio;
//            }
//        }

        // load data using the currently set VOI
        tf::uint8* loadData() ;

        friend class CViewer;

    signals:

        /*********************************************************************************
        * Send data (and metadata) to the listener throughout the loading process
        **********************************************************************************/
        void sendData(
                tf::uint8* data,                   // data (any dimension)
                tf::integer_array data_s,          // data start coordinates along X, Y, Z, C, t
                tf::integer_array data_c,          // data count along X, Y, Z, C, t
                QWidget* dest,                      // listener handle
                bool finished,                      // whether the loading operation is terminated
                tf::RuntimeException* ex = 0,      // exception (optional)
                qint64 elapsed_time = 0,            // elapsed time (optional)
                QString op_dsc="",                  // operation descriptor (optional)
                int step=0);                        // step number (optional)
};

#endif // CLOADSUBVOLUME_H
