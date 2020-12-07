#include <QFile>
#include <QProgressDialog>
#include <QApplication>
#include <QDir>
#include <QProcess>

#include "libIOctoMap.h"
#include <octomap/octomap.h>
#include <octomap/ColorOcTree.h>

libIOctoMap * libIOctoMap::mInstance = 0;

using namespace octomap;

bool libIOctoMap::getNewOriginForRayInBoundingBox(double xs, double ys, double zs, // start
                                                  double xe, double ye, double ze, // end
                                                  double minX, double minY,
                                                  double maxX, double maxY,
                                                  double& xo,
                                                  double& yo,
                                                  double& zo)
{
    xo=xs;yo=ys;zo=zs;
//    if(fabs(xs-minX)<LIBIOCTOMAP_LINEAR_TOLERANCE)
//    {
//        if(ys<minY)
//    }
    if((xs<minX&&xe<minX)
            ||(xs>maxX&&xe>maxX)
            ||(ys<minY&&ye<minY)
            ||(ys>maxY&&ye>maxY))
    {
        return(false);
    }
    double distanceStartToEnd=sqrt(pow((xe-xs),2.0)+pow((ye-ys),2.0));
    double slopeStartToEnd=(ze-zs)/distanceStartToEnd;
    for(int i=0;i<4;i++)
    {
        double x1,y1,x2,y2;
        if(i==0)
        {
            x1=minX;y1=minY;x2=minX;y2=maxY;
        }
        if(i==1)
        {
            x1=minX;y1=maxY;x2=maxX;y2=maxY;
        }
        if(i==2)
        {
            x1=maxX;y1=maxY;x2=maxX;y2=minY;
        }
        if(i==3)
        {
            x1=maxX;y1=minY;x2=minX;y2=minY;
        }
        float d = (x1 - x2) * (ys - ye) - (y1 - y2) * (xs - xe);
        // If d is zero, there is no intersection
        if (d == 0)
        {
            continue;
        }

        // Get the x and y
        float pre = (x1*y2 - y1*x2), post = (xs*ye - ys*xe);
        float x = ( pre * (xs - xe) - (x1 - x2) * post ) / d;
        float y = ( pre * (ys - ye) - (y1 - y2) * post ) / d;

        // Check if the x and y coordinates are within ray
        if(x<min(xs,xe)&&fabs(x-min(xs,xe))>LIBIOCTOMAP_LINEAR_TOLERANCE)
        {
            continue;
        }
        else if(x>max(xs,xe)&&fabs(max(xs,xe)-x)>LIBIOCTOMAP_LINEAR_TOLERANCE)
        {
            continue;
        }
        else if(y<min(ys,ye)&&fabs(y-min(ys,ye))>LIBIOCTOMAP_LINEAR_TOLERANCE)
        {
            continue;
        }
        else if(y>max(ys,ye)&&fabs(max(ys,ye)-y)>LIBIOCTOMAP_LINEAR_TOLERANCE)
        {
            continue;
        }
        // Check if the x and y coordinates are within bounding box
        if(i==0)
        {
            if(fabs(x-minX)>LIBIOCTOMAP_LINEAR_TOLERANCE
                    ||(y<minY&&fabs(y-minY)>LIBIOCTOMAP_LINEAR_TOLERANCE)
                    ||(y>maxY&&fabs(y-maxY)>LIBIOCTOMAP_LINEAR_TOLERANCE))
            {
                continue;
            }
        }
        else if(i==1)
        {
            if(fabs(y-maxY)>LIBIOCTOMAP_LINEAR_TOLERANCE
                    ||(x<minX&&fabs(x-minX)>LIBIOCTOMAP_LINEAR_TOLERANCE)
                    ||(x>maxX&&fabs(x-maxX)>LIBIOCTOMAP_LINEAR_TOLERANCE))
            {
                continue;
            }

        }
        else if(i==2)
        {
            if(fabs(x-maxX)>LIBIOCTOMAP_LINEAR_TOLERANCE
                    ||(y<minY&&fabs(y-minY)>LIBIOCTOMAP_LINEAR_TOLERANCE)
                    ||(y>maxY&&fabs(y-maxY)>LIBIOCTOMAP_LINEAR_TOLERANCE))
            {
                continue;
            }
        }
        else if(i==3)
        {
            if(fabs(y-minY)>LIBIOCTOMAP_LINEAR_TOLERANCE
                    ||(x<minX&&fabs(x-minX)>LIBIOCTOMAP_LINEAR_TOLERANCE)
                    ||(x>maxX&&fabs(x-maxX)>LIBIOCTOMAP_LINEAR_TOLERANCE))
            {
                continue;
            }

        }
        // Exists segments intersection
        xo=x;yo=y;
        double distanceStartToNewOrigin=sqrt(pow((xo-xs),2.0)+pow((yo-ys),2.0));
        zo=zs+slopeStartToEnd*distanceStartToNewOrigin;
        double disToNewOrigin=sqrt(pow(xo-xs,2.0)
                                   +pow(yo-ys,2.0)
                                   +pow(zo-zs,2.0));
        if(disToNewOrigin>LIBIOCTOMAP_LINEAR_TOLERANCE)
            return(true);
    }
    return(false);
}

bool libIOctoMap::createOctoMap(QVector<QVector<double> > &coordinates,
                                QVector<QVector<unsigned short> > &colors,
                                QVector<unsigned short> &sourceIds,
                                QMap<int, QVector<double> > &scanPositions,
                                float resolution,
                                bool computeFreeVoxels,
                                double &volume,
                                bool computeMetrics,
                                bool openViewer,
                                QString outputFileName,
                                QString &strError,
                                QWidget* ptrWidget)
{
    if(QFile::exists(outputFileName))
    {
        if(!QFile::remove(outputFileName))
        {
            strError="libIOctoMap::createOctoMap\n";
            strError+=QObject::tr("Error removing existing output file:\n%1")
                    .arg(outputFileName);
            return(false);
        }
    }
    // Comprobacion de las dimensiones de los contenedores
    if(coordinates.size()!=colors.size()
            ||colors.size()!=sourceIds.size())
    {
        strError="libIOctoMap::createOctoMap\n";
        strError+=QObject::tr("Different dimensions in input containers");
        return(false);
    }
    // Comprobacion de sourceIds
    if(scanPositions.size()>0)
    {
        for(int i=0;i<sourceIds.size();i++)
        {
            if(!scanPositions.contains(sourceIds[i]))
            {
                strError="libIOctoMap::createOctoMap\n";
                strError+=QObject::tr("Scan position not exists for: %1").arg(QString::number(sourceIds[i]));
                return(false);
            }
        }
    }
    QProgressDialog* ptrProgress=NULL;
    int pointsStep=LIBIOCTOMAP_NUMBER_OF_POINTS_TO_PROCESS_BY_STEP;
    int numberOfPoints=coordinates.size();
    int numberOfSteps=ceil((double)numberOfPoints/(double)pointsStep);
    if(ptrWidget!=NULL)
    {
        QString title=QObject::tr("OctoMap");
        QString msgGlobal=" Inserting rays for  ";
        msgGlobal+=QString::number(numberOfPoints,10);
        msgGlobal+=" points";
        ptrProgress=new QProgressDialog(title, "Abort",0,numberOfSteps, ptrWidget);
        ptrProgress->setWindowModality(Qt::WindowModal);
        ptrProgress->setLabelText(msgGlobal);
        ptrProgress->show();
        qApp->processEvents();
    }
    int step=0;
    int numberOfProcessedPoints=0;
    int numberOfProcessedPointsInStep=0;

    int numberOfPointsToProcessInStep=numberOfPoints;
    if(LIBIOCTOMAP_NUMBER_OF_POINTS_TO_PROCESS_BY_STEP<numberOfPointsToProcessInStep)
    {
        numberOfPointsToProcessInStep=LIBIOCTOMAP_NUMBER_OF_POINTS_TO_PROCESS_BY_STEP;
    }
    double minFc=10000000.0;
    double minSc=100000000.0;
    double minTc=10000.0;
    double maxFc=-10000000.0;
    double maxSc=-100000000.0;
    double maxTc=-10000.0;
    double centroidFc=0.0;
    double centroidSc=0.0;
    double centroidTc=0.0;
    for(int i=0;i<coordinates.size();i++)
    {
        if(coordinates[i][0]<minFc)
        {
            minFc=coordinates[i][0];
        }
        if(coordinates[i][1]<minSc)
        {
            minSc=coordinates[i][1];
        }
        if(coordinates[i][2]<minTc)
        {
            minTc=coordinates[i][2];
        }
        if(coordinates[i][0]>maxFc)
        {
            maxFc=coordinates[i][0];
        }
        if(coordinates[i][1]>maxSc)
        {
            maxSc=coordinates[i][1];
        }
        if(coordinates[i][2]>maxTc)
        {
            maxTc=coordinates[i][2];
        }
        centroidFc+=coordinates[i][0];
        centroidSc+=coordinates[i][1];
        centroidTc+=coordinates[i][2];
    }
    centroidFc/=((double)coordinates.size());
    centroidSc/=((double)coordinates.size());
    centroidTc/=((double)coordinates.size());
//    centroidFc-=minFc;
//    centroidSc-=minSc;
//    centroidTc-=minTc;
    ColorOcTree tree(resolution);
    for(int i=0;i<coordinates.size();i++)
    {
        numberOfProcessedPoints++;
        numberOfProcessedPointsInStep++;
        if(numberOfProcessedPointsInStep==numberOfPointsToProcessInStep)
        {
            step++;
            numberOfPointsToProcessInStep=numberOfPoints-numberOfProcessedPoints;
            if(LIBIOCTOMAP_NUMBER_OF_POINTS_TO_PROCESS_BY_STEP<numberOfPointsToProcessInStep)
            {
                numberOfPointsToProcessInStep=LIBIOCTOMAP_NUMBER_OF_POINTS_TO_PROCESS_BY_STEP;
            }
            numberOfProcessedPointsInStep=0;
            if(ptrWidget!=NULL)
            {
                ptrProgress->setValue(step);
                qApp->processEvents();
            }
        }
//        double endPointFc=coordinates[i][0]-minFc;
//        double endPointSc=coordinates[i][1]-minSc;
//        double endPointTc=coordinates[i][2]-minTc;
        double endPointFc=coordinates[i][0];
        double endPointSc=coordinates[i][1];
        double endPointTc=coordinates[i][2];
        int sourceId=sourceIds[i];
        double xo,yo,zo;
        bool insertRay=true;
        if(computeFreeVoxels)
        {
            double originPointFc=scanPositions[sourceId][0];
            double originPointSc=scanPositions[sourceId][1];
            double originPointTc=scanPositions[sourceId][2];
            if(fabs(endPointFc-minFc)<LIBIOCTOMAP_LINEAR_TOLERANCE
                    ||fabs(endPointFc-maxFc)<LIBIOCTOMAP_LINEAR_TOLERANCE
                    ||fabs(endPointSc-minSc)<LIBIOCTOMAP_LINEAR_TOLERANCE
                    ||fabs(endPointSc-maxSc)<LIBIOCTOMAP_LINEAR_TOLERANCE)
            {
                insertRay=false;
            }
            else
            {
                getNewOriginForRayInBoundingBox(originPointFc,originPointSc,originPointTc,
                                                endPointFc,endPointSc,endPointTc,
                                                minFc,minSc,maxFc,maxSc,
                                                xo,yo,zo);
                double disToNewOrigin=sqrt(pow(xo-originPointFc,2.0)
                                           +pow(yo-originPointSc,2.0)
                                           +pow(zo-originPointTc,2.0));
                if(disToNewOrigin<LIBIOCTOMAP_LINEAR_TOLERANCE
                        ||zo<minTc)
                {
                    QString strBoundingBox="POLYGON((";
                    strBoundingBox+=QString::number(minFc,'f',3);
                    strBoundingBox+=" ";
                    strBoundingBox+=QString::number(minSc,'f',3);
                    strBoundingBox+=",";

                    strBoundingBox+=QString::number(minFc,'f',3);
                    strBoundingBox+=" ";
                    strBoundingBox+=QString::number(maxSc,'f',3);
                    strBoundingBox+=",";

                    strBoundingBox+=QString::number(maxFc,'f',3);
                    strBoundingBox+=" ";
                    strBoundingBox+=QString::number(maxSc,'f',3);
                    strBoundingBox+=",";

                    strBoundingBox+=QString::number(maxFc,'f',3);
                    strBoundingBox+=" ";
                    strBoundingBox+=QString::number(minSc,'f',3);
                    strBoundingBox+=",";

                    strBoundingBox+=QString::number(minFc,'f',3);
                    strBoundingBox+=" ";
                    strBoundingBox+=QString::number(minSc,'f',3);

                    strBoundingBox+="))";
                    QString strRay="LINESTRING(";
                    strRay+=QString::number(originPointFc,'f',3);
                    strRay+=" ";
                    strRay+=QString::number(originPointSc,'f',3);
                    strRay+=",";

                    strRay+=QString::number(endPointFc,'f',3);
                    strRay+=" ";
                    strRay+=QString::number(endPointSc,'f',3);
                    strRay+=")";

                    QString strNewOrigin="POINT(";
                    strNewOrigin+=QString::number(xo,'f',3);
                    strNewOrigin+=" ";
                    strNewOrigin+=QString::number(yo,'f',3);
                    strNewOrigin+=")";

                    getNewOriginForRayInBoundingBox(originPointFc,originPointSc,originPointTc,
                                                    endPointFc,endPointSc,endPointTc,
                                                    minFc,minSc,maxFc,maxSc,
                                                    xo,yo,zo);

                }
            }
        }
        endPointFc-=centroidFc;
        endPointSc-=centroidSc;
        endPointTc-=centroidTc;
        int red=colors[i][0]/256;
        int green=colors[i][1]/256;
        int blue=colors[i][2]/256;
        point3d endPoint((float) endPointFc, (float) endPointSc, (float) endPointTc);
        if(computeFreeVoxels&&insertRay)
        {
            xo-=centroidFc;
            yo-=centroidSc;
            zo-=centroidTc;
            point3d originPoint((float) xo, (float) yo, (float) zo);
            tree.insertRay(originPoint,endPoint);
        }
        ColorOcTreeNode* n = tree.updateNode(endPoint, true);
        n->setColor(red,green,blue); // de 0 a 255
//        n->setColor(255,0,0); // set color to yellow
    }
    if(ptrWidget!=NULL)
    {
        ptrProgress->setValue(numberOfSteps);
        qApp->processEvents();
        ptrProgress->close();
        delete(ptrProgress);
    }
    tree.updateInnerOccupancy();
    // should already be pruned
    int nodesBeforePrune=tree.calcNumNodes();
//    EXPECT_EQ(tree.size(), tree.calcNumNodes());
//    const size_t initialSize = tree.size();
//    EXPECT_EQ(initialSize, 1034);
    tree.prune();
    int nodesAfterPrune=tree.calcNumNodes();

    int treeNumberOfLeafNodes=tree.getNumLeafNodes();
    tree.expand();
    int treeNumberOfLeafNodes2=tree.getNumLeafNodes();

//    std::vector<point3d> pcl;
//    for (ColorOcTree::iterator it = tree.begin(); it != tree.end(); ++it)
//    {
//        if(tree.isNodeOccupied(*it))
//        {
//            pcl.push_back(it.getCoordinate());
//        }
//    }

    double treeResolution=tree.getResolution();

    if(computeMetrics)
    {
        int numberOfOccupatedNodes=0;
        for (ColorOcTree::leaf_iterator it = tree.begin_leafs(),
            end = tree.end_leafs();  it != end; ++it)
        {
    //        OcTreeNode* n = tree2.search(it.getKey());
    //        if (!n)
    //        {
    //            OCTOMAP_ERROR("Could not find coordinate of 1st octree in 2nd octree\n");
    //        }
    //        else
    //        {
    //            // check occupancy prob:
    //            double p1 = it->getOccupancy();
    //        }
            if(tree.isNodeOccupied(*it))
            {
                numberOfOccupatedNodes+=1;
            }
            double p1 = it->getOccupancy();
    //        if (p1 < 0.0 || p1 > 1.0)
    //            OCTOMAP_ERROR("p1 wrong: %f", p1);

        }
        volume=((double)numberOfOccupatedNodes)*pow(resolution,3.0);
    }
//    EXPECT_EQ(tree.size(), tree.calcNumNodes());
//    EXPECT_EQ(initialSize, tree.size());
    // write color tree
    std::string stdOutputFileName=outputFileName.toStdString();
//    if(!tree.writeBinary(stdOutputFileName))
    if(!tree.write(stdOutputFileName))
    {
        strError="libIOctoMap::openViewerRGB\n";
        strError+=QObject::tr("Error writting output file:\n%1").arg(outputFileName);
        return(false);
    }
    if(!QFile::exists(outputFileName))
    {
        strError="libIOctoMap::openViewerRGB\n";
        strError+=QObject::tr("Error writting output file:\n%1").arg(outputFileName);
        return(false);
    }
    if(openViewer)
    {
        QString strAuxError;
        if(!openFileInViewer(outputFileName,
                             strAuxError))
        {
            strError="libIOctoMap::openViewerRGB\n";
            strError+=QObject::tr("Error opening viewer for output file:\n%1\nError:\n%2")
                    .arg(outputFileName).arg(strAuxError);
            return(false);
        }
    }
    return(true);
}

bool libIOctoMap::openFileInViewer(QString fileName,
                                   QString &strError)
{
    if(!QFile::exists(fileName))
    {
        strError="libIOctoMap::openFileInViewer\n";
        strError+=QObject::tr("Not exists file:\n%1").arg(fileName);
        return(false);
    }
    QDir programDir=qApp->applicationDirPath();
    QString openvisFileName=programDir.absolutePath();
    openvisFileName=openvisFileName+"/"+LIBIOCTOMAP_OPENVIS;
    if(!QFile::exists(openvisFileName))
    {
        strError=QObject::tr("libIOctoMap::openFileInViewer");
        strError+=QObject::tr("\nNot exists file:\n%1").arg(openvisFileName);
        return(false);
    }
    if(openvisFileName.contains(" "))
    {
        openvisFileName="\""+openvisFileName+"\"";
    }
    if(fileName.contains(" "))
    {
        fileName="\""+fileName+"\"";
    }
    QString command=openvisFileName+" "+fileName;

    QProcess proccess;
//	proccess.setStandardErrorFile(_fileNameStdError);
//	proccess.setWorkingDirectory(_path);
//  proccess.setStandardOutputFile(_fileNameStdOut);
//	proccess.setStandardOutputFile(metaInformationFileName);
    proccess.start(command);
    qApp->processEvents();
    if (!proccess.waitForFinished(-1))//&&!QFile::exists(metaInformationFileName))
    {
        strError=QObject::tr("libIOctoMap::openFileInViewer");
        strError+=QObject::tr("\nError in execution::\n%1").arg(command);
        proccess.close();
        return(false);
    }
    int success=proccess.exitCode();
    proccess.close();
    return(true);
}

double libIOctoMap::max(double v1, double v2)
{
    if(v1>=v2) return(v1);
    else return(v2);
}

double libIOctoMap::min(double v1, double v2)
{
    if(v1<=v2) return(v1);
    else return(v2);
}
